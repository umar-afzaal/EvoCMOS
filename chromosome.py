import networkx as nx
import matplotlib.pyplot as plt
from PIL import Image
import pydot
from io import BytesIO
import itertools
import sys
from subprocess import check_output
from random import randint, randrange, uniform
from os import remove
# system wide graphviz is required for rendering


class Chromosome():    
    def set_params(self, ni, no, nc, nr):
        '''
        - param: <int> ni: no. of graph inputs
        - param: <int> no: no. of graph outputs
        - param: <int> nc: no. of graph columns
        - param: <int> nr: no. of graph rows
        '''

        self.a = 2 # Arity
        self.Ln = nc * nr # Total nodes
        self.M = self.Ln + ni # Max no. of addresses on the graph
        self.Lg = (nc * nr*(self.a + 1)) + no # Integers in a chromosome
        self.ni = ni
        self.cases = 2**(ni-2) # -2 to ignore VCC and GND
        self.maxHamming = self.cases*no
        self.no = no
        self.nc = nc
        self.l = nc # levels-back
        self.nr = nr
        self.nn = self.a + 1 # Integers per node/ Size of each node
        self.module = 'ckt'
        self.inv = False
        self.nf = ['n', 'p', 'j', 'i']
        if 'm' in self.nf:
            self.mux = True
        else:
            self.mux = False


    def __init__(self):
        pass
        

    def encode(self, nc, nr, target):
        with open(target, 'r') as f:
            netlist = f.readlines()
        ni=0
        no=0
        outputs = []
        inputs = []
        for line in netlist:
            if line[0]!='/':
                if 'output ' in line:
                    no+=1
                    outputs.append( line.split('output ')[1].split(';')[0] )
                elif 'input ' in line:
                    ni+=1
                    inputs.append( line.split('input ')[1].split(';')[0] )
                elif 'supply' in line:
                    ni += 1

        intermediate = []
        for line in netlist:
            if line[0]!='/':          
                if 'mos' in line:
                    wires = line.split('(')[1].split(')')[0].split(',')
                    wires = [wire.strip() for wire in wires]
                    intermediate+=[wire for wire in wires if wire not in inputs and wire not in outputs and wire!='GND' and wire!='VDD']
        intermediate = list(set(intermediate))
        allWires = inputs + intermediate + outputs

        junctions = []
        for wire in intermediate+outputs:
            flag = 0
            for line in netlist:
                if line[0]!='/' and 'mos' in line:
                    wires = line.split('(')[1].split(')')[0].split(',')
                    wires = [wire.strip() for wire in wires]
                    if wires[0]==wire:
                        flag += 1
            if flag>1:
                junctions.append(wire)
        junctions = dict.fromkeys(junctions, [])

        # Assign addresses: start with inputs
        nf = ['n', 'p', 'j', 'i', 't', 'm']
        dicAddress = {'GND':0, 'VDD':1}
        index = 2
        for i in range(2, ni):
            dicAddress[inputs[i-2]] = i
            index += 1

        G = []
        blockedLines = []
        blockedJunctions = []
        addressPointer = ni

        def junction_recursion(wire, addressList, addressPointer):
            if len(addressList)==1:
                return (addressPointer)

            newAddressList = []
            if len(addressList)%2==0: # even
                for k in range(len(addressList)):
                    if k%2==0:
                        G.append([[i for i,x in enumerate(nf) if x == 'j'][0], addressList[k], addressList[k+1]])
                        dicAddress[wire] = addressPointer
                        newAddressList.append(addressPointer)
                        addressPointer += 1           
            else:
                for k in range(len(addressList)):
                    if k != len(addressList)-1 and k%2==0:
                        G.append([[i for i,x in enumerate(nf) if x == 'j'][0], addressList[k], addressList[k+1]])
                        dicAddress[wire] = addressPointer
                        newAddressList.append(addressPointer)
                        addressPointer += 1
                newAddressList.append(addressList[-1])

            addressPointer = junction_recursion(wire, newAddressList, addressPointer)
            return (addressPointer)

        while(True):   
            for line in netlist:
                if line[0]!='/' and 'mos' in line and line not in blockedLines:
                    func = line.split('mos')[0].strip()
                    wires = line.split('(')[1].split(')')[0].split(',')
                    wires = [wire.strip() for wire in wires]
                    try:
                        G.append([[i for i,x in enumerate(nf) if x == func][0], dicAddress[wires[1]], dicAddress[wires[2]]])
                        if wires[0] in junctions.keys():
                            junctions[wires[0]] = junctions[wires[0]] + [addressPointer]
                        else:
                            dicAddress[wires[0]] = addressPointer  
                        addressPointer += 1
                        blockedLines.append(line)
                    except:
                        pass
            for wire, addressList in junctions.items():
                if wire not in blockedJunctions:
                    if len(addressList)>1:
                        addressPointer = junction_recursion(wire, addressList, addressPointer)
                        if wire not in [set(netlist)-set(blockedLines)]:
                            blockedJunctions.append(wire)
            if len(blockedLines) == sum(1 for line in netlist if 'mos' in line):
                break

        Ln = nc*nr
        chromosome = list(G)
        if len(G)<Ln:
            toFill = Ln-len(G)
            while (len(chromosome)<Ln):
                for i in range(toFill):
                    chromosome.append( [0,0,0] )
                    if len(chromosome)==Ln:
                        break
                    if len(G)==(i+1):
                        break
        elif len(G)>Ln:
            print(f'Not enough nodes available. At least {len(G)} nodes are required.')
            sys.exit(0)

        chromosome = [item for ls in chromosome for item in ls]

        for k,v in dicAddress.items():
            if k in outputs:
                chromosome.append(v)

        # EOF
        self.chromosome = chromosome
        self.set_params(ni, no, nc, nr)
        self.outputs = [x for x in self.chromosome[self.Lg-self.no:]]
        self.inputs = [x for x in range(self.ni)]
        self.io = self.inputs + self.outputs
        #self.io = [int(x.replace('n','')) for x in io]
        self.mux = False



    def build_network(self):
        # Find active nodes
        nodes = self.chromosome[:self.Ln*self.nn]
        nodes = [nodes[x:x+self.nn] for x in range(0, len(nodes),self.nn)]
        activeNodesAddresses = []
        for output in self.outputs:
            activeNodesAddresses.append(output)
            trail = [output]
            while True:
                lesserTrail = []
                for address in trail:
                    node = nodes[address-self.ni]
                    if node[0]<0: # check if MUX
                        for gene in node:
                            if abs(gene) not in self.inputs:
                                lesserTrail.append(abs(gene))
                                activeNodesAddresses.append(abs(gene))
                    elif self.nf[node[0]]=='i': # check if inverter
                        if node[1] not in self.inputs:
                            lesserTrail.append(node[1])
                            activeNodesAddresses.append(node[1])
                    else: # its a regular 2-input gate
                        for gene in node[1:]:
                            if gene not in self.inputs:
                                lesserTrail.append(gene)
                                activeNodesAddresses.append(gene)
                if lesserTrail==[]:
                    break
                else:
                    trail = list(lesserTrail)
        activeNodesAddresses = list(set(activeNodesAddresses))
        activeNodesAddresses = self.inputs + activeNodesAddresses

        # Define a multilayered architecture
        extents = nx.utils.pairwise(itertools.accumulate( (0,) + tuple([self.ni] + self.nc*[self.nr]) ))
        layers = [range(start, end) for start, end in extents]
        G = nx.MultiDiGraph()

        for (layerIndex, layer) in enumerate(layers):
            for address in layer:
                if address in activeNodesAddresses: # active nodes
                    if layerIndex==0:
                        thisNode = (address, {'address':address, 'outgoing':[], 'color':'turquoise'})
                        edges = []
                    else:
                        if address not in self.outputs:
                            color = 'turquoise'
                        else:
                            color = 'lightcoral'
                        if nodes[address-self.ni][0]<0: # if mux
                            thisNode = (address, {'address':address, 'outgoing':[], 'function':'m', 'incoming':[abs(nodes[address-self.ni][0]), nodes[address-self.ni][1], nodes[address-self.ni][2]], 'color':color})
                            edges = [(thisNode[1]['incoming'][0], address), (thisNode[1]['incoming'][1], address), (thisNode[1]['incoming'][2], address)]
                        elif self.nf[nodes[address-self.ni][0]]=='i': # if inverter
                            thisNode = (address, {'address':address, 'outgoing':[], 'function':'i', 'incoming':[nodes[address-self.ni][1]], 'color':color})
                            edges = [(thisNode[1]['incoming'][0], address)]
                        else: # else regular 2-input
                            thisNode = (address, {'address':address, 'outgoing':[], 'function':self.nf[nodes[address-self.ni][0]], 'incoming':[nodes[address-self.ni][1], nodes[address-self.ni][2]], 'color':color})
                            edges = [(thisNode[1]['incoming'][0], address), (thisNode[1]['incoming'][1], address)]
                    G.add_nodes_from([thisNode], layer=layerIndex)
                    G.add_edges_from(edges)
                else: # inactive nodes
                    color = 'lavender'
                    if nodes[address-self.ni][0]<0: # if mux
                        thisNode = (address, {'address':address, 'outgoing':[], 'function':'m', 'incoming':[abs(nodes[address-self.ni][0]), nodes[address-self.ni][1], nodes[address-self.ni][2]], 'color':color})
                        edges = [(thisNode[1]['incoming'][0], address), (thisNode[1]['incoming'][1], address), (thisNode[1]['incoming'][2], address)]
                    elif self.nf[nodes[address-self.ni][0]]=='i': # if inverter
                        thisNode = (address, {'address':address, 'outgoing':[], 'function':'i', 'incoming':[nodes[address-self.ni][1]], 'color':color})
                        edges = [(thisNode[1]['incoming'][0], address)]
                    else: # else regular 2-input
                        thisNode = (address, {'address':address, 'outgoing':[], 'function':self.nf[nodes[address-self.ni][0]], 'incoming':[nodes[address-self.ni][1], nodes[address-self.ni][2]], 'color':color})
                        edges = [(thisNode[1]['incoming'][0], address), (thisNode[1]['incoming'][1], address)]
                    G.add_nodes_from([thisNode], layer=layerIndex)
                    G.add_edges_from([])      
        layers = tuple([self.ni] + self.nc*[self.nr])

        # Set outgoing attribute
        for node in range(self.M):
            if node in activeNodesAddresses:
                nx.set_node_attributes(G, {node:{'outgoing':[x[1] for x in list(G.edges([node]))]}})

        # EOF
        self.G = G
        self.nodes = nodes
        self.activeNodesAddresses = activeNodesAddresses
        



    def draw_nx(self):
        # Get all attributes: <dict> function = {nodeName:attrValue, nodeName:attrValue...}
        function = nx.get_node_attributes(self.G, 'function')
        incoming = nx.get_node_attributes(self.G, 'incoming')
        colorMap = nx.get_node_attributes(self.G, 'color').values()
        pos = nx.multipartite_layout(self.G, subset_key='layer')
        nx.draw(self.G, pos, node_color=colorMap, with_labels=True, arrows=True, width=0.1, connectionstyle='arc3, rad = 0.1') # connectionstyle='arc3,rad=-0.25'  |||  connectionstyle='arc3, rad = 0.1'

        # add text on plot
        for address in range(self.M):
            x,y = pos[address]
            if address==0:
                plt.text(x,y+0.05,s=f'GND', horizontalalignment='center', fontweight='bold')
            elif address==1:
                plt.text(x,y+0.05,s=f'VDD', horizontalalignment='center', fontweight='bold')
            elif address in range(2, self.ni):
                plt.text(x,y+0.05,s=f'in{address}', horizontalalignment='center', fontweight='bold')
            else:
                if function[address]=='m':
                    plt.text(x,y+0.05,s=f'm {incoming[address][0]} {incoming[address][1]} {incoming[address][2]}', horizontalalignment='center')
                elif function[address]=='i':
                    plt.text(x,y+0.05,s=f'i {incoming[address][0]}', horizontalalignment='center')
                else:
                    plt.text(x,y+0.05,s=f'{function[address]} {incoming[address][0]} {incoming[address][1]}', horizontalalignment='center')
        plt.show()


    def draw_dot(self):

        pydotG = nx.drawing.nx_pydot.to_pydot(self.G) # from networkx graph to pydot graph

        function = nx.get_node_attributes(self.G, 'function')
        incoming = nx.get_node_attributes(self.G, 'incoming')
        
        # add text on plot
        for address in range(self.M):
            if address==0:
                pydotG.add_node(pydot.Node(address, label=f'{address}\nGND'))
            elif address==1:
                pydotG.add_node(pydot.Node(address, label=f'{address}\nVDD'))
            elif address in range(2, self.ni):
                pydotG.add_node(pydot.Node(address, label=f'{address}\nin{address}'))
            else:
                if function[address]=='m':
                    pydotG.add_node(pydot.Node(address, label=f'{address}\nm {incoming[address][0]} {incoming[address][1]} {incoming[address][2]}'))
                elif function[address]=='i':
                    pydotG.add_node(pydot.Node(address, label=f'{address}\ni {incoming[address][0]}'))
                else:
                    pydotG.add_node(pydot.Node(address, label=f'{address}\n{function[address]} {incoming[address][0]} {incoming[address][1]}'))
                    
        Image.open(BytesIO(pydotG.create_png())).show()



    def decode(self):
        
        netlist = []
        self.outputs = []
        self.inputs = []

        # Outputs
        self.outputs = [x for x in self.chromosome[self.Lg-self.no:]]
        for item in self.outputs:
            netlist.append(f'output n{item};\n')

        # Inputs
        self.inputs = [x for x in range(self.ni)]
        for item in self.inputs[2:]:
            netlist.append(f'input n{item};\n')

        # Supply voltages
        netlist.append('supply0 GND;\n')
        netlist.append('supply1 VDD;\n')

        # Build a networkx graph of this chromosome
        self.build_network()

        # Find clusters of junctions
        junctions = [k for k,v in nx.get_node_attributes(self.G, 'function').items() if v=='j'] # get all junction addresses
        junctions = [junction for junction in junctions if self.G.nodes[junction]['color']!='lavender'] # filter for active junctions

        def related_junctions_recursion(outgoing, connected):
            for address in outgoing:
                if self.G.nodes[address]['function']=='j':
                    connected.append(address)
                    related_junctions_recursion(self.G.nodes[address]['outgoing'], connected)
            return (connected)

        blockedJunctions = []
        lesserConnectedJunctions = []
        for junction in junctions:
            if junction not in blockedJunctions:
                connected = [junction]
                connected = related_junctions_recursion(self.G.nodes[junction]['outgoing'], connected)
                connected = sorted(connected, key = lambda x:x)
                lesserConnectedJunctions.append(connected)
                blockedJunctions += connected
        lesserConnectedJunctions = sorted(lesserConnectedJunctions, key = lambda x:x[0])

        blockedCluster = []
        connectedJunctions = []
        for i in range(len(lesserConnectedJunctions)):
            if lesserConnectedJunctions[i] not in blockedCluster:
                connected = list(lesserConnectedJunctions[i])
                for j in range(i+1, len(lesserConnectedJunctions)):
                    if lesserConnectedJunctions[j] not in blockedCluster:
                        if len( set(connected) & set(lesserConnectedJunctions[j]) ) != 0:
                            connected += lesserConnectedJunctions[j]
                            blockedCluster.append(lesserConnectedJunctions[j])
                connectedJunctions.append(list(set(connected)))
            connected = []

        # Solve the problem of junctions
        dicReplace = {}
        tail = []
        for cluster in connectedJunctions:
            incoming = []
            outgoing = []
            for junction in cluster:
                incoming += self.G.nodes[junction]['incoming']
                outgoing += self.G.nodes[junction]['outgoing']
            incoming = list(set(incoming)-set(cluster))
            outgoing = list(set(outgoing)-set(cluster))
            
            poMembersOfThisCluster = list(set(self.outputs) & set(cluster))
            if poMembersOfThisCluster!=[]: # check if one or more junctions in this cluster are PO
                inputsIncomingToThisCluster = list(set(incoming) & set(self.inputs))
                if inputsIncomingToThisCluster!=[]: # Check if at least one of the incoming edges are a PI
                    successor = inputsIncomingToThisCluster[0]
                    if successor==0: successor = 'GND'
                    elif successor==1: successor = 'VDD'
                    else: successor = f'n{successor}'
                    for address in poMembersOfThisCluster: # Take care of outputs
                        tail.append(f'assign n{address} = {successor};\n') # Casting a wire
                        dicReplace[address] = successor
                    nonPoMembersOfThisCluster = list(set(cluster)-set(poMembersOfThisCluster)) # Rest of the cluster members
                    for address in nonPoMembersOfThisCluster:
                        dicReplace[address] = successor
                    incomingOtherThanInputs = list(set(incoming)-set(inputsIncomingToThisCluster)) # Now take care of incoming other than PIs
                    for address in incomingOtherThanInputs:
                        if address in self.outputs: # If this incoming is a PO
                            tail.append(f'assign n{address} = {successor};\n')
                            dicReplace[address] = successor
                        else:
                            dicReplace[address] = successor      
                else: # No edges are from PIs
                    successor = f'n{[x for x in incoming if x not in self.outputs][0]}'
                    for address in poMembersOfThisCluster: # Take care of outputs
                        tail.append(f'assign n{address} = {successor};\n')
                        dicReplace[address] = successor
                    nonPoMembersOfThisCluster = list(set(cluster)-set(poMembersOfThisCluster)) # Rest of the cluster members
                    for address in nonPoMembersOfThisCluster:
                        dicReplace[address] = successor
                    for address in incoming:
                        if address in self.outputs: # If this incoming is a PO
                            tail.append(f'assign n{address} = {successor};\n')
                            dicReplace[address] = successor
                        else:
                            dicReplace[address] = successor      
            else:
                inputsIncomingToThisCluster = list(set(incoming) & set(self.inputs))
                if inputsIncomingToThisCluster!=[]: # Check if at least one of the incoming edges are a PI
                    successor = inputsIncomingToThisCluster[0]
                    if successor==0: successor = 'GND'
                    elif successor==1: successor = 'VDD'
                    else: successor = f'n{successor}'
                    for address in cluster:
                        dicReplace[address] = successor
                    incomingOtherThanInputs = list(set(incoming)-set(inputsIncomingToThisCluster)) # Now take care of incoming other than PIs
                    for address in incomingOtherThanInputs:
                        if address in self.outputs: # If this incoming is a PO
                            tail.append(f'assign n{address} = {successor};\n')
                            dicReplace[address] = successor
                        else:
                            dicReplace[address] = successor
                else:
                    successor = f'n{[x for x in incoming if x not in self.outputs][0]}'
                    for address in cluster:
                        dicReplace[address] = successor
                    for address in incoming:
                        if address in self.outputs: # If this incoming is a PO
                            tail.append(f'assign n{address} = {successor};\n')
                            dicReplace[address] = successor
                        else:
                            dicReplace[address] = successor
                    
        nonJunctions = [k for k,v in nx.get_node_attributes(self.G, 'function').items() if v!='j'] # get all non-junction addresses
        nonJunctions = self.inputs + [junction for junction in nonJunctions if self.G.nodes[junction]['color']!='lavender'] # filter for active non-junctions
        nonJunctions = [x for x in nonJunctions if x not in dicReplace.keys()]
        for address in nonJunctions:
            if address==0:
                dicReplace[address] = 'GND'
            elif address==1:
                dicReplace[address] = 'VDD'
            else:
                dicReplace[address] = f'n{address}'

        # Decode to netlist
        for address, node in enumerate(self.nodes):
            address += self.ni
            if address in self.activeNodesAddresses:
                fn = node[0]
                A = node[1]
                B = node[2]
                if fn<0: # its a mux
                    fn = abs(fn)
                    netlist.append(f'pmos({dicReplace[address]}, {dicReplace[A]}, {dicReplace[fn]});\n')
                    netlist.append(f'nmos({dicReplace[address]}, {dicReplace[B]}, {dicReplace[fn]});\n')
                    
                    #netlist.append(f'pmos({dicReplace[address]}, {dicReplace[fn]}, {dicReplace[A]});\n')
                    #netlist.append(f'nmos({dicReplace[address]}, {dicReplace[fn]}, {dicReplace[B]});\n')
                    
                elif self.nf[fn]=='i': # its an inverter
                    netlist.append(f'pmos({dicReplace[address]}, VDD, {dicReplace[A]});\n')
                    netlist.append(f'nmos({dicReplace[address]}, GND, {dicReplace[A]});\n')
                elif self.nf[fn]=='n': # its nmos
                    netlist.append(f'nmos({dicReplace[address]}, {dicReplace[A]}, {dicReplace[B]});\n')
                elif self.nf[fn]=='p': # its pmos
                    netlist.append(f'pmos({dicReplace[address]}, {dicReplace[A]}, {dicReplace[B]});\n')
        netlist += tail
        netlist.append('endmodule')

        # Module definition should come at the top
        self.io = self.outputs + self.inputs[2:]
        preamble = []
        preamble.append(f'module {self.module}(\n')
        for i in range(len(self.io)):
            if (i == len(self.io)-1):
                preamble.append(f'n{self.io[i]}'+'\n);\n')
            else:
                preamble.append(f'n{self.io[i]}'+',\n')
        netlist = preamble + netlist

        # EOF
        self.netlist = list(dict.fromkeys(netlist))
        self.count_transistors()



    def read_truth_table(self, fileName):
        with open(fileName, 'r') as f:
            truthTable = f.readlines()
        self.simulation = {}
        for line in truthTable:
            if '.i' in line:
                self.ni = int(line.split('.i')[1].strip())
            elif '.o' in line:
                self.no = int(line.split('.o')[1].strip())
            elif '.e' not in line:
                self.simulation[ line.split(' ')[0] ] = line.split(' ')[1].strip()
        self.ni += 2
        self.cases = 2**(self.ni-2)
        self.maxHamming = self.cases * self.no
        self.calc_fitness(self.simulation)


    def count_transistors(self):
        # Gate counts
        self.nmos = 0
        self.pmos = 0
        for line in self.netlist:
            if 'nmos' in line:
                self.nmos += 1
            elif 'pmos' in line:
                self.pmos += 1                    
        self.gates = self.nmos + self.pmos


    def read_netlist(self, fileName):
        with open(fileName, 'r') as f:
            netlist = f.readlines()
        self.ni = 0
        self.no = 0
        inputs = []
        outputs = []
        for line in netlist:
            if not '/' in line and not '*' in line:
                if 'input' in line:
                    self.ni += 1
                    inputs.append( line.split('input ')[1].split(';')[0] )
                elif 'output' in line:
                    self.no += 1
                    outputs.append( line.split('output ')[1].split(';')[0] )
                elif 'supply' in line:
                    self.ni += 1

        ## Modify module name
        module = fileName.split('.')[0]
        for i in range(len(netlist)):
            if not '/' in netlist[i] and not '*' in netlist[i] and module in netlist[i]:
                netlist[i] = netlist[i].replace(module, 'ckt')
                break

        # EOF
        self.inputs = [0, 1] + [int(x.replace('n','')) for x in inputs]
        self.outputs = [int(x.replace('n','')) for x in outputs]
        self.io = self.inputs[2:] + self.outputs
        self.cases = 2**(self.ni-2)
        self.maxHamming = self.cases * self.no
        self.netlist = netlist



    def calc_fitness(self, goldenSim):
        '''
        - param: <dict> goldenSim: reference simulation data {'case':'output', 'case':'output', ...}
        '''
        
        def normalize(val,Rmin,Rmax,Tmin,Tmax):
            return (((val-Rmin)/(Rmax-Rmin)*(Tmax-Tmin))+Tmin)

        vectors = 0
        hamming = 0
        for key in goldenSim:
            if self.simulation[key] != goldenSim[key]:
                vectors += 1
            hamming += sum(c1 != c2 for c1, c2 in zip(self.simulation[key], goldenSim[key]))
                
        self.er = vectors/self.cases
        self.hd = hamming/self.maxHamming

        self.fitness = normalize(self.er,0.0,1.0,0.0,0.5) + normalize(self.hd,0.0,1.0,0.0,0.5)


    def simulate(self, delFiles=True):
        '''
        - Simulates netlist contained in netlist using iVerilog simulator
        - param: <bool> delFiles: delete simulation files y/n
        '''

        verilogFile = self.module+'.v'
        testbenchFile = 'tb_'+verilogFile
    
        tb = self.testbench()

        with open(verilogFile, 'w') as f:
            for line in self.netlist:
                f.write(line)

        with open(testbenchFile, 'w') as f:
            f.write(tb)

        # call icarus verilog
        check_output("iverilog -o dsn {} {}".format(verilogFile, testbenchFile), shell=True)
        console = check_output("vvp dsn", shell=True)
        remove("dsn")
        remove(self.module+".vcd")
        if delFiles==True:
            remove(verilogFile)
            remove(testbenchFile)

        # collect simulation output
        console = str(console.decode('utf-8')).split('\n')[1:-1]
        simulation = {}
        for line in console:
            simulation[line.split(' ')[0]] = line.split(' ')[1]

        self.simulation = dict(simulation)

        return (simulation)



    def testbench(self):
        '''
        - param: <str> tb: testbench
        '''
        
        tb = ''
        tb += f"module tb();\ninitial begin\n$dumpfile(\"{self.module}.vcd\");\n$dumpvars(0, {self.module});\nend\nreg "

        for i in range(len(self.inputs[2:])):
            if (i == len(self.inputs[2:])-1):
                tb += f'n{self.inputs[2:][i]}' + ';\nwire '
            else:
                tb += f'n{self.inputs[2:][i]}' + ', '

        for i in range(len(self.outputs)):
            if (i == len(self.outputs)-1):
                tb += f'n{self.outputs[i]}' + f";\ninteger int;\n{self.module} {self.module}(\n"
            else:
                tb += f'n{self.outputs[i]}' + ', '

        for i in range(len(self.io)):
            if (i == len(self.io)-1):
                tb += f".n{self.io[i]}(n{self.io[i]})\n);\n"
            else:
                tb += f".n{self.io[i]}(n{self.io[i]}),"

        tb += f"initial\nbegin\nfor (int=0; int<{2**len(self.inputs[2:])}; int=int+1) begin\n{{"

        for i in range(len(self.inputs[2:])):
            if (i == len(self.inputs[2:])-1):
                tb += f'n{self.inputs[2:][i]}'+"} = int;\n$write("
            else:
                tb += f'n{self.inputs[2:][i]}'+','

        for i in range(len(self.inputs[2:])):
            if (i == len(self.inputs[2:])-1):
                tb += f'n{self.inputs[2:][i]}'+', " ");\n#10;\n$write('
            else:
                tb += f'n{self.inputs[2:][i]}'+', '

        for i in range(len(self.outputs)):
            if (i == len(self.outputs)-1):
                tb += f'n{self.outputs[i]}'+", \"\\n\");\n"
            else:
                tb += f'n{self.outputs[i]}'+', '       
        tb += "end\n$finish;\nend\nendmodule"

        return (tb)


    def coin(self):
        if uniform(0,1) < uniform(0,1):
            return True
        else:
            return False



    def get_random_chromosome(self):
        '''
        - param <list> chromosome: circuit encoded as a list of integers
        - param: <list> netlist: decoded circuit as a list of strings
        '''

        while True:
            chromosome = []
            for j in range(self.nc):
                for k in range(self.nr):
                    # function bit
                    if self.coin() and self.mux: # its a MUX
                        if j>=self.l:
                            bit0 = randint(self.ni+(j-self.l)*self.nr, self.ni+j*self.nr-1)
                        else:
                            bit0 = randint(0, self.ni+j*self.nr-1)
                        bit0 = -bit0
                    else: # Just take from the gate library
                        bit0 = randint(0,len(self.nf)-1)
                    
                    # connection bit 1
                    if j>=self.l:
                        bit1 = randint(self.ni+(j-self.l)*self.nr, self.ni+j*self.nr-1)
                    else:
                        bit1 = randint(0, self.ni+j*self.nr-1)                          
                        
                    # connection bit 2
                    if j>=self.l:
                        bit2 = randint(self.ni+(j-self.l)*self.nr, self.ni+j*self.nr-1)
                    else:
                        bit2 = randint(0, self.ni+j*self.nr-1) 

                    # appending to chromosome
                    chromosome.append(bit0)
                    chromosome.append(bit1)
                    chromosome.append(bit2)

            for k in range(self.no):
                # output bit
                chromosome.append(randint(self.ni, self.ni+self.Ln-1))


            try:
                self.chromosome = chromosome
                self.decode()
                self.simulate()
                return None
            except:
                pass



    def point_mutation(self, chromosome, h=1):
        '''
        - param: <list> chromosome: circuit encoded as a list of integers
        - param: <list> netlist: decoded circuit as a list of strings
        '''
        
        child = list(chromosome)
        while True:
            for k in range(h):
                # generate a random index in range length(chromosome) or Lg
                iBit = randrange(self.Lg)

                if iBit >= (self.Lg-self.no):  # its an output bit
                    while (True):
                        toReplace = randint(self.ni, self.ni+self.Ln-1)
                        if not(chromosome[iBit] == toReplace):
                            chromosome[iBit] = toReplace
                            break

                elif iBit % 3 == 0 or iBit == '0': # its a function bit
                    if self.coin() and self.mux: # its a MUX
                        for index in range(self.nc):
                            if iBit in range(self.nr*self.nn*index, (self.nr*self.nn*(index+1))):
                                j = index

                        if j>=self.l:
                            while (True):
                                toReplace = randint(self.ni+(j-self.l)*self.nr, self.ni+j*self.nr-1)
                                if not(chromosome[iBit] == toReplace):
                                    chromosome[iBit] = -toReplace
                                    break
                        else:
                            while (True):
                                toReplace = randint(0, self.ni+j*self.nr-1)
                                if not(chromosome[iBit] == toReplace):
                                    chromosome[iBit] = -toReplace
                                    break
                    else:
                        while (True):
                            toReplace = randint(0,len(self.nf)-1)
                            if not(chromosome[iBit] == toReplace):
                                chromosome[iBit] = toReplace
                                break

                else: # else its a connection bit
                    for index in range(self.nc):
                        if iBit in range(self.nr*self.nn*index, (self.nr*self.nn*(index+1))):
                            j = index

                    if j>=self.l:
                        while (True):
                            toReplace = randint(self.ni+(j-self.l)*self.nr, self.ni+j*self.nr-1)
                            if not(chromosome[iBit] == toReplace):
                                chromosome[iBit] = toReplace
                                break
                    else:
                        while (True):
                            toReplace = randint(0, self.ni+j*self.nr-1)
                            if not(chromosome[iBit] == toReplace):
                                chromosome[iBit] = toReplace
                                break

            try:
                self.chromosome = chromosome
                self.decode()
                self.simulate()
                return None
            except:
                chromosome = list(child)




 



##if __name__ == '__main__':
##
##    nc = 30
##    nr = 2
##    target = 'compressor_4_2.v'
##    chromo = Chromosome()
##    chromo.encode(nc, nr, target)
##    #chromo.decode()
##    #chromo.draw_dot()
##    for x in range(1000):
##        chromo.randomChromosome()
##        #print(chromo.chromosome)
##    print('completed')







    
