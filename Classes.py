class Node:
    def __init__(self, node_id, x, y, bc = 0):
        self.node_id = node_id
        self.x = x
        self.y = y
        self.BC = bc

    def printNode(self):
        bc_info = f"BC {self.BC}" if self.BC > 0 else "Without BC"
        print(f"NodeID: {self.node_id}, x: {self.x}, y: {self.y}, {bc_info}")


class Element:
    def __init__(self, id):
        self.id = id
        self.connected_nodes = []

    def addNode(self, node: Node):
        self.connected_nodes.append(node)

    def printElement(self):
        print("Element ID: ", self.id, end="")
        print(", Nodes: ", end="")
        for node in self.connected_nodes:
            print(f"{node.node_id} ({node.x}, {node.y})", end=", ")
        print()


class Global:
    def __init__(self, simTime, simStepTime, conductivity, alfa, tot, initialTemp, density, specificHeat, nodesNo, elementsNo):
        self.simTime = simTime
        self.simStepTime = simStepTime
        self.conductivity = conductivity
        self.alfa = alfa
        self.tot = tot
        self.initialTemp = initialTemp
        self.density = density
        self.specificHeat = specificHeat
        self.nodesNo = nodesNo
        self.elementsNo = elementsNo

    def print_values(self):
        print("SimulationTime:", self.simTime)
        print("SimulationStepTime:", self.simStepTime)
        print("Conductivity:", self.conductivity)
        print("Alfa:", self.alfa)
        print("Tot:", self.tot)
        print("InitialTemp:", self.initialTemp)
        print("Density:", self.density)
        print("SpecificHeat:", self.specificHeat)
        print("Nodes number:", self.nodesNo)
        print("Elements number:", self.elementsNo)


class Grid:
    def __init__(self, nNodes, nElements):
        self.nNodes = nNodes
        self.nElements = nElements
        self.nodes = []
        self.elements = []

    def addNode(self, node):
        if len(self.nodes) > self.nNodes:
            print("Too many nodes")
        else:
            self.nodes.append(node)

    def addElement(self, element):
        if len(self.elements) > self.nElements:
            print("Too many elements")
        else:
            self.elements.append(element)

    def printGrid(self):
        print("Nodes:")
        count = 1
        for node in self.nodes:
            node.printNode()
            count += 1

        print("Elements:")
        for element in self.elements:
            element.printElement()

