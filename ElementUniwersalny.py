from Calkowanie_Gauss import GaussianIntegral
from Classes import Node


def n1_ksi(eta):
    return -(1 / 4) * (1 - eta)


def n2_ksi(eta):
    return (1 / 4) * (1 - eta)


def n3_ksi(eta):
    return (1 / 4) * (1 + eta)


def n4_ksi(eta):
    return -(1 / 4) * (1 + eta)


def n1_eta(ksi):
    return -(1/4) * (1-ksi)


def n2_eta(ksi):
    return -(1/4) * (1+ksi)


def n3_eta(ksi):
    return (1/4) * (1+ksi)


def n4_eta(ksi):
    return (1/4) * (1-ksi)


class UniversalElement:
    def __init__(self, no_int_nodes):
        self.no_int_nodes = no_int_nodes
        self.temp = GaussianIntegral(no_int_nodes)
        rows, cols = (4, no_int_nodes * no_int_nodes)
        self.ksi_derivatives = [[None for _ in range(cols)] for _ in range(rows)]
        self.eta_derivatives = [[None for _ in range(cols)] for _ in range(rows)]
        self.integration_points = []
        self.weights = []

        iterations = 1
        for i in range(no_int_nodes):
            for j in range(no_int_nodes):
                point = Node(iterations, self.temp.nodes[j], self.temp.nodes[i])
                weight = Node(iterations, self.temp.weights[j], self.temp.weights[i])
                self.integration_points.append(point)
                self.weights.append(weight)
                iterations += 1

        for i in range(cols):
            self.eta_derivatives[0][i] = n1_eta(self.integration_points[i].x)
            self.ksi_derivatives[0][i] = n1_ksi(self.integration_points[i].y)
            self.eta_derivatives[1][i] = n2_eta(self.integration_points[i].x)
            self.ksi_derivatives[1][i] = n2_ksi(self.integration_points[i].y)
            self.eta_derivatives[2][i] = n3_eta(self.integration_points[i].x)
            self.ksi_derivatives[2][i] = n3_ksi(self.integration_points[i].y)
            self.eta_derivatives[3][i] = n4_eta(self.integration_points[i].x)
            self.ksi_derivatives[3][i] = n4_ksi(self.integration_points[i].y)

    def print_ksi_array(self):
        print("Ksi:")
        for col in range(self.no_int_nodes ** 2):
            for row in range(4):
                value = self.ksi_derivatives[row][col]
                print(f"{value:.6f}", end="\t")
            print()

    def print_eta_array(self):
        print("Eta:")
        for col in range(self.no_int_nodes ** 2):
            for row in range(4):
                value = self.eta_derivatives[row][col]
                print(f"{value:.6f}", end="\t")
            print()

    def print_integration_points(self):
        print("Integration Points: ")
        for i in range(len(self.integration_points)):
            value = self.integration_points[i]
            value.printNode()
        print()

