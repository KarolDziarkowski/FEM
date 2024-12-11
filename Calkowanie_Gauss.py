import math


def functionX(x):
    return 5 * x ** 2 + 3 * x + 6


def functionXY(x, y):
    return 5 * x ** 2 * y ** 2 + 3 * x * y + 6


class GaussianIntegral2Nodes:
    def __init__(self):
        self.nodes = [-(1/math.sqrt(3)), 1/math.sqrt(3)]
        self.weight = 1

    def Integrate1D(self):
        result = functionX(self.nodes[0])*self.weight + functionX(self.nodes[1])*self.weight
        print(f"Integration result - 2Nodes (1D): {result}")
        return result

    def Integrate2D(self):
        result = 0
        for i in range(2):
            for j in range(2):
                result += functionXY(self.nodes[i], self.nodes[j]) * self.weight * self.weight

        print(f"Integration result - 2Nodes (2D): {result}")
        return result


class GaussianIntegral3Nodes:
    def __init__(self):
        self.nodes = [-(math.sqrt(3/5)), 0, math.sqrt(3/5)]
        self.weights = [5/9, 8/9, 5/9]

    def integrate1d(self):
        result = 0
        for i in range(3):
            result += functionX(self.nodes[i]) * self.weights[i]
        print(f"Integration result - 3Nodes (1D): {result}")
        return result

    def integrate2d(self):
        result = 0
        for i in range(3):
            for j in range(3):
                result += functionXY(self.nodes[i], self.nodes[j]) * self.weights[i] * self.weights[j]
        print(f"Integration result - 3Nodes (2D): {result}")
        return result


class GaussianIntegral:
    def __init__(self, no_nodes):
        self.no_nodes = no_nodes

        if no_nodes == 1:
            self.nodes = [0]
            self.weights = [2]

        elif no_nodes == 2:
            self.nodes = [-(1 / math.sqrt(3)), 1 / math.sqrt(3)]
            self.weights = [1, 1]

        elif no_nodes == 3:
            self.nodes = [-(math.sqrt(3 / 5)), 0, math.sqrt(3 / 5)]
            self.weights = [5 / 9, 8 / 9, 5 / 9]

        elif no_nodes == 4:
            self.nodes = [-(math.sqrt(3 / 7 + 2 / 7 * math.sqrt(6 / 5))),
                          -(math.sqrt(3 / 7 - 2 / 7 * math.sqrt(6 / 5))),
                          (math.sqrt(3 / 7 - 2 / 7 * math.sqrt(6 / 5))),
                          (math.sqrt(3 / 7 + 2 / 7 * math.sqrt(6 / 5)))]
            self.weights = [(18 - math.sqrt(30)) / 36, (18 + math.sqrt(30)) / 36,
                            (18 + math.sqrt(30)) / 36, (18 - math.sqrt(30)) / 36]

        elif no_nodes == 5:
            self.nodes = [-math.sqrt(5 + 2 * math.sqrt(10 / 7)) / 3,
                          -math.sqrt(5 - 2 * math.sqrt(10 / 7)) / 3,
                          0,
                          math.sqrt(5 - 2 * math.sqrt(10 / 7)) / 3,
                          math.sqrt(5 + 2 * math.sqrt(10 / 7)) / 3]

            self.weights = [(322 - 13 * math.sqrt(70)) / 900,
                            (322 + 13 * math.sqrt(70)) / 900,
                            128 / 225,
                            (322 + 13 * math.sqrt(70)) / 900,
                            (322 - 13 * math.sqrt(70)) / 900]

        else:
            raise ValueError("Unsupported number of nodes.")

    def integrate1d(self):
        result = 0
        for i in range(self.no_nodes):
            result += functionX(self.nodes[i]) * self.weights[i]
        print(f"Integration result - {self.no_nodes} Nodes (1D): {result}")
        return result

    def integrate2d(self):
        result = 0
        for i in range(self.no_nodes):
            for j in range(self.no_nodes):
                result += functionXY(self.nodes[i], self.nodes[j]) * self.weights[i] * self.weights[j]
        print(f"Integration result - {self.no_nodes} Nodes (2D): {result}")
        return result


# for i in range(1, 6):
#     quadrature = GaussianIntegral(i)
#     quadrature.integrate1d()
#     quadrature.integrate2d()

