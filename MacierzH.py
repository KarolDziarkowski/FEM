from Classes import Node
from ElementUniwersalny import UniversalElement
from Classes import Element
from tabulate import tabulate

no_integration_nodes = 2

n1 = Node(1, 0, 0)
n2 = Node(2, 0.025, 0)
n3 = Node(3, 0.025, 0.025)
n4 = Node(4, 0, 0.025)


elem = Element(1)
elem.addNode(n1)
elem.addNode(n2)
elem.addNode(n3)
elem.addNode(n4)

_universal = UniversalElement(no_integration_nodes)
_universal.print_ksi_array()
_universal.print_eta_array()

ksi = _universal.ksi_derivatives
eta = _universal.eta_derivatives

def dx_dksi(element: Element, no_nodes, integration_point):
    temp = UniversalElement(no_nodes)
    result = 0
    temporary = 0
    for i in range(len(element.connected_nodes)):
        temporary = element.connected_nodes[i].x * ksi[i][integration_point]
        # print("Dx_Dksi - x:", element.connected_nodes[i].x, "*", temp.ksi_derivatives[i][integration_point], " (ksi) =", temporary)
        # print("result-", result, "+", temporary)
        result += temporary

    #print(result)
    return result


def dy_dksi(element: Element, no_nodes, integration_point):
    temp = UniversalElement(no_nodes)
    result = 0
    for i in range(len(element.connected_nodes)):
        temporary = element.connected_nodes[i].y * ksi[i][integration_point]
        # print("Dy_Dksi - x:", element.connected_nodes[i].y, "*", temp.ksi_derivatives[i][integration_point], " (ksi) =",
        #       temporary)
        # print("result-", result, "+", temporary)
        result += temporary

    #    print(result)
    return result


def dx_deta(element: Element, no_nodes, integration_point):
    temp = UniversalElement(no_nodes)
    result = 0
    for i in range(len(element.connected_nodes)):
        temporary = element.connected_nodes[i].x * eta[i][integration_point]
        # print("Dx_Deta - x:", element.connected_nodes[i].x, "*", temp.eta_derivatives[i][integration_point], "(eta) =", temporary)
        # print("result-", result, "+", temporary)
        result += temporary

    #print(result)
    return result


def dy_deta(element: Element, no_nodes, integration_point):
    temp = UniversalElement(no_nodes)
    result = 0
    for i in range(len(element.connected_nodes)):
        temporary = element.connected_nodes[i].y * eta[i][integration_point]
        # print("Dx_eta - x:", element.connected_nodes[i].x, "*", temp.eta_derivatives[i][integration_point], " (ksi) =", temporary)
        # print("result-", result, "+", temporary)
        result += temporary
    #    print(result)
    return result


class JacobianMatrix:
    def __init__(self, element: Element, no_nodes, integration_point):
        self.matrix = [[0, 0], [0, 0]]

        self.matrix[0][0] = dx_dksi(element, no_nodes, integration_point)
        self.matrix[0][1] = dx_deta(element, no_nodes, integration_point)
        self.matrix[1][0] = dy_dksi(element, no_nodes, integration_point)
        self.matrix[1][1] = dy_deta(element, no_nodes, integration_point)

        self.detJ = self.matrix[0][0] * self.matrix[1][1] - self.matrix[0][1] * self.matrix[1][0]
        self.inverse_detJ = 1 / self.detJ

    def print_matrix(self):
        for row in self.matrix:
            print(row)
        print()
        print("det(J) - ", self.detJ)
        print("1/det(J) - ", self.inverse_detJ)

    def multiply_by_inverse(self):
        self.matrix[0][0] = self.matrix[0][0] * self.inverse_detJ
        self.matrix[0][1] = self.matrix[0][1] * self.inverse_detJ
        self.matrix[1][0] = self.matrix[1][0] * self.inverse_detJ
        self.matrix[1][1] = self.matrix[1][1] * self.inverse_detJ

    def get_matrix_ready_for_DNi(self):
        temp00 = self.matrix[0][0]
        temp01 = self.matrix[0][1]
        temp10 = self.matrix[1][0]
        temp11 = self.matrix[1][1]
        self.matrix[0][0] = temp11
        self.matrix[0][1] = -temp10
        self.matrix[1][0] = -temp01
        self.matrix[1][1] = temp00
        self.multiply_by_inverse()
        return self

    def get_matrix(self):
        return self.matrix


class dNi_dX:
    def __init__(self, element, no_nodes):
        self.no_nodes = no_nodes

        cols = no_nodes ** 2
        rows = 4
        self.j_matrices = [None for _ in range(cols)]
        self.j_matrices_ready = [None for _ in range(cols)]
        self.matrix = [[None for _ in range(cols)] for _ in range(rows)]

        for i in range(cols):
            temp = JacobianMatrix(element, no_nodes, i)
            self.j_matrices[i] = temp

        for i in range(cols):
            temp = self.j_matrices[i].get_matrix_ready_for_DNi()
            self.j_matrices_ready[i] = temp.get_matrix()

        for integration_point in range(cols):
            # self.print_j_matrix_ready(integration_point)
            jacobian_matrix = self.j_matrices_ready[integration_point]
            for shape_function in range(rows):
                result = (
                        jacobian_matrix[0][0] * ksi[shape_function][integration_point] +
                        jacobian_matrix[0][1] * eta[shape_function][integration_point]
                )
                self.matrix[shape_function][integration_point] = result

    def print_matrix(self):
        print("Matrix dN/dX:")
        for col in range(self.no_nodes ** 2):
            for row in range(4):
                value = self.matrix[row][col]
                print(f"{value:.6f}", end="\t")
            print()

    def print_j_matrix_ready(self, integration_point):
        print("Ready J Matrix", integration_point)
        for col in range(2):
            for row in range(2):
                value = self.j_matrices_ready[integration_point][row][col]
                print(f"{value:.6f}", end="\t")
            print()
        print("\n")


class dNi_dY:
    def __init__(self, element, no_nodes):
        self.no_nodes = no_nodes

        cols = no_nodes ** 2
        rows = 4
        self.j_matrices = []
        self.j_matrices_copy = [None for _ in range(cols)]
        self.j_matrices_ready = [None for _ in range(cols)]
        self.matrix = [[None for _ in range(cols)] for _ in range(rows)]

        for i in range(cols):
            temp1 = JacobianMatrix(element, no_nodes, i)
            temp2 = JacobianMatrix(element, no_nodes, i)
            self.j_matrices.append(temp1.get_matrix())
            self.j_matrices_copy[i] = temp2

        for i in range(cols):
            temp = self.j_matrices_copy[i].get_matrix_ready_for_DNi()
            self.j_matrices_ready[i] = temp.get_matrix()

        for integration_point in range(cols):
            jacobian_matrix = self.j_matrices_ready[integration_point]
            for shape_function in range(rows):
                result = (
                        jacobian_matrix[1][0] * ksi[shape_function][integration_point] +
                        jacobian_matrix[1][1] * eta[shape_function][integration_point]
                )
                self.matrix[shape_function][integration_point] = result

    def print_matrix(self):
        print("Matrix dN/dY:")
        for col in range(self.no_nodes ** 2):
            for row in range(4):
                value = self.matrix[row][col]
                print(f"{value:.6f}", end="\t")
            print()

    def print_j_matrix(self, integration_point):
        print("J Matrix", integration_point)
        for col in range(2):
            for row in range(2):
                value = self.j_matrices[integration_point][row][col]
                print(f"{value:.6f}", end="\t")
            print()
        print("\n")

    def print_j_matrix_ready(self, integration_point):
        print("Ready J Matrix", integration_point)
        for col in range(2):
            for row in range(2):
                value = self.j_matrices_ready[integration_point][row][col]
                print(f"{value:.6f}", end="\t")
            print()
        print("\n")


class TransposedMatrix:
    def __init__(self, elem_, no_nodes, k_value):
        self.temp_dx = dNi_dX(elem_, no_nodes)
        self.temp_dy = dNi_dY(elem_, no_nodes)
        self.no_nodes = no_nodes
        self.n_matrices = no_nodes ** 2
        self.rows = 4
        self.cols = no_nodes ** 2
        self.matricesX = []
        self.matricesY = []
        self.matricesSum = []
        self.jacobian_determinants = []
        # Loop over each integration point to calculate and store the determinant
        for i in range(no_nodes ** 2):
            temporary_jacobian = JacobianMatrix(elem_, no_nodes, i)
            self.jacobian_determinants.append(temporary_jacobian.detJ)

        for col in range(self.cols):
            matrix_x = self.calculateMatrixX(col)
            matrix_y = self.calculateMatrixY(col)
            self.matricesX.append(matrix_x)
            self.matricesY.append(matrix_y)
            self.matricesSum.append(self.add_matrices(matrix_x, matrix_y))

        for matrix, detJ in zip(self.matricesSum, self.jacobian_determinants):
            self.multiply_matrix(matrix, k_value * detJ)


    def calculateMatrixX(self, integration_point):
        temp_matrix = [[None for _ in range(4)] for _ in range(4)]
        for i in range(4):
            for j in range(4):
                temp_matrix[i][j] = self.temp_dx.matrix[i][integration_point] * self.temp_dx.matrix[j][integration_point]

        # print(f"Matrix {integration_point}")
        # for col in range(4):
        #     for row in range(4):
        #         value = temp_matrix[row][col]
        #         print(f"{value:.6f}", end="\t")
        #     print()

        return temp_matrix

    def calculateMatrixY(self, integration_point):
        temp_matrix = [[None for _ in range(4)] for _ in range(4)]
        for i in range(4):
            for j in range(4):
                temp_matrix[i][j] = self.temp_dy.matrix[i][integration_point] * self.temp_dy.matrix[j][integration_point]

        # print("Matrix")
        # for col in range(4):
        #     for row in range(4):
        #         value = temp_matrix[row][col]
        #         print(f"{value:.6f}", end="\t")
        #     print()

        return temp_matrix

    def print_matrices(self):
        print("MatricesX:")
        for matrix_index, matrix in enumerate(self.matricesX):
            print(f"MatrixX {matrix_index + 1}:")
            for row in matrix:
                for value in row:
                    print(f"{value:.6f}" if isinstance(value, (float, int)) else value, end="\t")
                print()

        print("\nMatricesY:")
        for matrix_index, matrix in enumerate(self.matricesY):
            print(f"MatrixY {matrix_index + 1}:")
            for row in matrix:
                for value in row:
                    print(f"{value:.6f}" if isinstance(value, (float, int)) else value, end="\t")
                print()

        print("\nMatricesSum:")
        for matrix_index, matrix in enumerate(self.matricesSum):
            print(f"MatrixSum {matrix_index + 1}:")
            for row in matrix:
                for value in row:
                    print(f"{value:.6f}" if isinstance(value, (float, int)) else value, end="\t")
                print()

    def add_matrices(self, matrix_a, matrix_b):
        result_matrix = [[None for _ in range(4)] for _ in range(4)]
        for i in range(4):
            for j in range(4):
                result_matrix[i][j] = matrix_a[i][j] + matrix_b[i][j]
        return result_matrix

    def multiply_matrix(self, matrix, factor):
        for i in range(4):
            for j in range(4):
                matrix[i][j] *= factor

class MatrixH:
    def __init__(self, _elem, no_nodes, k):
        self.element = _elem
        self.matrices = TransposedMatrix(_elem, no_nodes, k)
        self.k = k
        self.matrices_with_weights = []

        universal_el = UniversalElement(no_nodes)
        self.weights = []
        for weight in range(len(universal_el.weights)):
            tmp = universal_el.weights[weight]
            self.weights.append(tmp)

        # Print matrices before applying weights
        #print("MatricesSum before applying weights:")
        #for matrix_index, matrix in enumerate(self.matrices.matricesSum):
        #    print(f"MatrixSum {matrix_index + 1}:")
        #    for row in matrix:
        #        for value in row:
        #            print(f"{value:.6f}" if isinstance(value, (float, int)) else value, end="\t")
        #        print()

        # Multiply each matrix by its corresponding weight and add to matrices_with_weights
        for matrix, weight in zip(self.matrices.matricesSum, self.weights):
            weighted_matrix = [[element * weight.x * weight.y for element in row] for row in matrix]
            self.matrices_with_weights.append(weighted_matrix)

        self.total_matrix = self.calculate_total_matrix()

    def print_matrices_with_weights(self):
        print("Matrices with Weights:")
        for matrix_index, matrix in enumerate(self.matrices_with_weights):
            print(f"Matrix with Weight {matrix_index + 1}:")
            for row in matrix:
                for value in row:
                    print(f"{value:.6f}" if isinstance(value, (float, int)) else value, end="\t")
                print()

    def calculate_total_matrix(self):
        temp_matrix = [[0 for _ in range(4)] for _ in range(4)]

        # Sum up all matrices in self.matrices_with_weights
        for matrix in self.matrices_with_weights:
            for i in range(4):
                for j in range(4):
                    temp_matrix[i][j] += matrix[i][j]

        return temp_matrix

    def print_total_matrix(self):
        headers = [""] + list(range(1, len(self.total_matrix[0]) + 1))
        table = [[i + 1] + row for i, row in enumerate(self.total_matrix)]
        print("\nTotal Matrix:")
        print(tabulate(table, headers=headers, tablefmt="fancy_grid"))

    def get_matrix_h(self):
        return self.total_matrix

    def add_hbc_matrix(self, hbc_matrix):
        if len(hbc_matrix) != 4 or any(len(row) != 4 for row in hbc_matrix):
            raise ValueError("The provided matrix should be a 4x4 matrix")

        for i in range(4):
            for j in range(4):
                self.total_matrix[i][j] += hbc_matrix[i][j]

# a = TransposedMatrix(elem_, no_integration_nodes)
# # a.print_matrices()
#
# b = MatrixH(elem, no_integration_nodes, 30)
# b.matrices.print_matrices()
# b.print_matrixh()

# a = JacobianMatrix(elem, no_integration_nodes, 0)
# #a.multiply_by_inverse()
# a.print_matrix()
# print(a.determinant())

# a = dNi_dX(elem, no_integration_nodes)
# b = dNi_dY(elem, no_integration_nodes)
# a.print_matrix()
# b.print_matrix()

# a = MatrixH(elem, no_integration_nodes, 25)
# #a.print_matrices_with_weights()
# a.print_total_matrix()

