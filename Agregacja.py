from WektorP import MacierzHBC, WektorP
from MacierzH import no_integration_nodes
from tabulate import tabulate


class MacierzHGlobalna:
    def __init__(self, no_elements, no_nodes, h_matrices):
        self.no_nodes = no_nodes
        self.h_matrices = h_matrices
        self.elements = []
        self.element_IDs = [[0] * 4 for _ in range(no_elements)]
        self.h_matrix_global = [[0] * no_nodes for _ in range(no_nodes)]

        for matrix in self.h_matrices:
            self.elements.append(matrix.element)

        for i in range(no_elements):
            for j in range(4):
                self.element_IDs[i][j] = self.elements[i].connected_nodes[j].node_id

        for k in range(no_elements):
            h_matrix = self.h_matrices[k].total_matrix
            for i in range(4):
                for j in range(4):
                    global_row = self.element_IDs[k][i]
                    global_col = self.element_IDs[k][j]
                    self.h_matrix_global[global_row - 1][global_col - 1] += h_matrix[i][j]

    def print_global_matrix(self):
        headers = [""] + list(range(1, self.no_nodes + 1))
        table = [[i + 1] + row for i, row in enumerate(self.h_matrix_global)]
        print(tabulate(table, headers=headers, tablefmt="fancy_grid"))


class WektorPGlobalny:
    def __init__(self, no_elements, no_nodes, p_vectors):
        self.p_vectors = p_vectors
        self.p_vector_global = [0 for _ in range(no_nodes)]

        self.elements = []
        self.element_IDs = [[0] * 4 for _ in range(no_elements)]

        for vector in self.p_vectors:
            self.elements.append(vector.element)

        for i in range(no_elements):
            for j in range(4):
                self.element_IDs[i][j] = self.elements[i].connected_nodes[j].node_id

        for k in range(no_elements):
            p_vector = self.p_vectors[k].p_vector
            for i in range(4):
                global_index = self.element_IDs[k][i] - 1  # Adjusting for 0-based indexing
                self.p_vector_global[global_index] += p_vector[i]

    def print_global_vector(self):
        headers = ["Node ID", "P Vector Value"]
        table = [[i + 1, value] for i, value in enumerate(self.p_vector_global)]
        print(tabulate(table, headers=headers, tablefmt="fancy_grid"))


class MacierzCGlobalna:
    def __init__(self, no_elements, no_nodes, c_matrices):
        self.no_nodes = no_nodes
        self.c_matrices = c_matrices
        self.elements = []
        self.element_IDs = [[0] * 4 for _ in range(no_elements)]
        self.c_matrix_global = [[0] * no_nodes for _ in range(no_nodes)]

        for matrix in self.c_matrices:
            self.elements.append(matrix.element)

        for i in range(no_elements):
            for j in range(4):
                self.element_IDs[i][j] = self.elements[i].connected_nodes[j].node_id

        for k in range(no_elements):
            c_matrix = self.c_matrices[k].total_matrix
            for i in range(4):
                for j in range(4):
                    global_row = self.element_IDs[k][i]
                    global_col = self.element_IDs[k][j]
                    self.c_matrix_global[global_row - 1][global_col - 1] += c_matrix[i][j]

    def print_global_matrix(self):
        headers = [""] + list(range(1, self.no_nodes + 1))
        table = [[i + 1] + row for i, row in enumerate(self.c_matrix_global)]
        print(tabulate(table, headers=headers, tablefmt="fancy_grid"))

    def divide_matrix_by_dtau(self, dtau):
        for i in range(self.no_nodes):
            for j in range(self.no_nodes):
                self.c_matrix_global[i][j] /= dtau

    def multiply_matrix_by_vector(self, t0_vector):
        result_vector = [0] * self.no_nodes
        for i in range(self.no_nodes):
            for j in range(self.no_nodes):
                result_vector[i] += self.c_matrix_global[i][j] * t0_vector[j]
        return result_vector


def sum_matrices(c_matrix_total, h_matrix_total):
    # Assuming both matrices have the same size
    no_nodes = len(c_matrix_total)
    total_matrix_sum = [[0] * no_nodes for _ in range(no_nodes)]

    for i in range(no_nodes):
        for j in range(no_nodes):
            total_matrix_sum[i][j] = c_matrix_total[i][j] + h_matrix_total[i][j]

    return total_matrix_sum

def sum_vectors(c_multiplied, p_vector):
    # Assuming both vectors have the same size
    no_nodes = len(c_multiplied)
    total_vector_sum = [0] * no_nodes

    for i in range(no_nodes):
        total_vector_sum[i] = c_multiplied[i] + p_vector[i]

    return total_vector_sum