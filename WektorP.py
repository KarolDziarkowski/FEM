from Classes import Node, Element
from ElementUniwersalny import UniversalElement
from MacierzH import no_integration_nodes
from tabulate import tabulate


n1 = Node(1, 0, 0, 1)
n2 = Node(2, 0.025, 0, 1)
n3 = Node(3, 0.025, 0.025, 1)
n4 = Node(4, 0, 0.025, 1)

elem = Element(1)
elem.addNode(n1)
elem.addNode(n2)
elem.addNode(n3)
elem.addNode(n4)

_universal = UniversalElement(no_integration_nodes)
#_universal.print_integration_points()

# Function to print a 2D matrix
def print_matrix(matrix, name):
    print(f"{name}:")
    for row in matrix:
        print(row)
    print()

def N1(ksi, eta):
    result = 0.25 * (1-ksi) * (1-eta)
    return result

def N2(ksi, eta):
    result = 0.25 * (1+ksi) * (1-eta)
    return result

def N3(ksi, eta):
    result = 0.25 * (1+ksi) * (1+eta)
    return result

def N4(ksi, eta):
    result = 0.25 * (1-ksi) * (1+eta)
    return result


def calculate_distance(node1: Node, node2: Node) -> float:
    x1 = node1.x
    y1 = node1.y

    x2 = node2.x
    y2 = node2.y

    distance = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
    return distance


def powierzchnie_bc(element: Element):
    punkty_z_bc = [0,  # point 0 (First point in element ->  Bottom left)
                   0,  # point 1 (Node with MaxID - 1  ->  Bottom right)
                   0,  # point 2 (Node with MinID      ->  Top right)
                   0]  # point 3 (Node with MinID + 1  ->  Top left)

    powierzchnie_z_bc = [0,  # wall 0 (Between Node 0 and Node 1 ->  Bottom wall)
                         0,  # wall 1 (Between Node 1 and Node 2 ->  Right wall)
                         0,  # wall 2 (Between Node 2 and Node 3 ->  Top wall)
                         0]  # wall 3 (Between Node 3 and Node 0 -> Left wall)


    if element.connected_nodes[0].BC > 0:                        # point 0 bottom left
        punkty_z_bc[0] = element.connected_nodes[0].BC

    if element.connected_nodes[1].BC > 0:                        # point 1
        punkty_z_bc[1] = element.connected_nodes[1].BC

    if element.connected_nodes[2].BC > 0:                         # point 2
        punkty_z_bc[2] = element.connected_nodes[2].BC

    if element.connected_nodes[3].BC > 0:                         # point 3
        punkty_z_bc[3] = element.connected_nodes[3].BC

    if punkty_z_bc[0] > 0 and punkty_z_bc[0] == punkty_z_bc[1]:     # wall 0 *bottom* (node 0,1)
        powierzchnie_z_bc[0] = punkty_z_bc[0]

    if punkty_z_bc[1] > 0 and punkty_z_bc[1] == punkty_z_bc[2]:     # wall 1 *right* (node 1,2)
        powierzchnie_z_bc[1] = punkty_z_bc[1]

    if punkty_z_bc[2] > 0 and punkty_z_bc[2] == punkty_z_bc[3]:     # wall 2 *top* (node 2,3)
        powierzchnie_z_bc[2] = punkty_z_bc[2]

    if punkty_z_bc[3] > 0 and punkty_z_bc[3] == punkty_z_bc[0]:     # wall 3 *left* (node 3,0)
        powierzchnie_z_bc[3] = punkty_z_bc[3]

    #print(powierzchnie_z_bc)
    return powierzchnie_z_bc


class MacierzHBC:
    def __init__(self, element, no_int_nodes, alfa):
        self.element = element
        self.no_int_nodes = no_int_nodes
        self.sciany_z_bc = powierzchnie_bc(element)
        self.punkty_bc0 = []     #integration points - bottom wall
        self.punkty_bc1 = []     #integration points - right wall
        self.punkty_bc2 = []     #integration points - top wall
        self.punkty_bc3 = []     #integration points - left wall
        self.hbc_side_matrices = []
        self.hbc_matrix = [[0 for _ in range(4)] for _ in range(4)]

        self.weights = []
        for weight in range(no_int_nodes):
            tmp = _universal.weights[weight].x
            self.weights.append(tmp)

        # for weight in self.weights:
        #     print(weight)

        x_cords = []
        y_cords = []

        for i in range(no_int_nodes ** 2):
            y_cords.append(_universal.integration_points[i].y)
            x_cords.append(_universal.integration_points[i].x)

        x_cords = list(set(x_cords))        # unique values
        y_cords = list(set(y_cords))

        if self.sciany_z_bc[0] > 0:     #if bottom wall has bc
            for i in range(no_int_nodes):
                temp = Node(i, x_cords[i], -1, self.sciany_z_bc[0])
                self.punkty_bc0.append(temp)

        if self.sciany_z_bc[1] > 0:     #if right wall has bc
            for i in range(no_int_nodes):
                temp = Node(i, 1, y_cords[i], self.sciany_z_bc[1])
                self.punkty_bc1.append(temp)

        if self.sciany_z_bc[2] > 0:  # if top wall has bc
            for i in range(no_int_nodes):
                # Reverse the x_cords list for the top wall
                temp = Node(i, x_cords[no_int_nodes - 1 - i], 1, self.sciany_z_bc[2])
                self.punkty_bc2.append(temp)

        if self.sciany_z_bc[3] > 0:  # if left wall has bc
            for i in range(no_int_nodes):
                temp = Node(i, -1, y_cords[no_int_nodes - 1 - i], self.sciany_z_bc[3])
                self.punkty_bc3.append(temp)

        for i in range(4):
            self.hbc_side_matrices.append(self.calculate_surface(i, alfa))

        self.calculate_hbc()

    def print_matrix_hbc(self):
        for i in range(len(self.hbc_side_matrices)):
            print("Matrix HBC", i+1)
            for col in range(4):
                for row in range(4):
                    value = self.hbc_side_matrices[i][row][col]
                    print(f"{value:.6f}", end="\t")
                print()

    def print_integration_points(self):
        # Print points
        print("Points on the Bottom Wall (bc0):")
        for point in self.punkty_bc0:
            point.printNode()

        print("\nPoints on the Right Wall (bc1):")
        for point in self.punkty_bc1:
            point.printNode()

        print("\nPoints on the Top Wall (bc2):")
        for point in self.punkty_bc2:
            point.printNode()

        print("\nPoints on the Left Wall (bc3):")
        for point in self.punkty_bc3:
            point.printNode()

    def outer_product(self, vec1, vec2):
        # Calculate the outer product of two vectors
        return [[vec1[i] * vec2[j] for j in range(len(vec2))] for i in range(len(vec1))]

    def calculate_surface(self, nr_boku, conductivity):
        hbc = [[0 for j in range(4)] for i in range(self.no_int_nodes)]
        matrixHBC = [[0 for j in range(4)] for i in range(4)]

        if self.sciany_z_bc[nr_boku] > 0:
            if nr_boku == 0:
                detJ = calculate_distance(self.element.connected_nodes[0], self.element.connected_nodes[1]) * 0.5
                #print("DET J 0 - ", detJ)
                for i in range(len(hbc)):
                    hbc[i][0] = N1(self.punkty_bc0[i].x, self.punkty_bc0[i].y)
                    #print(self.punkty_bc0[i].x, self.punkty_bc0[i].y)
                    hbc[i][1] = N2(self.punkty_bc0[i].x, self.punkty_bc0[i].y)
                    hbc[i][2] = N3(self.punkty_bc0[i].x, self.punkty_bc0[i].y)
                    hbc[i][3] = N4(self.punkty_bc0[i].x, self.punkty_bc0[i].y)
                    #print(hbc[i])

            elif nr_boku == 1:
                detJ = calculate_distance(self.element.connected_nodes[1], self.element.connected_nodes[2]) * 0.5
                #print("DET J 1 - ", detJ)
                for i in range(len(hbc)):
                    #print(self.punkty_bc1[i].x, self.punkty_bc1[i].y)
                    hbc[i][0] = N1(self.punkty_bc1[i].x, self.punkty_bc1[i].y)
                    hbc[i][1] = N2(self.punkty_bc1[i].x, self.punkty_bc1[i].y)
                    hbc[i][2] = N3(self.punkty_bc1[i].x, self.punkty_bc1[i].y)
                    hbc[i][3] = N4(self.punkty_bc1[i].x, self.punkty_bc1[i].y)

            elif nr_boku == 2:
                detJ = calculate_distance(self.element.connected_nodes[2], self.element.connected_nodes[3]) * 0.5
                #print("DET J 2 - ", detJ)
                for i in range(len(hbc)):
                    #print(self.punkty_bc2[i].x, self.punkty_bc2[i].y)
                    hbc[i][0] = N1(self.punkty_bc2[i].x, self.punkty_bc2[i].y)
                    hbc[i][1] = N2(self.punkty_bc2[i].x, self.punkty_bc2[i].y)
                    hbc[i][2] = N3(self.punkty_bc2[i].x, self.punkty_bc2[i].y)
                    hbc[i][3] = N4(self.punkty_bc2[i].x, self.punkty_bc2[i].y)

            elif nr_boku == 3:
                detJ = calculate_distance(self.element.connected_nodes[3], self.element.connected_nodes[0]) * 0.5
                #print("DET J 3 - ", detJ)
                for i in range(len(hbc)):
                    #print(self.punkty_bc3[i].x, self.punkty_bc3[i].y)
                    hbc[i][0] = N1(self.punkty_bc3[i].x, self.punkty_bc3[i].y)
                    hbc[i][1] = N2(self.punkty_bc3[i].x, self.punkty_bc3[i].y)
                    hbc[i][2] = N3(self.punkty_bc3[i].x, self.punkty_bc3[i].y)
                    hbc[i][3] = N4(self.punkty_bc3[i].x, self.punkty_bc3[i].y)

            # print(f"\n=================== BOK {nr_boku} ======================== \n")
            # Calculate the outer product of each row with itself transposed
            for row_index, row in enumerate(hbc):
                multiplied_transposed = self.outer_product(row, row)

                # print_matrix(hbc, f"Original hbc Matrix (bok {nr_boku})")
                # print_matrix(multiplied_transposed, f"Step 1 - Outer Product (Row {row_index + 1})")

                # Scale the outer product
                for mt_row_index, mt_row in enumerate(multiplied_transposed):
                    for mt_value_index, mt_value in enumerate(mt_row):
                        scaled_value = conductivity * self.weights[row_index] * mt_value
                        # print(f"Scaling: {conductivity} * {self.weights[row_index]} * {mt_value} = {scaled_value}")
                        multiplied_transposed[mt_row_index][mt_value_index] = scaled_value

                # print_matrix(multiplied_transposed, f"Step 2 - Scaled Outer Product (Row {row_index + 1})")

                # Accumulate scaled outer product in matrixHBC
                for hbc_row_index, hbc_row in enumerate(matrixHBC):
                    for hbc_value_index, hbc_value in enumerate(hbc_row):
                        updated_value = matrixHBC[hbc_row_index][hbc_value_index] + \
                                        multiplied_transposed[hbc_row_index][hbc_value_index]
                        # print(
                        #     f"Accumulating in matrixHBC[{hbc_row_index}][{hbc_value_index}]: {matrixHBC[hbc_row_index][hbc_value_index]} + {multiplied_transposed[hbc_row_index][hbc_value_index]} = {updated_value}")
                        matrixHBC[hbc_row_index][hbc_value_index] = updated_value

                # print_matrix(matrixHBC, f"Step 3 - Updated matrixHBC (Row {row_index + 1})")
                # print("=" * 30)

            # Multiply matrixHBC by detJ after accumulation
            for hbc_row_index, hbc_row in enumerate(matrixHBC):
                for hbc_value_index, hbc_value in enumerate(hbc_row):
                    matrixHBC[hbc_row_index][hbc_value_index] *= detJ

            #print_matrix(matrixHBC, f"Final matrixHBC (After multiplying by detJ - {detJ})")

        return matrixHBC

    def calculate_hbc(self):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    self.hbc_matrix[i][j] += self.hbc_side_matrices[k][i][j]

    def print_total_matrix(self):
        print("\nTotal Matrix HBC:")
        for row in self.hbc_matrix:
            for value in row:
                print(f"{value:.6f}" if isinstance(value, (float, int)) else value, end="\t")
            print()


class WektorP:
    def __init__(self, element, no_int_nodes, alfa, ambient_temp):
        self.element = element
        self.no_int_nodes = no_int_nodes
        self.alfa = alfa
        self.ambient_temp = ambient_temp
        self.sciany_z_bc = powierzchnie_bc(element)
        self.punkty_bc0 = []     #integration points - bottom wall
        self.punkty_bc1 = []     #integration points - right wall
        self.punkty_bc2 = []     #integration points - top wall
        self.punkty_bc3 = []     #integration points - left wall
        self.p_side_vectors = []
        self.p_vector = [0 for _ in range(4)]

        self.weights = []
        for weight in range(no_int_nodes):
            tmp = _universal.weights[weight].x
            self.weights.append(tmp)

        x_cords = []
        y_cords = []

        for i in range(no_int_nodes ** 2):
            y_cords.append(_universal.integration_points[i].y)
            x_cords.append(_universal.integration_points[i].x)

        x_cords = list(set(x_cords))        # unique values
        y_cords = list(set(y_cords))

        if self.sciany_z_bc[0] > 0:     #if bottom wall has bc
            for i in range(no_int_nodes):
                temp = Node(i, x_cords[i], -1, self.sciany_z_bc[0])
                self.punkty_bc0.append(temp)

        if self.sciany_z_bc[1] > 0:     #if right wall has bc
            for i in range(no_int_nodes):
                temp = Node(i, 1, y_cords[i], self.sciany_z_bc[1])
                self.punkty_bc1.append(temp)

        if self.sciany_z_bc[2] > 0:  # if top wall has bc
            for i in range(no_int_nodes):
                # Reverse the x_cords list for the top wall
                temp = Node(i, x_cords[no_int_nodes - 1 - i], 1, self.sciany_z_bc[2])
                self.punkty_bc2.append(temp)

        if self.sciany_z_bc[3] > 0:  # if left wall has bc
            for i in range(no_int_nodes):
                temp = Node(i, -1, y_cords[no_int_nodes - 1 - i], self.sciany_z_bc[3])
                self.punkty_bc3.append(temp)

        for i in range(4):
            self.p_side_vectors.append(self.calculate_surface(i))

        self.calculate_p_vector()

    def print_p_vectors(self):
        for i in range(len(self.p_side_vectors)):
            print("Vector P", i+1)
            for row in range(4):
                value = self.p_side_vectors[i][row]
                print(f"{value:.6f}", end="\t")
            print()

    def print_integration_points(self):
        # Print points
        print("Points on the Bottom Wall (bc0):")
        for point in self.punkty_bc0:
            point.printNode()

        print("\nPoints on the Right Wall (bc1):")
        for point in self.punkty_bc1:
            point.printNode()

        print("\nPoints on the Top Wall (bc2):")
        for point in self.punkty_bc2:
            point.printNode()

        print("\nPoints on the Left Wall (bc3):")
        for point in self.punkty_bc3:
            point.printNode()

    def calculate_surface(self, nr_boku):
        n_funcs = [[0 for j in range(4)] for i in range(self.no_int_nodes)]
        p_vector = [0 for j in range(4)]

        if self.sciany_z_bc[nr_boku] > 0:
            if nr_boku == 0:
                detJ = calculate_distance(self.element.connected_nodes[0], self.element.connected_nodes[1]) * 0.5

                for i in range(len(n_funcs)):
                    n_funcs[i][0] = N1(self.punkty_bc0[i].x, self.punkty_bc0[i].y)
                    n_funcs[i][1] = N2(self.punkty_bc0[i].x, self.punkty_bc0[i].y)
                    n_funcs[i][2] = N3(self.punkty_bc0[i].x, self.punkty_bc0[i].y)
                    n_funcs[i][3] = N4(self.punkty_bc0[i].x, self.punkty_bc0[i].y)


            elif nr_boku == 1:
                detJ = calculate_distance(self.element.connected_nodes[1], self.element.connected_nodes[2]) * 0.5

                for i in range(len(n_funcs)):
                    n_funcs[i][0] = N1(self.punkty_bc1[i].x, self.punkty_bc1[i].y)
                    n_funcs[i][1] = N2(self.punkty_bc1[i].x, self.punkty_bc1[i].y)
                    n_funcs[i][2] = N3(self.punkty_bc1[i].x, self.punkty_bc1[i].y)
                    n_funcs[i][3] = N4(self.punkty_bc1[i].x, self.punkty_bc1[i].y)


            elif nr_boku == 2:
                detJ = calculate_distance(self.element.connected_nodes[2], self.element.connected_nodes[3]) * 0.5

                for i in range(len(n_funcs)):
                    n_funcs[i][0] = N1(self.punkty_bc2[i].x, self.punkty_bc2[i].y)
                    n_funcs[i][1] = N2(self.punkty_bc2[i].x, self.punkty_bc2[i].y)
                    n_funcs[i][2] = N3(self.punkty_bc2[i].x, self.punkty_bc2[i].y)
                    n_funcs[i][3] = N4(self.punkty_bc2[i].x, self.punkty_bc2[i].y)


            elif nr_boku == 3:
                detJ = calculate_distance(self.element.connected_nodes[3], self.element.connected_nodes[0]) * 0.5

                for i in range(len(n_funcs)):
                    n_funcs[i][0] = N1(self.punkty_bc3[i].x, self.punkty_bc3[i].y)
                    n_funcs[i][1] = N2(self.punkty_bc3[i].x, self.punkty_bc3[i].y)
                    n_funcs[i][2] = N3(self.punkty_bc3[i].x, self.punkty_bc3[i].y)
                    n_funcs[i][3] = N4(self.punkty_bc3[i].x, self.punkty_bc3[i].y)


            for i in range(self.no_int_nodes):
                for j in range(4):
                    temp = n_funcs[i][j] * self.weights[i] * self.ambient_temp
                    p_vector[j] += temp






            for i in range(4):
                temp = p_vector[i] * self.alfa * detJ
                p_vector[i] = temp


        return p_vector

    def calculate_p_vector(self):
        for i in range(4):
            for j in range(4):
                    self.p_vector[i] += self.p_side_vectors[j][i]

    def print_total_vector(self):
        print("\nTotal P Vector:")
        for value in self.p_vector:
            print(f"{value:.3f}" if isinstance(value, (float, int)) else value, end="\t")


a = MacierzHBC(elem, no_integration_nodes, 25)
#a.print_integration_points()
#a.print_matrix_hbc()
#a.print_total_matrix()

a = WektorP(elem, no_integration_nodes, 300, 1200)
#a.print_p_vectors()
#a.print_total_vector()