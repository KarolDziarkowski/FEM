def gaussian_elimination(matrix_c_h_summed, vectors_summed):
    # Combine the matrices and vectors into augmented matrices
    augmented_matrix = [row + [val] for row, val in zip(matrix_c_h_summed, vectors_summed)]

    # Ensure all rows in augmented_matrices have the same length
    max_row_length = max(len(row) for row in augmented_matrix)
    augmented_matrix = [row + [0] * (max_row_length - len(row)) for row in augmented_matrix]

    # Perform Gaussian elimination on the augmented matrix
    n = len(augmented_matrix)

    for i in range(n):
        # Find the pivot row
        pivot_row = max(range(i, n), key=lambda k: abs(augmented_matrix[k][i]) if isinstance(augmented_matrix[k][i], (int, float)) else 0)

        # Swap the current row with the pivot row
        augmented_matrix[i], augmented_matrix[pivot_row] = augmented_matrix[pivot_row], augmented_matrix[i]

        # Make the diagonal element 1
        pivot = augmented_matrix[i][i]
        augmented_matrix[i] = [x / pivot if isinstance(x, (int, float)) else x for x in augmented_matrix[i]]

        # Eliminate other rows
        for j in range(n):
            if i != j:
                factor = augmented_matrix[j][i]
                augmented_matrix[j] = [x - factor * y if isinstance(x, (int, float)) else x for x, y in zip(augmented_matrix[j], augmented_matrix[i])]

    # Extract the solution (temperature vector t)
    temp_solution = [row[-1] for row in augmented_matrix]

    return temp_solution
