import numpy as np

class Simplex:
    def __get_line(filepath, cast_type, sep):
        return np.array(list(map(cast_type,
            filepath.readline().strip().split(sep))))

    def read_file(filepath, minimization=True, sep=" "):
        """
            Argument description:
            1. minimization: (boolean) whether the problem is 
                minimizing (true) or not (false) the cost function
        """

        cost_func_coeffs = None
        coeff_matrix = None
        b_vector = None

        with open(filepath) as f:
            # Get cost function coefficients
            cost_func_coeffs = Simplex.__get_line(f, float, sep)

            # Get A matrix dimensions
            coeff_mat_dim = Simplex.__get_line(f, int, sep)

            # Add weight zero to the cost_function_coeffs associated
            # to the dummy variables
            cost_func_coeffs = np.concatenate((cost_func_coeffs, 
                [0] * coeff_mat_dim[0]))

            # Consider problem is in the form Ax <= b, then add space
            # for dummy variables already to set the problem to the 
            # "normal form" Ax = b, x >= 0.
            coeff_matrix = np.zeros(coeff_mat_dim +
                    [0, coeff_mat_dim[0]])
            b_vector = np.zeros((coeff_mat_dim[0], 1))
            
            for line_index, line in enumerate(f):
                dummy_coef = np.zeros(coeff_mat_dim[0])
                dummy_coef[line_index] = 1

                aux_coeffs = np.array(list(map(float, 
                    line.strip().split(sep))))
                
                # Assuming the last number of a line is in the
                # right side of the equation:
                #   x_0 * aux_coeffs[0] + 
                #   x_1 * aux_coeffs[1] + 
                #   ... + 
                #   x_(len_aux_coess)-2) * aux_coeffs[len(aux_coeffs) - 2]
                #   =
                #   aux_coeffs[len(aux_coeffs) - 1]
                coeff_matrix[line_index,:] =\
                        np.concatenate((aux_coeffs[:-1], dummy_coef))
                b_vector[line_index,:] = aux_coeffs[-1]

        return cost_func_coeffs, coeff_matrix, b_vector

    def initial_feasible_sol(coeffs_matrix, use_dummies=True):
        # Assuming all original equations is in the form A * x <= b
        m, n = coeffs_matrix.shape
        basis_indexes = np.array([i for i in range(m-1, n)]) # Expecting m variables
        non_basis_indexes = np.array([i for i in range(m)]) # Expecting m-n variables
        return basis_indexes, non_basis_indexes

    def non_basis_candidate(
            cost_func_coeffs, 
            coeffs_matrix, 
            non_basis_indexes,
            lambda_vector):

        # "k" is the index of a non-basis variable cadidate
        # to enter the basis (it has negative cost associated, 
        # supposing we're solving a MINIMIZATION problem)
        k_index = -1

        for non_basis_index in non_basis_indexes:
            cur_cost = cost_func_coeffs[non_basis_index] -\
                    np.dot(lambda_vector.T, coeffs_matrix[:, non_basis_index])

            if cur_cost < 0:
                k_index = non_basis_index
                break

        return k_index

    def solve(cost_func_coeffs, coeffs_matrix, b_vector):
        """
            Suppose A * x = b, x >= 0

            Where:
                dim(A) = (m, n)
                dim(x) = (n, 1)
                dim(b) = (m, 1)

            Parameters description:
            1. cost_funct_coeffs: 
                1.0 Description: coefficients of function to MINIMIZE
                1.1 Expected data type: np.array
                1.2 Expected dimension: (1, k), k <= n (n is the # of variables)

            2. coeffs_matrix: 
                1.0 Description: coefficients of matrix A
                1.1 Expected data type: np.array
                1.2 Expected dimension: (m, n)

            3. b_vector: 
                1.0 Description: "b" vector from A * x = b equation
                1.1 Expected data type: np.array
                1.2 Expected dimension: (m, 1)
        """
        solutions = []

        # 0. Obtain an initial Feasible Initial Solution
        basis_indexes, non_basis_indexes =\
                Simplex.initial_feasible_sol(coeffs_matrix)
        x_b = np.linalg.solve(coeffs_matrix[:, basis_indexes], b_vector)

        while True:
            # 1. Verify if it is optimal
            # B^T * lambda = cost_base_indexes
            lambda_vector = np.linalg.solve(
                    coeffs_matrix[:, basis_indexes].T, 
                    cost_func_coeffs[basis_indexes])

            k_index = Simplex.non_basis_candidate(
                    cost_func_coeffs, 
                    coeffs_matrix,
                    non_basis_indexes,
                    lambda_vector)

            if k_index < 0:
                # No non-basis variable candidate to enter the basis:
                # we're in the optimal solution.
                solutions.append(x_b)
                return solutions

            # 2. Calculate the "Simplex Direction"
            simplex_dir = np.linalg.solve(
                    coeffs_matrix[:, basis_indexes], 
                    coeffs_matrix[:, k_index])

            # 3. Calculate the "simplex step size" while checking for
            # unbounded solutions or inexistent alternative optimum 
            # solutions
            step_size = min([x_b[i]/simplex_dir[i] 
                for i in range(len(simplex_dir)) 
                if simplex_dir[i] > 0])

            if len(step_size) == 0:
                # There's no positive simples direction coordinates,
                # this is caused by an unbounded solution
                solutions.append(["Unbounded solution"])
                return solutions

            if step_size < 0:
                # In this case, we're probably trying to find
                # alternatives optima solutions and failed.
                return solutions

            # Update x_b
                 
            

        return solutions


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("usage:", sys.argv[0], "<filepath>")
        exit(1)

    
    """
    First supposing:
        - Always a minimization problem
        - Always in the form Ax <= b

    Input file format:
    """

    f, A, b = Simplex.read_file(sys.argv[1])

    print("Problem setup:", sep="\n", end=2*"\n")
    print("Function coeffs to minimize:", f, sep="\n", end=2*"\n")
    print("Coeffs matrix in normal form:", A, sep="\n", end=2*"\n")
    print("b vector coefficients (from Ax = b):", b, sep="\n")

    ans = Simplex.solve(f, A, b)

    print("Solution:", ans)
