import numpy as np

class Simplex:
    def __get_line(filepath, cast_type, sep):
        return np.array(list(map(cast_type,
            filepath.readline().strip().split(sep))))

    def read_file(filepath, minimization=True, sep=" ", bigm=1.0e+8):
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

            # Add weight "big-M" to the cost_function_coeffs associated
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

                coeffs_and_op = line.strip().split(sep)
                operation = coeffs_and_op[-1]
                aux_coeffs = np.array(list(map(float, 
                    coeffs_and_op[:-1])))

                if operation == ">":
                    # Change the sign of both A and b coeffs to
                    # make operation transform into "<"
                    aux_coeffs *= -1.0

                elif operation == "=":
                    # In this case, change the dummy variable
                    # cost in the objective function from 0 to big_m
                    # to evoid using it during the resolution
                    cost_func_coeffs[coeff_mat_dim[0] + line_index] = bigm
                
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

        if not minimization:
            # If maximization problem, just multiply the
            # cost function by -1.0
            cost_func_coeffs[:coeff_mat_dim[0]] *= -1.0

        return cost_func_coeffs, coeff_matrix, b_vector

    def initial_feasible_sol(coeffs_matrix, use_dummies=True):
        # Assuming all original equations is in the form A * x <= b
        m, n = coeffs_matrix.shape
        basis_indexes = np.array([i for i in range(m)]) # Expecting m variables
        non_basis_indexes = np.array([i for i in range(m-1, n)]) # Expecting m-n variables
        return basis_indexes, non_basis_indexes

    def non_basis_candidate(
            cost_func_coeffs, 
            coeffs_matrix, 
            non_basis_indexes,
            lambda_vector,
            epsilon=1.0e-8):

        # "k" is the index of a non-basis variable cadidate
        # to enter the basis (it has negative cost associated, 
        # supposing we're solving a MINIMIZATION problem)
        potential_alt_sol = -1
        for non_basis_index in non_basis_indexes:
            cur_cost = cost_func_coeffs[non_basis_index] -\
                    np.dot(lambda_vector.T, coeffs_matrix[:, non_basis_index])

            if cur_cost < 0:
                # Found an negative cost associated with an
                # variable, pick it to enter the new basis
                return False, non_basis_index

            elif potential_alt_sol < 0 and abs(cur_cost) <= epsilon:
                # Null costs may imply in alternative optimum
                # solutions. If no negative cost is found, then
                # try to find another optimal solution
                potential_alt_sol = non_basis_index

        return True, potential_alt_sol

    def build_answer(cur_solution_vector, basis_indexes, size):
        new_answer = np.zeros(size)

        for vec_index, basis_index in enumerate(basis_indexes):
            new_answer[basis_index] = cur_solution_vector[vec_index]

        return new_answer

    def solve(cost_func_coeffs, coeffs_matrix, b_vector, epsilon=1.0e-9):
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
        cur_solution_vector = np.linalg.solve(\
                coeffs_matrix[:, basis_indexes], b_vector)

        while True:
            # 1. Verify if it is optimal
            # B^T * lambda = cost_base_indexes
            lambda_vector = np.linalg.solve(
                    coeffs_matrix[:, basis_indexes].T, 
                    cost_func_coeffs[basis_indexes])

            alternative, basis_in_index =\
                    Simplex.non_basis_candidate(
                        cost_func_coeffs, 
                        coeffs_matrix,
                        non_basis_indexes,
                        lambda_vector)

            if basis_in_index < 0 or alternative:
                # No non-basis variable candidate to enter the basis:
                # we're in the optimal solution.
                new_solution = Simplex.build_answer(
                        cur_solution_vector, 
                        basis_indexes,
                        coeffs_matrix.shape[1])

                solutions.append(new_solution)

                # If we're trying to find an alternative solution,
                # proceed
                if not alternative:
                    return solutions

            # 2. Calculate the "Simplex Direction"
            simplex_dir = np.linalg.solve(
                    coeffs_matrix[:, basis_indexes], 
                    coeffs_matrix[:, basis_in_index])

            # 3. Calculate the "simplex step size" while checking for
            # unbounded solutions or inexistent alternative optimum 
            # solutions
            step_size = min([cur_solution_vector[i]/simplex_dir[i] 
                for i in range(len(simplex_dir)) 
                if simplex_dir[i] > 0])

            if len(step_size) == 0:
                # There's no positive simples direction coordinates,
                # this is caused by an unbounded solution
                solutions.append(["Unbounded solution"])
                return solutions

            if step_size <= 0:
                # In this case, we're probably trying to find
                # alternatives optima solutions and failed.
                return solutions

            # Update cur_solution_vector
            cur_solution_vector = cur_solution_vector -\
                    step_size * simplex_dir.reshape(cur_solution_vector.shape)
            
            # Use the Bland's Rule to find the index
            # which will leave the base (Bland's rule picks
            # up the smallest index associated with a null
            # value in cur_solution_vector)
            basis_vector_out_index = (abs(cur_solution_vector) <= epsilon).argmax()
            non_basis_vector_in_index = (non_basis_indexes == basis_in_index).argmax()
            basis_out_index = basis_indexes[basis_vector_out_index]

            if basis_in_index == basis_out_index or\
                    (basis_in_index in basis_indexes):
                return solutions

            # Swap base and non-basis indexes
            basis_indexes[basis_vector_out_index] = basis_in_index
            non_basis_indexes[non_basis_vector_in_index] = basis_out_index

        return solutions


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("usage: " + sys.argv[0] + "<filepath>",
                "[-sep separator, default is empty space]",
                "[-bigm big_m_value (float), default is 1.0e+9]",
                "[-max]", sep="\n\t")
        exit(1)
    
    try:
        big_m_val = float(sys.argv[1 + sys.argv.index("-bigm")])
    except:
        big_m_val = 1.0e+9

    try:
        sep = sys.argv[1 + sys.argv.index("-sep")]
    except:
        sep = " "

    minimize = "-max" not in sys.argv

    f, A, b = Simplex.read_file(sys.argv[1], 
            minimization=minimize,
            bigm=big_m_val)

    print("Problem setup:", sep="\n", end=2*"\n")

    if minimize:
        print("Function coeffs to minimize:", f, sep="\n", end=2*"\n")

    else:
        print("Function coeffs to maximize:", f * -1.0, sep="\n", end=2*"\n")

    print("Coeffs matrix in normal form:", A, sep="\n", end=2*"\n")
    print("b vector coefficients (from Ax = b):", b, sep="\n")

    ans = Simplex.solve(f, A, b)

    print("Solution:\n", ans)
