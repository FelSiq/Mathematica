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
