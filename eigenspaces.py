import sympy as sp
import matplotlib.pyplot as plt


def kernel(A, mod):
    """Return the basis for the image and the kernel of A (working modulo mod).
    Reduction on matrix A is done via column operations.

    Parameters:
    A -- a matrix of size m x n with nonempty domain of codomain (m, n > 0)
    mod -- positive integer

    Return:
    L, U -- m x n and n x n matrices such that it holds L = A * U and L is lower triangular
    first_zero_column -- index i of the first zero column in L (None if the kernel is trivial)
                         the basis for the kernel of A: columns of U from i to n
                         the basis for the image of A: columns of L from 1 to i-1
    """
    m = len(A)
    n = len(A[0])

    L = [[A[i][j] % mod for j in range(n)] for i in range(m)]
    U = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):  # U is identity matrix at the beginning.
        U[i][i] = 1

    curr_line = 0  # Row which we currently want to zero out.
    first_zero_column = None

    for j in range(n):
        curr_col = j  # Becomes the index of first nonzero element in the current row (greater than j only).
        pivot = False  # True if we found a nonzero element.

        while True:
            if curr_col == n:
                curr_col = j
                curr_line += 1
                if curr_line == m:
                    break

            if L[curr_line][curr_col] != 0:
                pivot = True
                break
            else:
                curr_col += 1

        if not pivot:  # All the remaining elementsare zero.
            first_zero_column = j
            break
        else:
            if curr_col != j:
                # Swap columns.
                for i in range(curr_line, m):
                    L[i][j], L[i][curr_col] = L[i][curr_col], L[i][j]
                for i in range(n):
                    U[i][j], U[i][curr_col] = U[i][curr_col], U[i][j]

            inv = sp.mod_inverse(L[curr_line][j], mod)
            # Eliminate all element if the current row.
            for k in range(j + 1, n):
                coeff = L[curr_line][k] * inv
                for i in range(curr_line, m):
                    L[i][k] = (L[i][k] - coeff * L[i][j]) % mod
                for i in range(n):
                    U[i][k] = (U[i][k] - coeff * U[i][j]) % mod

    return L, U, first_zero_column


def mat_mult(A, B, mod):
    """Return the matrix C = A * B modulo mod where A and B are matrices with nonempty domains and codomains.
    Dimensions: A -- m x n, B -- n x p.
    """
    m = len(A)
    n = len(B)
    p = len(B[0])

    C = [[0 for _ in range(p)] for _ in range(m)]
    for i in range(m):
        for j in range(p):
            for k in range(n):
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % mod

    return C


def vect_mult(A, vect, mod):
    """Return the vector C = A * vect modulo mod where A is a matrix and vect is a vector (as a list).
    Dimensions: A -- m x n, vect -- length n.
    """
    n = len(vect)
    m = len(A)

    C = [0 for _ in range(m)]
    for i in range(m):
        for j in range(n):
            C[i] = (C[i] + A[i][j] * vect[j]) % mod

    return C


def mat_sum(A, B, t, mod):
    """Return the matrix A - t * B modulo mod where A and B are matrices of the same dimensions and t is a number.
    A and B have nonempty domains and codomains.
    """
    m = len(A)
    n = len(A[0])

    return [[(A[i][j] - t * B[i][j]) % mod for j in range(n)] for i in range(m)]


def intersection(A, B, mod):
    """Return the basis for the intersection of subspaces. Bases for subspaces are given as A and B which are
    lists (nonempty) of n and p vectors. All vectors have length m. Basis for the intersection is given
    as a list of vectors.
    """
    m = len(A[0])
    n = len(A)
    p = len(B)

    # Matrix of size m x (n + p), C = [A^T B^T].
    C = [[A[j][i] % mod for j in range(n)] + [B[j][i] % mod for j in range(p)] for i in range(m)]
    _, U, col = kernel(C, mod)

    intersect_basis = []
    if col is not None:
        A_mat = [line[:n] for line in C]
        for j in range(col, n + p):
            vect = [U[i][j] for i in range(n)]
            intersect_basis.append(vect_mult(A_mat, vect, mod))

    return intersect_basis


def kernel_basis(A, mod):
    """Return the basis for the kernel of A as a list of basis vectors. This function calls
    the function kernel(A, mod).
    """
    _, U, col = kernel(A, mod)
    basis = []
    if col is not None:
        for j in range(col, len(U)):
            vect = [U[i][j] for i in range(len(U))]
            basis.append(vect)

    return basis


def reduce_vect(top_list, vect, mod):
    """Return the reduced vector and its coefficients with respect to the basis in top_list, working modulo mod.

    Parameters:
    top_list -- a list of basis elements as returned in the function complete_basis
    vect -- a vector to reduce against the basis in top_list

    Return:
    vect -- reduced form of the vector
    top -- top nonzero index in vect (reduced form), if vect reduces to the zero vector then top is None
    coeffs -- a dictionary: key = index of basis vector in top_list
                            value = coefficient at this basis vector when expanding vect in terms of top_list
    """
    vect = vect.copy()
    top = None
    coeffs = dict()
    n = len(vect)
    for i in range(n):
        if (vect[i] % mod) == 0:
            continue

        if top_list[i][1] is not None:
            base = top_list[i][0]
            inv = sp.mod_inverse(base[i], mod)
            coeff = (inv * vect[i]) % mod
            coeffs[i] = coeff
            for j in range(i, n):
                vect[j] = (vect[j] - coeff * base[j]) % mod
        else:
            top = i
            break

    return vect, top, coeffs


def complete_basis(A, B, mod):
    """Return the reduced basis for space B (nonempty) and also indices of basis elements which form the basis
    for the quotient space B/A (A is a subspace of B). First we reduced the basis for A, then the basis for B
    against the reduced basis for A.

    Parameters:
    A, B -- bases for subspaces given as lists of basis vectors (of equal length p)

    Return:
    top_non_zero -- a list of basis vectors for B, top_non_zero[i] contains the pair
                   (basis vector with the top nonzero coefficient at index i, string which tells us if it belongs to A),
                    if no such vector exists then top_non_zero[i] is ([], None)
    quotient basis -- indices of vectors in top_non_zero which form the basis for B/A
    """
    p = len(B[0])
    top_non_zero = [([], None) for _ in range(p)]
    quotient_basis = []

    for vect in A:
        vect_reduced, top, _ = reduce_vect(top_non_zero, vect, mod)
        if top is not None:
            top_non_zero[top] = (vect_reduced, 'intersection')

    for i in range(len(B)):
        vect_reduced, top, _ = reduce_vect(top_non_zero, B[i], mod)
        if top is not None:  # Was not reduced to the zero vector.
            top_non_zero[top] = (vect_reduced, 'eigenspace')
            quotient_basis.append(top)

    return top_non_zero, quotient_basis


def quotient_coeffs(top_list, quotient_basis, vect, mod):
    """Return the list of coefficients when expanding vector vect against the quotient basis from top_list.

    Parameters:
    top_list -- as top_non_zero in the function complete_basis
    quotient_basis -- as in the function complete_basis
    """
    _, _, coeffs = reduce_vect(top_list, vect, mod)
    return [coeffs[i] if i in coeffs else 0 for i in quotient_basis]


def basis_conversion(vect_one, basis_one, basis_two):
    """Return vect_one represented in basis_two instead of basis_one.

    Parameters:
    vect_one -- a list of length n
    basis_one, basis_two -- lists of length n and m, each contain indices of basis elements (from some bigger space)
                            lists are ordered by increasing indices

    Return:
    vect_two -- a list of length m, vect_one[i] = vect_two[j] iff basis_one[i] = basis_two[j]
    """
    n = len(basis_one)
    m = len(basis_two)
    vect_two = [0 for _ in range(m)]

    j = 0
    for i in range(m):
        current = basis_two[i]
        while j < n and basis_one[j] < current:
            j += 1

        if j == n:
            break
        elif basis_one[j] == current:
            vect_two[i] = vect_one[j]

    return vect_two


def eigenspace_tower(inclusion_tower, mapped_tower, domain_filt_basis, max_filt_index, eigenvalue, mod):
    """Return the tower of eigenspaces for the given eigenvalue, working modulo mod.

    Parameters:
    inclusion_tower, mapped_tower -- as returned by the function towers.tower_of_pairs
    domain_filt_basis -- as returned by the function towers.filtration_basis
    max_filt_index -- as in the function rips.rips_complex
    eigenvalue -- an integer
    mod -- a positive integer

    Return:
    eigen_tower_quotient_basis -- a list of length max_filt_index, at index i is a list of basis elements
                                  for E_i which are given as lists of coefficients in terms of domain_filt_basis[i]
    eigen_tower -- a list of length max_filt_index - 1, at index i is the matrix E_i -> E_{i+1}
                   if either E_i or E_{i+1} is trivial, then eigen_tower[i] is None
    """
    # At index i is a tuple (top_non_zero, quotient_basis) as returned by the function complete_basis.
    # If E_i is trivial, then we store None instead.
    eigen_tower_basis = [None for _ in range(max_filt_index)]
    eigen_tower_quotient_basis = [[] for _ in range(max_filt_index)]
    eigen_tower = [None for _ in range(max_filt_index - 1)]

    for i in range(max_filt_index):
        A = mapped_tower[i]
        B = inclusion_tower[i]

        if len(A) == 0:  # Dimension of codomain is 0.
            dim = len(domain_filt_basis[i])  # Dimension of domain.
            if dim == 0:  # Dimension of domain is 0. There are no basis elements.
                continue
            else:  # All basis elements are in the intersection of the kernels.
                top_non_zero = [([], None) for _ in range(dim)]
                for d in range(dim):
                    vect = [0 for _ in range(dim)]
                    vect[d] = 1
                    top_non_zero[d] = (vect, 'intersection')

                eigen_tower_basis[i] = (top_non_zero, [])
        elif len(A[0]) == 0:  # Dimension of domain is 0. There are no basis elements.
            continue
        else:
            C = mat_sum(A, B, eigenvalue, mod)
            eigen_space_basis = kernel_basis(C, mod)
            if len(eigen_space_basis) == 0:
                continue

            ker_basisA = kernel_basis(A, mod)
            ker_basisB = kernel_basis(B, mod)

            if len(ker_basisA) == 0 or len(ker_basisB) == 0:
                intersect_basis = []
            else:
                intersect_basis = intersection(ker_basisA, ker_basisB, mod)

            basis_list, quotient_basis = complete_basis(intersect_basis, eigen_space_basis, mod)
            eigen_tower_basis[i] = (basis_list, quotient_basis)

            for index in quotient_basis:
                eigen_tower_quotient_basis[i].append(basis_list[index][0])

    for i in range(max_filt_index - 1):
        if eigen_tower_basis[i] is None or eigen_tower_basis[i + 1] is None:  # Domain or codomain is zero.
            continue
        else:
            indices_list_first, basis_first = eigen_tower_basis[i]
            indices_list_second, basis_second = eigen_tower_basis[i + 1]
            m = len(basis_second)
            n = len(basis_first)
            if n == 0 or m == 0:  # Quotient in domain or codomain is zero.
                continue

            eigen_matrix = [[None for _ in range(n)] for _ in range(m)]  # Matrix of dimension m x n.

            for j in range(n):
                vect_first = indices_list_first[basis_first[j]][0]  # Coefficients with respect to basis of domain.
                vect_second = basis_conversion(vect_first, domain_filt_basis[i], domain_filt_basis[i + 1])
                coeffs = quotient_coeffs(indices_list_second, basis_second, vect_second, mod)

                for k in range(m):
                    eigen_matrix[k][j] = coeffs[k]

            eigen_tower[i] = eigen_matrix

    return eigen_tower, eigen_tower_quotient_basis


def tower_normal_form(tower, tower_basis, max_filt_index, mod):
    """Transform a given tower and its basis into its normal form. Each matrix in the tower has to become a matrix
    of a matching. The function changes tower and tower_basis. We work modulo mod.

    Parameters:
    tower, tower_basis -- as returned by the function eigenspace_tower
    max_filt_index -- length of tower_basis
    mod -- a positive integer
    """
    for ind in range(max_filt_index - 2, -1, -1):
        right = tower[ind]  # Matrix of dimensions m x n.
        if right is None:  # Domain or codomain is zero.
            continue

        basis = tower_basis[ind]  # Basis for domain of matrix right and codomain of matrix left.
        left = None
        if ind > 0:  # There is a matrix preceding right.
            left = tower[ind - 1]

        m = len(right)
        n = len(basis)
        pivot = 0  # Current column.

        for i in range(m):
            column = pivot  # This will be the column with the first nonzero element in the row i.
            while column < n and right[i][column] % mod == 0:
                column += 1
            if column == n:  # Zero row.
                continue

            # Swapping columns.
            basis[pivot], basis[column] = basis[column], basis[pivot]
            for row in range(m):
                right[row][pivot], right[row][column] = right[row][column], right[row][pivot]
            # Swapping rows.
            if left is not None:
                left[pivot], left[column] = left[column], left[pivot]

            # Set the pivot element to 1 and correct both matrices and basis accordingly.
            pivot_el = right[i][pivot]
            inv = sp.mod_inverse(pivot_el, mod)
            for row in range(i, m):
                right[row][pivot] = (right[row][pivot] * inv) % mod
            for k in range(len(basis[pivot])):
                basis[pivot][k] = (basis[pivot][k] * inv) % mod
            if left is not None:
                for k in range(len(left[pivot])):
                    left[pivot][k] = (left[pivot][k] * pivot_el) % mod

            # Eliminate all other elements in the row.
            for j in range(pivot + 1, n):
                coeff = right[i][j]
                for row in range(i, m):
                    right[row][j] = (right[row][j] - coeff * right[row][pivot]) % mod
                for k in range(len(basis[j])):
                    basis[j][k] = (basis[j][k] - coeff * basis[pivot][k]) % mod
                if left is not None:
                    for k in range(len(left[pivot])):
                        left[pivot][k] = (left[pivot][k] + coeff * left[j][k]) % mod

            pivot += 1

    for ind in range(max_filt_index - 1):
        left = tower[ind]  # Matrix of dimension m x n.
        if left is None:
            continue

        basis = tower_basis[ind + 1]
        right = None
        if ind < max_filt_index - 2:
            right = tower[ind + 1]

        m = len(left)
        n = len(left[0])
        pivot = 0  # Index of the row with the current pivot element.

        for i in range(n):
            while pivot < m and left[pivot][i] % mod == 0:
                pivot += 1
            if pivot == m:  # We eliminated all the rows already.
                break

            # Eliminate all the other elements in the column.
            for j in range(pivot + 1, m):
                coeff = left[j][i]
                left[j][i] = 0
                for k in range(len(basis[pivot])):
                    basis[pivot][k] = (basis[pivot][k] + coeff * basis[j][k]) % mod
                if right is not None:
                    for k in range(len(right)):
                        right[k][pivot] = (right[k][pivot] + coeff * right[k][j]) % mod


def tower_persistence(tower, tower_basis, max_filt_index):
    """Return the persistent intervals and their persistent generators from a tower in normal form.

    Parameters:
    tower, tower_basis -- as transformed by the function tower_normal_form
    max_filt_index -- length of the tower

    Return:
    intervals -- a list of persistent intervals given as pairs (birth, length of interval)
    interval_generators -- a list of persistent generators given as lists of coefficients in terms of domain basis
    """
    intervals = []
    interval_generators = []
    # A dictionary of current active intervals.
    # Key = index of current basis element for the interval, value = index of corresponding interval in intervals.
    active_intervals = dict()

    for ind in range(max_filt_index - 1):
        M = tower[ind]
        basis = tower_basis[ind]
        if M is None:
            for i in range(len(basis)):
                if i not in active_intervals:  # This is a new interval.
                    intervals.append([ind, 0])  # The interval ends immediately (we map to or from a zero vector space).
                    interval_generators.append(basis[i])

            active_intervals = dict()
        else:
            pivot = 0
            new_active_intervals = dict()  # Update active intervals.
            for i in range(len(basis)):  # Columns of matrix M.
                while pivot < len(M) and M[pivot][i] == 0:
                    pivot += 1

                if pivot == len(M):  # All other columns are zero.
                    for j in range(i, len(basis)):
                        # The remaining columns end an existing interval or create a new interval (ends immediately).
                        if j not in active_intervals:
                            intervals.append([ind, 0])
                            interval_generators.append(basis[j])

                    break
                else:
                    if i in active_intervals:  # Interval has to be extended.
                        position = active_intervals[i]
                        intervals[position][1] += 1
                        new_active_intervals[pivot] = position
                    else:  # We create a new interval which doesn't end immediately.
                        position = len(intervals)
                        intervals.append([ind, 1])
                        interval_generators.append(basis[i])
                        new_active_intervals[pivot] = position

            active_intervals = new_active_intervals

    return intervals, interval_generators


def transform_intervals(intervals, max_filt_index):
    """Return a list of transformed intervals.

    Parameters:
    intervals -- as returned by the function tower_persistence
    max_filt_index -- maximal death time of an interval

    Return:
    intervals in the form (birth, death), where death is None instead of max_filt_index
    """
    return [(birth, birth + time) if birth + time < max_filt_index - 1 else (birth, None) for birth, time in intervals]


def generators_to_cycles(intervals, interval_generators, domain_filt_basis, domain_basis, domain_index_to_simp, mod):
    """Return persistent generators as cycles.

    Parameters:
    intervals -- as returned by the function transform_intervals
    interval_generators -- as returned by the function tower_persistence
    domain_filt_basis -- as returned by the function towers.filtration_basis
    domain_basis -- as returned by the function towers.persistence_intervals
    domein_index_to_simp -- as returned by the function towers.sparse_boundary_matrix
    mod -- positive integer

    Return:
    cycles -- a list of persistent generators, given as dictionaries:
              key = ordered simplex, value = coefficient of this simplex in the cycle
    """
    cycles = []
    for i in range(len(intervals)):
        birth = intervals[i][0]
        basis = domain_filt_basis[birth]  # A list of indices of basis elements in domain_basis.
        coeffs = interval_generators[i]  # Coefficients of the generator in this basis.
        cycle = dict()

        for j in range(len(basis)):
            base_vect = domain_basis[basis[j]]  # A dictionary: key = index of a simplex, value = its coefficient
            base_coeff = coeffs[j]  # Can be zero.

            for simp_ind in base_vect:
                simp = domain_index_to_simp[simp_ind]
                value = (base_coeff * base_vect[simp_ind]) % mod

                if simp not in cycle:
                    if value != 0:
                        cycle[simp] = value
                else:
                    value = (value + cycle[simp]) % mod
                    if value != 0:
                        cycle[simp] = value
                    else:
                        del cycle[simp]

        cycles.append(cycle)

    return cycles


def stretched_intervals(intervals, epsilon):
    """Return the intervals in the proper scale given by epsilon.

    Parameters:
    intervals -- as returned by the function transform_intervals
    epsilon -- a list of inreasing threshold values, as in the function rips.rips_complex

    Return:
    Intervals which are scaled according to epsilon.
    """
    return [(epsilon[i], j) if j is None else (epsilon[i], epsilon[j]) for i, j in intervals]


def visualize_barcode(intervals, final_value, set_ticks=False):
    """Plot the barcode of the persistent intervals. If death = None, it is replaced by final_value and plotted in red.
    Intervals are as returned by the function stretched_intervals.
    """
    inter = intervals.copy()
    for i, (birth, death) in enumerate(inter):
        if death is None:
            inter[i] = (birth, final_value)

    inter.sort(key=lambda x: (x[0], -x[1]))

    for i in range(len(inter)):
        birth, death = inter[i]
        if death == final_value:
            plt.plot(inter[i], (i + 1, i + 1), color='r', linewidth=2)
        else:
            plt.plot(inter[i], (i + 1, i + 1), color='b', linewidth=2)
        plt.scatter(final_value, i + 1, color='r', s=5)

    if set_ticks:
        plt.yticks(list(range(1, len(inter) + 1)), list(range(1, len(inter) + 1)))


def visualize_persistent_diagram(intervals, final_value):
    """Plot the persistent diagram from the list of intervals. If death = None, it is replaced by final_value and
    plotted in red.
    """
    plt.axes().set_aspect('equal')

    for birth, death in intervals:
        if death is None:
            plt.scatter(birth, final_value, color='r', s=10)
        else:
            plt.scatter(birth, death, color='b', s=10)

    plt.plot([0, final_value], [0, final_value], linewidth=1, color='black')


if __name__ == '__main__':

    A = [[2, 4, 3, 0, 0], [8, 9, 6, 1, 2], [8, 0, 10, 7, 1], [0, 0, 0, 0, 6]]
    B = [[4, 1, 3], [2, -1, 3], [2, 1, 1], [1, 1, 2]]

    U = [[1, 0, 0], [0, 1, 0], [0, 0, 3]]
    V = [[2, 1, 0], [0, 0, 1]]

    tower = [[[4, 2, 3], [4, 1, 3], [5, 6, 10]], [[1, 3, 1], [4, 6, 10], [0, 0, 0]], [[2, 7, 3]]]
    tower_basis = [[[1], [2], [1]], [[4, 3], [3, 5], [2, 1]], [[1], [6], [10]], [[1, 2, 0, 1, 4]]]
