import rips
import sympy as sp
import copy


def add_columns(col1, col2, d, mod):
    """Add d-times col2 to col1, return new column as a result. All operations are done modulo mod.

    Parameters:
    col1, col2 -- columns, given as lists of pairs (row index, value in row), ordered by decreasing row index
    d, mod -- integers

    Return:
    new_col -- column ordered by decreasing row index
    """
    new_col = []
    i = 0
    j = 0

    while i < len(col1) and j < len(col2):
        ind1, v1 = col1[i]
        ind2, v2 = col2[j]
        if ind1 > ind2:
            new_col.append((ind1, v1))
            i += 1
        elif ind2 > ind1:
            new_col.append((ind2, (d * v2) % mod))
            j += 1
        else:
            new_val = (v1 + d * v2) % mod
            i += 1
            j += 1
            if new_val != 0:
                new_col.append((ind1, new_val))

    if i < len(col1):
        new_col += col1[i:]
    else:
        remainder = col2[j:]
        mult_column(remainder, d, mod)
        new_col += remainder
        
    return new_col


def mult_column(col, d, mod):
    """Multiply column col by d (working modulo mod). The function changes column col.

    Parameters:
    col -- columns, given as a list of pairs (row index, value in row), ordered by decreasing row index
    d, mod -- integers
    """
    for i in range(len(col)):
        ind, value = col[i]
        col[i] = (ind, (d * value) % mod)


def add_chains(chain1, chain2, d, mod):
    """Add d-times chain2 to chain1 (working modulo mod). The function changes chain1.

    Parameters:
    chain1, chain2 -- chains, given as dictionaries: key = index of a simplex (with non-zero coefficient),
                      value = coefficient of the simplex in the chain
    d, mod -- integers
    """
    for key in chain2:
        if key not in chain1:
            chain1[key] = (d * chain2[key]) % mod
        else:
            x = (chain1[key] + d * chain2[key]) % mod
            if x == 0:
                del chain1[key]
            else:
                chain1[key] = x


def mult_chain(chain, d, mod):
    """Multiply chain by d (working modulo mod). The function changes chain.

    Parameters:
    chain -- given as dictionary: key = index of a simplex (with non-zero coefficient),
             value = coefficient of the simplex in the chain
    d, mod -- integers
    """
    for key in chain:
        chain[key] = (d * chain[key]) % mod


def sparse_boundary_matrix(filtration):
    """Return sparse boundary matrix of the filtration. All the simplices in the filtration are at most 3D.

    Parameters:
    filtration -- a list of pairs (simplex, minimal filtration index) ordered by increasing minimal filtration index

    Return:
    matrix -- sparse boundary matrix, each column is given as a list of pairs (row index, value in row) for non-zero
              rows only, ordered by decreasing row index
    simp_to_index_dict -- a dictionary: key = simplex,
                          value = index of the column of matrix which belongs to the simplex
    index_to_simp_dict -- a dictionary: key = index of a column of matrix
                          value = simplex whose boundary is in this column
    """
    simp_to_index_dict = dict()
    index_to_simp_dict = dict()
    matrix = []
    
    for i in range(len(filtration)):
        _, simp = filtration[i]
        simp_to_index_dict[simp] = i
        index_to_simp_dict[i] = simp

        if len(simp) == 1:
            matrix.append([])
        elif len(simp) == 2:
            column = []
            column.append((simp_to_index_dict[(simp[1],)], 1))
            column.append((simp_to_index_dict[(simp[0],)], -1))
            matrix.append(sorted(column, reverse=True))
        elif len(simp) == 3:
            column = []
            column.append((simp_to_index_dict[(simp[0], simp[1])], 1))
            column.append((simp_to_index_dict[(simp[0], simp[2])], -1))
            column.append((simp_to_index_dict[(simp[1], simp[2])], 1))
            matrix.append(sorted(column, reverse=True))
        else:
            column = []
            column.append((simp_to_index_dict[(simp[0], simp[1], simp[2])], -1))
            column.append((simp_to_index_dict[(simp[0], simp[1], simp[3])], 1))
            column.append((simp_to_index_dict[(simp[0], simp[2], simp[3])], -1))
            column.append((simp_to_index_dict[(simp[1], simp[2], simp[3])], 1))
            matrix.append(sorted(column, reverse=True))
            
    return matrix, simp_to_index_dict, index_to_simp_dict


def reduction_sparse_matrix(boundary_matrix, mod):
    """Return reduced boundary matrix in column echelon form, basis elements for each column and indices of all
    zero columns. All operations are done modulo mod.

    Parameters:
    boundary matrix -- as in the function sparse_boundary_matrix
    mod -- positive integer

    Return:
    reduced_matrix -- reduced matrix in echelon form in the sparse format, given as a list of pairs which are initially
                      all equal to ([], None)
                      reduced_matrix[i] = (col, j) iff j-th column has lowest nonzero entry in i-th row
                      column col is given in ordered sparse format with lowest entry equal to 1
    zero_columns -- a list of indices of zero columns in reduced matrix
    columns -- a list of dictionaries, columns[i] represents the basis element (chain) for the i-th column
               columns[i] is a dictionary: key = indeks of a simplex in the total ordering of the filtration
                                           value = coefficient (only nonzero) for the simplex in the chain
    """
    m = len(boundary_matrix)
    reduced_matrix = [([], None) for _ in range(m)]

#    print_sparse(boundary_matrix)
#    print('\n')

    zero_columns = []
    columns = [{i: 1} for i in range(m)]  # We start with the standard basis.
    
    for j in range(m):
        column = boundary_matrix[j]
        while len(column) > 0 and len(reduced_matrix[column[0][0]][0]) > 0:  # Reduce the column.
            add_chains(columns[j], columns[reduced_matrix[column[0][0]][1]], - column[0][1], mod)
            column = add_columns(column, reduced_matrix[column[0][0]][0], - column[0][1], mod)

        if len(column) > 0:  # New pivot column.
            mult_column(column, sp.mod_inverse(column[0][1], mod), mod)
            mult_chain(columns[j], sp.mod_inverse(column[0][1], mod), mod)
            reduced_matrix[column[0][0]] = (column, j)
        else:  # New zero column.
            zero_columns.append(j)

#        test_matrix = reduced_matrix[:j+1] + boundary_matrix[j+1:]
#        print_sparse(test_matrix)
#        print('\n')

#    print('boundary_matrix')
#    print_sparse1(reduced_matrix)

    return reduced_matrix, zero_columns, columns


def to_dict(column):
    """Return column, given as a list, to column_dict, given as a dictionary."""
    column_dict = dict()
    for indeks, koeficient in column:
        column_dict[indeks] = koeficient
    return column_dict


def persistence_intervals(reduced_matrix, zero_columns, columns, index_to_simp, simp_dict, dim2=False):
    """Return persistent intervals and their persistent generators in dimensions 0 and 1,
    if dim2 = True also return them in dimension 2.

    Parameters:
    reduced_matrix -- as in the function reduction_sparse_matrix
    zero columns -- as in the function reduction_sparse_matrix
    columns -- as in the function reduction_sparse_matrix
    index_to_simp -- as in the function sparse_boundary_matrix
    simp_dict -- as in the function rips.rips_complex (or as domain_dict in the function rips.domain)
    dim2 -- True of False if we want to compute also 2D persistent homology

    Return:
    homology_intervals -- a list, at i-th index is a list of persistent intervals in dimension i
                          intervals are given as pairs (birth, death) or (birth, None) if death is at infinity
    homology_basis -- a list, at i-th index is a list of persistent generators in dimension i
                      generators are given as chains (dictionaries as in the function reduction_sparse_matrix
    """
    homology_intervals = [[], []]
    if dim2:
        homology_intervals.append([])

    homology_basis = [[], []]
    if dim2:
        homology_basis.append([])

    # Degrees of homology to compute.
    degree = 4 if dim2 else 3

    for col in zero_columns:
        new_interval = False  # True, if we created a new interval.
        simp_i = index_to_simp[col]
        birth = simp_dict[simp_i]  # Filtration index when the cycle is born.
        
        if len(reduced_matrix[col][0]) == 0:
            interval = (birth, None)  # The cycle never dies.
            new_interval = True
        else:
            simp_j = index_to_simp[reduced_matrix[col][1]]
            death = simp_dict[simp_j]
            if death > birth:  # The cycle doesn't die immediately.
                interval = (birth, death - 1)
                new_interval = True

        if new_interval and len(simp_i) < degree:
            d = len(simp_i) - 1
            homology_intervals[d].append(interval)
            birth, death = interval
            
            if death is None:  # Basis element is already stored in the corresponding column of the reduced matrix.
                homology_basis[d].append(columns[col])
            else:  # Basis element is the image of the current basis element in the column (its boundary).
                homology_basis[d].append(to_dict(reduced_matrix[col][0]))
            
    return homology_intervals, homology_basis


def filtration_basis(intervals, max_filt_index):
    """Return a collection of bases for each index in the filtration.

    Parameters:
    intervals -- a list of persistent intervals for one dimension only (as in the function persistence_intervals)
    max_filt_index -- positive integer,  maximal filtration index

    Return:
    filt_basis -- a list of length max_filt_index,
                  at index i is a list of indices of persistent intervals which contain i (their generators
                  form the basis at index i in the filtration)
    """
    filt_basis = [[] for _ in range(max_filt_index)]

    for x in range(len(intervals)):
        birth, death = intervals[x]
        if death is None:
            for k in range(birth - 1, max_filt_index):
                filt_basis[k].append(x)
        else:
            for k in range(birth - 1, death):
                filt_basis[k].append(x)

    return filt_basis


def mapped_basis(dom_basis, dom_index_to_simp, simp_to_index, mod, mapped_simp=False):
    """Return the mapped persistent generators for the filtration of the domain by the mapping given by mapped_simp
    (if mapped_simp = False, we assume that the map is the inclusion of the domain). We are working modulo mod.

    Parameters:
    dom_basis -- persistent generators in one dimension only (as given in homology_basis in the function
                 persistent_intervals)
    dom_index_to_simp -- a dictionary as in the function sparse_boundary_matrix (for the filtration of domains)
    simp_to_index -- a dictionary as in the function sparse_boundary_matric
    mod -- positive integer
    mapped simp -- False (in the case of inclusion) or a dictionary as in the function rips.domain

    Return:
    mapped_basis -- a list of mapped basis elements from dom_basis, given as dictionaries representing their chains
                    (in the same convention as dom_basis)
    """
    mapped_basis = []
    if mapped_simp == False:  # We map with inclusion.
        for cycle in dom_basis:
            mapped_cycle = dict()

            for index in cycle:
                simp = dom_index_to_simp[index]
                mapped_index = simp_to_index[simp]
                mapped_cycle[mapped_index] = cycle[index]

            mapped_basis.append(mapped_cycle)

    else:
        for cycle in dom_basis:
            mapped_cycle = dict()

            for index in cycle:
                simp = dom_index_to_simp[index]
                image_simp, signature, collapse = mapped_simp[simp]

                if not collapse:  # Mapped simplex has the same dimension.
                    mapped_index = simp_to_index[image_simp]
                    if mapped_index not in mapped_cycle:
                        mapped_cycle[mapped_index] = (cycle[index] * signature) % mod
                    else:
                        value = (mapped_cycle[mapped_index] + cycle[index] * signature) % mod
                        if value != 0:
                            mapped_cycle[mapped_index] = value
                        else:
                            del mapped_cycle[mapped_index]

            mapped_basis.append(mapped_cycle)

    return mapped_basis


def to_list(column):
    """Return a column, given as a chain dictionary, as a list of pairs (row index, value in row)
    ordered by decreasing row index.
    """
    return [(index, column[index]) for index in sorted(column, reverse=True)]


def basis_coefficients(map_basis, basis, reduced_matrix, mod):
    """Return the matrix of coefficients for the 'map_basis' represented in 'basis', working modulo mod.

    Parameters:
    map_basis -- as in the function mapped_basis
    basis -- as an element of homology_basis from the function persistent_intervals
    reduced matrix -- as in the function reduction_sparse_matrix
    mod -- positive integer

    Return:
    coefficients -- a transition matrix from map_basis to basis
    """
    coefficients = [[0 for _ in range(len(basis))] for _ in range(len(map_basis))]

    # Convert basis elements from dictionaries to lists.
    map_basis = [to_list(b) for b in map_basis]
    basis = [to_list(b) for b in basis]

    # We add persistent generators of infinite intervals from basis, since their corresponding columns in
    # reduced_matrix are zero.
    alt_reduced_matrix = copy.deepcopy(reduced_matrix)
    for cycle in basis:
        if len(alt_reduced_matrix[cycle[0][0]][0]) == 0:
            alt_reduced_matrix[cycle[0][0]] = (cycle, cycle[0][0])

    for i in range(len(map_basis)):
        cycle = map_basis[i]
        coeffs = dict()  # Keeps track of all the reductions.
        while len(cycle) > 0:  # Reducing map_basis via reduced_matrix.
            assert len(alt_reduced_matrix[cycle[0][0]][0]) > 0
            coeffs[alt_reduced_matrix[cycle[0][0]][1]] = cycle[0][1]
            cycle = add_columns(cycle, alt_reduced_matrix[cycle[0][0]][0], - cycle[0][1], mod)

        for j in range(len(basis)):
            col = alt_reduced_matrix[basis[j][0][0]][1]
            if col in coeffs:  # Store only coeffs from elements of basis.
                coefficients[i][j] = coeffs[col]

    return coefficients


def tower_of_pairs(dom_filt_basis, filt_basis, coeffs, max_filt_index):
    """Return a tower of linear maps (matrices).

    Parameters:
    dom_filt_basis, filt_basis -- basis for filtrations of the domain and the image,
                                  given as in the function filtration_basis
    coeffs -- matrix of coefficients as in the function basis_coefficients
    max_filt_index -- positive integer as in the function rips.rips_complex

    Return:
    tower -- a list of length max_filt_index, i-th entry equals the coeffs matrix restricted to the columns and rows
             given by dom_filt_basis[i] and filt_basis[i] (if any of them are empty, tower[i] = None)
    """
    tower = [None for _ in range(max_filt_index)]

    for ind in range(max_filt_index):
        dom_basis = dom_filt_basis[ind]
        basis = filt_basis[ind]
        matrix = [[None for _ in range(len(dom_basis))] for _ in range(len(basis))]

        for i in range(len(basis)):
            for j in range(len(dom_basis)):
                matrix[i][j] = coeffs[dom_basis[j]][basis[i]]

        tower[ind] = matrix

    return tower


def print_sparse(matrix):
    """Print a sparse matrix.

    Parameters:
    matrix -- a list of columns, each column is given as a list of pairs (row index, value in row)
    """
    m = len(matrix)
    dense_matrix = [[0 for _ in range(m)] for _ in range(m)]
    
    for i in range(m):
        if isinstance(matrix[i], tuple):
            col = matrix[i][0]
        else:
            col = matrix[i]
        for j, v in col:
            dense_matrix[j][i] = v

    for line in dense_matrix:
        print(line)


def print_sparse1(matrix):
    """Print a sparse matrix.

    Parameters:
    matrix -- a list of columns, each column is given as a tuple (index of column, list of pairs (row index, value))
    """
    m = len(matrix)
    dense_matrix = [[0 for _ in range(m)] for _ in range(m)]
    
    for i in range(m):
        data, ind = matrix[i]
        if len(data) > 0:
            for j, val in data:
                dense_matrix[j][ind] = val

    for line in dense_matrix:
        print(line)


if __name__ == '__main__':
    # test_filtracija = [(1, (1,)), (2, (2,)), (3, (3,)), (4, (1, 2)),
    #                    (5, (2, 3)), (6, (1, 3)), (7, (1, 2, 3))]
    # dict1 = {0: (1,), 1: (2,), 2: (3,), 3: (1, 2), 4: (2, 3), 5: (1, 3), 6: (1, 2, 3)}
    # dict2 = {(1,): 1, (2,): 2, (3,): 3, (1, 2): 4, (2, 3): 5, (1, 3): 6, (1, 2, 3): 7}
    # sbm, _, _ = sparse_boundary_matrix(test_filtracija)
    # rm, zero, cols = reduction_sparse_matrix(sbm, 2)
    # h = persistence_intervals(rm, zero, cols, dict1, dict2)
    #
    # test1_filtracija = [(1, (0,)), (1, (1,)), (2, (2,)), (2, (3,)), (2, (0, 1)),
    #                     (2, (1, 2)), (3, (0, 3)), (3, (2, 3)), (4, (0, 2)),
    #                     (5, (0, 1, 2)), (6, (0, 2, 3))]
    # dict3 = {0: (0,), 1: (1,), 2: (2,), 3: (3,), 4: (0, 1), 5: (1, 2), 6: (0, 3),
    #          7: (2, 3), 8: (0, 2), 9: (0, 1, 2), 10: (0, 2, 3)}
    # dict4 = {(0,): 1, (1,): 1, (2,): 2, (3,): 2, (0, 1): 2, (1, 2): 2, (0, 3): 3,
    #          (2, 3): 3, (0, 2): 4, (0, 1, 2): 5, (0, 2, 3): 6}
    # sbm1, _, _ = sparse_boundary_matrix(test1_filtracija)
    # rm1, zero1, cols1 = reduction_sparse_matrix(sbm1, 11)
    # h1 = persistence_intervals(rm1, zero1, cols1, dict3, dict4)

    points = rips.points_torus_polar(10, 5, 10, 5)
    images = rips.rotate_torus_polar(points, 3, 2)
    mapped = rips.mapped_points(points, images, rips.distance_torus_polar)
    
    rips_comp, simp_dict, max_filt_index, epsilon = rips.rips_complex(points, rips.distance_torus_polar)
    domain_filt, mapped_simp, domain_dict = rips.domain(rips_comp, simp_dict, mapped)

    # For filtration of K_i.
    boundary_matrix, simp_to_index, index_to_simp = sparse_boundary_matrix(rips_comp)
    reduced_matrix, zero_columns, columns = reduction_sparse_matrix(boundary_matrix, 29)
    [hom0, hom1], [hom0_basis, hom1_basis] = \
        persistence_intervals(reduced_matrix, zero_columns, columns, index_to_simp, simp_dict)
    hom0_filt_basis = filtration_basis(hom0, max_filt_index)
    hom1_filt_basis = filtration_basis(hom1, max_filt_index)
    # hom2_filt_basis = filtration_basis(hom2, max_filt_index)

    # For filtration of domains.
    dom_boundary_matrix, dom_simp_to_index, dom_index_to_simp = sparse_boundary_matrix(domain_filt)
    dom_reduced_matrix, dom_zero_columns, dom_columns = reduction_sparse_matrix(dom_boundary_matrix, 29)
    [dom_hom0, dom_hom1], [dom_hom0_basis, dom_hom1_basis] = \
        persistence_intervals(dom_reduced_matrix, dom_zero_columns, dom_columns, dom_index_to_simp, domain_dict)
    dom_hom0_filt_basis = filtration_basis(dom_hom0, max_filt_index)
    dom_hom1_filt_basis = filtration_basis(dom_hom1, max_filt_index)
    # dom_hom2_filt_basis = filtration_basis(dom_hom2, max_filt_index)

    # Map the basis.
    inclusion_basis0 = mapped_basis(dom_hom0_basis, dom_index_to_simp, simp_to_index, 29)
    inclusion_basis1 = mapped_basis(dom_hom1_basis, dom_index_to_simp, simp_to_index, 29)
    # inclusion_basis2 = mapped_basis(dom_hom2_basis, dom_index_to_simp, simp_to_index, 29)
    mapped_basis0 = mapped_basis(dom_hom0_basis, dom_index_to_simp, simp_to_index, 29, mapped_simp)
    mapped_basis1 = mapped_basis(dom_hom1_basis, dom_index_to_simp, simp_to_index, 29, mapped_simp)
    # mapped_basis2 = mapped_basis(dom_hom2_basis, dom_index_to_simp, simp_to_index, 29, mapped_simp)

    # Write mapped basis elements in terms of the original basis elements.
    inclusion_coeffs0 = basis_coefficients(inclusion_basis0, hom0_basis, reduced_matrix, 29)
    inclusion_coeffs1 = basis_coefficients(inclusion_basis1, hom1_basis, reduced_matrix, 29)
    # inclusion_coeffs2 = basis_coefficients(inclusion_basis2, hom2_basis, reduced_matrix, 29)
    mapped_coeffs0 = basis_coefficients(mapped_basis0, hom0_basis, reduced_matrix, 29)
    mapped_coeffs1 = basis_coefficients(mapped_basis1, hom1_basis, reduced_matrix, 29)
    # mapped_coeffs2 = basis_coefficients(mapped_basis2, hom2_basis, reduced_matrix, 29)

    # Creating a tower of pairs.
    inclusion_tower0 = tower_of_pairs(dom_hom0_filt_basis, hom0_filt_basis, inclusion_coeffs0, max_filt_index)
    inclusion_tower1 = tower_of_pairs(dom_hom1_filt_basis, hom1_filt_basis, inclusion_coeffs1, max_filt_index)
    # inclusion_tower2 = tower_of_pairs(dom_hom2_filt_basis, hom2_filt_basis, inclusion_coeffs2, max_filt_index)
    mapped_tower0 = tower_of_pairs(dom_hom0_filt_basis, hom0_filt_basis, mapped_coeffs0, max_filt_index)
    mapped_tower1 = tower_of_pairs(dom_hom1_filt_basis, hom1_filt_basis, mapped_coeffs1, max_filt_index)
    # mapped_tower2 = tower_of_pairs(dom_hom2_filt_basis, hom2_filt_basis, mapped_coeffs2, max_filt_index)
