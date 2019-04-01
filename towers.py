import rips
import sympy as sp
import copy


def add_columns(col1, col2, d, mod):
    #stolpcu col1 pristeje d-kratnik stolpca col2 modulo m
    #stolpca imata elemente urejene po padajocem indeksu in oblike (indeks vrstice, vrednost)
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
        remainer = col2[j:]
        mult_column(remainer, d, mod)
        new_col += remainer
        
    return new_col

def mult_column(col, d, mod):
    for i in range(len(col)):
        ind, value = col[i]
        col[i] = (ind, (d * value) % mod)

def add_chains(chain1, chain2, d, mod):
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
    for key in chain:
        chain[key] = (d * chain[key]) % mod


def sparse_boundary_matrix(filtration):
    simp_to_index_dict = dict() #kljuc: simplex, vrednost:indeks stolpca, ki pripada simplexu
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
        else:
            column = []
            column.append((simp_to_index_dict[(simp[0], simp[1])], 1))
            column.append((simp_to_index_dict[(simp[0], simp[2])], -1))
            column.append((simp_to_index_dict[(simp[1], simp[2])], 1))
            matrix.append(sorted(column, reverse=True))
            
    return matrix, simp_to_index_dict, index_to_simp_dict

def reduction_sparse_matrix(boundary_matrix, mod):
    m = len(boundary_matrix)
    reduced_matrix = [([], None) for _ in range(m)]
    #j-ti element (prva koordinata) je stolpec v R, ki ima najbolj spodnji nenicelen element v j-ti vrstici
    #zapisan je kot padajoc seznam parov (indeks_vrstice, koeficient pri tem simpleksu)
    #vedno poskrbimo, da je najbolj spodnji nenicelen element 1
    #dejansko mesto stolpca je zapisano v drugi koordinati

##    print_sparse(boundary_matrix)
##    print('\n')

    zero_columns = []
    #stolpci, kjer se rodijo novi cikli (zacetki persistance intervalov)
    columns = [{i:1} for i in range(m)]
    #hrani verige, ki se trenutno nahajajo v stolpcih
    #verigo predstavlja slovar: kljuc: indeks simplexa, vrednost: koeficient v verigi
    
    for j in range(m):
        column = boundary_matrix[j]
        while len(column) > 0 and len(reduced_matrix[column[0][0]][0]) > 0:
            add_chains(columns[j], columns[reduced_matrix[column[0][0]][1]], - column[0][1], mod)
            column = add_columns(column, reduced_matrix[column[0][0]][0], - column[0][1], mod)

        if len(column) > 0:
            mult_column(column, sp.mod_inverse(column[0][1], mod), mod)
            mult_chain(columns[j], sp.mod_inverse(column[0][1], mod), mod)
            reduced_matrix[column[0][0]] = (column, j)
        else:
            zero_columns.append(j)

##        test_matrix = reduced_matrix[:j+1] + boundary_matrix[j+1:]
##        print_sparse(test_matrix)
##        print('\n')

##    print('boundary_matrix')
##    print_sparse1(reduced_matrix)

    return reduced_matrix, zero_columns, columns

def to_dict(column):
    column_dict = dict()
    for indeks, koeficient in column:
        column_dict[indeks] = koeficient
    return column_dict

def persistence_intervals(reduced_matrix, zero_columns, columns,
                          index_to_simp, simp_dict):
    homology_0_intervals = []
    homology_1_intervals = []

    #za 1 zamaknjeno nazaj: indeks i pomeni bazo za K_(i+1)
    homology_0_basis = []
    homology_1_basis = []

    for col in zero_columns:
        new_interval = False #True, ce smo ustvarili nov interval
        simp_i = index_to_simp[col]
        birth = simp_dict[simp_i] #index filtracije, ko se cikel rodi
        
        if len(reduced_matrix[col][0]) == 0:
            interval = (birth, None) #cikel nikoli ne umre
            new_interval = True
        else:
            simp_j = index_to_simp[reduced_matrix[col][1]]
            death = simp_dict[simp_j]
            if death > birth: #prezivi vsaj do naslednje razsiritve filtracije
                interval = (birth, death - 1)
                new_interval = True

        if new_interval and len(simp_i) == 1:
            homology_0_intervals.append(interval)
            birth, death = interval
            
            if death is None: #dodamo nicelni stolpec v bazo (verigo katere boundary je 0) in njegov indeks v reduced_matrix
                #v tem primeru bo na zapisanem indeksu prazen seznam, ker smo dodali nicelni stolpec v bazo
                homology_0_basis.append(columns[col])
            else: #dodamo nenicelni stolpec (boundary verige v stolpcu)
                homology_0_basis.append(to_dict(reduced_matrix[col][0]))
                    
        if new_interval and len(simp_i) == 2:
            homology_1_intervals.append(interval)
            birth, death = interval
            
            if death is None:
                homology_1_basis.append(columns[col])
            else:
                homology_1_basis.append(to_dict(reduced_matrix[col][0]))
            
    return homology_0_intervals, homology_1_intervals, homology_0_basis, homology_1_basis

def filtration_basis(intervals, basis, max_filt_index):
    filt_basis = [[] for _ in range(max_filt_index)]

    for x in range(len(intervals)):
        birth, death = intervals[x]
        if death is None:
            for k in range(birth - 1, max_filt_index):
                filt_basis[k].append(basis[x])
        else:
            for k in range(birth - 1, death):
                filt_basis[k].append(basis[x])

    return filt_basis

def mapped_basis(dom_basis, dom_index_to_simp, simp_to_index, mod, mapped_simp = False):
    mapped_basis = []
    if mapped_simp == False: #gledamo inkluzijo
        for cycle in dom_basis:
            mapped_cycle = dict()

            for index in cycle: #pogledamo kateri index ustreza temu simpleksu v celotni filtraciji
                simp = dom_index_to_simp[index]
                mapped_index = simp_to_index[simp]
                mapped_cycle[mapped_index] = cycle[index]

            mapped_basis.append(mapped_cycle)

    else: #mapped_simp je slovar preslikanih simpleksev
        for cycle in dom_basis:
            mapped_cycle = dict()

            for index in cycle:
                simp = dom_index_to_simp[index]
                image_simp, signature, collapse = mapped_simp[simp]

                if not collapse:
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
    return [(index, column[index]) for index in sorted(column, reverse=True)]


def basis_coefficients(map_basis, basis, reduced_matrix, mod):
    coefficients = [[0 for _ in range(len(basis))] for _ in range(len(map_basis))]

    #uredimo cikle spet po padajocem indeksu simpleksov
    map_basis = [to_list(b) for b in map_basis]
    basis = [to_list(b) for b in basis]

    #dodamo v matriko tiste bazne elemente, ki pridejo od intervalov z death = None
    #na najnizjem mestu bodo imeli 1 in sicer bo to ravno na diagonali
    alt_reduced_matrix = copy.deepcopy(reduced_matrix)
    for cycle in basis:
        if len(alt_reduced_matrix[cycle[0][0]][0]) == 0:
            alt_reduced_matrix[cycle[0][0]] = (cycle, cycle[0][0])

    for i in range(len(map_basis)):
        cycle = map_basis[i]
        coeffs = dict()
        while len(cycle) > 0:
            #vsak stolpec bo najvec enkrat na vrsti ker se indeks najbolj spodnjega nenicelnega strogo manjsa
            assert len(alt_reduced_matrix[cycle[0][0]][0]) > 0
            coeffs[alt_reduced_matrix[cycle[0][0]][1]] = cycle[0][1]
            cycle = add_columns(cycle, alt_reduced_matrix[cycle[0][0]][0], - cycle[0][1], mod)

        for j in range(len(basis)):
            col = alt_reduced_matrix[basis[j][0][0]][1]
            if col in coeffs:
                coefficients[i][j] = coeffs[col]

    return coefficients




def print_sparse(matrix):
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

##    test_filtracija = [(1, (1,)), (2, (2,)), (3, (3,)), (4, (1, 2)),
##                       (5, (2, 3)), (6, (1, 3)), (7, (1, 2, 3))]
##    dict1 = {0:(1,), 1:(2,), 2:(3,), 3:(1,2), 4:(2,3), 5:(1,3), 6:(1,2,3)}
##    dict2 = {(1,):1, (2,):2, (3,):3, (1,2):4, (2,3):5, (1,3):6, (1,2,3):7}
##    sbm, _, _ = sparse_boundary_matrix(test_filtracija)
##    rm, zero, cols = reduction_sparse_matrix(sbm, 2)
##    h0, h1, h0b, h1b = persistence_intervals(rm, zero, cols, dict1, dict2, 7)

##    test1_filtracija = [(1,(0,)), (1,(1,)), (2,(2,)), (2,(3,)), (2,(0,1)),
##                        (2, (1,2)), (3, (0,3)), (3, (2,3)), (4, (0,2)),
##                        (5, (0,1,2)), (6, (0,2,3))]
##    dict3 = {0:(0,), 1:(1,), 2:(2,), 3:(3,), 4:(0,1), 5:(1,2), 6:(0,3),
##             7:(2,3), 8:(0,2), 9:(0,1,2), 10:(0,2,3)}
##    dict4 = {(0,):1, (1,):1, (2,):2, (3,):2, (0,1):2, (1,2):2, (0,3):3,
##             (2,3):3, (0,2):4, (0,1,2):5, (0,2,3):6}
##    sbm1, _, _ = sparse_boundary_matrix(test1_filtracija)
##    rm1, zero1, cols1 = reduction_sparse_matrix(sbm1, 11)
##    hh0, hh1, hh0b, hh1b = persistence_intervals(rm1, zero1, cols1, dict3, dict4, 6)

    
    points = rips.points_circle_polar(1, 50, 0.2)
    images = rips.power_polar(points, 2)
    mapped = rips.mapped_points(points, images, rips.distance_polar)
    
    rips_comp, simp_dict, max_filt_index = rips.rips_complex(points, rips.distance_polar)
    domain_filt, mapped_simp, domain_dict = rips.domain(rips_comp, simp_dict, mapped)

    #za K_i
    boundary_matrix, simp_to_index, index_to_simp = sparse_boundary_matrix(rips_comp)
    reduced_matrix, zero_columns, columns = reduction_sparse_matrix(boundary_matrix, 29)
    hom0, hom1, hom0_basis, hom1_basis = persistence_intervals(reduced_matrix, zero_columns, columns,
                                                               index_to_simp, simp_dict)
    hom0_filt_basis = filtration_basis(hom0, hom0_basis, max_filt_index)
    hom1_filt_basis = filtration_basis(hom1, hom1_basis, max_filt_index)

    #za dom_i
    dom_boundary_matrix, dom_simp_to_index, dom_index_to_simp = sparse_boundary_matrix(domain_filt)
    dom_reduced_matrix, dom_zero_columns, dom_columns = reduction_sparse_matrix(dom_boundary_matrix, 29)
    dom_hom0, dom_hom1, dom_hom0_basis, dom_hom1_basis = \
        persistence_intervals(dom_reduced_matrix, dom_zero_columns, dom_columns, dom_index_to_simp, domain_dict)
    dom_hom0_filt_basis = filtration_basis(dom_hom0, dom_hom0_basis, max_filt_index)
    dom_hom1_filt_basis = filtration_basis(dom_hom1, dom_hom1_basis, max_filt_index)

    #preslikamo bazo
    inclusion_basis0 = mapped_basis(dom_hom0_basis, dom_index_to_simp, simp_to_index, 29)
    inclusion_basis1 = mapped_basis(dom_hom1_basis, dom_index_to_simp, simp_to_index, 29)
    mapped_basis0 = mapped_basis(dom_hom0_basis, dom_index_to_simp, simp_to_index, 29, mapped_simp)
    mapped_basis1 = mapped_basis(dom_hom1_basis, dom_index_to_simp, simp_to_index, 29, mapped_simp)

    #razvijemo preslikano bazo po originalni bazi
    inclusion_coeffs0 = basis_coefficients(inclusion_basis0, hom0_basis, reduced_matrix, 29)
    inclusion_coeffs1 = basis_coefficients(inclusion_basis1, hom1_basis, reduced_matrix, 29)
    mapped_coeffs0 = basis_coefficients(mapped_basis0, hom0_basis, reduced_matrix, 29)
    mapped_coeffs1 = basis_coefficients(mapped_basis1, hom1_basis, reduced_matrix, 29)
    

