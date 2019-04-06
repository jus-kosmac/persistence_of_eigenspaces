import sympy as sp

def kernel(A, mod):
    #A je m x n matrika, NEPRAZNI DOMENA IN KODOMENA!
    #izracuna jedro via column operations
    m = len(A)
    n = len(A[0])

    L = [[A[i][j] % mod for j in range(n)] for i in range(m)] #spodnje trikotna na koncu
    U = [[0 for _ in range(n)] for _ in range(n)] #velja L = A * U, U je lahko polna matrika?
    for i in range(n):
        U[i][i] = 1

    curr_line = 0 #vrstica, kjer trenutno želimo pridobiti nicle
    first_zero_column = None #None: matrika ima trivialno jedro, sicer indeks prvega stolpca, ki je v jedru

    for j in range(n):
        curr_col = j #trenutni stolpec
        pivot = False

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

        if not pivot: #vsi preostali elementi so nicelni
            first_zero_column = j
            break
        else:
            if curr_col != j:
                #zamenjamo stolpca
                for i in range(curr_line, m):
                    L[i][j], L[i][curr_col] = L[i][curr_col], L[i][j]
                for i in range(n):
                    U[i][j], U[i][curr_col] = U[i][curr_col], U[i][j]

            inv = sp.mod_inverse(L[curr_line][j], mod)
            #eliminiramo ostale elemente v vrstici
            for k in range(j + 1, n):
                coeff = L[curr_line][k] * inv
                for i in range(curr_line, m):
                    L[i][k] = (L[i][k] - coeff * L[i][j]) % mod
                for i in range(n):
                    U[i][k] = (U[i][k] - coeff * U[i][j]) % mod

    return L, U, first_zero_column

def mat_mult(A, B, mod):
    #privzamemo da sta matriki dimenzij m x n in n x p NEPRAZNI!
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
    #vect je vektor predstavljen kot seznam (vrstica)
    #vrne vrstico kot rezultat
    n = len(vect)
    m = len(A)

    C = [0 for _ in range(m)]
    for i in range(m):
        for j in range(n):
            C[i] = (C[i] + A[i][j] * vect[j]) % mod

    return C

def mat_sum(A, B, t, mod):
    #vrne A - t * B, obe velikosti m x n NEPRAZNI!
    m = len(A)
    n = len(A[0])

    return [[(A[i][j] - t * B[i][j]) % mod for j in range(n)] for i in range(m)]

def intersection(A, B, mod):
    #A in B sta BAZI za 2 podprostora, predstavljeni kot seznam vrstic (vektorjev) NEPRAZNI!
    #moci baz sta n in p, dolzina vektorjev je m
    #vrne bazo kot seznam vektorjev
    m = len(A[0])
    n = len(A)
    p = len(B)

    C = [[A[j][i] % mod for j in range(n)] + [B[j][i] % mod for j in range(p)] for i in range(m)] #m x (n + p)
    _, U, col = kernel(C, mod)

    intersect_basis = []
    if col is not None:
        A_mat = [line[:n] for line in C]
        for j in range(col, n + p):
            vect = [U[i][j] for i in range(n)]
            intersect_basis.append(vect_mult(A_mat, vect, mod))

    return intersect_basis

def kernel_basis(A, mod):
    #vrne seznam vrstic
    _, U, col = kernel(A, mod)
    basis = []
    if col is not None:
        for j in range(col, len(U)):
            vect = [U[i][j] for i in range(len(U))]
            basis.append(vect)

    return basis

def reduce_vect(top_list, vect, mod):
    vect = vect.copy()
    top = None #indeks, kjer je prvi nenicelni element v reducirani obliki vect
    coeffs = dict() #koeficienti v razvoju po bazi v top_list
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
    #A je podprostor B, dopolni bazo za A (jo popravi do spodnje trikotne) do baze za B
    #pri tem skrbi, da je baza vedno spodaj trikotna, da bomo lazje razvijali po tej bazi kasneje
    #B NEPRAZEN!
    p = len(B[0])
    #na indeksu i je bazni vektor, ki ima na i-tem mestu od zgoraj prvi nenicelni element
    top_non_zero = [([], None) for _ in range(p)]
    quotient_basis = [] #indeksi, na katerih so elementi, ki smo jih dobili iz B

    for vect in A:
        vect_reduced, top, _ = reduce_vect(top_non_zero, vect, mod)
        if top is not None: #se ni zreduciral do nicelnega vektorja
            top_non_zero[top] = (vect_reduced, 'intersection')

    for i in range(len(B)):
        vect_reduced, top, _ = reduce_vect(top_non_zero, B[i], mod)
        if top is not None:
            top_non_zero[top] = (vect_reduced, 'eigenspace')
            quotient_basis.append(top)

    return top_non_zero, quotient_basis

def quotient_coeffs(top_list, quotient_basis, vect, mod):
    _, _, coeffs = reduce_vect(top_list, vect, mod)
    return [coeffs[i] if i in coeffs else 0 for i in quotient_basis]

def basis_conversion(vect_one, basis_one, basis_two):
    # basis_two ima vedno dodane le vecje indekse od basis_one, vsi indeksi so urejeni
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
    # seznam tuplov oblike (basis list glede na top nenicelni index, indeksi kjer so bazni elementi kvocienta)
    # ce je nicelen prostor je vrednost None
    eigen_tower_basis = [None for _ in range(max_filt_index)]
    # bazni vektorji za kvocient, koeficienti pripadajo baznim vektorjem v domain_filt_basis
    eigen_tower_quotient_basis = [[] for _ in range(max_filt_index)]
    # matrike, ki tvorijo stolp
    # ce je preslikava iz nicelnega ali v nicelni prostor je vrednost None
    eigen_tower = [None for _ in range(max_filt_index - 1)]

    for i in range(max_filt_index):
        A = mapped_tower[i]
        B = inclusion_tower[i]

        if len(A) == 0: #homologija K_i ima dimenzijo 0
            dim = len(domain_filt_basis[i]) #dimenzija homologije dom_i
            if dim == 0: #trivialen vektorski prostor = ni baze
                continue
            else: #vsi bazni elementi so v preseku jeder
                top_non_zero = [([], None) for _ in range(dim)]
                for d in range(dim):
                    vect = [0 for _ in range(dim)]
                    vect[d] = 1
                    top_non_zero[d] = (vect, 'intersection')

                eigen_tower_basis[i] = (top_non_zero, [])
        elif len(A[0]) == 0: #nicelen prostor
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
        if eigen_tower_basis[i] is None or eigen_tower_basis[i + 1] is None: #celoten prostor je trivialen
            continue
        else:
            indices_list_first, basis_first = eigen_tower_basis[i]
            indices_list_second, basis_second = eigen_tower_basis[i + 1]
            m = len(basis_second)
            n = len(basis_first)
            if n == 0 or m == 0: #kvocient je trivialen
                continue

            eigen_matrix = [[None for _ in range(n)] for _ in range(m)] #m x n matrika

            for j in range(n):
                vect_first = indices_list_first[basis_first[j]][0] #koeficienti po prvotni bazi celega prostora
                vect_second = basis_conversion(vect_first, domain_filt_basis[i], domain_filt_basis[i + 1])
                coeffs = quotient_coeffs(indices_list_second, basis_second, vect_second, mod)

                for k in range(m):
                    eigen_matrix[k][j] = coeffs[k]

            eigen_tower[i] = eigen_matrix

    return eigen_tower, eigen_tower_quotient_basis

def tower_normal_form(tower, tower_basis, max_filt_index, mod):
    #SPREMINJA TOWER IN TOWER_BASIS!
    for ind in range(max_filt_index - 2, -1, -1):
        right = tower[ind] #m x n matrix
        if right is None: #slikamo v prazen prostor ali iz praznega prostora
            continue

        basis = tower_basis[ind] #baza domene right in kodomene left
        left = None #ce nimamo leve matrike ali pa je levi prostor niceln
        if ind > 0: #imamo se levo matriko
            left = tower[ind - 1]

        m = len(right)
        n = len(basis)
        pivot = 0 #trenutni stolpec

        for i in range(m):
            column = pivot #bo postal stolpec, ki ima nenicelni element v tej vrstici
            while column < n and right[i][column] % mod == 0:
                column += 1
            if column == n: #nicelna vrstica
                continue

            #zamenjamo stolpca
            basis[pivot], basis[column] = basis[column], basis[pivot]
            for row in range(m):
                right[row][pivot], right[row][column] = right[row][column], right[row][pivot]
            #zamenjamo vrstici
            if left is not None:
                left[pivot], left[column] = left[column], left[pivot]

            #pivotni element nastavimo na 1 in popravimo bazo ter obe matriki
            pivot_el = right[i][pivot]
            inv = sp.mod_inverse(pivot_el, mod)
            for row in range(i, m):
                right[row][pivot] = (right[row][pivot] * inv) % mod
            for k in range(len(basis[pivot])):
                basis[pivot][k] = (basis[pivot][k] * inv) % mod
            if left is not None:
                for k in range(len(left[pivot])):
                    left[pivot][k] = (left[pivot][k] * pivot_el) % mod

            #eliminiramo ostale elemente v vrstici in popravimo matrike
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
        left = tower[ind] #m x n matrika
        if left is None:
            continue

        basis = tower_basis[ind + 1]
        right = None
        if ind < max_filt_index - 2:
            right = tower[ind + 1]

        m = len(left)
        n = len(left[0])
        pivot = 0 #indeks vrstice, kjer je trenutno pivotni element

        for i in range(n):
            while pivot < m and left[pivot][i] % mod == 0:
                pivot += 1
            if pivot == m: #vse vrstice smo ze unicili
                break

            #eliminiramo ostale elemente v stolpcu in popravimo matrike
            for j in range(pivot + 1, m):
                coeff = left[j][i]
                left[j][i] = 0
                for k in range(len(basis[pivot])):
                    basis[pivot][k] = (basis[pivot][k] + coeff * basis[j][k]) % mod
                if right is not None:
                    for k in range(len(right)):
                        right[k][pivot] = (right[k][pivot] + coeff * right[k][j]) % mod

def tower_persistence(tower, tower_basis, max_filt_index):
    intervals = [] #zaenkrat podani kot [zacetni indeks, dolzina]
    interval_generators = [] #bazni element, ki predstavlja interval, dodan ob zacetku intervala
    active_intervals = dict() #intervali, ki so trenutno zivi
    #Kljuc: indeks baznega elementa, vrednost: indeks pripadajocega intervala v intervals
    for ind in range(max_filt_index - 1):
        M = tower[ind]
        basis = tower_basis[ind]
        if M is None:
            for i in range(len(basis)):
                if i not in active_intervals: #ta bazni element ni nadaljevanje intervala od prej
                    intervals.append([ind, 0]) #interval se takoj zakljuci, ker slikamo iz nicelnega ali v nicelni prostor
                    interval_generators.append(basis[i])

            active_intervals = dict()
        else:
            pivot = 0
            new_active_intervals = dict()
            for i in range(len(basis)): #stolpci matrike M
                while pivot < len(M) and M[pivot][i] == 0:
                    pivot += 1

                if pivot == len(M): #vsi preostali stolpci so nicelni
                    for j in range(i, len(basis)):
                        #preostali stolpci zakljucijo obstojec interval ali tvorijo novega trivialnega
                        if j not in active_intervals:
                            intervals.append([ind, 0])
                            interval_generators.append(basis[j])

                    break
                else:
                    if i in active_intervals: #podaljsamo interval
                        position = active_intervals[i]
                        intervals[position][1] += 1
                        new_active_intervals[pivot] = position
                    else: #ustvarimo nov netrivialen interval
                        position = len(intervals)
                        intervals.append([ind, 1])
                        interval_generators.append(basis[i])
                        new_active_intervals[pivot] = position

            active_intervals = new_active_intervals

    return intervals, interval_generators

def transform_intervals(intervals, max_filt_index):
    return [(birth, birth + time) if birth + time < max_filt_index - 1 else (birth, None) for birth, time in intervals]

def generators_to_cycles(intervals, interval_generators, domain_filt_basis, domain_basis, domain_index_to_simp, mod):
    cycles = [] #cikli so shranjeni kot slovarji, KLJUC: simplex, VREDNOST: koeficient pri simplexu v ciklu
    for i in range(len(intervals)):
        birth = intervals[i][0]
        basis = domain_filt_basis[birth] #seznam indeksov elementov v domain_basis, ki tvorijo bazo te filtracije
        coeffs = interval_generators[i] #koeficienti generatorji po bazi basis
        cycle = dict()

        for j in range(len(basis)):
            base_vect = domain_basis[basis[j]] #slovar, KLJUC: indeks simplexa, VREDNOST: koeficient
            base_coeff = coeffs[j] #lahko je tudi 0

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



if __name__ == '__main__':

    A = [[2, 4, 3, 0, 0], [8, 9, 6, 1, 2], [8, 0, 10, 7, 1], [0, 0, 0, 0, 6]]
    B = [[4, 1, 3], [2, -1, 3], [2, 1, 1], [1, 1, 2]]

    U = [[1, 0, 0], [0, 1, 0], [0, 0, 3]]
    V = [[2, 1, 0], [0, 0, 1]]

    tower = [[[4, 2, 3], [4, 1, 3], [5, 6, 10]], [[1, 3, 1], [4, 6, 10], [0, 0, 0]], [[2, 7, 3]]]
    tower_basis = [[[1], [2], [1]], [[4, 3], [3, 5], [2, 1]], [[1], [6], [10]], [[1, 2, 0, 1, 4]]]

