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
    #vrne A + t * B, obe velikosti m x n NEPRAZNI!
    m = len(A)
    n = len(A[0])

    return [[(A[i][j] + t * B[i][j]) % mod for j in range(n)] for i in range(m)]

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

def complete_basis(A, B, mod, empty = False):
    #A je podprostor B, dopolni bazo za A (jo popravi do spodnje trikotne) do baze za B
    #pri tem skrbi, da je baza vedno spodaj trikotna, da bomo lazje razvijali po tej bazi kasneje
    #empty pove, ce je slucajno A trivialen podprostor, B NEPRAZEN!
    p = len(B[0])
    #na indeksu i je bazni vektor, ki ima na i-tem mestu od zgoraj prvi nenicelni element
    top_non_zero = [([], None) for _ in range(p)]
    quotient_basis = [] #indeksi, na katerih so elementi, ki smo jih dobili iz B

    if not empty:
        for vect in A:
            vect_reduced, top, _ = reduce_vect(top_non_zero, vect, mod)
            if top is not None: #se ni zreduciral do nicelnega vektorja
                top_non_zero[top] = (vect_reduced, 'A')

    for i in range(len(B)):
        vect_reduced, top, _ = reduce_vect(top_non_zero, B[i], mod)
        if top is not None:
            top_non_zero[top] = (vect_reduced, 'B')
            quotient_basis.append(top)

    return top_non_zero, quotient_basis

def quotient_coeffs(top_list, quotient_basis, vect, mod):
    _, _, coeffs = reduce_vect(top_list, vect, mod)
    return [coeffs[i] for i in quotient_basis]


if __name__ == '__main__':

    A = [[2, 4, 3, 0, 0], [8, 9, 6, 1, 2], [8, 0, 10, 7, 1], [0, 0, 0, 0, 6]]
    B = [[4, 1, 3], [2, -1, 3], [2, 1, 1], [1, 1, 2]]

    U = [[1, 0, 0], [0, 1, 0], [0, 0, 3]]
    V = [[2, 1, 0], [0, 0, 1]]

