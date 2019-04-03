import sympy as sp

def kernel(A, mod):
    #A je m x n matrika, NEPRAZNI DOMENA IN KODOMENA!
    #izracuna jedro via column operations
    m = len(A)
    n = len(A[0])

    L = [[A[i][j] % mod for j in range(n)] for i in range(m)] #spodnje trikotna na koncu
    U = [[0 for _ in range(n)] for _ in range(n)] #velja L = A * U, U je lahko polna matrika
    for i in range(n):
        U[i][i] = 1

    curr_line = 0 #vrstica, kjer trenutno želimo pridobiti nicle
    first_zero_column = None

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

def mult(A, B, mod):
    m = len(A)
    n = len(B)
    p = len(B[0])

    C = [[0 for _ in range(p)] for _ in range(m)]
    for i in range(m):
        for j in range(p):
            for k in range(n):
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % mod

    return C

#TODO: intersection basis
#TODO: razvij po koeficientih baze


A = [[2, 4, 3, 0, 0], [8, 9, 6, 1, 2], [8, 0, 10, 7, 1], [0, 0, 0, 0, 6]]
B = [[4, 1, 3], [2, -1, 3], [2, 1, 1], [1, 1, 2]]
