import numpy as np

def construir_tensor_superdiagonales(vals):
    """
    Construye un tensor de 4to orden (3x3x3x3) con 3 simetrías.
    
    Orden de entrada vals(1:21):
    - Diagonales: (11, 22, 33, 44, 55, 66) -> indices 0:6
    - Super-diag 1: (12, 23, 34, 45, 56)  -> indices 6:11
    - Super-diag 2: (13, 24, 35, 46)      -> indices 11:15
    - Super-diag 3: (14, 25, 36)          -> indices 15:18
    - Super-diag 4: (15, 26)              -> indices 18:20
    - Super-diag 5: (16)                  -> index 20
    """
    if len(vals) != 21:
        raise ValueError("Se necesitan exactamente 21 valores.")

    # 1. Creamos la matriz 6x6
    M = np.zeros((6, 6))
    
    # Llenar la diagonal principal
    for i in range(6):
        M[i, i] = vals[i]
        
    # Llenar las super-diagonales según tu orden específico
    idx_actual = 6
    for k in range(1, 6): # k es el offset de la diagonal
        for i in range(6 - k):
            M[i, i + k] = vals[idx_actual]
            idx_actual += 1
            
    # Hacerla simétrica (Simetría Mayor)
    M = M + M.T - np.diag(M.diagonal())

    # 2. Mapeo Voigt estándar para pasar de 6x6 a 3x3x3x3
    # 1->xx(0,0), 2->yy(1,1), 3->zz(2,2), 4->xy(0,1), 5->yz(1,2), 6->xz(0,2)
    voigt_map = {
        0: (0, 0), 1: (1, 1), 2: (2, 2),
        3: (0, 1), 4: (1, 2), 5: (0, 2)
    }

    # 3. Expandir al tensor de 4to orden
    C = np.zeros((3, 3, 3, 3))
    for I in range(6):
        for J in range(6):
            i, j = voigt_map[I]
            k, l = voigt_map[J]
            val = M[I, J]
            
            # Aplicar simetrías menores
            C[i, j, k, l] = val
            C[j, i, k, l] = val
            C[i, j, l, k] = val
            C[j, i, l, k] = val
            
    return C

# --- Ejemplo de Test ---
# Valores del 1 al 21 para identificar cada posición fácilmente
datos_test = np.arange(1, 22, dtype=float)
C_tensor = construir_tensor_superdiagonales(datos_test)

test1 = np.zeros((6,6))
test1[0,0] = C_tensor[0,0,0,0]
test1[0,1] = C_tensor[0,0,1,1]
test1[0,2] = C_tensor[0,0,2,2]
test1[0,3] = C_tensor[0,0,0,1]
test1[0,4] = C_tensor[0,0,1,2]
test1[0,5] = C_tensor[0,0,0,2]

test1[1,0] = C_tensor[1,1,0,0]
test1[1,1] = C_tensor[1,1,1,1]
test1[1,2] = C_tensor[1,1,2,2]
test1[1,3] = C_tensor[1,1,0,1]
test1[1,4] = C_tensor[1,1,1,2]
test1[1,5] = C_tensor[1,1,0,2]

test1[2,0] = C_tensor[2,2,0,0]
test1[2,1] = C_tensor[2,2,1,1]
test1[2,2] = C_tensor[2,2,2,2]
test1[2,3] = C_tensor[2,2,0,1]
test1[2,4] = C_tensor[2,2,1,2]
test1[2,5] = C_tensor[2,2,0,2]

test1[3,0] = C_tensor[0,1,0,0]
test1[3,1] = C_tensor[0,1,1,1]
test1[3,2] = C_tensor[0,1,2,2]
test1[3,3] = C_tensor[0,1,0,1]
test1[3,4] = C_tensor[0,1,1,2]
test1[3,5] = C_tensor[0,1,0,2]

test1[4,0] = C_tensor[1,2,0,0]
test1[4,1] = C_tensor[1,2,1,1]
test1[4,2] = C_tensor[1,2,2,2]
test1[4,3] = C_tensor[1,2,0,1]
test1[4,4] = C_tensor[1,2,1,2]
test1[4,5] = C_tensor[1,2,0,2]

test1[5,0] = C_tensor[0,2,0,0]
test1[5,1] = C_tensor[0,2,1,1]
test1[5,2] = C_tensor[0,2,2,2]
test1[5,3] = C_tensor[0,2,0,1]
test1[5,4] = C_tensor[0,2,1,2]
test1[5,5] = C_tensor[0,2,0,2]
print(test1)


test2 = np.zeros((6,6))
test2[0,0] = C_tensor[0,0,0,0]
test2[0,1] = C_tensor[0,0,1,1]
test2[0,2] = C_tensor[0,0,2,2]
test2[0,3] = C_tensor[0,0,1,0]
test2[0,4] = C_tensor[0,0,2,1]
test2[0,5] = C_tensor[0,0,2,0]

test2[1,0] = C_tensor[1,1,0,0]
test2[1,1] = C_tensor[1,1,1,1]
test2[1,2] = C_tensor[1,1,2,2]
test2[1,3] = C_tensor[1,1,1,0]
test2[1,4] = C_tensor[1,1,2,1]
test2[1,5] = C_tensor[1,1,2,0]

test2[2,0] = C_tensor[2,2,0,0]
test2[2,1] = C_tensor[2,2,1,1]
test2[2,2] = C_tensor[2,2,2,2]
test2[2,3] = C_tensor[2,2,1,0]
test2[2,4] = C_tensor[2,2,2,1]
test2[2,5] = C_tensor[2,2,2,0]

test2[3,0] = C_tensor[0,1,0,0]
test2[3,1] = C_tensor[0,1,1,1]
test2[3,2] = C_tensor[0,1,2,2]
test2[3,3] = C_tensor[0,1,1,0]
test2[3,4] = C_tensor[0,1,2,1]
test2[3,5] = C_tensor[0,1,2,0]

test2[4,0] = C_tensor[1,2,0,0]
test2[4,1] = C_tensor[1,2,1,1]
test2[4,2] = C_tensor[1,2,2,2]
test2[4,3] = C_tensor[1,2,1,0]
test2[4,4] = C_tensor[1,2,2,1]
test2[4,5] = C_tensor[1,2,2,0]

test2[5,0] = C_tensor[0,2,0,0]
test2[5,1] = C_tensor[0,2,1,1]
test2[5,2] = C_tensor[0,2,2,2]
test2[5,3] = C_tensor[0,2,1,0]
test2[5,4] = C_tensor[0,2,2,1]
test2[5,5] = C_tensor[0,2,2,0]
print(test2)


test3 = np.zeros((6,6))
test3[0,0] = C_tensor[0,0,0,0]
test3[0,1] = C_tensor[0,0,1,1]
test3[0,2] = C_tensor[0,0,2,2]
test3[0,3] = C_tensor[0,0,1,0]
test3[0,4] = C_tensor[0,0,2,1]
test3[0,5] = C_tensor[0,0,2,0]

test3[1,0] = C_tensor[1,1,0,0]
test3[1,1] = C_tensor[1,1,1,1]
test3[1,2] = C_tensor[1,1,2,2]
test3[1,3] = C_tensor[1,1,1,0]
test3[1,4] = C_tensor[1,1,2,1]
test3[1,5] = C_tensor[1,1,2,0]

test3[2,0] = C_tensor[2,2,0,0]
test3[2,1] = C_tensor[2,2,1,1]
test3[2,2] = C_tensor[2,2,2,2]
test3[2,3] = C_tensor[2,2,1,0]
test3[2,4] = C_tensor[2,2,2,1]
test3[2,5] = C_tensor[2,2,2,0]

test3[3,0] = C_tensor[1,0,0,0]
test3[3,1] = C_tensor[1,0,1,1]
test3[3,2] = C_tensor[1,0,2,2]
test3[3,3] = C_tensor[1,0,1,0]
test3[3,4] = C_tensor[1,0,2,1]
test3[3,5] = C_tensor[1,0,2,0]

test3[4,0] = C_tensor[2,1,0,0]
test3[4,1] = C_tensor[2,1,1,1]
test3[4,2] = C_tensor[2,1,2,2]
test3[4,3] = C_tensor[2,1,1,0]
test3[4,4] = C_tensor[2,1,2,1]
test3[4,5] = C_tensor[2,1,2,0]

test3[5,0] = C_tensor[2,0,0,0]
test3[5,1] = C_tensor[2,0,1,1]
test3[5,2] = C_tensor[2,0,2,2]
test3[5,3] = C_tensor[2,0,1,0]
test3[5,4] = C_tensor[2,0,2,1]
test3[5,5] = C_tensor[2,0,2,0]
print(test3)


test4 = np.zeros((6,6))
test4[0,0] = C_tensor[0,0,0,0]
test4[0,1] = C_tensor[0,0,1,1]
test4[0,2] = C_tensor[0,0,2,2]
test4[0,3] = C_tensor[0,0,0,1]
test4[0,4] = C_tensor[0,0,1,2]
test4[0,5] = C_tensor[0,0,0,2]

test4[1,0] = C_tensor[1,1,0,0]
test4[1,1] = C_tensor[1,1,1,1]
test4[1,2] = C_tensor[1,1,2,2]
test4[1,3] = C_tensor[1,1,0,1]
test4[1,4] = C_tensor[1,1,1,2]
test4[1,5] = C_tensor[1,1,0,2]

test4[2,0] = C_tensor[2,2,0,0]
test4[2,1] = C_tensor[2,2,1,1]
test4[2,2] = C_tensor[2,2,2,2]
test4[2,3] = C_tensor[2,2,0,1]
test4[2,4] = C_tensor[2,2,1,2]
test4[2,5] = C_tensor[2,2,0,2]

test4[3,0] = C_tensor[1,0,0,0]
test4[3,1] = C_tensor[1,0,1,1]
test4[3,2] = C_tensor[1,0,2,2]
test4[3,3] = C_tensor[1,0,0,1]
test4[3,4] = C_tensor[1,0,1,2]
test4[3,5] = C_tensor[1,0,0,2]

test4[4,0] = C_tensor[2,1,0,0]
test4[4,1] = C_tensor[2,1,1,1]
test4[4,2] = C_tensor[2,1,2,2]
test4[4,3] = C_tensor[2,1,0,1]
test4[4,4] = C_tensor[2,1,1,2]
test4[4,5] = C_tensor[2,1,0,2]

test4[5,0] = C_tensor[2,0,0,0]
test4[5,1] = C_tensor[2,0,1,1]
test4[5,2] = C_tensor[2,0,2,2]
test4[5,3] = C_tensor[2,0,0,1]
test4[5,4] = C_tensor[2,0,1,2]
test4[5,5] = C_tensor[2,0,0,2]
print(test4)
print("¿Son iguales test1, test2, test3 y test4? Deberían ser iguales debido a las simetrías.")
print(test1-test2)
print(test1-test3)
print(test1-test4)


print("Tensor de 4to orden con simetrías construido a partir de vals(1:21):")
# # Verificación de una posición: vals(7) es la posición (1,2) de la matriz 6x6
# # En términos tensoriales Voigt, eso es M[0,1] -> C[0,0,1,1]
# print(f"Valor esperado en C[0,0,1,1] (vals 7): {C_tensor[0,0,1,1]}")



print("*************************************")

A = np.arange(1, 22, dtype=float)
print(A)
A = construir_tensor_superdiagonales(A)

B = np.arange(21, 0, -1, dtype=float)
print(B)
B = construir_tensor_superdiagonales(B)

C = np.einsum('ijkl,klmn->ijmn', A, B)
print("Resultado de la contracción A : B (C[i,j,m,n] = A[i,j,k,l] * B[k,l,m,n]):")

C_mat = np.zeros((6,6))
C_mat[0,0] = C[0,0,0,0]
C_mat[0,1] = C[0,0,1,1]
C_mat[0,2] = C[0,0,2,2]
C_mat[0,3] = C[0,0,0,1]
C_mat[0,4] = C[0,0,1,2]
C_mat[0,5] = C[0,0,0,2]

C_mat[1,0] = C[1,1,0,0]
C_mat[1,1] = C[1,1,1,1]
C_mat[1,2] = C[1,1,2,2]
C_mat[1,3] = C[1,1,0,1]
C_mat[1,4] = C[1,1,1,2]
C_mat[1,5] = C[1,1,0,2]

C_mat[2,0] = C[2,2,0,0]
C_mat[2,1] = C[2,2,1,1]
C_mat[2,2] = C[2,2,2,2]
C_mat[2,3] = C[2,2,0,1]
C_mat[2,4] = C[2,2,1,2]
C_mat[2,5] = C[2,2,0,2]

C_mat[3,0] = C[0,1,0,0]
C_mat[3,1] = C[0,1,1,1]
C_mat[3,2] = C[0,1,2,2]
C_mat[3,3] = C[0,1,0,1]
C_mat[3,4] = C[0,1,1,2]
C_mat[3,5] = C[0,1,0,2]

C_mat[4,0] = C[1,2,0,0]
C_mat[4,1] = C[1,2,1,1]
C_mat[4,2] = C[1,2,2,2]
C_mat[4,3] = C[1,2,0,1]
C_mat[4,4] = C[1,2,1,2]
C_mat[4,5] = C[1,2,0,2]

C_mat[5,0] = C[0,2,0,0]
C_mat[5,1] = C[0,2,1,1]
C_mat[5,2] = C[0,2,2,2]
C_mat[5,3] = C[0,2,0,1]
C_mat[5,4] = C[0,2,1,2]
C_mat[5,5] = C[0,2,0,2]
print(C_mat)