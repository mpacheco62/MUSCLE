import numpy as np
import numdifftools as nd

def voigt_to_matrix(v):
    """Convierte vector Voigt a matriz simétrica 3x3"""
    S = np.zeros((3, 3))
    S[0, 0] = v[0] # xx
    S[1, 1] = v[1] # yy
    S[2, 2] = v[2] # zz
    S[0, 1] = S[1, 0] = v[3] # xy
    S[1, 2] = S[2, 1] = v[4] # yz
    S[0, 2] = S[2, 0] = v[5] # xz
    return S

# 1. TENSOR DE PRUEBA (Componentes Voigt: xx, yy, zz, xy, yz, xz)
Sigma_test = np.array([1.0, 1E-2, -1E-2, 0.0, 0.0, 0.0])

theta = np.radians(15)  # Ángulo de rotación
omega = np.radians(60)  # Ángulo de rotación

S = voigt_to_matrix(Sigma_test)  # Tensor inicial en formato matriz
# R1 = np.array([[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
# R2 = np.array([[1, 0, 0], [0, np.cos(omega), np.sin(omega)], [0, -np.sin(omega), np.cos(omega)]])
# R = R2 @ R1  # Rotación compuesta
# S = R @ S @ R.T  # Rotación del tensor
Sigma_test = np.array([S[0,0], S[1,1], S[2,2], S[0,1], S[1,2], S[0,2]])
print("Tensor inicial (matriz 3x3):\n", S)


def get_eigenvalues(v):
    """Retorna los 3 autovalores como un vector (salida vectorial)"""
    S = voigt_to_matrix(Sigma_test)
    v2 = v.copy()
    v2[3:] /= 2.0  # Ajuste para componentes de corte
    S2 = S + voigt_to_matrix(v2)  # Perturbación del tensor
    return np.linalg.eigvalsh(S2)

def get_first_derivative(v):
    S = voigt_to_matrix(Sigma_test)
    v2 = v.copy()
    v2[3:] /= 2.0  # Ajuste para componentes de corte
    S2 = S + voigt_to_matrix(v2)  # Perturbación del tensor
    I1 = S2[0,0] + S2[1,1] + S2[2,2]
    I2 = 0.5 * (I1**2 - np.einsum('ij,ij->', S2, S2))
    eig_values = np.linalg.eigvalsh(S2)

    num_sing = 0
    good_idx = 1
    sing_idx = np.zeros((3,), dtype=int)
    res = np.zeros((3, 3, 3))
    for i in range(3):
        Sk = eig_values[i]
        D_scal = 3.0 * Sk**2 - 2*I1*Sk + I2
        if (abs(D_scal) > 1e-12):
            F = np.einsum('ij,jk->ik', S2, S2) - (I1-Sk)*S2 + (Sk**2 - I1*Sk + I2)*np.eye(3)
            res[i] = F / D_scal
            good_idx = i
        else:
            sing_idx[num_sing] = i
            num_sing += 1

    if num_sing == 2:
        F = 0.5*(np.eye(3) - res[good_idx])
        res[sing_idx[0]] = F
        res[sing_idx[1]] = F
    elif num_sing == 3:
        F = (1/3)*np.eye(3)
        res[0] = F
        res[1] = F
        res[2] = F

    res2 = np.zeros((3, 6))
    for i in range(3):
        res2[i, 0] = res[i, 0, 0] # d_lambda/d_xx
        res2[i, 1] = res[i, 1, 1] # d_lambda/d_yy
        res2[i, 2] = res[i, 2, 2] # d_lambda/d_zz
        res2[i, 3] = res[i, 0, 1] # d_lambda/d_xy
        res2[i, 4] = res[i, 1, 2] # d_lambda/d_yz
        res2[i, 5] = res[i, 0, 2] # d_lambda/d_xz
    return res2

# =====================================================================
# CÁLCULO AUTOMÁTICO CON NUMDIFFTOOLS
# =====================================================================

# 1. PRIMERA DERIVADA (Gradiente de 3x6):
# Como 'get_eigenvalues' devuelve 3 valores, su primera derivada es el Jacobiano.
S_num = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # Tensor de perturbación en formato Voigt
jacobian_fun = nd.Jacobian(get_eigenvalues)
grad_automatico = jacobian_fun(S_num)  # Devuelve una matriz de 3x6
grad_analitico = np.zeros((3, 6))
grad_analitico = get_first_derivative(S_num) # Cálculo analítico para comparación


# 2. SEGUNDA DERIVADA (Hessiana de 6x6 para cada autovalor):
# Como la Hessiana requiere una función que devuelva un solo escalar,
# creamos una función para cada autovalor (a = 0, 1, 2).
hessian_automatico = []
hessian_automatico_der = []
for a in range(3):
    # Función que retorna solo el autovalor 'a'
    f_scalar = lambda v: get_eigenvalues(v)[a]
    f_scalar2 = lambda v: get_first_derivative(v)[a]
    
    hessian_fun = nd.Hessian(f_scalar)
    hessian_automatico.append(hessian_fun(S_num))  # Guarda la matriz de 6x6

    jacobian_fun = nd.Jacobian(f_scalar2)
    hessian_automatico_der.append(jacobian_fun(S_num))  # Guarda la matriz de 6x6 de la derivada del gradiente


# =====================================================================
# IMPRESIÓN DE RESULTADOS
# =====================================================================
np.set_printoptions(precision=6, suppress=True)
components = ["xx", "yy", "zz", "xy", "yz", "xz"]

print("=================================================================")
print("          RESULTADOS AUTOMÁTICOS CON NUMDIFFTOOLS")
print("=================================================================\n")

eigenvalues = get_eigenvalues(S_num)

for a in range(3):
    print(f"★★★ AUTOVALOR {a+1} {eigenvalues[a]:.6f} ★★★")
    
    # Gradiente
    print("\n  -> GRADIENTE:")
    for i in range(6):
        print(f"     d_lambda / d_{components[i]:<2} = {grad_automatico[a, i]:.8f}")

    print("\n  -> GRADIENTE ANALÍTICO:")
    for i in range(6):
        print(f"     d_lambda / d_{components[i]:<2} = {grad_analitico[a, i]:.8f}")

    # Diferencia con analítico
    print("\n  -> DIFERENCIA CON ANALÍTICO:")
    for i in range(6):
        diff = grad_automatico[a, i] - grad_analitico[a, i]
        print(f"     d_lambda / d_{components[i]:<2} = {diff: .2e}")
        
    # Hessiana
    print("\n  -> HESSIANA:")
    print("           xx                  yy                  zz                  xy                  yz                  xz")
    for i in range(6):
        row_str = "  ".join(f"{val:18.12f}" for val in hessian_automatico[a][i, :])
        print(f" {components[i]:<2} [{row_str}]")

    print("\n  -> HESSIANA2:")
    print("           xx                  yy                  zz                  xy                  yz                  xz")
    for i in range(6):
        row_str = "  ".join(f"{val:18.12f}" for val in hessian_automatico_der[a][i, :])
        print(f" {components[i]:<2} [{row_str}]")
    print("\n" + "="*65)


print(" Sum of all Hessians (should be close to zero if consistent):")
hessian_sum = np.zeros((6, 6))
hessian_sum2 = np.zeros((6, 6))
for a in range(3):
    hessian_sum += hessian_automatico[a]
    hessian_sum2 += hessian_automatico_der[a]
print("Hessian sum from eigenvalues:")
print("           xx                  yy                  zz                  xy                  yz                  xz")
for i in range(6):
    row_str = "  ".join(f"{val:18.12f}" for val in hessian_sum[i, :])
    print(f" {components[i]:<2} [{row_str}]")

print("Hessian sum from first derivatives:")
print("           xx                  yy                  zz                  xy                  yz                  xz")
for i in range(6):
    row_str = "  ".join(f"{val:18.12f}" for val in hessian_sum2[i, :])
    print(f" {components[i]:<2} [{row_str}]")