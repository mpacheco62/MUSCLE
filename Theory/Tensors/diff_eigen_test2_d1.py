from mpmath import mp
import pyperclip

# =====================================================================
# CONFIGURACIÓN DE PRECISIÓN
# =====================================================================
mp.dps = 150 

# El tensor que deseas evaluar (Componentes Voigt: xx, yy, zz, xy, yz, xz)
v1 = [mp.mpf('1.0'), mp.mpf('1.0'), mp.mpf('1.0')-mp.mpf('1E-1'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')]
# v1 = [mp.mpf('1.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')]
# S = mp.matrix(3, 3)
# S[0,0] = v1[0]; S[1,1] = v1[1]; S[2,2] = v1[2]
# S[0,1] = S[1,0] = v1[3]
# S[1,2] = S[2,1] = v1[4]
# S[0,2] = S[2,0] = v1[5]
# theta = mp.radians(15)
# omega = mp.radians(60)
# R1 = mp.matrix([[mp.cos(theta), mp.sin(theta), 0], [-mp.sin(theta), mp.cos(theta), 0], [0, 0, 1]])
# R2 = mp.matrix([[1, 0, 0], [0, mp.cos(omega), mp.sin(omega)], [0, -mp.sin(omega), mp.cos(omega)]])
# R = R2 @ R1
# S = R @ S @ R.T
# v1 = [S[0,0], S[1,1], S[2,2], S[0,1], S[1,2], S[0,2]]
# print("Tensor inicial (matriz 3x3):\n", S)
# print("Tensor inicial (matriz 3x3):\n", v1)
# print(f"{float(v1[0]):21.16f}D0  {float(v1[1]):21.16f}D0  {float(v1[2]):21.16f}D0  {float(v1[3]):21.16f}D0  {float(v1[4]):21.16f}D0  {float(v1[5]):21.16f}D0")


# Caso singular (descomentar para probar):
# v1 = [mp.mpf('1.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')]

def get_eigenvalues_mp(v):
    """
    Retorna los 3 autovalores calculados con precisión arbitraria.
    v es una lista de 6 componentes de Voigt.
    """
    S = mp.matrix(3, 3)
    S[0,0] = v[0]; S[1,1] = v[1]; S[2,2] = v[2]
    
    # Al perturbar el vector de Voigt, se altera la matriz simétricamente.
    # Esto es exactamente lo que hace ten_3D2Osym en Fortran.
    S[0,1] = S[1,0] = v[3]/mp.mpf('2.0')
    S[1,2] = S[2,1] = v[4]/mp.mpf('2.0')
    S[0,2] = S[2,0] = v[5]/mp.mpf('2.0')

    S2 = mp.matrix(3, 3)
    S2[0,0] = v1[0]; S2[1,1] = v1[1]; S2[2,2] = v1[2]
    S2[0,1] = S2[1,0] = v1[3]
    S2[1,2] = S2[2,1] = v1[4]
    S2[0,2] = S2[2,0] = v1[5]

    S = S2 + S  # Perturbación del tensor original
    
    # mp.eigsy calcula autovalores de matrices simétricas con la precisión definida
    E, _ = mp.eigsy(S)
    
    return [E[0], E[1], E[2]]

def compute_gradient_mp(v_list, a, eps_val='1e-20'):
    """
    Calcula el gradiente (1ra derivada) del autovalor 'a' usando diferencias finitas
    con un paso hiper-pequeño gracias a la precisión arbitraria.
    """
    eps = mp.mpf(eps_val)
    grad = mp.matrix(1, 6) # Tensor de 1x6 (Voigt)
    
    for i in range(6):
        # Derivada central
        v_p = list(v_list); v_p[i] += eps
        v_m = list(v_list); v_m[i] -= eps
        
        f_p = get_eigenvalues_mp(v_p)[a]
        f_m = get_eigenvalues_mp(v_m)[a]
        
        grad[0, i] = (f_p - f_m) / (2 * eps)
        
    return grad

def print_fortran_grad(g):
    """Imprime el vector g en formato Fortran (ten_3D2Osym) para copiar y pegar."""
    text =  f'    call expected%init(xx= {float(g[0]):20.14f}D0, yy= {float(g[1]):20.14f}D0, zz= {float(g[2]):20.14f}D0, &\n'
    text += f'                       xy= {float(g[3]):20.14f}D0, yz= {float(g[4]):20.14f}D0, xz= {float(g[5]):20.14f}D0  &\n'
    text += f'                       )'
    pyperclip.copy(text)
    print(text)

# =====================================================================
# EJECUCIÓN DEL TEST
# =====================================================================

print(f"Calculando PRIMERA DERIVADA con {mp.dps} decimales de precisión...")
print("Paso de diferencias finitas: epsilon = 10^-20")
print("-" * 65)

v = [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')]
eig = get_eigenvalues_mp(v)

# Calculamos los 3 gradientes
gradients = []
for a in range(3):
    g_a = compute_gradient_mp(v, a, eps_val='1e-20')
    gradients.append(g_a)

components = ["xx", "yy", "zz", "xy", "yz", "xz"]

for a in range(3):
    print(f"\n★★★ GRADIENTE DEL AUTOVALOR {a+1}: {eig[a]} ★★★")
    # Imprimimos redondeando a 6 decimales para que sea legible en consola
    row_str = "  ".join(f"{components[i]}: {float(gradients[a][0, i]):10.6f}" for i in range(6))
    print(f" [{row_str}]")

# Suma de Gradientes
grad_sum = mp.matrix(1, 6)
for a in range(3):
    grad_sum += gradients[a]

print("\n★★★ SUMA DE LOS TRES GRADIENTES (Debe ser la Identidad) ★★★")
row_str = "  ".join(f"{components[i]}: {float(grad_sum[0, i]):10.6f}" for i in range(6))
print(f" [{row_str}]")


# =====================================================================
# GENERACIÓN DE CÓDIGO FORTRAN
# =====================================================================
print("\n★★★ CÓDIGO FORTRAN (Copia automática al portapapeles) ★★★")

print("\n--- GRADIENTE DEL AUTOVALOR 3 ---")
print_fortran_grad(gradients[2])
input("\nPresiona Enter para copiar el gradiente del autovalor 2 al portapapeles...")

print("\n--- GRADIENTE DEL AUTOVALOR 2 ---")
print_fortran_grad(gradients[1])
input("\nPresiona Enter para copiar el gradiente del autovalor 1 al portapapeles...")

print("\n--- GRADIENTE DEL AUTOVALOR 1 ---")
print_fortran_grad(gradients[0])
input("\nPresiona Enter para copiar el PROMEDIO de los autovalores 1 y 2 al portapapeles...")

print("\n--- PROMEDIO DEL AUTOVALOR 1 y 2 (Útil para parche de raíces repetidas) ---")
print_fortran_grad((gradients[0] + gradients[1]) / 2)