from mpmath import mp
import pyperclip

# =====================================================================
# CONFIGURACIÓN DE PRECISIÓN
# =====================================================================
# Definimos 50 dígitos decimales de precisión (mucho mayor a 128 bits)
mp.dps = 150 
e_val = mp.mpf('2e-0')
v1 = [mp.mpf('1.0'), e_val, -e_val, mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')]
v1 = [mp.mpf('1.0'), mp.mpf('2.0'), mp.mpf('3.0'), mp.mpf('4.0'), mp.mpf('5.0'), mp.mpf('6.0')]

def get_eigenvalues_mp(v):
    """
    Retorna los 3 autovalores calculados con precisión arbitraria.
    v es una lista de 6 componentes de Voigt.
    """
    S1 = mp.matrix(3, 3)
    S1[0,0] = v1[0]; S1[1,1] = v1[1]; S1[2,2] = v1[2]
    S1[0,1] = S1[1,0] = v1[3]
    S1[1,2] = S1[2,1] = v1[4]
    S1[0,2] = S1[2,0] = v1[5]

    # Construimos la matriz 3x3 en mpmath
    S2 = mp.matrix(3, 3)
    S2[0,0] = v[0]; S2[1,1] = v[1]; S2[2,2] = v[2]
    
    # Ajuste Voigt para componentes de corte (dividiendo por 2)
    S2[0,1] = S2[1,0] = v[3] / mp.mpf('2.0')
    S2[1,2] = S2[2,1] = v[4] / mp.mpf('2.0')
    S2[0,2] = S2[2,0] = v[5] / mp.mpf('2.0')

    S = S1 + S2
    
    # mp.eigsy calcula autovalores de matrices simétricas con la precisión definida
    E, _ = mp.eigsy(S)
    
    return [E[0], E[1], E[2]]

def compute_hessian_mp(v_list, a, eps_val='1e-20'):
    """
    Calcula la Hessiana del autovalor 'a' usando diferencias finitas
    con un paso hiper-pequeño gracias a la precisión arbitraria.
    """
    eps = mp.mpf(eps_val)
    H = mp.matrix(6, 6)
    
    for i in range(6):
        for j in range(6):
            if i == j:
                # Derivada segunda diagonal
                v_p = list(v_list); v_p[i] += eps
                v_m = list(v_list); v_m[i] -= eps
                
                f_p = get_eigenvalues_mp(v_p)[a]
                f_c = get_eigenvalues_mp(v_list)[a]
                f_m = get_eigenvalues_mp(v_m)[a]
                
                H[i,i] = (f_p - 2*f_c + f_m) / (eps**2)
            elif i < j:
                # Derivada segunda cruzada
                v_pp = list(v_list); v_pp[i] += eps; v_pp[j] += eps
                v_pm = list(v_list); v_pm[i] += eps; v_pm[j] -= eps
                v_mp = list(v_list); v_mp[i] -= eps; v_mp[j] += eps
                v_mm = list(v_list); v_mm[i] -= eps; v_mm[j] -= eps
                
                f_pp = get_eigenvalues_mp(v_pp)[a]
                f_pm = get_eigenvalues_mp(v_pm)[a]
                f_mp = get_eigenvalues_mp(v_mp)[a]
                f_mm = get_eigenvalues_mp(v_mm)[a]
                
                val = (f_pp - f_pm - f_mp + f_mm) / (4 * eps**2)
                H[i,j] = val
                H[j,i] = val
    return H

def print_fortran(m):
    """Imprime la matriz m en formato Fortran para copiar y pegar fácilmente."""
    text = f'    call expected%init(xxxx= {float(m[0,0]):20.14f}D0, yyyy= {float(m[1,1]):20.14f}D0, zzzz= {float(m[2,2]):20.14f}D0, &\n'
    text += f'                       xyxy= {float(m[3,3]):20.14f}D0, yzyz= {float(m[4,4]):20.14f}D0, xzxz= {float(m[5,5]):20.14f}D0, &\n'
    text += f'                       xxyy= {float(m[0,1]):20.14f}D0, yyzz= {float(m[1,2]):20.14f}D0,                               &\n'
    text += f'                       zzxy= {float(m[2,3]):20.14f}D0, xyyz= {float(m[3,4]):20.14f}D0, yzxz= {float(m[4,5]):20.14f}D0, &\n'
    text += f'                       xxzz= {float(m[0,2]):20.14f}D0,                                                             &\n'
    text += f'                       yyxy= {float(m[1,3]):20.14f}D0, zzyz= {float(m[2,4]):20.14f}D0, xyxz= {float(m[3,5]):20.14f}D0, &\n'
    text += f'                       xxxy= {float(m[0,3]):20.14f}D0, yyyz= {float(m[1,4]):20.14f}D0, zzxz= {float(m[2,5]):20.14f}D0, &\n'
    text += f'                       xxyz= {float(m[0,4]):20.14f}D0, yyxz= {float(m[1,5]):20.14f}D0, xxxz= {float(m[0,5]):20.14f}D0  &\n'
    text += f'                       )'
    pyperclip.copy(text)
    print(text)

# =====================================================================
# EJECUCIÓN DEL TEST
# =====================================================================
# Nuestro tensor con una perturbación minúscula pero exacta
e_val = mp.mpf('1e-3')
Sigma_test = [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')]

print(f"Calculando con {mp.dps} decimales de precisión...")
print("Paso de diferencias finitas: epsilon = 10^-20")
print("-" * 65)
eig = get_eigenvalues_mp(Sigma_test)

# Calculamos las 3 Hessianas
hessians = []
for a in range(3):
    H_a = compute_hessian_mp(Sigma_test, a, eps_val='1e-20')
    hessians.append(H_a)

components = ["xx", "yy", "zz", "xy", "yz", "xz"]

for a in range(3):
    print(f"\n★★★ HESSIANA DEL AUTOVALOR {a+1}: {eig[a]} ★★★")
    print("           xx            yy            zz            xy            yz            xz")
    for i in range(6):
        # Imprimimos redondeando a 6 decimales para que sea legible
        row_str = "  ".join(f"{float(val):12.6f}" for val in hessians[a][i, :])
        print(f" {components[i]:<2} [{row_str}]")

# Suma de Hessianas
hessian_sum = mp.matrix(6, 6)
for a in range(3):
    hessian_sum += hessians[a]

print("\n★★★ SUMA DE LAS TRES HESSIANAS (Precisión Extrema) ★★★")
print("           xx            yy            zz            xy            yz            xz")
for i in range(6):
    # Imprimimos en formato científico para ver los ceros absolutos
    row_str = "  ".join(f"{float(val):12.2e}" for val in hessian_sum[i, :])
    print(f" {components[i]:<2} [{row_str}]")

print(f"\n★★★ HESSIANA DEL AUTOVALOR: ★★★")
print("           xx            yy            zz            xy            yz            xz")
for i in range(6):
    # Imprimimos redondeando a 6 decimales para que sea legible
    to_print = hessians[2] + hessians[1] + hessians[0]
    to_print /= 3
    row_str = "  ".join(f"{float(val):12.6f}" for val in to_print[i, :])
    print(f" {components[i]:<2} [{row_str}]")

print("\n★★★ HESSIANA DEL AUTOVALOR 3 ★★★")
print_fortran(hessians[2])
input("\nPresiona Enter para copiar la Hessiana del autovalor 2 al portapapeles...")

print("\n★★★ HESSIANA DEL AUTOVALOR 2 ★★★")
print_fortran(hessians[1])
input("\nPresiona Enter para copiar la Hessiana del autovalor 1 al portapapeles...")

print("\n★★★ HESSIANA DEL AUTOVALOR 1 ★★★")
print_fortran(hessians[0])
input("\nPresiona Enter para copiar la Hessiana el promedio del autovalor 1 y 2 al portapapeles...")

print("\n★★★ HESSIANA DEL AUTOVALOR 1 y 2 ★★★")
print_fortran((hessians[0] + hessians[1]) / 2)

