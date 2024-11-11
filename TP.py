# Constantes y condiciones iniciales
G = 6.67430e-11  # Constante gravitacional en m^3 kg^-1 s^-2
M = 5.972e24  # Masa de la Tierra (en kg)
R_SAT = 2820000000  # Posición radial inicial (en metros)
pos_radial_asteroide = 4122500000  # Posición radial del asteroide (en metros)
pos_angular_asteroide = 3.025  # Posición angular del asteroide (en radianes)

# Condiciones iniciales
s_0 = 0  # Velocidad radial inicial

def ds_dt(s, omega, r):
    # Derivada de la velocidad radial (aceleración radial)
    return r * omega**2 - G * M / r**2  # Aceleración centrífuga menos gravitacional

def domega_dt(s, omega, r):
    # Derivada de omega
    return -2 * s * omega / r  # Ecuación de movimiento para la velocidad angular

# Método de Runge-Kutta de 2º orden
def runge_kutta_2(s, omega, r, delta_t):
    # Calcular k1 y l1
    k1_s = delta_t * ds_dt(s, omega, r)
    k1_omega = delta_t * domega_dt(s, omega, r)

    # Calcular k2 y l2
    k2_s = delta_t * ds_dt(s + k1_s / 2, omega + k1_omega / 2, r)
    k2_omega = delta_t * domega_dt(s + k1_s / 2, omega + k1_omega / 2, r)

    # Actualizar s y omega usando promedio de k1, k2 y l1, l2
    s_next = s + (k1_s + k2_s) / 2
    omega_next = omega + (k1_omega + k2_omega) / 2

    return s_next, omega_next

# Método de Adams-Bashforth de 3er orden
def adams_bashforth_3(s_0, omega_0, r_0, delta_t, t_end):
    # Inicialización
    s_vals = [s_0]
    omega_vals = [omega_0]
    r_vals = [r_0]

    # Paso preliminar con Runge-Kutta 2º orden (3 pasos para iniciar Adams-Bashforth)
    for _ in range(3):
        s_next, omega_next = runge_kutta_2(s_vals[-1], omega_vals[-1], r_vals[-1], delta_t)
        s_vals.append(s_next)
        omega_vals.append(omega_next)
        r_vals.append(r_vals[-1] + s_next * delta_t)  # Actualizar r

    # Usamos Adams-Bashforth 3er orden para el resto de los pasos
    t = 3
    while t * delta_t <= t_end:
        # Siguiente valor de s usando Adams-Bashforth
        s_next = (s_vals[t] + delta_t / 12 * (
                23 * ds_dt(s_vals[t], omega_vals[t], r_vals[t]) -
                16 * ds_dt(s_vals[t - 1], omega_vals[t - 1], r_vals[t - 1]) +
                5 * ds_dt(s_vals[t - 2], omega_vals[t - 2], r_vals[t - 2])
        ))

        # Siguiente valor de omega usando Adams-Bashforth
        omega_next = (omega_vals[t] + delta_t / 12 * (
                23 * domega_dt(s_vals[t], omega_vals[t], r_vals[t]) -
                16 * domega_dt(s_vals[t - 1], omega_vals[t - 1], r_vals[t - 1]) +
                5 * domega_dt(s_vals[t - 2], omega_vals[t - 2], r_vals[t - 2])
        ))

        # Agregar resultados al historial
        s_vals.append(s_next)
        omega_vals.append(omega_next)
        r_vals.append(r_vals[-1] + s_next * delta_t)  # Actualizar r

        t += 1

    return s_vals, omega_vals, r_vals

# Función f(v0) para evaluar la diferencia entre la posición final del proyectil y el asteroide
def f(v0):
    delta_t = 10  # Paso de tiempo en segundos
    t_end = 1e6  # Tiempo total de simulación en segundos

    # Condición inicial omega
    omega_0 = v0 / R_SAT

    # Llamamos a la función Adams-Bashforth para simular la trayectoria
    s_vals, omega_vals, r_vals = adams_bashforth_3(s_0, omega_0, R_SAT, delta_t, t_end)

    # Calculamos la diferencia radial
    final_s = r_vals[-1]

    return final_s - pos_radial_asteroide

# Derivada de f(v0) para el método de Newton-Raphson
def df(v0):
    epsilon = 1e-6  # Un pequeño incremento para aproximar la derivada numéricamente
    return (f(v0 + epsilon) - f(v0)) / epsilon

# Método de Newton-Raphson para encontrar la velocidad inicial óptima
def newton_raphson(f, df, v0_init, tol=1e-6, max_iter=1000):
    v0 = v0_init
    for _ in range(max_iter):
        v0_next = v0 - f(v0) / df(v0)
        if abs(v0_next - v0) < tol:
            return v0_next
        v0 = v0_next
    raise ValueError("No se encontró una solución dentro del número máximo de iteraciones")

# Valor inicial para el método de Newton-Raphson
v0_init = 1000  # Valor inicial para la velocidad tentativa

# Ejecutamos el método de Newton-Raphson
try:
    v0_optimo = newton_raphson(f, df, v0_init)
    print(f"Velocidad inicial óptima para impactar el asteroide: {v0_optimo:.2f} m/s")
except ValueError as e:
    print(e)
