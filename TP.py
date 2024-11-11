import numpy as np

# Constantes y condiciones iniciales
G = 6.67430e-11  # Constante gravitacional en m^3 kg^-1 s^-2
M = 7.342e22  # Masa de la luna (en kg)
R_SAT = 2.82e10  # Posición radial inicial (en metros)
pos_radial_asteroide = 4.1225e10  # Posición radial del asteroide (en metros)
pos_angular_asteroide = 3.025  # Posición angular del asteroide (en radianes)

# Condiciones iniciales
s_0 = 0  # Velocidad radial inicial

def ds_dt(s, omega, r):
    # Derivada de s
    omega_prime = -2 * s * omega / r
    return r * omega_prime - G * M / r ** 2

def domega_dt(s, omega, r):
    # Derivada de omega
    return -2 * s * omega / r

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
    t_end = 1e7  # Tiempo total de simulación en segundos

    # Condición inicial omega
    omega_0 = v0 / R_SAT

    # Llamamos a la función Adams-Bashforth para simular la trayectoria
    s_vals, omega_vals, r_vals = adams_bashforth_3(s_0, omega_0, R_SAT, delta_t, t_end)

    # Calculamos la diferencia radial
    final_s = r_vals[-1]

    return final_s - pos_radial_asteroide

# Método de Bisección para encontrar la velocidad inicial óptima
def bisection(f, a, b, tol=1e-6, max_iter=1000):
    if f(a) * f(b) >= 0:
        raise ValueError("La función debe cambiar de signo en el intervalo [a, b].")

    iter_count = 0
    while (b - a) / 2 > tol and iter_count < max_iter:
        c = (a + b) / 2
        if f(c) == 0:
            return c  # Encontramos la raíz exacta
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
        iter_count += 1

    return (a + b) / 2

# Intervalo inicial para Bisección
a = 500  # Límite inferior de la velocidad inicial tentativa
b = 500000  # Límite superior de la velocidad inicial tentativa

# Ejecutamos el método de Bisección
try:
    v0_optimo = bisection(f, a, b)
    print(f"Velocidad inicial óptima para impactar el asteroide: {v0_optimo:.2f} m/s")
except ValueError as e:
    print(e)
