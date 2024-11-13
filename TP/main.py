import numpy as np
import matplotlib.pyplot as plt

# Constantes
G = 6.67430e-11  # Constante gravitacional (m^3 kg^-1 s^-2)
M = 5.972e24  # Masa de la Tierra (kg)
R_SAT = 2.82e9  # Posición radial inicial del satélite (metros)
pos_radial_asteroide = 4.1225e9  # Posición radial del asteroide (metros)
pos_angular_asteroide = 0  # Posición angular del asteroide (alineado con el eje x a 0 radianes)
initial_angular_position_satellite = 3.025  # Posición angular inicial del satélite (en radianes)
s_0 = 0  # Velocidad radial inicial
v0_optimo = 409.85  # Valor inicial para la velocidad óptima

# Ecuaciones diferenciales
def ds_dt(s, w, r):
    return r * w ** 2 - G * M / r ** 2

def dw_dt(s, w, r):
    return -2 * s * w / r

# Método de Adams-Bashforth
def adams_bashforth_3(s_0, w_0, r_0, theta_0, delta_t, t_end):
    s_vals = [s_0]
    w_vals = [w_0]
    r_vals = [r_0]
    theta_vals = [theta_0]

    for t in range(3):
        k1_s = delta_t * ds_dt(s_vals[-1], w_vals[-1], r_vals[-1])
        k1_w = delta_t * dw_dt(s_vals[-1], w_vals[-1], r_vals[-1])
        k2_s = delta_t * ds_dt(s_vals[-1] + k1_s / 2, w_vals[-1] + k1_w / 2, r_vals[-1])
        k2_w = delta_t * dw_dt(s_vals[-1] + k1_s / 2, w_vals[-1] + k1_w / 2, r_vals[-1])
        s_next = s_vals[-1] + (k1_s + k2_s) / 2
        w_next = w_vals[-1] + (k1_w + k2_w) / 2
        r_next = r_vals[-1] + s_next * delta_t
        theta_next = theta_vals[-1] + w_next * delta_t

        s_vals.append(s_next)
        w_vals.append(w_next)
        r_vals.append(r_next)
        theta_vals.append(theta_next)

    t = 3
    while t * delta_t <= t_end:
        s_next = (s_vals[t] + delta_t / 12 * (
                23 * ds_dt(s_vals[t], w_vals[t], r_vals[t]) -
                16 * ds_dt(s_vals[t - 1], w_vals[t - 1], r_vals[t - 1]) +
                5 * ds_dt(s_vals[t - 2], w_vals[t - 2], r_vals[t - 2])
        ))
        w_next = (w_vals[t] + delta_t / 12 * (
                23 * dw_dt(s_vals[t], w_vals[t], r_vals[t]) -
                16 * dw_dt(s_vals[t - 1], w_vals[t - 1], r_vals[t - 1]) +
                5 * dw_dt(s_vals[t - 2], w_vals[t - 2], r_vals[t - 2])
        ))
        r_next = r_vals[-1] + s_next * delta_t
        theta_next = theta_vals[-1] + w_next * delta_t

        s_vals.append(s_next)
        w_vals.append(w_next)
        r_vals.append(r_next)
        theta_vals.append(theta_next)
        t += 1

    return r_vals, theta_vals

def f(v0):
    # Función objetivo para encontrar la velocidad inicial que permita impactar el asteroide
    delta_t = 50  # Paso temporal en segundos
    t_end = 3.5e7 # Tiempo total de la simulación en segundos
    w_0 = v0 / R_SAT  # Calculamos la velocidad angular inicial
    r_vals, w_vals = adams_bashforth_3(s_0, w_0, R_SAT, initial_angular_position_satellite, delta_t, t_end)
    final_s = r_vals[-1]  # Posición radial final
    return final_s - pos_radial_asteroide  # Comparamos con la posición del asteroide

def df(v0):
    # Derivada numérica de la función f(v0) usando un valor pequeño de epsilon
    epsilon = 1e-6  # Diferencia pequeña para la aproximación numérica
    return (f(v0 + epsilon) - f(v0)) / epsilon  # Derivada numérica de f(v0)

def newton_raphson(f, df, v0_init, tol=1e-6, max_iter=1000):
    # Método de Newton-Raphson para encontrar la raíz de una función
    v0 = v0_init
    for _ in range(max_iter):
        v0_next = v0 - f(v0) / df(v0)  # Fórmula de Newton-Raphson
        if abs(v0_next - v0) < tol:  # Si la diferencia es menor que la tolerancia, hemos encontrado la solución
            return v0_next
        v0 = v0_next
    raise ValueError("No se encontró una solución dentro del número máximo de iteraciones")

# Encontrar la velocidad inicial óptima
v0_init = 1000  # Valor inicial para la velocidad
try:
    v0_optimo = newton_raphson(f, df, v0_init)  # Encontramos la velocidad óptima usando el método de Newton-Raphson
    print(f"Velocidad inicial óptima para impactar el asteroide: {v0_optimo:.2f} m/s")
except ValueError as e:
    print(e)

# Simular la trayectoria con la velocidad óptima
delta_t = 50  # Paso de tiempo en segundos
t_end = 3.4e7  # Tiempo final de la simulación en segundos
w_0 = v0_optimo / R_SAT  # Velocidad angular inicial
r_vals, theta_vals = adams_bashforth_3(s_0, w_0, R_SAT, initial_angular_position_satellite, delta_t, t_end)

# Convertir coordenadas polares a cartesianas para la gráfica
x_vals = [r * np.cos(theta) for r, theta in zip(r_vals, theta_vals)]
y_vals = [r * np.sin(theta) for r, theta in zip(r_vals, theta_vals)]

# Posición del asteroide en coordenadas cartesianas
asteroid_x = pos_radial_asteroide * np.cos(pos_angular_asteroide)
asteroid_y = pos_radial_asteroide * np.sin(pos_angular_asteroide)

# Graficar la trayectoria
plt.figure(figsize=(8, 8))
plt.plot(0, 0, 'bo', markersize=10, label="Tierra")  # Tierra en el centro
plt.plot(x_vals, y_vals, label="Trayectoria del satélite")
plt.plot(x_vals[0], y_vals[0], 'go', markersize=5, label="Posición inicial del satélite")
plt.plot(asteroid_x, asteroid_y, 'ro', markersize=5, label="Posición del asteroide")
plt.xlabel("Posición X (m)")
plt.ylabel("Posición Y (m)")
plt.legend()
plt.title("Trayectoria del satélite desde la Tierra al asteroide")
plt.grid()
plt.axis('equal')
plt.show()

# Análisis de error en función del paso temporal delta_t
delta_t_vals = [1, 10, 50, 100, 200, 500, 1000]  # Valores específicos de delta_t
errors = []

for delta_t in delta_t_vals:
    w_0 = v0_optimo / R_SAT  # Velocidad angular inicial para cada delta_t
    r_vals, _ = adams_bashforth_3(s_0, w_0, R_SAT, initial_angular_position_satellite, delta_t, t_end)
    final_s = r_vals[-1]
    error = abs(final_s - pos_radial_asteroide)  # El error es la diferencia en la posición final
    errors.append(error)

# Graficar el error en función de delta_t
plt.figure(figsize=(8, 6))
plt.plot(delta_t_vals, errors, marker='o', color='b')
plt.xlabel("Paso de tiempo (delta_t) [segundos]")
plt.ylabel("Error en la posición radial final [metros]")
plt.title("Error vs Paso de tiempo en la trayectoria del satélite")
plt.grid(True)
plt.show()
