import matplotlib.pyplot as plt

# Constantes y condiciones iniciales
G = 6.67430e-11  # Constante gravitacional en m^3 kg^-1 s^-2
M = 5.972e24  # Masa de la Tierra (en kg)
R_SAT = 2820000000  # Posición radial inicial (en metros)
pos_radial_asteroide = 4122500000  # Posición radial del asteroide (en metros)
pos_angular_asteroide = 3.025  # Posición angular del asteroide (en radianes)
s_0 = 0  # Velocidad radial inicial


# Modificar ds_dt para aceptar tres argumentos
def ds_dt(s, w, r):
    # Aceleración centrífuga menos la aceleración gravitacional
    return r * w ** 2 - G * M / r ** 2


def dw_dt(s, w, r):
    # Ecuación de movimiento para la velocidad angular (derivada de la segunda ley de Newton)
    return -2 * s * w / r


def runge_kutta_2(s, w, r, delta_t):
    # Método de Runge-Kutta de segundo orden para resolver las ecuaciones diferenciales.
    # Este método se utiliza para calcular la posición y velocidad en el siguiente paso de tiempo.
    k1_s = delta_t * ds_dt(s, w, r)
    k1_w = delta_t * dw_dt(s, w, r)
    k2_s = delta_t * ds_dt(s + k1_s / 2, w + k1_w / 2, r)
    k2_w = delta_t * dw_dt(s + k1_s / 2, w + k1_w / 2, r)

    # Calculamos los siguientes valores de s y w usando un promedio ponderado de los incrementos
    s_next = s + (k1_s + k2_s) / 2
    w_next = w + (k1_w + k2_w) / 2
    return s_next, w_next


def adams_bashforth_3(s_0, w_0, r_0, delta_t, t_end):
    # Método de Adams-Bashforth de tercer orden para resolver las ecuaciones diferenciales.
    # Este método usa los tres primeros pasos de Runge-Kutta y luego aplica el método de Adams-Bashforth
    # para los pasos posteriores.
    s_vals = [s_0]
    w_vals = [w_0]
    r_vals = [r_0]

    # Inicializamos con tres pasos de Runge-Kutta para obtener los primeros tres valores
    for _ in range(3):
        s_next, w_next = runge_kutta_2(s_vals[-1], w_vals[-1], r_vals[-1], delta_t)
        s_vals.append(s_next)
        w_vals.append(w_next)
        r_vals.append(r_vals[-1] + s_next * delta_t)

    t = 3
    # Usamos el método de Adams-Bashforth para los siguientes pasos
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
        s_vals.append(s_next)
        w_vals.append(w_next)
        r_vals.append(r_vals[-1] + s_next * delta_t)
        t += 1

    return s_vals, w_vals, r_vals


def f(v0):
    # Función objetivo para encontrar la velocidad inicial que permita impactar el asteroide.
    # La función evalúa la diferencia entre la posición radial final y la posición radial del asteroide.
    delta_t = 1  # Paso de tiempo en segundos
    t_end = 1e6  # Tiempo total de simulación en segundos
    w_0 = v0 / R_SAT  # Calculamos la velocidad angular inicial
    s_vals, w_vals, r_vals = adams_bashforth_3(s_0, w_0, R_SAT, delta_t, t_end)
    final_s = r_vals[-1]  # Posición radial final
    return final_s - pos_radial_asteroide  # Comparamos con la posición del asteroide


def df(v0):
    # Derivada numérica de la función f(v0) usando un valor pequeño de epsilon.
    # Esto es útil para el método de Newton-Raphson.
    epsilon = 1e-6  # Diferencia pequeña para aproximación numérica
    return (f(v0 + epsilon) - f(v0)) / epsilon  # Derivada numérica de f(v0)


def newton_raphson(f, df, v0_init, tol=1e-6, max_iter=1000):
    # Método de Newton-Raphson para encontrar la raíz de una función.
    # Busca la velocidad inicial que permita impactar el asteroide utilizando una aproximación iterativa.
    v0 = v0_init
    for _ in range(max_iter):
        v0_next = v0 - f(v0) / df(v0)  # Fórmula de Newton-Raphson
        if abs(v0_next - v0) < tol:  # Si la diferencia es menor que la tolerancia, hemos encontrado la solución
            return v0_next
        v0 = v0_next
    raise ValueError("No se encontró una solución dentro del número máximo de iteraciones")


# Evaluación del error con diferentes pasos de tiempo
def evaluar_error(v0_optimo):
    # Función para evaluar el error en la posición radial final con diferentes pasos de tiempo.
    # Esto nos permite ver cómo el paso de tiempo afecta la precisión de la simulación.
    pasos = [1, 5, 10, 50, 100]  # Diferentes pasos de tiempo a evaluar
    errores = []
    for delta_t in pasos:
        w_0 = v0_optimo / R_SAT  # Calculamos la velocidad angular inicial
        s_vals, w_vals, r_vals = adams_bashforth_3(s_0, w_0, R_SAT, delta_t, t_end=1e6)
        error = abs(r_vals[-1] - pos_radial_asteroide)  # Error en la posición radial final
        errores.append(error)

    # Gráfica del error en función del paso de tiempo
    plt.figure(figsize=(10, 6))
    plt.plot(pasos, errores, marker='o')
    plt.xlabel("Paso de tiempo (s)")
    plt.ylabel("Error radial (m)")
    plt.title("Error en función del paso de tiempo")
    plt.grid(True)
    plt.show()


# Ejecutamos el método de Newton-Raphson
v0_init = 1000  # Valor inicial para la velocidad
try:
    v0_optimo = newton_raphson(f, df, v0_init)  # Encontramos la velocidad óptima usando el método de Newton-Raphson
    print(f"Velocidad inicial óptima para impactar el asteroide: {v0_optimo:.2f} m/s")
    evaluar_error(v0_optimo)  # Evaluamos el error para diferentes pasos de tiempo
except ValueError as e:
    print(e)
