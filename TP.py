from scipy.optimize import bisect

pos_radial_proyectil = 28.2e10
pos_radial_asteroide = 4.1225e10
pos_angular_asteroide = 3.025


def aceleracion_asteroide(r):
  G = 6.667430e-11
  M = 1000
  return -G*M/r**2

def adam_bashforth_3(y_valores,f_valores,delta_t):
  return y_valores[-1] + (delta_t / 12) * (23 * f_valores[-1] - 16 * f_valores[-2] + 5 * f_valores[-3])

def f(v0):
  r = pos_radial_proyectil
  velocidad_angular = v0/r
  theta = 0

  r_valores = [r]
  r_deriv = [v0]
  theta_valores = [theta]
  theta_deriv = [velocidad_angular]

  for i in range(2):
    a_r = aceleracion_asteroide(r)
    r += r_deriv[-1] * delta_t
    theta += theta_deriv[-1] * delta_t
    r_deriv.append(a_r)
    velocidad_angular = v0 / r

    r_valores.append(r)
    theta_valores.append(theta)
    theta_deriv.append(velocidad_angular)

  while theta < pos_angular_asteroide:
    a_r = aceleracion_asteroide(r_valores[-1])
    r_next = adams_bashforth_3(r_valores, r_deriv, delta_t)
    theta_next = adams_bashforth_3(theta_vals, theta_deriv, delta_t)

    r_valores.append(r_next)
    r_deriv.append(a_r)

    theta_valores.append(theta_next)
    theta_deriv.append(v0 / r_next)

    r = r_next
    theta = theta_next

  return r - pos_radial_asteroide

  v0_min = 500
  v0_max = 1e5
  v0_optimo = bisect(f, v0_min, v0_max, xtol=1e-5)
  print(f"La velocidad inicial optima para impactar el asteroide es: {v0_optimo}m/s")

