pos_radial_proyectil = 28.2e10
pos_radial_asteroide = 4.1225e10
pos_angular_asteroide = 3.025


def aceleracion_asteroide(r):
  G = 6.667430e-11
  M = masa_asteroide
  return -G*M/r**2

def adam_bashforth_3(y_valores,f_valores,delta_t):
  return y_valores[-1] + (delta_t / 12) * (23 * f_valores[-1] - 16 * f_valores[-2] + 5 * f_valores[-3])

def f(v0):
  r = pos_radial_proyectil
  velocidad_ang = v0/r
  theta = 0

  r_valores = [r]
  r_deriv = [v0]
  theta_valores = [theta]
  theta_deriv = [velocidad_ang]

  for i in range(2)
    a_r = aceleracion_asteroide(r_vals[-1])
    v_r = r_deriv[-1] + a_r * delta_t
    r_next = r_vals[-1] + v_r * delta_t

    r_valores.append(r_next)
    r_deriv.append(v_r)
    theta_next = theta_valores[-1] + theta_deriv[-1] * delta_t
    theta_valores.append(theta_next)
    theta_deriv.append(v0 / r_next)
