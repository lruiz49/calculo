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
  velocidad_angular = v0/r
  theta = 0

  r_valores = [r]
  r_deriv = [v0]
  theta_valores = [theta]
  theta_deriv = [velocidad_ang]

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
    
    
