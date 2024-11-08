pos_radial_proyectil = 28.2e10
pos_radial_asteroide = 4.1225e10
pos_angular_asteroide = 3.025


def aceleracion_asteroide(r):
  G = 6.667430e-11
  masa_asteroide = M
  return -G*M/r**2

def adam_bashforth_3(y,f,delta_t):
  return y[-1] + (delta_t / 12) * (23 * f[-1] - 16 * f[-2] + 5 * f[-3])

def f(v0):
  r = pos_radial_proyectil
  velocidad_ang = v0/r

  valores_r = [r]
  r_deriv = [v0]
  valores_theta = [theta]
  theta_deriv = [velocidad_ang]

