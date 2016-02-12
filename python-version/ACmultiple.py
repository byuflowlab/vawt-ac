import numpy as np
from math import pi, sin, cos, sqrt, atan2

ntheta = 6
dtheta = 2*pi/ntheta
theta = np.arange(dtheta/2, 2*pi, dtheta)

xJ = -np.sin(theta)
yJ = np.cos(theta)
xJc = 0.0
yJc = 0.0

xI = 100 - np.sin(theta)
yI = 100 + np.cos(theta)

Ay = np.zeros((ntheta, ntheta))
Ay2 = np.zeros((ntheta, ntheta))

for i in range(ntheta):
    for j in range(ntheta):

        xstar = xI[i] - xJc
        ystar = yI[i] - yJc

        xstar2 = xI[i] - xJ[j]
        ystar2 = yI[i] - yJ[j]

        phi_j = np.linspace(theta[j] - dtheta/2, theta[j] + dtheta/2, 100)

        v1 = xstar + np.sin(phi_j)
        v2 = ystar - np.cos(phi_j)

        v1_2 = xstar2 - sin(theta[j]) + np.sin(phi_j)
        v2_2 = ystar2 + cos(theta[j]) - np.cos(phi_j)

        value_y = (v1*np.cos(phi_j) + v2*np.sin(phi_j))/(v1**2 + v2**2)
        value_y2 = (v1_2*np.cos(phi_j) + v2_2*np.sin(phi_j))/(v1_2**2 + v2_2**2)

        # value_y2 = (xstar2*np.cos(phi_j) + ystar2*np.sin(phi_j))/(xstar2**2 + ystar2**2)

        Ay[i, j] = np.trapz(value_y, phi_j)/(2*pi)
        Ay2[i, j] = np.trapz(value_y2, phi_j)/(2*pi)


print np.amax(np.abs(Ay - Ay2))
print Ay
# print Ay2

# def generate_ac_matrix(theta):

#     # dtheta = 2*pi/ntheta
#     # theta = np.arange(dtheta/2, 2*pi, dtheta)

#     x = -np.sin(theta)
#     y = np.cos(theta)

#     Ry = np.zeros((ntheta, ntheta))

#     for i in range(ntheta):
#         for j in range(ntheta):

#             phi_j = np.linspace(theta[j] - dtheta/2, theta[j] + dtheta/2, 100)

#             v1 = x[i] + np.sin(phi_j)
#             v2 = y[i] - np.cos(phi_j)

#             value_y = (v1*np.cos(phi_j) + v2*np.sin(phi_j))/(v1**2 + v2**2)

#             Ry[i, j] = np.trapz(value_y, phi_j)


#     Rx = dtheta/2*np.ones((ntheta, ntheta))

#     for i in range(ntheta):
#         if i < ntheta/2:
#             Rx[i, i] = pi*(-1 + 1.0/ntheta)
#         else:
#             Rx[i, i] = pi*(1 + 1.0/ntheta)


#     Wx = np.zeros((ntheta, ntheta))
#     for i in range(ntheta/2, ntheta):
#         Wx[i, ntheta-(i+1)] = -1

#     Ax = Rx/(2*pi) + Wx
#     Ay = Ry/(2*pi)



# def velocities(theta, Qr, x, y, Uinf, R):
#     """x, y relative to center of this VAWT"""

#     x /= R
#     y /= R
#     # Qr is already normalized

#     v1 = x + np.sin(theta)
#     v2 = y - np.cos(theta)

#     integrand = (v1*np.sin(theta) - v2*np.cos(theta))/(v1**2 + v2**2)
#     u = 1.0/(2*pi)*np.trapz(Qr*integrand, theta)

#     # check if in wake
#     if y > -1 and y < 1 and x > 1:
#         u += -np.interp(np.arccos(y), theta, Qr) + np.interp(2*pi-np.arccos(y), theta, Qr)

#     integrand = (v1*np.cos(theta) + v2*np.sin(theta))/(v1**2 + v2**2)
#     v = 1.0/(2*pi)*np.trapz(Qr*integrand, theta)

#     # unnormlize velocities
#     u *= Uinf
#     v *= Uinf

#     return u, v


# def loading(Uinf, Omega, u, v, r, chord, twist, af, B, theta, rho, mu):

#     Vt = Omega*r + Uinf*(1 + u)*cos(theta) + Uinf*v*sin(theta)
#     Vn = Uinf*(1 + u)*sin(theta) - Uinf*v*cos(theta)
#     W = sqrt(Vt**2 + Vn**2)
#     phi = atan2(Vn, Vt)
#     alpha = phi - twist
#     Re = rho*chord*Uinf/mu

#     cl, cd = af.evaluate(alpha, Re)

#     # rotate force coefficients
#     cn = cl*cos(phi) + cd*sin(phi)
#     ct = cl*sin(phi) - cd*cos(phi)

#     # pressure
#     sigma = B*chord/(2*r)
#     Qr = sigma/(2*pi)*cn*(W/Uinf)**2

#     return Qr


# if __name__ == '__main__':

#     from ccblade import CCAirfoil

#     D = 5.0
#     R = D/2.0
#     tsr = 4.0
#     Omega = 127.0*pi/30.0  # rad/s
#     Uinf = Omega*R/tsr

#     u = 0.0
#     v = 0.0
#     r = R
#     chord = 0.1524
#     twist = 0.0
#     af = CCAirfoil.initFromAerodynFile('airfoils/NACA_0015_extrap.dat')
#     B = 3
#     ntheta = 20
#     theta = np.linspace(0, 2*pi, 20)
#     Qr = np.zeros(ntheta)
#     rho = 1.225
#     mu = 1.7894e-5

#     for i in range(ntheta):
#         Qr[i] = loading(Uinf, Omega, u, v, r, chord, twist, af, B, theta[i], rho, mu)

#     print Qr

#     x = 1.0*R
#     y = 1.0*R

#     u, v = velocities(theta, Qr, x, y, Uinf, r)

#     print u, v
