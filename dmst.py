#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Andrew Ning on 2013-02-07.
Copyright (c) NREL. All rights reserved.
"""

from math import cos, pi, fabs
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import UnivariateSpline, RectBivariateSpline

import _vbemroutines


class VAWT:


    def __init__(self, z, r, dtheta, chord, twist, af, B=2, rho=1.225, mu=1.7894e-5, nt=60):

        self.z = z
        self.r = r
        self.dtheta = dtheta
        self.chord = chord
        self.twist = twist
        self.af = af

        self.B = B
        self.rho = rho
        self.mu = mu
        self.nt = nt



    def __checkCrossing(self, args, afixed, arange):

        # initialized
        ffixed = self.bem_residual(afixed, *args)
        found = False
        afound = None

        for a in arange:
            fother = self.bem_residual(a, *args)

            # check for change in sign
            if ffixed*fother < 0:
                afound = a
                found = True
                break

        return found, afound



    def forces(self, Uinf, Omega):

        nz = len(self.z)

        # fit spline (z-r)
        spline = UnivariateSpline(self.z, self.r)

        # compute swept area
        S = 2.0*spline.integral(self.z[0], self.z[-1])

        # compute slope of blade (z-r plane)
        m = np.zeros(nz)
        for i in range(nz):
            m[i] = spline.derivatives(self.z[i])[1]
        delta = np.arctan(-m)


        # # fit blade with spline for sweep (theta-z)
        # if (np.max(self.dtheta) - np.min(self.dtheta)) < 1e-6:
        #     Lambda = np.zeros(nz)
        # else:
        #     spline = UnivariateSpline(self.dtheta, self.z)

        #     m = np.zeros(nz)
        #     for i in range(nz):
        #         m[i] = spline.derivatives(self.dtheta[i])[1]
        #     Lambda = np.arctan(self.r/m)


        # dynamic pressure
        q = 0.5*self.rho*Uinf**2

        # discretize in theta
        nt = self.nt
        thetaV = np.linspace(-pi/2, 3*pi/2, nt)

        # compute forces on a grid in height and theta
        Np = np.zeros((nz, nt))
        Tp = np.zeros((nz, nt))
        Zp = np.zeros((nz, nt))

        R = np.max(self.r)


        # loop across height
        for i in range(nz):

            # save induction factors from first disk
            afirst = np.zeros(nt)

            # loop across theta
            for j in range(nt):

                # incoming velocity
                if thetaV[j] <= pi/2:  # first actuator disk
                    V = Uinf

                else:  # second actuator disk
                    # TODO: may want a better interpolation method here
                    a = np.interp(pi-thetaV[j], thetaV, afirst)
                    V = Uinf*(1 - 2*a)

                # save arguments
                args = (self.r[i], self.chord[i], self.twist[i], thetaV[j], delta[i], self.af[i], V, Omega)


                # special cases
                if self.r[i]/R < 1e-3:
                    Np[i, j] = 0.0
                    Tp[i, j] = 0.0
                    Zp[i, j] = 0.0
                    afirst[j] = 0.0
                    continue
                elif fabs(cos(thetaV[j])) < 1e-6:
                    a = 0.0  # (TODO: check that this is right)
                    residual, Np[i, j], Tp[i, j], Zp[i, j] = self.bem(a, *args)
                    afirst[j] = a
                    continue


                # find lower and upper limits
                amin = 0.0
                found, amax = self.__checkCrossing(args, amin, np.arange(0.09, 1.0, 0.1))
                if not found:
                    amax = 0.0
                    found, amin = self.__checkCrossing(args, amax, np.arange(-0.09, -1.0, -0.1))
                if not found:
                    print 'solution not found'
                    Np[i, j] = 0.0
                    Tp[i, j] = 0.0
                    afirst[j] = 0.0
                else:
                    a = brentq(self.bem_residual, amin, amax, args=args)
                    residual, Np[i, j], Tp[i, j], Zp[i, j] = self.bem(a, *args)
                    afirst[j] = a


        # integrate across one blade (forces relative to blade center)
        N_theta_oneblade = np.zeros(nt)
        T_theta_oneblade = np.zeros(nt)
        Z_theta_oneblade = np.zeros(nt)
        Q_theta_oneblade = np.zeros(nt)

        # interpolation
        fN = RectBivariateSpline(self.z, thetaV, Np, kx=1, ky=1)
        fT = RectBivariateSpline(self.z, thetaV, Tp, kx=1, ky=1)
        fZ = RectBivariateSpline(self.z, thetaV, Zp, kx=1, ky=1)

        for j in range(nt):
            theta_intersect = thetaV[j] + self.dtheta
            theta_intersect[theta_intersect < -pi/2] += 2*pi
            theta_intersect[theta_intersect > 3*pi/2] -= 2*pi
            Np_along_blade = fN.ev(self.z, theta_intersect)
            Tp_along_blade = fT.ev(self.z, theta_intersect)
            Zp_along_blade = fZ.ev(self.z, theta_intersect)

            # integrate to get total force from blade
            Np_rotated = Np_along_blade*np.cos(self.dtheta) - Tp_along_blade*np.sin(self.dtheta)
            Tp_rotated = Np_along_blade*np.sin(self.dtheta) + Tp_along_blade*np.cos(self.dtheta)
            N_theta_oneblade[j] = np.trapz(Np_rotated, self.z)
            T_theta_oneblade[j] = np.trapz(Tp_rotated, self.z)
            Z_theta_oneblade[j] = np.trapz(Zp_along_blade, self.z)
            Q_theta_oneblade[j] = np.trapz(self.r*Tp_along_blade, self.z)


        # add other blades
        X_theta = np.zeros(nt)
        Y_theta = np.zeros(nt)
        Z_theta = np.zeros(nt)

        for j in range(self.B):
            delta_theta = 2*pi*j/float(self.B)
            other_theta = thetaV + delta_theta
            other_theta[other_theta > 3*pi/2] -= 2*pi
            # TODO: maybe use other interpolation scheme
            N_other_blade = np.interp(other_theta, thetaV, N_theta_oneblade)
            T_other_blade = np.interp(other_theta, thetaV, T_theta_oneblade)
            Z_other_blade = np.interp(other_theta, thetaV, Z_theta_oneblade)

            X_theta += -N_other_blade*np.cos(other_theta) + T_other_blade*np.sin(other_theta)
            Y_theta += N_other_blade*np.sin(other_theta) + T_other_blade*np.cos(other_theta)
            Z_theta += Z_other_blade



        # average across theta for torque and power
        Qbar = self.B/(2*pi) * np.trapz(Q_theta_oneblade, thetaV)
        P = Omega*Qbar

        # normalize
        CN_theta_oneblade = N_theta_oneblade / (q*S)
        CT_theta_oneblade = T_theta_oneblade / (q*S)
        CX_theta = X_theta / (q*S)
        CY_theta = Y_theta / (q*S)
        CZ_theta = Z_theta / (q*S)
        CP = P / (q * Uinf * S)

        return thetaV, CN_theta_oneblade, CT_theta_oneblade, CX_theta, CY_theta, CZ_theta, CP


    def bem(self, a, r, chord, twist, theta, delta, af, Uinf, Omega):

        if (r < 1e-6):
            return 0.0, 0.0, 0.0

        alpha, Re, phi, W = _vbemroutines.velocity(a, r, chord, twist, theta, delta,
                                                   Uinf, Omega, self.rho, self.mu)

        # airfoil lift and drag
        cl, cd = af.evaluate(alpha, Re)

        residual, Np, Tp, Zp = _vbemroutines.vbem(a, r, chord, theta, delta,
                                                  cl, cd, phi, W, Uinf, self.rho, self.B)

        return residual, Np, Tp, Zp



    def bem_residual(self, a, r, chord, twist, theta, delta, af, Uinf, Omega):

        residual, Np, Tp, Zp = self.bem(a, r, chord, twist, theta, delta, af, Uinf, Omega)

        return residual


def run():

    D = 5.0
    # H = 1.02*D
    H = 1.5*D

    hh = H/2
    R = D/2
    B = 3

    r0 = 0.5*R

    z = np.linspace(-hh, hh, 10)
    r = r0 + (R-r0) * (1 - (z/hh)**2)
    # chord = 0.1524 * np.ones_like(r)
    chord = 0.2 * np.ones_like(r)
    twist = 0.0*np.ones_like(r)
    # sigma = 0.15
    # chord = sigma * R / B * np.ones_like(r)
    dtheta = .3*np.sign(z)*np.sqrt(np.abs(z))
    # dtheta = np.zeros_like(z)

    from twister.rotoraero.ccblade import CCAirfoil

    af = [CCAirfoil.initFromAerodynFile('NACA_0012.dat')]*len(r)


    v = VAWT(z, r, dtheta, chord, twist, af, B)








    Uinf = 28.36  # 5.0
    tsr = 5.0  # 3.0
    Omega = tsr*Uinf/R

    theta, CN, CT, CX, CY, CZ, CP = v.forces(Uinf, Omega)

    print CP
    import matplotlib.pyplot as plt
    plt.plot(theta*180/pi, CN, label='$C_N$', color='r')
    plt.plot(theta*180/pi, CT, label='$C_T$', color='k')
    plt.xlabel('$\\theta$')
    plt.legend(loc='lower right')

    # plt.savefig('/Users/sning/Dropbox/NREL/VAWT/memo/loads.pdf')
    # plt.show()


    plt.figure()
    plt.plot(theta*180/pi, CX, label='$C_X$')
    plt.plot(theta*180/pi, CY, label='$C_Y$')
    plt.plot(theta*180/pi, CZ, label='$C_Z$')
    plt.legend()
    plt.show()

    n = 20
    tsr_vec = np.linspace(3, 10, n)
    cp_vec = np.zeros(n)
    for i in range(n):
        Omega = tsr_vec[i]*Uinf/R
        theta, CN, CT, CX, CY, CZ, cp_vec[i] = v.forces(Uinf, Omega)


    plt.plot(tsr_vec, cp_vec, 'k')
    plt.xlabel('$\\lambda$')
    plt.ylabel('$C_p$')
    # plt.savefig('/Users/sning/Dropbox/NREL/VAWT/memo/cp.pdf')
    plt.show()


    # plot geometry


    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.view_init(10, 315)
    # ax = fig.gca(projection='3d')
    # theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
    # z = np.linspace(-2, 2, 100)
    # r = z**2 + 1
    x = -r * np.cos(dtheta)
    y = r * np.sin(dtheta)
    ax.plot(x, y, z, 'r')
    x = -r * np.cos(dtheta+2.0*pi/3)
    y = r * np.sin(dtheta+2.0*pi/3)
    ax.plot(x, y, z, 'r')
    x = -r * np.cos(dtheta+4*pi/3)
    y = r * np.sin(dtheta+4*pi/3)
    ax.plot(x, y, z, 'r')
    theta_circle = np.linspace(0, 2*pi, 50)
    x = -r[0]*np.cos(theta_circle)
    y = r[0]*np.sin(theta_circle)
    ax.plot(x, y, z[0]*np.ones_like(x), 'k--')
    x = -r[-1]*np.cos(theta_circle)
    y = r[-1]*np.sin(theta_circle)
    ax.plot(x, y, z[-1]*np.ones_like(x), 'k--')
    x = -r[len(r)/2]*np.cos(theta_circle)
    y = r[len(r)/2]*np.sin(theta_circle)
    ax.plot(x, y, z[len(r)/2]*np.ones_like(x), 'k--')

    # ax.set_xlabel('x (m)')
    # ax.set_ylabel('y (m)')
    # ax.set_zlabel('z (m)')
    ax.grid(False)

    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(x.max()+x.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(y.max()+y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(z.max()+z.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w')

    # plt.axis('equal')
    # plt.axes().set_aspect('equal')
    # ax.xlim([-2.5, 2.5])
    # ax.ylim([-2.5, 2.5])
    # ax.zlim([-2.5, 2.5])
    # ax.legend()

    # plt.savefig('/Users/sning/Dropbox/NREL/VAWT/memo/geometry.pdf')
    plt.show()




if __name__ == '__main__':
    run()

    # import cProfile
    # cProfile.run('run()', 'deleteme')
    # import pstats
    # p = pstats.Stats('deleteme')
    # p.sort_stats('cumulative').print_stats(10)

