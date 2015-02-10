#!/usr/bin/env python
# encoding: utf-8
"""
vawt.py

Created by Andrew Ning on 2013-02-04.
Copyright (c) NREL. All rights reserved.
"""

import cPickle as pickle
from math import cos, pi, tan, fabs
from scipy.optimize import brentq, root
from scipy.interpolate import UnivariateSpline, RectBivariateSpline
import numpy as np

from _ac import actuatorcylinder as _ac
import _vbemroutines


class NACA0015_linear:

    def __init__(self):

        self.alpha_data = np.array([-20, -19, -18, -17, -16, -15, -14, -13,
                                    -12, -11, -10, -9, -8, -7, -6, -5, -4,
                                    -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8,
                                    9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                                    19, 20])*pi/180.0
        self.cl_data = np.array([-2.4235, -2.3069, -2.1896, -2.0716, -1.9531, -1.8339,
                                 -1.7141, -1.5939, -1.4731, -1.352, -1.2304, -1.1084,
                                 -0.9861, -0.8635, -0.7406, -0.6175, -0.4942, -0.3708,
                                 -0.2473, -0.1236, 0, 0.1237, 0.2473, 0.3708, 0.4943,
                                 0.6175, 0.7406, 0.8635, 0.9861, 1.1084, 1.2304, 1.352,
                                 1.4732, 1.5939, 1.7142, 1.8339, 1.9531, 2.0717, 2.1896,
                                 2.3069, 2.4235])

    def evaluate(self, alpha, Re):

        alpha = np.array(alpha)
        cl = np.interp(alpha, self.alpha_data, self.cl_data)
        cd = 0.023

        if alpha.size > 1:
            cd = cd*np.ones_like(alpha)

        return cl, cd





def generate_ac_matrix(ntheta=36):

    dtheta = 2*pi/ntheta
    theta = np.arange(dtheta/2, 2*pi, dtheta)

    x = -np.sin(theta)
    y = np.cos(theta)

    Ry = np.zeros((ntheta, ntheta))

    for i in range(ntheta):
        for j in range(ntheta):

            phi_j = np.linspace(theta[j] - dtheta/2, theta[j] + dtheta/2, 100)

            v1 = x[i] + np.sin(phi_j)
            v2 = y[i] - np.cos(phi_j)

            value_y = (v1*np.cos(phi_j) + v2*np.sin(phi_j))/(v1**2 + v2**2)

            Ry[i, j] = np.trapz(value_y, phi_j)


    Rx = dtheta/2*np.ones((ntheta, ntheta))

    for i in range(ntheta):
        if i < ntheta/2:
            Rx[i, i] = pi*(-1 + 1.0/ntheta)
        else:
            Rx[i, i] = pi*(1 + 1.0/ntheta)


    Wx = np.zeros((ntheta, ntheta))
    for i in range(ntheta/2, ntheta):
        Wx[i, ntheta-(i+1)] = -1

    Ax = Rx/(2*pi) + Wx
    Ay = Ry/(2*pi)

    A = np.vstack((Ax, Ay))

    pickle.dump((theta, A), open('coeff.p', 'wb'))



class VAWT:
    """docstring for VAWT"""

    def __init__(self, z, r, chord, twist, dtheta, af, zbase, zref, shear, B=3, rho=1.225,
                 mu=1.7894e-5, method='ac'):

        self.z = z
        self.r = r
        self.chord = chord
        self.twist = twist
        self.dtheta = dtheta
        self.af = af

        self.zbase = zbase
        self.zref = zref
        self.shear = shear

        self.B = B
        self.rho = rho
        self.mu = mu
        self.method = method


        # grid size
        self.nz = len(self.z)

        # fit spline (z-r)
        spline = UnivariateSpline(self.z, self.r)

        # compute swept area
        self.S = 2.0*spline.integral(self.z[0], self.z[-1])

        # compute slope of blade (z-r plane)
        m = np.zeros(self.nz)
        for i in range(self.nz):
            m[i] = spline.derivatives(self.z[i])[1]
        self.delta = np.arctan(-m)

        # length of blade
        self.L = np.trapz(np.sqrt(1 + m**2), self.z)

        # max radius
        self.R = np.max(self.r)

        if method == 'ac':
            self._ac_init()
        elif method == 'dmst':
            self._dmst_init()
        else:
            print 'method does not exist'
            exit()


    def _ac_init(self):

        # influence coefficients for a specific discretization
        _ac.theta, _ac.am = pickle.load(open('coeff.p', 'rb'))
        self.theta = _ac.theta
        self.nt = len(self.theta)


    def _dmst_init(self):

        # discretize in theta
        self.nt = 36
        # self.theta = np.linspace(0, 2*pi, self.nt)
        dtheta = 2*pi/self.nt
        self.theta = np.arange(dtheta/2.0, 2*pi, dtheta)


    def _ac_eval(self, w, r, chord, twist, af, B, delta, Uinf, Omega, rho, mu):

        W, phi, alpha, Re = _ac.velocity(w, r, chord, twist, Uinf, Omega,
                                         rho, mu)
        cl, cd = af.evaluate(alpha, Re)

        cn, ct, error_array = _ac.error(w, r, chord, Uinf, delta, B, W,
                                        phi, cl, cd, modlin=True)

        return cn, ct, W, error_array


    def _ac_error(self, w, r, chord, twist, af, B, delta, Uinf, Omega, rho, mu):

        cn, ct, W, error_array = self._ac_eval(w, r, chord, twist, af, B, delta, Uinf, Omega, rho, mu)

        return error_array


    def _ac_forces(self, Uinf, Omega):

        nz = self.nz
        nt = self.nt

        # compute forces on a grid in height and theta
        Rp = np.zeros((nz, nt))
        Tp = np.zeros((nz, nt))
        Zp = np.zeros((nz, nt))

        # loop across height
        for i in range(nz):

            # account for shear
            zturbine = self.z - self.z[0]
            V = Uinf*((zturbine[i] + self.zbase)/(self.zref + self.zbase))**self.shear

            args = (self.r[i], self.chord[i], self.twist[i], self.af[i], self.B,
                    self.delta[i], V, Omega, self.rho, self.mu)

            # special cases
            if self.r[i]/self.R < 1e-3:
                Rp[i, :] = np.zeros(nt)
                Tp[i, :] = np.zeros(nt)
                Zp[i, :] = np.zeros(nt)
                continue

            w = np.zeros(2*nt)
            sol = root(self._ac_error, w, args=(args), method='hybr')
            w = sol.x

            cn, ct, W, error = self._ac_eval(w, *args)
            q = 0.5*self.rho*W**2
            Rp[i, :] = -cn*q*self.chord[i]
            Tp[i, :] = ct*q*self.chord[i]/cos(self.delta[i])
            Zp[i, :] = -cn*q*self.chord[i]*tan(self.delta[i])

        return Rp, Tp, Zp




    def __checkCrossing(self, args, afixed, arange):

        # initialized
        ffixed = self._dmst_error(afixed, *args)
        found = False
        afound = None

        for a in arange:
            fother = self._dmst_error(a, *args)

            # check for change in sign
            if ffixed*fother < 0:
                afound = a
                found = True
                break

        return found, afound


    def _dmst_eval(self, a, r, chord, twist, theta, delta, af, Uinf, Omega):

        if (r < 1e-6):
            return 0.0, 0.0, 0.0

        alpha, Re, phi, W = _vbemroutines.velocity(a, r, chord, twist, theta, delta,
                                                   Uinf, Omega, self.rho, self.mu)

        # airfoil lift and drag
        cl, cd = af.evaluate(alpha, Re)

        residual, Np, Tp, Zp = _vbemroutines.vbem(a, r, chord, theta, delta,
                                                  cl, cd, phi, W, Uinf, self.rho, self.B)

        return residual, Np, Tp, Zp



    def _dmst_error(self, a, r, chord, twist, theta, delta, af, Uinf, Omega):

        residual, Np, Tp, Zp = self._dmst_eval(a, r, chord, twist, theta, delta, af, Uinf, Omega)

        return residual


    def _dmst_forces(self, Uinf, Omega):

        nz = self.nz
        nt = self.nt

        # compute forces on a grid in height and theta
        Rp = np.zeros((nz, nt))
        Tp = np.zeros((nz, nt))
        Zp = np.zeros((nz, nt))

        # loop across height
        for i in range(nz):

            # save induction factors from first disk
            afirst = np.zeros(nt)

            # account for shear
            zturbine = self.z - self.z[0]
            Vshear = Uinf*((zturbine[i] + self.zbase)/(self.zref + self.zbase))**self.shear

            # loop across theta
            for j in range(nt):

                # incoming velocity
                if self.theta[j] <= pi:  # first actuator disk
                    V = Vshear

                else:  # second actuator disk
                    # TODO: may want a better interpolation method here
                    a = np.interp(pi-self.theta[j], self.theta, afirst)
                    V = Vshear*(1 - 2*a)

                # save arguments
                args = (self.r[i], self.chord[i], self.twist[i], self.theta[j], self.delta[i], self.af[i], V, Omega)


                # special cases
                if self.r[i]/self.R < 1e-3:
                    Rp[i, j] = 0.0
                    Tp[i, j] = 0.0
                    Zp[i, j] = 0.0
                    afirst[j] = 0.0
                    continue
                elif fabs(cos(self.theta[j])) < 1e-6:
                    a = 0.0  # (TODO: check that this is right)
                    residual, Rp[i, j], Tp[i, j], Zp[i, j] = self._dmst_eval(a, *args)
                    afirst[j] = a
                    continue


                # find lower and upper limits
                amin = 0.0
                found, amax = self.__checkCrossing(args, amin, np.arange(0.09, 1.0, 0.1))
                if not found:
                    amax = 0.0
                    found, amin = self.__checkCrossing(args, amax, np.arange(-0.09, -1.0, -0.1))
                if not found:
                    print 'dmst solution not found at this station'
                    Rp[i, j] = 0.0
                    Tp[i, j] = 0.0
                    afirst[j] = 0.0
                else:
                    a = brentq(self._dmst_error, amin, amax, args=args)
                    residual, Rp[i, j], Tp[i, j], Zp[i, j] = self._dmst_eval(a, *args)
                    afirst[j] = a

        return Rp, Tp, Zp


    def evaluate(self, Uinf, Omega):

        if self.method == 'ac':
            Rp, Tp, Zp = self._ac_forces(Uinf, Omega)
        elif self.method == 'dmst':
            Rp, Tp, Zp = self._dmst_forces(Uinf, Omega)
        else:
            ValueError('method does not exist')



        nt = self.nt

        # integrate across one blade (forces relative to blade center)
        R_theta_oneblade = np.zeros(nt)
        T_theta_oneblade = np.zeros(nt)
        Z_theta_oneblade = np.zeros(nt)
        Q_theta_oneblade = np.zeros(nt)

        # interpolation
        fR = RectBivariateSpline(self.z, self.theta, Rp, kx=1, ky=1)
        fT = RectBivariateSpline(self.z, self.theta, Tp, kx=1, ky=1)
        fZ = RectBivariateSpline(self.z, self.theta, Zp, kx=1, ky=1)

        for j in range(nt):
            theta_intersect = self.theta[j] + self.dtheta
            theta_intersect[theta_intersect < 0] += 2*pi
            theta_intersect[theta_intersect > 2*pi] -= 2*pi
            Rp_along_blade = fR.ev(self.z, theta_intersect)
            Tp_along_blade = fT.ev(self.z, theta_intersect)
            Zp_along_blade = fZ.ev(self.z, theta_intersect)

            # integrate to get total force from blade
            Rp_rotated = Rp_along_blade*np.cos(self.dtheta) - Tp_along_blade*np.sin(self.dtheta)
            Tp_rotated = Rp_along_blade*np.sin(self.dtheta) + Tp_along_blade*np.cos(self.dtheta)
            R_theta_oneblade[j] = np.trapz(Rp_rotated, self.z)
            T_theta_oneblade[j] = np.trapz(Tp_rotated, self.z)
            Z_theta_oneblade[j] = np.trapz(Zp_along_blade, self.z)
            Q_theta_oneblade[j] = np.trapz(self.r*Tp_along_blade, self.z)


        # add other blades
        X_theta = np.zeros(nt)
        Y_theta = np.zeros(nt)
        Z_theta = np.zeros(nt)

        for j in range(self.B):
            delta_theta = 2*pi*j/float(self.B)
            other_theta = self.theta + delta_theta
            other_theta[other_theta > 2*pi] -= 2*pi
            # TODO: maybe use other interpolation scheme
            R_other_blade = np.interp(other_theta, self.theta, R_theta_oneblade)
            T_other_blade = np.interp(other_theta, self.theta, T_theta_oneblade)
            Z_other_blade = np.interp(other_theta, self.theta, Z_theta_oneblade)

            X_theta += -R_other_blade*np.sin(other_theta) - T_other_blade*np.cos(other_theta)
            Y_theta += R_other_blade*np.cos(other_theta) - T_other_blade*np.sin(other_theta)
            Z_theta += Z_other_blade



        # average across theta for torque and power
        Qbar = self.B/(2*pi) * np.trapz(Q_theta_oneblade, self.theta)
        P = Omega*Qbar

        # normalize
        q = 0.5*self.rho*Uinf**2
        CR = R_theta_oneblade / (q*self.S)  # for one blade
        CT = T_theta_oneblade / (q*self.S)
        CX = X_theta / (q*self.S)  # for all blades
        CY = Y_theta / (q*self.S)
        CZ = Z_theta / (q*self.S)
        CP = P / (q * Uinf * self.S)

        return self.theta, CP, CR, CT, CX, CY, CZ, self.S



if __name__ == '__main__':


    # from wisdem.rotor.ccblade import CCAirfoil

    D = 5.0
    H = 1.02*D
    B = 3
    chord0 = 0.1524

    hh = H/2
    R = D/2
    r0 = 0.0*R
    z = np.linspace(-hh, hh, 10)
    r = r0 + (R-r0) * (1 - (z/hh)**2)
    chord = chord0 * np.ones_like(r)
    twist = np.zeros_like(r)
    dtheta = np.zeros_like(z)
    af = [NACA0015_linear()]*len(r)
    # af = [CCAirfoil.initFromAerodynFile('NACA_0015_extrap.dat')]*len(r)

    rho = 0.82*1.225
    zbase = H/2
    zref = H
    shear = 0.2




    vawt_ac = VAWT(z, r, chord, twist, dtheta, af, zbase, zref, shear, B, rho, method='ac')
    vawt_dmst = VAWT(z, r, chord, twist, dtheta, af, zbase, zref, shear, B, rho, method='dmst')


    tsr = np.linspace(1, 8, 20)
    ntsr = len(tsr)
    CT = np.zeros(ntsr)
    CP = np.zeros(ntsr)
    Kp = np.zeros(ntsr)
    CTdmst = np.zeros(ntsr)
    CPdmst = np.zeros(ntsr)
    Kpdmst = np.zeros(ntsr)


    for i in range(ntsr):
        Omega = 150.0*pi/30.0  # rad/s
        Uinf = Omega*R/tsr[i]

        theta, CP[i], CcR, CcT, CX, CY, CZ, Sref = vawt_ac.evaluate(Uinf, Omega)
        CT[i] = 1/(2*pi)*np.trapz(CX, theta)
        Kp[i] = CP[i]*Uinf**3/(Omega*R)**3

        # CP[i] *= Sref/(D*H)

        theta, CPdmst[i], CcR, CcT, CX, CY, CZ, Sref = vawt_dmst.evaluate(Uinf, Omega)
        CTdmst[i] = 1/(2*pi)*np.trapz(CX, theta)
        Kpdmst[i] = CPdmst[i]*Uinf**3/(Omega*R)**3

        # CPdmst[i] *= Sref/(D*H)




    import matplotlib.pyplot as plt
    plt.plot(tsr, CP, '-ko')
    plt.plot(tsr, CPdmst, '-ro')
    plt.xlabel('$\\lambda$')
    plt.ylim([0, 0.6])
    plt.grid()

    plt.figure()
    plt.plot(tsr, CT, '-ko')
    plt.plot(tsr, CTdmst, '-ro')


    plt.show()



    # # Madsen test case

    # r = 1.0
    # chord = 0.2/3
    # twist = 0.0
    # af = NACA0015_linear()
    # B = 3
    # Uinf = 9.0
    # rho = 1.0
    # mu = 1.0

    # z = [0.0, 2.5, 5.0, 7.5, 10.0]
    # dtheta = 0.0
    # n = len(z)
    # vawt = VAWT(z, [r]*n, [chord]*n, [twist]*n,
    #             [dtheta]*n, [af]*n, B, rho, mu, method='ac')
    # vawt2 = VAWT(z, [r]*n, [chord]*n, [twist]*n,
    #              [dtheta]*n, [af]*n, B, rho, mu, method='dmst')


    # lambda_r = [2, 2.25, 2.5, 3, 4, 5]
    # ntsr = len(lambda_r)
    # CT = np.zeros(ntsr)
    # CP = np.zeros(ntsr)
    # CTdmst = np.zeros(ntsr)
    # CPdmst = np.zeros(ntsr)

    # for i in range(ntsr):
    #     Omega = Uinf*lambda_r[i]/r

    #     theta, CP[i], CcR, CcT, CX, CY, CZ, Sref = vawt.evaluate(Uinf, Omega)
    #     CT[i] = 1/(2*pi)*np.trapz(CX, theta)

    #     theta, CPdmst[i], CcR, CcT, CX, CY, CZ, Sref = vawt2.evaluate(Uinf, Omega)
    #     CTdmst[i] = 1/(2*pi)*np.trapz(CX, theta)



    # import matplotlib.pyplot as plt
    # plt.plot(lambda_r, CT, '-ko')
    # plt.plot(lambda_r, CP, '-ko')
    # plt.plot(lambda_r, CTdmst, '-ro')
    # plt.plot(lambda_r, CPdmst, '-ro')
    # plt.xlabel('$\\lambda$')
    # # plt.legend(loc='upper left')
    # plt.xlim([0, 7])
    # plt.ylim([0, 1.2])
    # plt.grid()
    # plt.show()

    # exit()
