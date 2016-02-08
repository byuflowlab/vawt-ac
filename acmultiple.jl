using Root  # using my wrapper to hybrd a nD root finder
using AirfoilRead  # my airfoil file reader
# using PyPlot

# --- Influence Coefficients ---

"""
applies for both Ay and Rx depending on which function ifunc(x, y, phi)
is passed in
"""
function panelIntegration(xvec, yvec, thetavec, ifunc)

    ntheta = length(thetavec)
    dtheta = thetavec[2] - thetavec[1]  # assumes equally spaced
    A = zeros(ntheta, ntheta)

    for i = 1:ntheta
        for j = 1:ntheta
            # redefine function so it has one parameter for use in quadgk
            integrand(phi) = ifunc(xvec[i], yvec[i], phi)

            # an Adaptive Gauss-Kronrod quadrature integration.  May want something simpler like trapz later.
            A[i, j], error = quadgk(integrand, thetavec[j]-dtheta/2, thetavec[j]+dtheta/2, abstol=1e-10)
        end
    end

    A /= (2*pi)
    return A
end

"""
version for VAWTs that are very far apart
"""
function panelIntegrationFar(xc, yc, thetavec, ifunc)

    ntheta = length(thetavec)
    dtheta = thetavec[2] - thetavec[1]  # assumes equally spaced
    A = zeros(ntheta, ntheta)

    for j = 1:ntheta
        # redefine function so it has one parameter for use in quadgk
        integrand(phi) = ifunc(xc, yc, phi)

        # an Adaptive Gauss-Kronrod quadrature integration.  May want something simpler like trapz later.
        A[:, j], error = quadgk(integrand, thetavec[j]-dtheta/2, thetavec[j]+dtheta/2, abstol=1e-10)
    end

    A /= (2*pi)
    return A
end

"""
integrand used for computing Rx
"""
function Rxintegrand(x, y, phi)
    v1 = x + sin(phi)
    v2 = y - cos(phi)
    # v1 and v2 should never both be zero b.c. we never integrate self.  RxII handles that case.
    return (v1*sin(phi) - v2*cos(phi))/(v1^2 + v2^2)
end


"""
integrand used for computing Ay
"""
function Ayintegrand(x, y, phi)
    v1 = x + sin(phi)
    v2 = y - cos(phi)
    if abs(v1) < 1e-12 && abs(v2) < 1e-12  # occurs when integrating self, function symmetric around singularity, should integrate to zero
        return 0.0
    end
    return (v1*cos(phi) + v2*sin(phi))/(v1^2 + v2^2)
end


function AyIJ(xvec, yvec, thetavec)
    return panelIntegration(xvec, yvec, thetavec, Ayintegrand)
end

function RxIJ(xvec, yvec, thetavec)
    return panelIntegration(xvec, yvec, thetavec, Rxintegrand)
end

function AyIJFar(xc, yc, thetavec)
    return panelIntegrationFar(xc, yc, thetavec, Ayintegrand)
end

function RxIJFar(xc, yc, thetavec)
    return panelIntegrationFar(xc, yc, thetavec, Rxintegrand)
end


function WxIJ(xvec, yvec, thetavec)

    ntheta = length(thetavec)
    dtheta = thetavec[2] - thetavec[1]  # assumes equally spaced
    Wx = zeros(ntheta, ntheta)

    for i = 1:ntheta
        # if yvec[i] >= -1 && yvec[i] <= 1 && xvec[i] >= 1
        if yvec[i] >= -1 && yvec[i] <= 1 && xvec[i] >= 0 && xvec[i]^2 + yvec[i]^2 >= 1
            thetak = acos(yvec[i])
            k = findfirst(thetavec + dtheta/2 .> thetak)  # index of intersection
            Wx[i, k] = -1
            Wx[i, ntheta-k+1] = 1
        end
    end

    return Wx
end

function RxII(thetavec)

    ntheta = length(thetavec)
    dtheta = thetavec[2] - thetavec[1]  # assumes equally spaced
    Rx = dtheta/(4*pi)*ones(ntheta, ntheta)

    for i = 1:ntheta
        if i <= ntheta/2
            Rx[i, i] = (-1 + 1/ntheta)/2
        else
            Rx[i, i] = (1 + 1/ntheta)/2
        end
    end

    return Rx
end

function WxII(thetavec)

    ntheta = length(thetavec)
    Wx = zeros(ntheta, ntheta)

    for i = ntheta/2+1:ntheta
        Wx[i, ntheta+1-i] += -1
    end

    return Wx
end

function matrixAssemble(centerX, centerY, radii, ntheta, approxfar=false)

    # setup discretization
    dtheta = 2*pi/ntheta
    theta = [dtheta/2:dtheta:2*pi]

    # initialize matrices
    nturbines = length(radii)
    Rx = zeros(nturbines*ntheta, nturbines*ntheta)
    Wx = zeros(nturbines*ntheta, nturbines*ntheta)
    Ay = zeros(nturbines*ntheta, nturbines*ntheta)

    # precompute (TODO: may want to do this even outside of this function)
    Rxself = RxII(theta)
    Wxself = WxII(theta)
    Ayself = AyIJ(-sin(theta), cos(theta), theta)

    # iterate through turbines
    for I = 1:nturbines
        for J = 1:nturbines

            # find i locations relative to center of turbine J
            x = (centerX[I]-radii[I]*sin(theta) - centerX[J])/radii[J]
            y = (centerY[I]+radii[I]*cos(theta) - centerY[J])/radii[J]

            # self-influence is precomputed
            if I == J
                Rxsub = Rxself
                Wxsub = Wxself
                Aysub = Ayself

            # pairs can be mapped for same radius
            elseif J < I && radii[I] == radii[J]

                # grab cross-diagonal I,J -> J,I matrix
                Rxsub = Rx[(J-1)*ntheta+1:J*ntheta, (I-1)*ntheta+1:I*ntheta]
                Aysub = Ay[(J-1)*ntheta+1:J*ntheta, (I-1)*ntheta+1:I*ntheta]

                # mapping index for coefficients that are the same
                idx = [[ntheta/2+1:ntheta], [1:ntheta/2]]

                # directly map over
                Rxsub = Rxsub[idx, idx]
                Aysub = Aysub[idx, idx]

                # wake term must be recomptued
                Wxsub = WxIJ(x, y, theta)

            # if VAWTs are very far apart we can approximate some of the influence coefficients
            elseif approxfar && sqrt((centerX[I]-centerX[J])^2 + (centerY[I]-centerY[J])^2) > 10*radii[I]
                println("far apart")
                xc = (centerX[I] - centerX[J])/radii[J]
                yc = (centerY[I] - centerY[J])/radii[J]

                Rxsub = RxIJFar(xc, yc, theta)
                Wxsub = zeros(ntheta, ntheta)  # should have negligible wake impact
                Aysub = AyIJFar(xc, yc, theta)

            else
                Rxsub = RxIJ(x, y, theta)
                Wxsub = WxIJ(x, y, theta)
                Aysub = AyIJ(x, y, theta)
            end

            # assemble into global matrix
            # IJidx =
            Rx[(I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta] = Rxsub
            Wx[(I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta] = Wxsub
            Ay[(I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta] = Aysub

        end
    end
    Ax = Rx + Wx

    return Ax, Ay, theta
end
# ---------------------------


# ------- Force Coefficients ---------

function radialforce(uvec, vvec, thetavec, r, chord, twist, delta, af, B, Vinf, Omega, rho, mu, rotation)
    # u, v, theta - arrays of size ntheta
    # r, chord, twist, Ving, Omega, rho, mu - scalars
    # rotation - 1.0 for ccw or -1 for cw

    # velocity components and angles
    # Vn = Vinf*(1 + uvec).*sin(thetavec) - Vinf*vvec.*cos(thetavec)
    # Vt = Vinf*(1 + uvec).*cos(thetavec) + Vinf*vvec.*sin(thetavec) + Omega*r
    Vn = Vinf*(1 + uvec).*sin(thetavec) - Vinf*vvec.*cos(thetavec)
    Vt = rotation*(Vinf*(1 + uvec).*cos(thetavec) + Vinf*vvec.*sin(thetavec)) + Omega*r
    W = sqrt(Vn.^2 + Vt.^2)
    phi = atan2(Vn, Vt)
    alpha = phi - twist
    Re = rho*chord*Vinf/mu

    # airfoil
    cl, cd = afinterp(af, alpha, Re)

    # rotate force coefficients
    cn = cl.*cos(phi) + cd.*sin(phi)
    ct = cl.*sin(phi) - cd.*cos(phi)

    # radial force
    sigma = B*chord/r
    q = sigma/(4*pi)*cn.*(W/Vinf).^2

    # instantaneous forces
    qdyn = 0.5*rho*W.^2
    Rp = -cn.*qdyn*chord
    Tp = ct.*qdyn*chord/cos(delta)
    Zp = -cn.*qdyn*chord*tan(delta)

    # nonlinear correction factor
    # integrand = (W/Vinf).^2 .* (cn.*sin(thetavec) - ct.*cos(thetavec)/cos(delta))
    integrand = (W/Vinf).^2 .* (cn.*sin(thetavec) - rotation*ct.*cos(thetavec)/cos(delta))
    CT = sigma/(4*pi) * pInt(thetavec, integrand)
    if CT > 2
        a = 0.5*(1 + sqrt(1 + CT))
        ka = 1.0 / (a-1)

    elseif CT > 0.96
        a = 1.0/7*(1 + 3*sqrt(7.0/2*CT - 3))
        ka = 18.0*a / (7*a^2 - 2*a + 4)

    else
        a = 0.5*(1 - sqrt(1 - CT))
        ka = 1.0 / (1-a)
    end

    q *= ka

    # power coefficient
    H = 1.0  # per unit depth
    Sref = (2*r*H)
    Q = r.*Tp
    # println(Q)
    P = Omega*B/(2*pi)*pInt(thetavec, Q)
    CP = P / (0.5*rho*Vinf^3 * Sref)

    return q, CT, CP, Rp, Tp, Zp
end


# --- Solution

function residual(w, A, theta, r, chord, twist, delta, af, B, Vinf, Omega, rho, mu, rotation)

    ntheta = length(theta)
    nturbines = int(length(w)/2/ntheta)
    p = zeros(ntheta*nturbines)

    for i = 1:nturbines
        idx = [(i-1)*ntheta+1:i*ntheta]

        u = w[idx]
        v = w[ntheta*nturbines + idx]
        p[idx], _, _ = radialforce(u, v, theta, r, chord, twist, delta, af, B, Vinf, Omega, rho, mu, rotation[i])
    end

    # TODO: add k
    return A*p - w
end


function actuatorcylinder(centerX, centerY, radii, ntheta, r, chord, twist,
                          delta, af, B, Vinf, Omega, rho, mu, rotation)


    Ax, Ay, theta = matrixAssemble(centerX, centerY, radii, ntheta)
    A = [Ax; Ay]
    args = (A, theta, r, chord, twist, delta, af, B, Vinf, Omega, rho, mu, rotation)
    tol = 1e-6
    nturbines = length(centerX)
    w0 = zeros(nturbines*ntheta*2)
    w, zz, info = hybrd(residual, w0, args, tol)

    if info != 1
        println("hybrd terminated prematurely. info = ", info)
    end

    CT = zeros(nturbines)
    CP = zeros(nturbines)
    Rp = zeros(ntheta, nturbines)
    Tp = zeros(ntheta, nturbines)
    Zp = zeros(ntheta, nturbines)

    for i = 1:nturbines
        idx = [(i-1)*ntheta+1:i*ntheta]

        u = w[idx]
        v = w[ntheta*nturbines + idx]
        q, CT[i], CP[i], Rp[:, i], Tp[:, i], Zp[:, i] = radialforce(u, v, theta, r, chord, twist, delta, af, B, Vinf, Omega, rho, mu, rotation[i])


        Vn = Vinf*(1 + u).*sin(theta) - Vinf*v.*cos(theta)
        Vt = rotation[i]*(Vinf*(1 + u).*cos(theta) + Vinf*v.*sin(theta)) + Omega*r
        W = sqrt(Vn.^2 + Vt.^2)
        # println("Vn =", Vn)
        # println("Vt =", Vt)
        # println("W =", W)
        # println("theta =", theta)
        # figure()
        # plot(theta, u)
        # plot(theta, v)
        # figure()
        # plot(theta, q)
        #

        # figure()
        # plot(theta, Rp[:, i])
        # plot(theta, Tp[:, i])
        # plot(theta, Zp[:, i])
        # show()

        # Fn = -Rp[:, i]
        # Ft = Tp[:, i]*cos(delta)
        # figure()
        # plot(theta*180/pi, Fn)
        # plot(theta*180/pi, Ft)
        # figure()
        # plot(theta*180/pi, sqrt((1 + u).^2 + v.^2))
        # figure()
        # plot(theta*180/pi, atan(v./(1 + u))*180/pi)
        #
        # show()

    end


    # println(CT)
    # println(CP)
    # show()

    return CT, CP
end


function linear_interp(x::Float64, xp::Array{Float64,1}, yp::Array{Float64,1})  # Not sure why something like this doesn't exist in core Julia.  Saw Grid, but this is not necessarily equally spaced.

    if x < xp[1] || x > xp[end]
        error("out of interp bounds")
    end

    idx = findfirst(x .<= xp) - 1  # findlast(x .> xp).  findlast doesn't exist in 0.3.  arggh.  jumped on julia bandwagon too early
    if idx == 0  # edge-case if right at first index, equivalent to interpolating in above range.
        idx = 1
    end

    yinterp = yp[idx] + (yp[idx+1] - yp[idx])*(x - xp[idx])/(xp[idx+1] - xp[idx])

    return yinterp
end

function linear_interp(xvec::Array{Float64,1}, xp::Array{Float64,1}, yp::Array{Float64,1})

    yvec = zeros(xvec)
    for i = 1:length(xvec)  # eachindex(xvec)  - not in julia 0.3??
        yvec[i] = linear_interp(xvec[i], xp, yp)
    end

    return yvec
end

# arrg. this too.
function trapz(x::Array{Float64,1}, y::Array{Float64,1})  # integrate y w.r.t. x

    integral = 0
    for i = 1:length(x)-1
        integral += (x[i+1]-x[i])*0.5*(y[i] + y[i+1])
    end
    return integral
end

# integration for a periodic function (trapezoidal method)
function pInt(theta, f)

    integral = trapz(theta, f)

    # add end points
    dtheta = 2*theta[1]  # assumes equally spaced, starts at 0
    integral += dtheta * 0.5*(f[1] + f[end])

    return integral
end

# println(LOAD_PATH)

function afinterp(data::AirfoilData, alpha, Re)

    cl = linear_interp(alpha, data.alpha, data.cl)
    cd = linear_interp(alpha, data.alpha, data.cd)

    return cl, cd
end

function madsen_test()
    alpha = linspace(-pi, pi, 400)
    cl = 2*pi*alpha
    cd = 0.023*ones(400)

    af = AirfoilData(alpha, cl, cd)

    sigma = 0.1
    B = 2
    R = 1.0
    c = sigma*2*R/B
    tsr = 4.0
    twist = 0.0
    delta = 0.0

    centerX = [0.0]
    centerY = [0.0]
    radii = [R]
    rotation = [1.0]

    ntheta = 36

    Vinf = 1.0
    Omega = Vinf*tsr/R
    rho = 1.0
    mu = 1.0

    tsrvec = linspace(2, 5, 10)
    CTvec = zeros(10)
    CPvec = zeros(10)
    for i = 1:10
        Omega = Vinf*tsrvec[i]/R
        CT, CP = actuatorcylinder(centerX, centerY, radii, ntheta, R, c, twist,
                                  delta, af, B, Vinf, Omega, rho, mu, rotation)
        CTvec[i] = CT[1]
        CPvec[i] = CP[1]
    end


    return tsrvec, CTvec, CPvec


end

function dagan_compare()

    ntheta = 36

    r = 3.0
    chord = 0.125  # also solidity if r=3.0
    twist = 0.0
    delta = 0.0
    af = readaerodyn("airfoils/NACA_0021.dat")
    B = 3

    Vinf = 1.0
    tsr = 3.25
    Omega = Vinf*tsr/r
    rho = 1.225
    mu = 1.7894e-5

    centerX = [0.0]
    centerY = [0.0]
    radii = [r]
    rotation = [1.0]

    tsrvec = linspace(2, 7, 10)
    CTvec = zeros(10)
    CPvec = zeros(10)
    for i = 1:10
        Omega = Vinf*tsrvec[i]/r
        CT, CP = actuatorcylinder(centerX, centerY, radii, ntheta, r, chord, twist,
                                  delta, af, B, Vinf, Omega, rho, mu, rotation)
        CTvec[i] = CT[1]
        CPvec[i] = CP[1]
    end

    tsr_cfd = [2, 2.5, 2.75, 3, 3.25, 3.5, 4]
    cp_cfd = [0.111496, 0.256455, 0.31859575, 0.35238, 0.40143025, 0.413539, 0.401584]
    figure()
    plot(tsrvec, CPvec, "-o")
    plot(tsr_cfd, cp_cfd, "-o")
    show()
    sleep(60)


    nD = 1
    nT = 20
    Theta = linspace(-90, 90, nT)*pi/180.0
    Diameters = [2.0]  # linspace(1.0, 6.0, nD)

    D = r*2

    radii = [1.0*r, 1.0*r]
    rotation = [1.0, -1.0]

    # CP1 = zeros(nD, nT)
    # CP2 = zeros(nD, nT)
    CP1 = zeros(nT)
    CP2 = zeros(nT)
    X2 = zeros(nD, nT)
    Y2 = zeros(nD, nT)

    for i = 1:nD
        for j = 1:nT
            centerX = [0.0, Diameters[i]*D*cos(Theta[j])]
            centerY = [0.0, Diameters[i]*D*sin(Theta[j])]

            CT, CP = actuatorcylinder(centerX, centerY, radii, ntheta, r, chord, twist,
                             delta, af, B, Vinf, Omega, rho, mu, rotation)

            X2[i, j] = centerX[2]
            Y2[i, j] = centerY[2]
            CP1[j] = CP[1]
            CP2[j] = CP[2]
        end
    end


    theta_cfd = [90, 80, 70, 60, 50, 40, 30, 20, 10, 0, -10, -20, -30, -40, -50, -60, -70, -80, -90]
    cp1_cfd = [0.439335, 0.43134, 0.42237, 0.4112225, 0.4029025, 0.3952325, 0.393315, 0.4039425, 0.40248, 0.4064775, 0.4064125, 0.4117425, 0.4113525, 0.4173, 0.4242875, 0.428155, 0.431795, 0.44109, 0.46319]
    cp2_cfd = [0.43762225, 0.4517435, 0.4626765, 0.46826, 0.4741295, 0.4800185, 0.46843225, 0.31384275, 0.1020955, 0.0588185, 0.2234895, 0.43908475, 0.47478925, 0.46431125, 0.4580745, 0.453778, 0.45521125, 0.45118775, 0.4619875]


    figure()
    plot(Theta*180.0/pi, CP1)
    plot(theta_cfd, cp1_cfd, "o")

    figure()
    plot(Theta*180.0/pi, CP2)
    plot(theta_cfd, cp2_cfd, "o")

    show()
    sleep(60)

end

# # --- dagan compare ---
# using PyPlot
# dagan_compare()
# exit()
# # ----------------------

# tsr, CT, CP = madsen_test()
# using PyPlot
# figure()
# plot(tsr, CP, "-o", label="C_P")
# plot(tsr, CT, "-o", label="C_T")
# ylim([0, 1.0])
# xlim([0, 7.0])
# show()
# exit()

# function run()
# run script
ntheta = 36

# centerX = [0.0]
# centerY = [0.0]
# radii = [1.0]

r = 3.0
chord = 0.25  # solidity = 0.25
twist = 0.0
delta = 0.0
af = readaerodyn("airfoils/NACA_0012_mod.dat")
B = 3

Vinf = 1.0
# Omega = 127.0*pi/30.0  # rad/s
# tsr = 4.0
# tsr = 2.0
tsr = 3.0
Omega = Vinf*tsr/r
rho = 1.225
mu = 1.7894e-5

centerX = [0.0]
centerY = [0.0]
radii = [r]
rotation = [1.0]

CT0, CP0 = actuatorcylinder(centerX, centerY, radii, ntheta, r, chord, twist,
                          delta, af, B, Vinf, Omega, rho, mu, rotation)

CP0 = CP0[1]

println(CP0)
exit()
# n = 100
# tsrvec = linspace(1.0, 7.0, n)
# cpvec = zeros(n)
# for i = 1:n
#     tsr = tsrvec[i]
#     Omega = Vinf*tsr/r
#     CT0, CP0 = actuatorcylinder(centerX, centerY, radii, ntheta, r, chord, twist,
#                               delta, af, B, Vinf, Omega, rho, mu, rotation)
#     cpvec[i] = CP0[1]
# end
#
# println(tsrvec)
# println(cpvec)
# exit()


# print(CT0, ' ', CP0)
# 0.9241, 0.4479
# exit()

# radii = [1.0*r, 1.0*r]
# rotation = [1.0, 1.0]
# centerX = [0.0, 3*2*r*cos(20*pi/180)]
# centerY = [0.0, 3*2*r*cos(20*pi/180)]
#
# CT1, CP1 = actuatorcylinder(centerX, centerY, radii, ntheta, r, chord, twist,
#                           delta, af, B, Vinf, Omega, rho, mu, rotation)

# print(CP0)
# print(CP1)
#
#
# exit()

# println(CT)
# println(CP)
# end
#
# run()


# @profile run()
# print("here")
# Profile.print()

# [0.9241022702618533]
# [0.4479423198294058]
# [Finished in 2.746s]

# Profile.init(delay=0.01)
# run()
# Profile.clear()
# @profile run()
# using ProfileView
# ProfileView.view()
# using PyPlot
# show()

# exit()

# n = 10
# pos = linspace(2.0, 10.0, n)
# pos = linspace(-5.0, 5.0, n)

# nT = 50
# nD = 18


nT = 20
nD = 10
Theta = linspace(-90, 90, nT)*pi/180.0
Diameters = linspace(1.0, 6.0, nD)

D = r*2

radii = [1.0*r, 1.0*r]
rotation = [1.0, -1.0]

# nT = length(Theta)
# nD = length(Diameters)
CP1 = zeros(nD, nT)
CP2 = zeros(nD, nT)
X2 = zeros(nD, nT)
Y2 = zeros(nD, nT)

for i = 1:nD
    for j = 1:nT
        centerX = [0.0, Diameters[i]*D*cos(Theta[j])]
        centerY = [0.0, Diameters[i]*D*sin(Theta[j])]

        CT, CP = actuatorcylinder(centerX, centerY, radii, ntheta, r, chord, twist,
                         delta, af, B, Vinf, Omega, rho, mu, rotation)

        X2[i, j] = centerX[2]
        Y2[i, j] = centerY[2]
        CP1[i, j] = CP[1]
        CP2[i, j] = CP[2]
    end
    println(i)
end

println("X2 = ", X2)
println("Y2 = ", Y2)
println("CP0 = ", CP0)
println("CP1 = ", CP1)
println("CP2 = ", CP2)


using PyPlot

# figure()
# contourf(X2./D, Y2./D, CP1./CP0, linspace(0.85, 1.1, 50))
# contourf(X2./D, Y2./D, CP1./CP0, linspace(0.85, 1.1, 50))
# xlabel("x/D")
# ylabel("y/D")
# colorbar()
#
# figure()
# contourf(X2./D, Y2./D, CP2./CP0, linspace(0.85, 1.1, 50))
# contourf(X2./D, Y2./D, CP2./CP0, linspace(0.85, 1.1, 50))
# xlabel("x/D")
# ylabel("y/D")
# colorbar()
#
# figure()
# contourf(X2./D, Y2./D, (CP1+CP2)./(2*CP0), linspace(0.3, 1.0, 50))
# contourf(X2./D, Y2./D, (CP1+CP2)./(2*CP0), linspace(0.3, 1.0, 50))
# xlabel("x/D")
# ylabel("y/D")
# colorbar()

scale = linspace(0.875, 1.125, 50)
ticks = 0.875:0.025:1.126

cmap = get_cmap("coolwarm")

figure()
contourf(X2./D, Y2./D, CP1./CP0, scale, cmap=cmap)
contourf(X2./D, Y2./D, CP1./CP0, scale, cmap=cmap)
xlabel("x/D")
ylabel("y/D")
colorbar(ticks=ticks)
axis("equal")

figure()
contourf(X2./D, Y2./D, CP2./CP0, scale, cmap=cmap)
contourf(X2./D, Y2./D, CP2./CP0, scale, cmap=cmap)
xlabel("x/D")
ylabel("y/D")
colorbar(ticks=ticks)
axis("equal")


scale = linspace(0.98, 1.02, 50)
ticks = 0.98:0.005:1.02

figure()
contourf(X2./D, Y2./D, (CP1+CP2)./(2*CP0), scale, cmap=cmap)
contourf(X2./D, Y2./D, (CP1+CP2)./(2*CP0), scale, cmap=cmap)
xlabel("x/D")
ylabel("y/D")
colorbar(ticks=ticks)
axis("equal")


show()
sleep(60)


# for i = 1:length(pos)
#
#     centerX = [0.0*D, 3.0*D]
#     centerY = [0.0*D, pos[i]*D]
#     radii = [1.0*r, 1.0*r]
#
#     CT, CP = actuatorcylinder(centerX, centerY, radii, ntheta, r, chord, twist,
#                      delta, af, B, Vinf, Omega, rho, mu)
#
#     P1[i] = CP[1]
#     P2[i] = CP[2]
# end

# using PyPlot
# plot(pos, P1)
# plot(pos, P2, "--")
# show()
#
# Ax, Ay, theta = matrixAssemble(centerX, centerY, radii, ntheta)
#
# approxfar = true
# Ax2, Ay2, theta2 = matrixAssemble(centerX, centerY, radii, ntheta, approxfar)
#
#
# println("Ax = ", maximum(abs(Ax - Ax2)))
# println("Ay = ", maximum(abs(Ay - Ay2)))
# println(Ax)
# println(Ax2)
