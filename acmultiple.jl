include("airfoilread.jl")  # my airfoil file reader
push!(LOAD_PATH, "minpack")
using Root  # mywrapper to minpack
using HDF5


# --- Influence Coefficients ---

"""
applies for both Ay and Rx depending on which function ifunc(x, y, phi)
is passed in
"""
function panelIntegration(xvec::Array{Float64,1}, yvec::Array{Float64,1}, thetavec::Array{Float64,1}, ifunc)

    # initialize
    nx = length(xvec)
    ntheta = length(thetavec)
    dtheta = thetavec[2] - thetavec[1]  # assumes equally spaced
    A = zeros(nx, ntheta)

    for i in eachindex(xvec)
        # redefine function so it has one parameter for use in quadgk
        integrand(phi) = ifunc(xvec[i], yvec[i], phi)

        for j in eachindex(thetavec)
            # an Adaptive Gauss-Kronrod quadrature integration.  Tried trapz but this was faster.
            A[i, j], error = quadgk(integrand, thetavec[j]-dtheta/2.0, thetavec[j]+dtheta/2.0, abstol=1e-10)
        end

    end

    return A
end


"""
integrand used for computing Dx
"""
function Dxintegrand(x::Float64, y::Float64, phi::Float64)
    v1 = x + sin(phi)
    v2 = y - cos(phi)
    # v1 and v2 should never both be zero b.c. we never integrate self.  RxII handles that case.
    return (v1*sin(phi) - v2*cos(phi))/(2*pi*(v1*v1 + v2*v2))
end


"""
integrand used for computing Ay
"""
function Ayintegrand(x::Float64, y::Float64, phi::Float64)
    v1 = x + sin(phi)
    v2 = y - cos(phi)
    if abs(v1) < 1e-12 && abs(v2) < 1e-12  # occurs when integrating self, function symmetric around singularity, should integrate to zero
        return 0.0
    end
    return (v1*cos(phi) + v2*sin(phi))/(2*pi*(v1*v1 + v2*v2))
end


function AyIJ(xvec::Array{Float64,1}, yvec::Array{Float64,1}, thetavec::Array{Float64,1})
    return panelIntegration(xvec, yvec, thetavec, Ayintegrand)
end

function DxIJ(xvec::Array{Float64,1}, yvec::Array{Float64,1}, thetavec::Array{Float64,1})
    return panelIntegration(xvec, yvec, thetavec, Dxintegrand)
end

function WxIJ(xvec::Array{Float64,1}, yvec::Array{Float64,1}, thetavec::Array{Float64,1})

    # initialize
    nx = length(xvec)
    ntheta = length(thetavec)
    dtheta = thetavec[2] - thetavec[1]  # assumes equally spaced
    Wx = zeros(nx, ntheta)

    for i in eachindex(xvec)
        if yvec[i] >= -1.0 && yvec[i] <= 1.0 && xvec[i] >= 0.0 && xvec[i]^2 + yvec[i]^2 >= 1.0
        # if yvec[i] >= -1.0 && yvec[i] <= 1.0 && (xvec[i] >= 0.0 || (xvec[i] >= -1 && xvec[i]^2 + yvec[i]^2 <= 1.0))
            thetak = acos(yvec[i])
            k = findfirst(thetavec + dtheta/2 .> thetak)  # index of intersection
            Wx[i, k] = -1.0
            Wx[i, ntheta-k+1] = 1.0
        end
    end

    return Wx
end

function DxII(thetavec::Array{Float64,1})

    # initialize
    ntheta = length(thetavec)
    dtheta = thetavec[2] - thetavec[1]  # assumes equally spaced
    Rx = dtheta/(4*pi)*ones(ntheta, ntheta)

    for i in eachindex(thetavec)
        if i <= ntheta/2
            Rx[i, i] = (-1 + 1.0/ntheta)/2.0
        else
            Rx[i, i] = (1 + 1.0/ntheta)/2.0
        end
    end

    return Rx
end

function WxII(thetavec::Array{Float64,1})

    # initialize
    ntheta = length(thetavec)
    Wx = zeros(ntheta, ntheta)

    for i = div(ntheta,2)+1:ntheta
        Wx[i, ntheta+1-i] = -1
    end

    return Wx
end

function precomputeMatrices(ntheta)

    # precompute self influence matrices

    # setup discretization (all the same, and uniformly spaced in theta)
    dtheta = 2*pi/ntheta
    theta = collect(dtheta/2:dtheta:2*pi)

    Dxself = DxII(theta)
    Wxself = WxII(theta)
    Ayself = AyIJ(-sin(theta), cos(theta), theta)

    # write to file
    h5open("theta-$ntheta.h5", "w") do file
        write(file, "theta", theta)
        write(file, "Dx", Dxself)
        write(file, "Wx", Wxself)
        write(file, "Ay", Ayself)
    end

end


function matrixAssemble(centerX::Array{Float64,1}, centerY::Array{Float64,1}, radii::Array{Float64,1}, ntheta::Int64)
    """
    centerX, centerY: array of x,y coordinates for centers of the VAWTs in the farm
    radii: corresponding array of their radii
    """

    file = "theta-$ntheta.h5"
    if !isfile(file)
        precomputeMatrices(ntheta)
    end

    theta = h5read(file, "theta")
    Dxself = h5read(file, "Dx")
    Wxself = h5read(file, "Wx")
    Ayself = h5read(file, "Ay")

    # initialize global matrices
    nturbines = length(radii)
    Dx = zeros(nturbines*ntheta, nturbines*ntheta)
    Wx = zeros(nturbines*ntheta, nturbines*ntheta)
    Ay = zeros(nturbines*ntheta, nturbines*ntheta)

    # iterate through turbines
    for I in eachindex(radii)
        for J in eachindex(radii)

            # find normalized i locations relative to center of turbine J
            x = (centerX[I]-radii[I]*sin(theta) - centerX[J])/radii[J]
            y = (centerY[I]+radii[I]*cos(theta) - centerY[J])/radii[J]

            # self-influence is precomputed
            if I == J
                Dxsub = Dxself
                Wxsub = Wxself
                Aysub = Ayself

            # pairs can be mapped for same radius
            elseif J < I && radii[I] == radii[J]

                # grab cross-diagonal I,J -> J,I matrix
                Dxsub = Dx[(J-1)*ntheta+1:J*ntheta, (I-1)*ntheta+1:I*ntheta]
                Aysub = Ay[(J-1)*ntheta+1:J*ntheta, (I-1)*ntheta+1:I*ntheta]

                # mapping index for coefficients that are the same
                idx = [div(ntheta,2)+1:ntheta; 1:div(ntheta,2)]

                # directly map over
                Dxsub = Dxsub[idx, idx]
                Aysub = Aysub[idx, idx]

                # wake term must be recomptued
                Wxsub = WxIJ(x, y, theta)

            # # if VAWTs are very far apart we can approximate some of the influence coefficients
            # elseif approxfar && sqrt((centerX[I]-centerX[J])^2 + (centerY[I]-centerY[J])^2) > 10*radii[I]
            #     println("far apart")
            #     xc = (centerX[I] - centerX[J])/radii[J]
            #     yc = (centerY[I] - centerY[J])/radii[J]

            #     Rxsub = RxIJFar(xc, yc, theta)
            #     Wxsub = zeros(ntheta, ntheta)  # should have negligible wake impact
            #     Aysub = AyIJFar(xc, yc, theta)

            else
                Dxsub = DxIJ(x, y, theta)
                Wxsub = WxIJ(x, y, theta)
                Aysub = AyIJ(x, y, theta)
            end

            # assemble into global matrix
            Dx[(I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta] = Dxsub
            Wx[(I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta] = Wxsub
            Ay[(I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta] = Aysub

        end
    end
    Ax = Dx + Wx

    return Ax, Ay, theta
end
# ---------------------------


# ------- Force Coefficients ---------

type Turbine
    r::Float64
    chord::Float64
    twist::Float64
    delta::Float64
    B::Int64
    af::AirfoilData
    Omega::Float64
    centerX::Float64
    centerY::Float64
end

type Environment
    Vinf::Float64
    rho::Float64
    mu::Float64
end


function radialforce(uvec::Array{Float64,1}, vvec::Array{Float64,1}, thetavec::Array{Float64,1},
    turbine::Turbine, env::Environment)
    # u, v, theta - arrays of size ntheta
    # r, chord, twist, Vinf, Omega, rho, mu - scalars

    # unpack
    r = turbine.r
    chord = turbine.chord
    twist = turbine.twist
    delta = turbine.delta
    B = turbine.B
    Omega = turbine.Omega

    Vinf = env.Vinf
    rho = env.rho

    # set the rotation direction
    rotation = sign(Omega)

    # velocity components and angles
    Vn = Vinf*(1.0 + uvec).*sin(thetavec) - Vinf*vvec.*cos(thetavec)
    Vt = rotation*(Vinf*(1.0 + uvec).*cos(thetavec) + Vinf*vvec.*sin(thetavec)) + abs(Omega)*r
    W = sqrt(Vn.^2 + Vt.^2)
    phi = atan2(Vn, Vt)
    alpha = phi - twist
    # Re = rho*W*chord/mu  # currently no Re dependence

    # airfoil
    cl = turbine.af.cl[alpha]
    cd = turbine.af.cd[alpha]

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
    integrand = (W/Vinf).^2 .* (cn.*sin(thetavec) - rotation*ct.*cos(thetavec)/cos(delta))
    CT = sigma/(4*pi) * pInt(thetavec, integrand)
    if CT > 2.0
        a = 0.5*(1.0 + sqrt(1.0 + CT))
        ka = 1.0 / (a-1)

    elseif CT > 0.96
        a = 1.0/7*(1 + 3.0*sqrt(7.0/2*CT - 3))
        ka = 18.0*a / (7*a^2 - 2*a + 4)

    else
        a = 0.5*(1 - sqrt(1.0 - CT))
        ka = 1.0 / (1-a)
    end

    # power coefficient
    H = 1.0  # per unit height
    Sref = 2*r*H
    Q = r*Tp
    P = abs(Omega)*B/(2*pi)*pInt(thetavec, Q)
    CP = P / (0.5*rho*Vinf^3 * Sref)

    return q, ka, CT, CP, Rp, Tp, Zp
end

# -----------------------------------------



# ------ Solve System --------------

function residual(w::Array{Float64,1}, A::Array{Float64,2}, theta::Array{Float64,1},
    k::Array{Float64,1}, turbines::Array{Turbine,1}, env::Environment)

    # setup
    ntheta = length(theta)
    nturbines = length(turbines)  #  int(length(w)/2/ntheta)
    q = zeros(ntheta*nturbines)
    ka = 0.0

    for i in eachindex(turbines)
        idx = (i-1)*ntheta+1:i*ntheta

        u = w[idx]
        v = w[ntheta*nturbines + idx]

        q[idx], ka, _, _, _, _, _ = radialforce(u, v, theta, turbines[i], env)
    end

    if nturbines == 1  # if only one turbine use the k from the analysis
        k = [ka]
    end  # otherwise, use k that was input to this function

    # reformat to multiply in correct locations
    kmult = repeat(k, inner=[ntheta])
    kmult = [kmult; kmult]

    return (A*q).*kmult - w
end


function actuatorcylinder(turbines::Array{Turbine,1}, env::Environment, ntheta::Int64)

    # list comprehensions
    centerX = [turbine.centerX for turbine in turbines]
    centerY = [turbine.centerY for turbine in turbines]
    radii = [turbine.r for turbine in turbines]

    # assemble global matrices
    Ax, Ay, theta = matrixAssemble(centerX, centerY, radii, ntheta)

    # setup
    ntheta = length(theta)
    nturbines = length(turbines)
    tol = 1e-6
    CT = zeros(nturbines)
    CP = zeros(nturbines)
    Rp = zeros(ntheta, nturbines)
    Tp = zeros(ntheta, nturbines)
    Zp = zeros(ntheta, nturbines)

    q = zeros(ntheta)

    # compute nonlinear correction factors (each turbine individaully)
    k = zeros(nturbines)

    for i in eachindex(turbines)
        w0 = zeros(ntheta*2)

        idx = (i-1)*ntheta+1:i*ntheta
        args = ([Ax[idx, idx]; Ay[idx, idx]], theta, [1.0], [turbines[i]], env)

        w, zz, info = hybrd(residual, w0, args, tol)

        if info != 1
            println("hybrd terminated prematurely. info = ", info)
        end

        idx = 1:ntheta
        u = w[idx]
        v = w[ntheta + idx]
        q, k[i], CT[i], CP[i], Rp[:, i], Tp[:, i], Zp[:, i] = radialforce(u, v, theta, turbines[i], env)

    end


    if nturbines == 1
        return CT, CP, Rp, Tp, Zp, theta
    end

    # Solve coupled system
    w0 = zeros(nturbines*ntheta*2)
    args = ([Ax; Ay], theta, k, turbines, env)

    w, zz, info = hybrd(residual, w0, args, tol)

    if info != 1
        println("hybrd terminated prematurely. info = ", info)
    end

    for i in eachindex(turbines)
        idx = (i-1)*ntheta+1:i*ntheta

        u = w[idx]
        v = w[ntheta*nturbines + idx]
        _, _, CT[i], CP[i], Rp[:, i], Tp[:, i], Zp[:, i] = radialforce(u, v, theta, turbines[i], env)

    end

    return CT, CP, Rp, Tp, Zp, theta
end

# -----------------------------------------

# ---------- helper methods --------------

# trapezoidal integration
function trapz(x::Array{Float64,1}, y::Array{Float64,1})  # integrate y w.r.t. x

    integral = 0.0
    for i = 1:length(x)-1
        integral += (x[i+1]-x[i])*0.5*(y[i] + y[i+1])
    end
    return integral
end

# integration for a periodic function where end points don't reach ends (uses trapezoidal method)
function pInt(theta::Array{Float64,1}, f::Array{Float64,1})

    integral = trapz(theta, f)

    # add end points
    dtheta = 2*theta[1]  # assumes equally spaced, starts at 0
    integral += dtheta * 0.5*(f[1] + f[end])

    return integral
end

# define fixed conditions geometry

# ntheta = 36
# precomputeMatrices(ntheta)

# ntheta = 36
#
# r = 3.0
# twist = 0.0
# delta = 0.0
# af = readaerodyn("airfoils/NACA_0012_mod.dat")
# B = 3
# solidity = 0.25
# chord = solidity*r/B  # solidity = 0.25
#
# Vinf = 1.0
# rho = 1.225
# mu = 1.7894e-5
# tsr = 3.5
# Omega = Vinf*tsr/r
#
# # ccw1 = true
#
# turbines = Array{Turbine}(2)
# turbines[1] = Turbine(r, chord, twist, delta, B, af, Omega, 0.0, 0.0)
# turbines[2] = Turbine(r, chord, twist, delta, B, af, -Omega, 0.0, 2*r)
#
# env = Environment(Vinf, rho, mu)
#
# using PyPlot
# close("all")
#
# CT, CP, Rp, Tp, Zp, theta = actuatorcylinder(turbines, env, ntheta)
# println(CP)
#
#
#
#
# # # println(CT)
# # # println(CP/0.47318637495058935 -1)
# #
# figure(1)
# plot(theta, r*Tp)
#
# figure(2)
# plot(theta, -Rp)
#
#
# turbines[1] = Turbine(r, chord, twist, delta, B, af, -Omega, 0.0, 0.0)
# CT, CP, Rp, Tp, Zp, theta = actuatorcylinder(turbines, env, ntheta)
#
# figure(1)
# plot(theta, r*Tp)
#
# figure(2)
# plot(theta, -Rp)
#
# turbines = Array{Turbine}(2)
# turbines[1] = Turbine(r, chord, twist, delta, B, af, Omega, 0.0, 0.0)
# turbines[2] = Turbine(r, chord, twist, delta, B, af, -Omega, 0.0, -2*r)
# CT, CP, Rp, Tp, Zp, theta = actuatorcylinder(turbines, env, ntheta)
#
# figure(1)
# plot(theta, r*Tp)
#
# figure(2)
# plot(theta, -Rp)






# ntheta = 6
# dtheta = 2*pi/ntheta
# theta2 = collect(dtheta/2:dtheta:2*pi)
#
# x = collect(linspace(-1.25, 1.25, 4))
# y = 1.25*ones(4)
# Dxsub = DxIJ(x, y, theta2)
# # println(Dxsub)
# Dxsub = DxII(theta2)
# nothing

# #
#
# function mytest(n)
#     for i = 1:n
#         CT, CP, Rp, Tp, Zp = actuatorcylinder(turbines, env, ntheta)
#         println(CT, CP)
#     end
# end
#
# mytest(1)
# #
# @time mytest(20)
# nothing

# using ProfileView
# Profile.clear()
# @profile mytest(20)
# # # Profile.print()

# # nothing
# ProfileView.view()
# ProfileView.svgwrite("profile_results.svg")
# show()

# turbines[2] = Turbine(r, chord, twist, delta, B, af, ccw2, 0.0, 0.0)

# # ---- testing ---
# x = [-0.087155742747658166, -0.25881904510252074, -0.42261826174069939, -0.57357643635104594, -0.70710678118654746, -0.81915204428899169, -0.90630778703664983, -0.9659258262890682, -0.99619469809174555, -0.99619469809174555, -0.96592582628906842, -0.90630778703665027, -0.81915204428899224, -0.70710678118654791, -0.57357643635104671, -0.42261826174069989, -0.25881904510252141, -0.087155742747659082, 0.0871557427476575, 0.25881904510251991, 0.42261826174069889, 0.57357643635104538, 0.7071067811865468, 0.81915204428899102, 0.90630778703664938, 0.96592582628906809, 0.99619469809174543, 0.99619469809174566, 0.96592582628906865, 0.9063077870366506, 0.81915204428899235, 0.70710678118654835, 0.57357643635104716, 0.42261826174070077, 0.2588190451025224, 0.087155742747659207]
# y = [0.99619469809174555, 0.96592582628906831, 0.90630778703665005, 0.81915204428899191, 0.70710678118654768, 0.57357643635104627, 0.42261826174069983, 0.25881904510252118, 0.087155742747658568, -0.087155742747657791, -0.25881904510252041, -0.42261826174069894, -0.57357643635104538, -0.70710678118654713, -0.81915204428899135, -0.90630778703664971, -0.96592582628906809, -0.99619469809174543, -0.99619469809174555, -0.96592582628906853, -0.90630778703665027, -0.81915204428899224, -0.70710678118654824, -0.57357643635104716, -0.42261826174070072, -0.25881904510252146, -0.087155742747659137, 0.087155742747657014, 0.25881904510251946, 0.422618261740698, 0.57357643635104538, 0.70710678118654668, 0.81915204428899102, 0.90630778703664938, 0.96592582628906787, 0.99619469809174543]
# theta = [0.087266462599716474, 0.26179938779914941, 0.43633231299858233, 0.61086523819801519, 0.78539816339744817, 0.95993108859688114, 1.1344640137963138, 1.3089969389957468, 1.4835298641951797, 1.6580627893946127, 1.8325957145940457, 2.0071286397934784, 2.1816615649929112, 2.3561944901923444, 2.5307274153917771, 2.7052603405912103, 2.8797932657906431, 3.0543261909900759, 3.228859116189509, 3.4033920413889418, 3.577924966588375, 3.7524578917878078, 3.9269908169872405, 4.1015237421866733, 4.276056667386106, 4.4505895925855397, 4.6251225177849724, 4.7996554429844052, 4.9741883681838379, 5.1487212933832707, 5.3232542185827043, 5.4977871437821371, 5.6723200689815698, 5.8468529941810026, 6.0213859193804353, 6.195918844579869]

# @time AyIJ(x, y, theta)
# @time AyIJ(x, y, theta)
# # Ay = AyIJ(x, y, theta)
# # println(Ay[1, 1])
# # println(Ay[1, 4])
# # println(Ay[1, 6])
# # println(Ay[3, 7])
# # println(Ay[33, 23])
# # println(Ay[15, 16])
# # println(Ay[3, 23])
# # println(Ay[end, end])

# # 0.012839 seconds (186.89 k allocations: 3.564 MB, 30.58% gc time)
# # -1.1079654327854169e-17
# # 0.0523334537067746
# # 0.029891290254785752
# # 0.03836839306486461
# # -0.011666784681921115
# # 0.17444530820614323
# # -0.002450589617136626
# # -1.5343167042433112e-16

# r = 3.0
# ntheta = 36

# centerX = [0.0, 6.0*r]
# centerY = [0.0, 2.0*r]
# radii = [1.0*r, 1.0*r]

# Ax, Ay, theta = matrixAssemble(centerX, centerY, radii, ntheta)

# af = readaerodyn("airfoils/NACA_0012_mod.dat")
# alpha = collect(0:.01:.1)
# println(af.cl[alpha])
# println(af.cd[alpha])


# interpolate(af.alpha, afl.cl, Gridded(Linear()))


# ntheta = 36

# r = 3.0
# chord = 0.25  # solidity = 0.25
# twist = 0.0
# delta = 0.0
# af = readaerodyn("airfoils/NACA_0012_mod.dat")
# B = 3

# Vinf = 1.0
# tsr = 3.0
# Omega = Vinf*tsr/r
# rho = 1.225
# mu = 1.7894e-5

# RSEP = 6.0
# THETA = 30.0*pi/180

# turbines = Array{Turbine}(1)

# turbines[1] = Turbine(r, chord, twist, delta, B, af, true, 0.0, 0.0)
# # turbines[2] = Turbine(r, chord, twist, delta, B, af, true, RSEP*r*cos(THETA), RSEP*r*sin(THETA))

# env = Environment(Vinf, Omega, rho, mu)


# CT0, CP0, Rp, Tp, Zp = actuatorcylinder(turbines, env, ntheta)
# println(CT0, CP0)
# # @time actuatorcylinder(turbines, env, ntheta)
# # @time actuatorcylinder(turbines, env, ntheta)
# # println("here")

# # Profile.clear()
# # @profile actuatorcylinder(turbines, env, ntheta)
# # Profile.print()

# # using ProfileView
# # ProfileView.view()




# [0.7071796277049306][0.421509550862392]
#   0.161361 seconds (220.61 k allocations: 8.910 MB)
