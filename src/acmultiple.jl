import NLsolve
import HDF5
import QuadGK


# --- Influence Coefficients ---

"""
applies for both Ay and Rx depending on which function ifunc(x, y, phi)
is passed in
"""
function panelIntegration(xvec, yvec, thetavec, ifunc)

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
            A[i, j], error = QuadGK.quadgk(integrand, thetavec[j]-dtheta/2.0, thetavec[j]+dtheta/2.0, atol=1e-10)
        end

    end

    return A
end


"""
integrand used for computing Dx
"""
function Dxintegrand(x, y, phi)
    v1 = x + sin(phi)
    v2 = y - cos(phi)
    # v1 and v2 should never both be zero b.c. we never integrate self.  RxII handles that case.
    return (v1*sin(phi) - v2*cos(phi))/(2*pi*(v1*v1 + v2*v2))
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
    return (v1*cos(phi) + v2*sin(phi))/(2*pi*(v1*v1 + v2*v2))
end


function AyIJ(xvec, yvec, thetavec)
    return panelIntegration(xvec, yvec, thetavec, Ayintegrand)
end

function DxIJ(xvec, yvec, thetavec)
    return panelIntegration(xvec, yvec, thetavec, Dxintegrand)
end

function WxIJ(xvec, yvec, thetavec)

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

function DxII(thetavec)

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

function WxII(thetavec)

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
    Ayself = AyIJ(-sin.(theta), cos.(theta), theta)

    # write to file
    HDF5.h5open("$modulepath/theta-$ntheta.h5", "w") do file
        HDF5.write(file, "theta", theta)
        HDF5.write(file, "Dx", Dxself)
        HDF5.write(file, "Wx", Wxself)
        HDF5.write(file, "Ay", Ayself)
    end

end


function matrixAssemble(centerX, centerY, radii, ntheta)
    """
    centerX, centerY: array of x,y coordinates for centers of the VAWTs in the farm
    radii: corresponding array of their radii
    """

    file = "$modulepath/theta-$ntheta.h5"
    if !isfile(file)
        precomputeMatrices(ntheta)
    end

    theta = HDF5.h5read(file, "theta")
    Dxself = HDF5.h5read(file, "Dx")
    Wxself = HDF5.h5read(file, "Wx")
    Ayself = HDF5.h5read(file, "Ay")

    # initialize global matrices
    nturbines = length(radii)
    Dx = zeros(nturbines*ntheta, nturbines*ntheta)
    Wx = zeros(nturbines*ntheta, nturbines*ntheta)
    Ay = zeros(nturbines*ntheta, nturbines*ntheta)

    # iterate through turbines
    for I in eachindex(radii)
        for J in eachindex(radii)

            # find normalized i locations relative to center of turbine J
            x = (centerX[I].-radii[I]*sin.(theta) .- centerX[J])/radii[J]
            y = (centerY[I].+radii[I]*cos.(theta) .- centerY[J])/radii[J]

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

mutable struct Turbine{TF,TI,TFN}
    r::TF
    chord::TF
    twist::TF
    delta::TF
    B::TI
    af::TFN
    Omega::TF
    centerX::TF
    centerY::TF
end

mutable struct Environment{TF}
    Vinf::TF
    rho::TF
    mu::TF
end

function radialforce(uvec, vvec, thetavec,
    turbine::Turbine, env)
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
    Vn = Vinf*(1.0 .+ uvec).*sin.(thetavec) .- Vinf*vvec.*cos.(thetavec)
    Vt = rotation*(Vinf*(1.0 .+ uvec).*cos.(thetavec) .+ Vinf*vvec.*sin.(thetavec)) .+ abs(Omega)*r
    W = sqrt.(Vn.^2 .+ Vt.^2)
    phi = atan.(Vn, Vt)
    alpha = phi .- twist
    # Re = rho*W*chord/mu  # currently no Re dependence

    # airfoil
    # cl = turbine.af.cl[alpha]
    # cd = turbine.af.cd[alpha]
    cl, cd = turbine.af(alpha)

    # rotate force coefficients
    cn = cl.*cos.(phi) .+ cd.*sin.(phi)
    ct = cl.*sin.(phi) .- cd.*cos.(phi)

    # radial force
    sigma = B*chord/r
    q = sigma/(4*pi)*cn.*(W/Vinf).^2

    # instantaneous forces
    qdyn = 0.5*rho*W.^2
    Rp = -cn.*qdyn*chord
    Tp = ct.*qdyn*chord/cos(delta)
    Zp = -cn.*qdyn*chord*tan(delta)

    # nonlinear correction factor
    integrand = (W/Vinf).^2 .* (cn.*sin.(thetavec) .- rotation*ct.*cos.(thetavec)/cos(delta))
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

function residual(w, A, theta,
    k, turbines, env)

    # setup
    ntheta = length(theta)
    nturbines = Int(length(w)/2/ntheta)
    q = zeros(ntheta*nturbines)
    ka = 0.0

    for i = 1:nturbines
        idx = collect((i-1)*ntheta+1:i*ntheta)

        u = w[idx]
        v = w[ntheta*nturbines .+ idx]

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


function actuatorcylinder(turbines, env, ntheta)

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

    for i = 1:nturbines
        w0 = zeros(ntheta*2)

        idx = collect((i-1)*ntheta+1:i*ntheta)

        resid_single(x) = residual(x,[Ax[idx, idx]; Ay[idx, idx]], theta, [1.0], turbines, env)
        result = NLsolve.nlsolve(resid_single, w0, ftol=tol)
        w = result.zero
        if !NLsolve.converged(result)
            println("NLsolve terminated prematurely. info = ", info)
        end

        idx = collect(1:ntheta)
        u = w[idx]
        v = w[ntheta .+ idx]
        q, k[i], CT[i], CP[i], Rp[:, i], Tp[:, i], Zp[:, i] = radialforce(u, v, theta, turbines[i], env)

    end


    if nturbines == 1
        return CT, CP, Rp, Tp, Zp, theta
    end

    # Solve coupled system
    w0 = zeros(nturbines*ntheta*2)

    resid_multiple(x) = residual(x,[Ax; Ay], theta, k, turbines, env)
    result = NLsolve.nlsolve(resid_multiple, w0, ftol=tol)
    w = result.zero
    if !NLsolve.converged(result)
        println("NLsolve terminated prematurely. info = ", info)
    end

    for i = 1:nturbines
        idx = collect((i-1)*ntheta+1:i*ntheta)

        u = w[idx]
        v = w[ntheta*nturbines .+ idx]
        _, _, CT[i], CP[i], Rp[:, i], Tp[:, i], Zp[:, i] = radialforce(u, v, theta, turbines[i], env)

    end

    return CT, CP, Rp, Tp, Zp, theta
end

# -----------------------------------------

# ---------- helper methods --------------

# trapezoidal integration
function trapz(x, y)  # integrate y w.r.t. x

    integral = 0.0
    for i = 1:length(x)-1
        integral += (x[i+1]-x[i])*0.5*(y[i] + y[i+1])
    end
    return integral
end

# integration for a periodic function where end points don't reach ends (uses trapezoidal method)
function pInt(theta, f)

    integral = trapz(theta, f)

    # add end points
    dtheta = 2*theta[1]  # assumes equally spaced, starts at 0
    integral += dtheta * 0.5*(f[1] + f[end])

    return integral
end
