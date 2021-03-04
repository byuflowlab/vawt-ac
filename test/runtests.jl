# using VAWTAC
using Test
using PyPlot
close("all")

const path = splitdir(@__FILE__)[1]
include("$path/../src/VAWTAC.jl")


# @testset "VAWTAC.jl one turbine" begin

    # --------- one turbine -------------
    ntheta = 36

    r = 3.0
    twist = 0.0
    delta = 0.0
    af = VAWTAC.readaerodyn("$path/airfoils/NACA_0012_mod.dat")
    B = 3
    solidity = 0.25
    chord = solidity*r/B

    Vinf = 1.0
    rho = 1.225
    mu = 1.7894e-5

    Omega = 0.0

    # baseline
    turbines = Array{VAWTAC.Turbine}(undef,1)
    turbines[1] = VAWTAC.Turbine(r, chord, twist, delta, B, af, Omega, 0.0, 0.0)
    env = VAWTAC.Environment(Vinf, rho, mu)


    n = 20
    tsrvec = LinRange(1, 7, n)
    cpvec = zeros(n)
    for i in eachindex(tsrvec)
        tsr = tsrvec[i]
        turbines[1].Omega = Vinf*tsr/r
        CT, CP, Rp, Tp, Zp, theta = VAWTAC.actuatorcylinder(turbines, env, ntheta)
        cpvec[i] = CP[1]
    end

    figure()
    plot(tsrvec, cpvec)
    xlabel("\$\\lambda\$")
    ylabel("\$C_p\$")
# end
#
# @testset "VAWTAC.jl two turbines" begin
    # -----------  two turbines -----------------
    tsr = 3.5
    Omega = Vinf*tsr/r

    turbines = Array{VAWTAC.Turbine}(undef,2)
    turbines[1] = VAWTAC.Turbine(r, chord, twist, delta, B, af, Omega, 0.0, 0.0)
    turbines[2] = VAWTAC.Turbine(r, chord, twist, delta, B, af, -Omega, 0.0, 2*r)
    CT, CP, Rp, Tp, Zp, theta = VAWTAC.actuatorcylinder(turbines, env, ntheta)

    figure()
    plot(theta, r*Tp[:, 1])  # r *T = Q
    plot(theta, r*Tp[:, 2])
    xlabel("\$\\theta\$")
    ylabel("Q (torque)")
    xlim([0, 2*pi])

# end
