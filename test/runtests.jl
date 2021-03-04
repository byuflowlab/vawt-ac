using VAWTAC
using Test
using HDF5
# using PyPlot
# close("all")

const path = splitdir(@__FILE__)[1]
# include("$path/../src/VAWTAC.jl")

atol = 1e-6

#Common variables
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


@testset "VAWTAC.jl one turbine" begin

    # --------- one turbine -------------

    n = 20
    tsrvec = LinRange(1, 7, n)
    CPvec = zeros(n)
    CTvec = zeros(n)
    Rpvec = zeros(n)
    Tpvec = zeros(n)
    Zpvec = zeros(n)
    thetavec = zeros(n,ntheta)
    for i in eachindex(tsrvec)
        tsr = tsrvec[i]
        turbines[1].Omega = Vinf*tsr/r
        CT, CP, Rp, Tp, Zp, theta = VAWTAC.actuatorcylinder(turbines, env, ntheta)
        CPvec[i] = CP[1]
        CTvec[i] = CT[1]
        Rpvec[i] = Rp[1]
        Tpvec[i] = Tp[1]
        Zpvec[i] = Zp[1]
        thetavec[i,:] = theta
    end

    #write to file
    filename = "$path/single_unit_test_data.h5"
    # HDF5.h5open(filename, "w") do file
    #     HDF5.write(file,"CPvec_old",CPvec)
    #     HDF5.write(file,"CTvec_old",CTvec)
    #     HDF5.write(file,"Rpvec_old",Rpvec)
    #     HDF5.write(file,"Tpvec_old",Tpvec)
    #     HDF5.write(file,"Zpvec_old",Zpvec)
    #     HDF5.write(file,"thetavec_old",thetavec)
    # end

    #read from file
    CPvec_old = HDF5.h5read(filename,"CPvec_old")
    CTvec_old = HDF5.h5read(filename,"CTvec_old")
    Rpvec_old = HDF5.h5read(filename,"Rpvec_old")
    Tpvec_old = HDF5.h5read(filename,"Tpvec_old")
    Zpvec_old = HDF5.h5read(filename,"Zpvec_old")
    thetavec_old = HDF5.h5read(filename,"thetavec_old")

    #test
    @test isapprox(CPvec,CPvec_old;atol)
    @test isapprox(CTvec,CTvec_old;atol)
    @test isapprox(Rpvec,Rpvec_old;atol)
    @test isapprox(Tpvec,Tpvec_old;atol)
    @test isapprox(Zpvec,Zpvec_old;atol)
    @test isapprox(thetavec,thetavec_old;atol)

    # figure()
    # plot(tsrvec, CPvec)
    # xlabel("\$\\lambda\$")
    # ylabel("\$C_p\$")
end

@testset "VAWTAC.jl two turbines" begin
    # -----------  two turbines -----------------
    tsr = 3.5
    Omega = Vinf*tsr/r

    turbines = Array{VAWTAC.Turbine}(undef,2)
    turbines[1] = VAWTAC.Turbine(r, chord, twist, delta, B, af, Omega, 0.0, 0.0)
    turbines[2] = VAWTAC.Turbine(r, chord, twist, delta, B, af, -Omega, 0.0, 2*r)
    CT, CP, Rp, Tp, Zp, theta = VAWTAC.actuatorcylinder(turbines, env, ntheta)

    #write to file
    filename = "$path/dual_unit_test_data.h5"
    # HDF5.h5open(filename, "w") do file
    #     HDF5.write(file,"CP_old",CP)
    #     HDF5.write(file,"CT_old",CT)
    #     HDF5.write(file,"Rp_old",Rp)
    #     HDF5.write(file,"Tp_old",Tp)
    #     HDF5.write(file,"Zp_old",Zp)
    #     HDF5.write(file,"theta_old",theta)
    # end

    #read from file
    CP_old = HDF5.h5read(filename,"CP_old")
    CT_old = HDF5.h5read(filename,"CT_old")
    Rp_old = HDF5.h5read(filename,"Rp_old")
    Tp_old = HDF5.h5read(filename,"Tp_old")
    Zp_old = HDF5.h5read(filename,"Zp_old")
    theta_old = HDF5.h5read(filename,"theta_old")

    #test
    @test isapprox(CP,CP_old;atol)
    @test isapprox(CT,CT_old;atol)
    @test isapprox(Rp,Rp_old;atol)
    @test isapprox(Tp,Tp_old;atol)
    @test isapprox(Zp,Zp_old;atol)
    @test isapprox(theta,theta_old;atol)

    # figure()
    # plot(theta, r*Tp[:, 1])  # r *T = Q
    # plot(theta, r*Tp[:, 2])
    # xlabel("\$\\theta\$")
    # ylabel("Q (torque)")
    # xlim([0, 2*pi])

end
