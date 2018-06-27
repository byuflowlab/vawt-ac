using PyPlot
using vawtac
close("all")
path,_ = splitdir(@__FILE__)
LP="$path/../data/"

using Base.Test

@testset "1 Turbine" begin
# --------- one turbine -------------
ntheta = 36

r = 3.0
twist = 0.0
delta = 0.0
af = vawtac.readaerodyn(LP*"airfoils/NACA_0012_mod.dat")
B = 3
solidity = 0.25
chord = solidity*r/B

Vinf = 1.0
rho = 1.225
mu = 1.7894e-5

Omega = 0.0

# baseline
turbines = Array{vawtac.Turbine}(1)
turbines[1] = vawtac.Turbine(r, chord, twist, delta, B, af, Omega, 0.0, 0.0)
env = vawtac.Environment(Vinf, rho, mu)


n = 20
tsrvec = linspace(1, 7, n)
cpvec = zeros(n)
for i in eachindex(tsrvec)
    tsr = tsrvec[i]
    turbines[1].Omega = Vinf*tsr/r
    CT, CP, Rp, Tp, Zp, theta = vawtac.actuatorcylinder(turbines, env, ntheta)
    cpvec[i] = CP[1]
end
res=[0.0162131, 0.0353388, 0.0660244, 0.111485, 0.18195, 0.284533, 0.401221, 0.464719, 0.472799, 0.459879, 0.437125, 0.409037, 0.375927, 0.338411, 0.296753, 0.250839, 0.200475, 0.145378, 0.0852209, 0.0196674]
@testset "CP Vector, Position $i" for i in 1:20
    @test isapprox(cpvec[i],res[i],atol=1e-6)
end

end #1 turbine test


@testset "2 Turbines" begin
# -----------  two turbines -----------------

ntheta = 36

r = 3.0
twist = 0.0
delta = 0.0
af = vawtac.readaerodyn(LP*"airfoils/NACA_0012_mod.dat")
B = 3
solidity = 0.25
chord = solidity*r/B

Vinf = 1.0
rho = 1.225
mu = 1.7894e-5

Omega = 0.0


tsr = 3.5
Omega = Vinf*tsr/r

turbines = Array{vawtac.Turbine}(2)
turbines[1] = vawtac.Turbine(r, chord, twist, delta, B, af, Omega, 0.0, 0.0)
turbines[2] = vawtac.Turbine(r, chord, twist, delta, B, af, -Omega, 0.0, 2*r)
env = vawtac.Environment(Vinf, rho, mu)


CT, CP, Rp, Tp, Zp, theta = vawtac.actuatorcylinder(turbines, env, ntheta)


tpc1=r*Tp[:,1]
tpc2=r*Tp[:,2]


tp1=[-0.083021, 0.0478864, 0.223784, 0.430954, 0.653567, 0.866291, 1.03769, 1.14349, 1.17436, 1.13524, 1.03798, 0.892611, 0.70451, 0.481848, 0.248877, 0.0503489, -0.0325917, 0.0631221, 0.273867, 0.469695, 0.607922, 0.693477, 0.737407, 0.746173, 0.727156, 0.685551, 0.627816, 0.563094, 0.501119, 0.448174, 0.403128, 0.357709, 0.294945, 0.190764, 0.0395836, -0.0818098]
tp2=[0.0631221, -0.0325917, 0.0503489, 0.248877, 0.481848, 0.70451, 0.892611, 1.03798, 1.13524, 1.17436, 1.14349, 1.03769, 0.866291, 0.653567, 0.430954, 0.223784, 0.0478864, -0.083021, -0.0818098, 0.0395836, 0.190764, 0.294945, 0.357709, 0.403128, 0.448174, 0.501119, 0.563094, 0.627816, 0.685551, 0.727156, 0.746173, 0.737407, 0.693477, 0.607922, 0.469695, 0.273867]


@testset "Torque Vector, Potion $i" for i in 1:36
    @test isapprox(tpc1[i],tp1[i],atol=9e-5)
    @test isapprox(tpc2[i],tp2[i],atol=9e-5)
end

end #2 turbines
