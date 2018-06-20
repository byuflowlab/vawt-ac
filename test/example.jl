include("acmultiple.jl")
using PyPlot
close("all")


# --------- one turbine -------------
ntheta = 36

r = 3.0
twist = 0.0
delta = 0.0
af = readaerodyn("airfoils/NACA_0012_mod.dat")
B = 3
solidity = 0.25
chord = solidity*r/B

Vinf = 1.0
rho = 1.225
mu = 1.7894e-5

Omega = 0.0

# baseline
turbines = Array{Turbine}(1)
turbines[1] = Turbine(r, chord, twist, delta, B, af, Omega, 0.0, 0.0)
env = Environment(Vinf, rho, mu)


n = 20
tsrvec = linspace(1, 7, n)
cpvec = zeros(n)
for i in eachindex(tsrvec)
    tsr = tsrvec[i]
    turbines[1].Omega = Vinf*tsr/r
    CT, CP, Rp, Tp, Zp, theta = actuatorcylinder(turbines, env, ntheta)
    cpvec[i] = CP[1]
end

figure()
plot(tsrvec, cpvec)
xlabel("\$\\lambda\$")
ylabel("\$C_p\$")


# -----------  two turbines -----------------
tsr = 3.5
Omega = Vinf*tsr/r

turbines = Array{Turbine}(2)
turbines[1] = Turbine(r, chord, twist, delta, B, af, Omega, 0.0, 0.0)
turbines[2] = Turbine(r, chord, twist, delta, B, af, -Omega, 0.0, 2*r)
CT, CP, Rp, Tp, Zp, theta = actuatorcylinder(turbines, env, ntheta)

figure()
plot(theta, r*Tp[:, 1])  # r *T = Q
plot(theta, r*Tp[:, 2])
xlabel("\$\\theta\$")
ylabel("Q (torque)")
xlim([0, 2*pi])
