# module AirfoilRead

# export AirfoilData, readaerodyn

import Interpolations: interpolate, Gridded, Linear, GriddedInterpolation

# type AirfoilData
#     alpha::Array{Float64,1}
#     cl::Array{Float64,1}
#     cd::Array{Float64,1}
# end

type AirfoilData
    cl::GriddedInterpolation
    cd::GriddedInterpolation
end

function readaerodyn(filename)
    """currently only reads one Reynolds number if multiple exist"""

    alpha = Float64[]
    cl = Float64[]
    cd = Float64[]

    open(filename) do f

        # skip header
        for i = 1:13
            readline(f)
        end

        # read until EOT
        while true
            line = readline(f)
            if contains(line, "EOT")
                break
            end
            parts = split(line)
            push!(alpha, float(parts[1]))
            push!(cl, float(parts[2]))
            push!(cd, float(parts[3]))
        end
    end

    # af = AirfoilData(alpha*pi/180.0, cl, cd)

    # 1D interpolations for now.  ignoring Re dependence (which is very minor)
    afcl = interpolate((alpha*pi/180.0,), cl, Gridded(Linear()))
    afcd = interpolate((alpha*pi/180.0,), cd, Gridded(Linear()))
    af = AirfoilData(afcl, afcd)

    return af
end

# end

# filename = "airfoils/NACA_0012_mod.dat"
# alpha, cl, cd = readaerodyn(filename)
# using PyPlot
# plot(alpha, cl)
# plot(alpha, cd)
# show()
