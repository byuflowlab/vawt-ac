module VAWTAC

const modulepath = splitdir(@__FILE__)[1]

include("$modulepath/acmultiple.jl")
include("$modulepath/airfoilread.jl")

end
