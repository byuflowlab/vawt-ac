# VAWT Multiple Actuator Cylinder Model

[![DOI](https://zenodo.org/badge/51314790.svg)](https://zenodo.org/badge/latestdoi/51314790)

Author: Andrew Ning

A modified version of actuator cylinder theory to account for multiple vertical axis wind turbines.  Can run with just a single turbine.  See examples in example.jl.

Written in Julia, but an older Python version of the single actuator cylinder theory is also available in a branch.

See paper for theory details: Ning, A., "Actuator Cylinder Theory for Multiple Vertical Axis Wind Turbines," Wind Energy Science, Nov 2016, (accepted). doi:10.5194/wes-2016-19 

To install use (changes name of package and updates to module branch)

``` julia
Pkg.clone("https://github.com/byuflowlab/vawt-ac.git","vawtac")
Pkg.checkout("vawtac","module")
```

 To run the code, the fortran libraries will need to be compiled. On OSX this is accomplished by navigating to the vawtac package in your julia instalation and in the data/minpack diretory running the following code:
 
```
 gfortran -shared -O2 *.f -o libhybrd.dylib -fPIC
```
 
 You can test the package to ensure that this compilation worked.
