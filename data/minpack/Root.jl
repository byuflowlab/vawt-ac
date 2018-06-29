module Root
export hybrd
LP,_ = splitdir(@__FILE__)
# function residual_wrapper(nw::Ref{Cint}, xw::Ref{Cdouble}, fw::Ref{Cdouble}, iflag::Ref{Cint})
function residual_wrapper(nw::Int32, xw::Ptr{Cdouble}, fw::Ptr{Cdouble}, iflag::Int32)
    """
    Callback function from hybrd.  Unpacks the Fortran pointers, calls a
    residual function written in Julia, then stores results in another Fortran
    pointer.
    """

    # initialize
    # nvec = unsafe_load(nw)  # length of arrays
    xvec = zeros(nw)  # initialize copy

    # copy x values into Julia
    for i = 1:nw
        xvec[i] = unsafe_load(xw, i)
    end

    # call residual function
    # note that res and res_args are global b.c. hybrd does not provide a
    # hook to pass arbitrary data through.  same approach used in scipy.optimize
    fvec = res(xvec, res_args...)

    # copy f values into C pointer
    for i = 1:nw
        unsafe_store!(fw, fvec[i], i)
    end

    return
end


function hybrd(residual, x0::Array{Float64,1}, args,tol=1e-6)
    """
    Inputs
    ------
    residual : function
        function handle to a function of the form
        out = residual(x, arg1, arg2, ...)
        where out is an array that should equal all zeros
        x is the evaluation point and
        arg1, ... are other necessary parameters

    x0 : array
        an initial estimate of the solution vector

    args : tuple of other arguments that residual takes (arg1, arg2, ...)

    tol : float
        nonnegative input variable. termination occurs
        when the algorithm estimates that the relative error
        between x and the solution is at most tol.

    Outputs
    -------
    x : array of length(x0)
        the final estimate of the solution vector.

    f : array of length(x0)
        contains the functions evaluated at the output x.

    info : integer
        if the user has terminated execution, info is set to the (negative)
        value of iflag.otherwise, info is set as follows.

        info = 0   improper input parameters.

        info = 1   algorithm estimates that the relative error
                   between x and the solution is at most tol.

        info = 2   number of calls to fcn has reached or exceeded
                   200*(n+1).

        info = 3   tol is too small. no further improvement in
                   the approximate solution x is possible.

        info = 4   iteration is not making good progress.

    """

    # initialize
    info_in = Cint[0]  # need to set as C pointer to get data back
    x = copy(x0)  # copy so we don't modify input
    n = length(x)
    f = zeros(x)
    lwa = round(Int, (n*(3*n+13))/2)  # cast to int
    wa = zeros(lwa)

    # unfortunately using global vars is necessary because minpack does not
    # provide hooks to pass arbitrary pointers as many modern C/fortran
    # functions do. Note that scipy.optimize uses global variables for this same
    # purpose.

    global res = residual
    global res_args = args

    # define the callback function
    const res_func = cfunction(residual_wrapper, Void, (Ref{Cint}, Ptr{Cdouble},
        Ptr{Cdouble}, Ref{Cint}))

    # call hybrd.  must pass by reference for Fortran
    # compilation command I used (OS X with gfortran):
    # gfortran -shared -O2 *.f -o libhybrd.dylib -fPIC
    ccall( (:hybrd1_, LP*"../data/minpack/libhybrd"), Void, (Ptr{Void}, Ptr{Cint}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}),
        res_func, &n, x, f, &tol, info_in, wa, &lwa)
    info = info_in[1]  # extract termination info

    return x, f, info
end


function test()
    """
    just an example
    """

    function residual(xvec, arg1, arg2)

        fvec = zeros(xvec)

        fvec[1] = xvec[1] - arg1
        fvec[2] = xvec[1] + xvec[2] - arg2

        return fvec
    end

    x0 = [1.0, 0.0]
    args = (3.0, 1.0)
    tol = 1e-4
    x, f, info = hybrd(residual, x0, args, tol)
    println(x)
    println(f)
    println(info)
end


end
