
module actuatorcylinder

    implicit none

    real(8), allocatable, dimension(:) :: theta
    real(8), allocatable, dimension(:, :) :: Am
!     real(8), allocatable, dimension(:, :) :: Rx, Ry
    real(8) :: pi = 3.141592654


    contains



    subroutine velocity(ntheta, ww, r, chord, twist, Uinf, Omega, rho, mu, &
        W, phi, alpha, Re)

        implicit none

        ! in
        integer, intent(in) :: ntheta
        real(8), dimension(2*ntheta), intent(in) :: ww
        real(8), intent(in) :: r, chord, twist, Uinf, Omega, rho, mu

        ! out
        real(8), dimension(ntheta), intent(out) :: W, phi, alpha, Re

        ! local
        real(8), dimension(ntheta) :: wx, wy, Vt, Vn


        wx = ww(:ntheta)
        wy = ww(ntheta+1:)

        ! velocity components and angle
        Vt = Omega*r + Uinf*(1 + wx)*cos(theta) + Uinf*wy*sin(theta)
        Vn = Uinf*(1 + wx)*sin(theta) - Uinf*wy*cos(theta)
        W = sqrt(Vt**2 + Vn**2)
        phi = atan2(Vn, Vt)
        alpha = phi - twist
        Re = rho*chord*Uinf/mu


    end subroutine velocity


    subroutine error(ntheta, ww, r, chord, Uinf, delta, B, W, phi, cl, cd, &
        modlin, cn, ct, error_array)

        implicit none

        ! in
        integer, intent(in) :: ntheta
        real(8), dimension(2*ntheta), intent(in) :: ww
        real(8), intent(in) :: r, chord, Uinf, delta
        real(8), dimension(ntheta), intent(in) :: W, phi, cl, cd
        integer, intent(in) :: B
        logical :: modlin

        ! out
        real(8), dimension(ntheta), intent(out) :: cn, ct
        real(8), dimension(2*ntheta), intent(out) :: error_array

        ! local
        real(8) :: sigma, CThrust, a, ka
        real(8), dimension(ntheta) :: integrand, Pn

        ! rotate force coefficients
        cn = cl*cos(phi) + cd*sin(phi)
        ct = cl*sin(phi) - cd*cos(phi)

        ! pressure
        sigma = B*chord/(2*r)
        Pn = sigma/(2*pi)*cn*(W/Uinf)**2

        ! correction factor based on thrust coefficient
        if (modlin) then
            integrand = (W/Uinf)**2 * (cn*sin(theta) - ct*cos(theta)/cos(delta))
            CThrust = sigma/(2*pi) * pInt(ntheta, integrand)
            if (CThrust > 2) then
                a = 0.5*(1 + sqrt(1 + CThrust))
                ka = 1.0 / (a-1)

            else if (CThrust > 0.96) then
                a = 1.0/7*(1 + 3*sqrt(7.0/2*CThrust - 3))
                ka = 18.0*a / (7*a**2 - 2*a + 4)

            else
                a = 0.5*(1 - sqrt(1 - CThrust))
                ka = 1 / (1-a)
            end if

        else
            ka = 1.0
        end if

        error_array = ka*matmul(Am, Pn) - ww

    end subroutine error


    ! integration for a periodic function
    function pInt(ntheta, f)
        implicit none

        integer, intent(in) :: ntheta
        real(8), dimension(ntheta), intent(in) :: f
        real(8) :: pInt, dtheta
        integer :: i

        pInt = 0.0
        do i = 1, ntheta-1
            pInt = pInt + (theta(i+1)-theta(i)) * (f(i+1) + f(i)) / 2.0
        end do

        ! add end points
        dtheta = 2*theta(1)  ! equally spaced, starts at 0
        pInt = pInt + dtheta * (f(1) + f(ntheta)) /2.0

    end function pInt


end module actuatorcylinder


