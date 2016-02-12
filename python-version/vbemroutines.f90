! module spline_mod

! implicit none
! private

! public :: spline

! integer, parameter :: dp = selected_real_kind(15, 307)


! type spline
!     private

!     integer :: nx, ny, kx, ky
!     real(dp), dimension(:), allocatable :: tx, ty, c

!     contains
!     private
!     procedure :: init
!     procedure, public :: eval

! end type spline


! interface spline
! procedure constructor
! end interface


! contains


! function constructor(mx, x, my, y, z, s)
!     type(spline) :: constructor

!     integer :: mx, my
!     real(dp), dimension(mx) :: x
!     real(dp), dimension(my) :: y
!     real(dp), dimension(mx*my) :: z
!     real(dp) :: s

!     call constructor%init(mx, x, my, y, z, s)

! end function


! subroutine init(self, mx, x, my, y, z, s)  !, &
! !     kx, ky, nx, tx, ny, ty, c)

!     implicit none

!     ! inputs
!     class(spline), intent(inout) :: self
!     integer, intent(in) :: mx, my
!     real(dp), dimension(mx), intent(in) :: x
!     real(dp), dimension(my), intent(in) :: y
!     real(dp), dimension(mx*my), intent(in) :: z
!     real(dp), intent(in) :: s

! !     ! outputs
! !     integer, intent(out) :: nx, ny, kx, ky
! !     real(dp), dimension(:), allocatable, intent(out) :: tx, ty, c

!     ! local
!     integer :: iopt, nxest, nyest
!     real(dp) :: xb, xe, yb, ye
!     real(dp) :: fp
!     real(dp), dimension(:), allocatable :: wrk
!     integer, dimension(:), allocatable :: iwrk
!     integer :: lwrk, kwrk, ier



! !         alpha = (/0.0, 1.0, 2.0, 3.0, 4.0/)
! !         Re =  [1.0e6, 2.0e6]
! !         !     cl = (/ 0.0, 0.1, 0.2, 0.3, 0.4, 0.01, 0.11, 0.21, 0.31, 0.41/)
! !         cl = [0.0, 0.01, 0.1, 0.11, 0.2, 0.21, 0.3, 0.31, 0.4, 0.41]

!     ! use smoothing
!     iopt = 0

!     ! set bounds
!     xb = x(1)
!     xe = x(mx)
!     yb = y(1)
!     ye = y(my)

!     ! degree of spline
!     self%kx = min(3, mx-1)
!     self%ky = min(3, my-1)

!     ! max number of knots
!     nxest = mx + self%kx + 1
!     nyest = my + self%ky + 1
!     allocate(self%tx(2*mx))
!     allocate(self%ty(2*my))
!     allocate(self%c(mx*my))

!     ! work arrays
!     lwrk = 4 + nxest*(my+2*self%kx+5) + nyest*(2*self%ky+5) + &
!         mx*(self%kx+1)+my*(self%ky+1) + max(my, nxest)
!     kwrk = 3 + mx + my + nxest + nyest
!     allocate(wrk(lwrk))
!     allocate(iwrk(kwrk))


!     call regrid(iopt, mx, x, my, y, z, xb, xe, yb, ye, &
!         self%kx, self%ky, s, nxest, nyest, self%nx, self%tx, &
!         self%ny, self%ty, self%c, fp, wrk, lwrk, iwrk, kwrk, ier)


!     ! TODO: error handling

! end subroutine init




! subroutine eval(self, mx, x, my, y, z)  !, &
! !     kx, ky, nx, tx, ny, ty, c)

!     implicit none

!     ! inputs
!     class(spline), intent(inout) :: self
!     integer, intent(in) :: mx, my
!     real(dp), dimension(mx), intent(in) :: x
!     real(dp), dimension(my), intent(in) :: y
! !     integer, intent(in) :: nx, ny, kx, ky
! !     real(dp), dimension(nx), intent(in) :: tx
! !     real(dp), dimension(ny), intent(in) :: ty
! !     real(dp), dimension((nx-kx-1)*(ny-ky-1)), intent(in) :: c

!     ! outputs
!     real(dp), dimension(mx*my), intent(out) :: z

!     ! local
!     real(dp), dimension(:), allocatable :: wrk
!     integer, dimension(:), allocatable :: iwrk
!     integer :: lwrk, kwrk, ier

!     ! create work arrays
!     lwrk = mx*(self%kx+1)+my*(self%ky+1)
!     kwrk = mx+my
!     allocate(wrk(lwrk))
!     allocate(iwrk(kwrk))


!     ! evaluate
!     call bispev(self%tx, self%nx, self%ty, self%ny, self%c, &
!         self%kx, self%ky, x, mx, y, my, z, wrk, lwrk, iwrk, kwrk, ier)


! end subroutine eval


! end module spline_mod



! program vbem

!     use spline_mod
!     implicit none

!     integer, parameter :: dp = selected_real_kind(15, 307)



!     real(dp), dimension(5) :: alpha
!     real(dp), dimension(2) :: Re
!     real(dp), dimension(10) :: cl
!     real(dp) :: s

!     real(dp), dimension(3) :: alpha_pt
!     real(dp), dimension(1) :: Re_pt
!     real(dp), dimension(3) :: cl_pt

!     type(spline) :: sp

! !     integer :: iopt, mx, my, kx, ky, nxest, nyest, nx, ny
! !     real(dp) :: xb, xe, yb, ye, s
! !     real(dp), dimension(20) :: tx
! !     real(dp), dimension(20) :: ty
! !     real(dp), dimension(20) :: c
! !     real(dp) :: fp
! !     real(dp), dimension(1000) :: wrk
! !     integer :: lwrk
! !     integer, dimension(1000) :: iwrk
! !     integer :: kwrk, ier

! !     lwrk = 1000
! !     kwrk = 1000

!     alpha = (/0.0, 1.0, 2.0, 3.0, 4.0/)
!     Re =  [1.0e6, 2.0e6]
!     cl = [0.0, 0.01, 0.1, 0.11, 0.2, 0.21, 0.3, 0.31, 0.4, 0.41]
!     s = 0.1

!     sp = spline(5, alpha, 2, Re, cl, s)

! !     iopt = 0
! !     mx = 5
! !     my = 2
! !     kx = 3
! !     ky = 1
! !     xb = alpha(1)
! !     xe = alpha(5)
! !     yb = Re(1)
! !     ye = Re(2)
! !     s = 0.1
! !     nxest = mx + kx + 1
! !     nyest = my + ky + 1

! !     print *, alpha

! !     call regrid(iopt,mx,alpha,my,Re,cl,xb,xe,yb,ye,kx,ky,s, &
! !         nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)

! !     print *, ier

!     alpha_pt = [0.5, 1.5, 2.5]
!     Re_pt = [2.0e6]

! !     call bispev(tx,nx,ty,ny,c,kx,ky,alpha_pt,mx,Re_pt,my,cl_pt,wrk,lwrk, &
! !         iwrk,kwrk,ier)

!     call sp%eval(3, alpha_pt, 1, Re_pt, cl_pt)

!     print *, cl_pt


! end program vbem





subroutine velocity(a, r, chord, twist, theta, delta, &
    Uinf, Omega, rho, mu, alpha, Re, phi, W)

    implicit none

    ! in
    real(8), intent(in) :: a, r, chord, twist, theta, delta
    real(8), intent(in) :: Uinf, Omega, rho, mu

    ! out
    real(8), intent(out) :: alpha, Re, phi, W

    ! local
    real(8) :: V, Vt, Vn

    ! apply induction factor
    V = Uinf*(1-a)

    ! angle of attack, local velocity
    Vn = V*sin(theta)*cos(delta)
    Vt = Omega*r + V*cos(theta)

    phi = atan2(Vn, Vt)
    W = sqrt(Vn**2 + Vt**2)

    alpha = phi - twist
    Re = rho*W*chord/mu

end subroutine velocity



subroutine vbem(a, r, chord, theta, delta, &
    cl, cd, phi, W, Uinf, rho, B, residual, Np, Tp, Zp)

    implicit none

    ! in
    real(8), intent(in) :: a, r, chord, theta, delta
    real(8), intent(in) :: cl, cd, phi, Uinf, W, rho
    integer, intent(in) :: B

    ! out
    real(8), intent(out) :: residual, Np, Tp, Zp

    ! local
    real(8) :: cn, ct, force, CT_bem, CT_mom, q, pi

    pi = 3.141592654

    ! normal and tangential forces
    cn = cl*cos(phi) + cd*sin(phi)
    ct = cl*sin(phi) - cd*cos(phi)

    ! forces
    q = 0.5*rho*W**2
    Np = -cn * q * chord
    Tp = ct * q * chord / cos(delta)
    Zp = -cn * q * chord * tan(delta)

    ! residual
    force = (cn*sin(theta) - ct*cos(theta)/cos(delta)) / abs(sin(theta))
    CT_bem = B*chord/(2*pi*r)*(W/Uinf)**2*force

    if (a < 0.4) then
        CT_mom = 4*a*(1-a)
    else ! Glauert
        CT_mom = 2.0/9*(7*a**2 - 2*a + 4)
    end if

    residual = CT_mom - CT_bem



end subroutine vbem