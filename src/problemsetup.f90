module problemsetup
  use type_defs
  implicit none
  integer, parameter :: nvar = 1
  integer, parameter :: q = 9
  integer, parameter :: nint = q+3

  !number of elements along each direction
  integer, parameter :: nelemx = 10, nelemy = 10
  !the x and y coordinate of the source
  integer, parameter :: source_x = 4, source_y = 4
  real(kind = dp), parameter :: CFL = 0.001d0, tend = 1
  ! real(kind = dp), parameter :: tend = 1.d0
  real(kind = dp), parameter :: nu = 1.0_dp, strgth = 15_dp
  real(kind = dp) :: bc(10:99,nvar)

  integer, parameter :: nplot = 2, plot_freq = 100
  logical, parameter :: upwind = .true.
  logical, parameter :: plot = .true.

contains

  subroutine set_bc
    implicit none
    ! This routine is used to set boundary conditions
    ! on boundary curve xx
    !
    ! Default is Dirichlet for all boundaries.
    bc = 0.0_dp
    !Here we will eventually wish to implement radiating
    !boundary conditions
  end subroutine set_bc

  real(kind = dp) function init_u(x,y)
  !This subroutine returns the desired
  !initial condition/data for the problem
  !we wish to solve.
  !Inputs:
  !   - x  : x coordinate in physical space 
  !   - y  : y coordinate in physical space 
    use type_defs
    implicit none
    real(kind = dp) :: x,y
    real(kind = dp), parameter :: pi = acos(-1.d0)
    ! init_u = 0.0_dp
    ! init_u = EXP(-36.0_dp*(x**2.0_dp + y**2.0_dp))
    init_u = SIN(2.0_dp*pi*y)*SIN(2.0_dp*pi*x)
    ! init_u = (x**2.0_dp + y**2.0_dp)
    ! init_u = x**2.0_dp
    ! init_u = 5.0_dp*(y**2.0_dp)*(x**2.0_dp)
    ! init_u = 3.d0*x**2.d0*y**4.d0
    return
  end function init_u

  subroutine pis(xy,s,xy_start,xy_end,curve_type)
  !This subroutine computes the metric for a specified
  !geometry. If curvilinear elements are used, 
  !the Gordon-Hall mapping is used to create the metric.
  !Inputs:
  !   - s           : parameterization variable 
  !   - xy_start    : physical space starting coordinate 
  !   - xy_end      : physical space ending coordinate 
  !   - curve_type  : Integer specifying the curve type
  !                   (see below) 
  !Outputs:
  !   - xy          : x coordinate in physical space 


    use type_defs
    implicit none
    real(kind=dp) :: xy(2),xy_start(2),xy_end(2),s
    real(kind=dp) :: theta_s, theta_e
    integer :: curve_type
    if (curve_type .eq. 10 ) then
     ! Straight lines.
     xy = xy_start + 0.5d0*(1.d0 + s)*(xy_end - xy_start)
    elseif (curve_type .eq. 12) then
     ! Circle of radius 0.5 and center in (0,0).
     theta_s = atan2(xy_start(2),xy_start(1))
     theta_e = atan2(xy_end(2),xy_end(1))
     xy(1) = 0.5d0*cos(theta_s + 0.5d0*(1.d0 + s)*(theta_e - theta_s))
     xy(2) = 0.5d0*sin(theta_s + 0.5d0*(1.d0 + s)*(theta_e - theta_s))
    elseif (curve_type .eq. 14) then
     ! Circle of radius 1 and center in (0,0).
     theta_s = atan2(xy_start(2),xy_start(1))
     theta_e = atan2(xy_end(2),xy_end(1))
     xy(1) = cos(theta_s + 0.5d0*(1.d0 + s)*(theta_e - theta_s))
     xy(2) = sin(theta_s + 0.5d0*(1.d0 + s)*(theta_e - theta_s))
    elseif (curve_type .eq. 100 ) then
     ! Straight line for internal boundaries, used for the mapping of
     ! curved elements.
     xy = xy_start + 0.5d0*(1.d0 + s)*(xy_end - xy_start)
    end if
  end subroutine pis

end module problemsetup
