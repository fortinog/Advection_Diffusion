module problem_burger_2d
  !
  ! This sets parameters for the heat equation in 1D
  !
  use type_defs
  implicit none
  integer, parameter :: nvar = 1 ! No. variables
  integer, parameter :: q = 5  ! Degree of polynomials used
  integer, parameter :: nref = 20
  integer, parameter :: nelemx = 10 ! Number of elements
  integer, parameter :: nelemy = 10 ! Number of elements
  real(kind = dp), parameter :: xl = -pi
  real(kind = dp), parameter :: xr =  pi
  real(kind = dp), parameter :: yb = -pi
  real(kind = dp), parameter :: yt =  pi
  integer, parameter :: nplot = 20
  integer, parameter :: nint = q+2
  real(kind=dp), parameter :: CFL = 0.1_dp
  real(kind=dp), parameter :: nu  = 0.01_dp
  real(kind=dp), parameter :: tend = 10.0_dp
  real(kind=dp), parameter :: gamma = 1.4_dp
  integer, parameter :: npl_mod = 10

contains

  real(dp) function init_data(x,y,ivar)
    use type_defs
    implicit none
    real(dp) :: x,y,urand
    integer :: ivar
    call random_number(urand)

    init_data = urand-0.5_dp
    
  end function init_data

end module problem_burger_2d
