module legfuns
  
contains

  real(dp) function legendre(x,n)
    use type_defs
    implicit none
    integer :: n,k
    real(dp) :: x,p(0:n)
    p(0) = 1.0_dp
    if(n .gt. 0) then
     p(1) = x
     ! k*P_k(x)=(2*k-1)*x*P_{k-1}(x)-(k-1)*P_{k-2}(x)
     do k=2,n
      p(k)=(real(2*k-1,dp)*x*p(k-1)-real(k-1,dp)*p(k-2) )/real(k,dp)
     end do
    end if
    legendre = p(n)
    return
  end function legendre
  
  real(dp) function legendre_der(x,n)
    use type_defs
    implicit none
    integer :: n,k
    real(dp) :: x,p(0:n),px(0:n)

    p(0)  = 1.0_dp
    px(0) = 0.0_dp
    if(n .gt. 0) then
     p(1) = x
     ! k*P_k(x)=(2*k-1)*x*P_{k-1}(x)-(k-1)*P_{k-2}(x)
     do k=2,n
      p(k)=(real(2*k-1,dp)*x*p(k-1)-real(k-1,dp)*p(k-2) )/real(k,dp)
     end do
     px(1)=1.0_dp
     do k=2,n
      px(k) = px(k-2)+real(2*k-1,dp)*p(k-1)
     end do
    end if
    legendre_der = px(n)
    return
  end function legendre_der
  
end module legfuns
