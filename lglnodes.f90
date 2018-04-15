subroutine lglnodes(x,w,n)
  !
  ! 
  ! lglnodes.m
  !
  ! Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde 
  ! matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
  ! integration and spectral methods. 
  !
  ! Reference on LGL nodes and weights:  
  !   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
  !   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
  !
  ! Written by Greg von Winckel - 04/17/2004
  ! Contact: gregvw@chtm.unm.edu
  !
  !
  use type_defs
  implicit none
  integer :: n,n1
  real(kind=dp) :: w(0:n),x(0:n),xold(0:n)
  integer :: i,k
  real(kind=dp) :: P(1:n+1,1:n+1),eps
  ! Truncation + 1
  N1=N+1
  eps = 2.2204d-16

  ! Use the Chebyshev-Gauss-Lobatto nodes as the first guess
  do i = 0,n
     x(i) = -cos(pi*dble(i)/dble(N))
  end do
  
  ! The Legendre Vandermonde Matrix
  !  P=zeros(N1,N1);
  
  ! Compute P_(N) using the recursion relation
  ! Compute its first and second derivatives and 
  ! update x using the Newton-Raphson method.
  
  xold = 2.d0
  
  !do while (maxval(abs(x-xold)) .gt. eps) 
  
  do i = 1,100   
     xold = x
     
     P(:,1) = 1.d0
     P(:,2) = x
     
     do  k=2,n
        P(:,k+1)=( dble(2*k-1)*x*P(:,k)-dble(k-1)*P(:,k-1) )/dble(k);
     end do
     x = xold-( x*P(:,N1)-P(:,N) )/( dble(N1)*P(:,N1) )
     if (maxval(abs(x-xold)).lt. eps ) exit
  end do
  
  w=2.d0/(dble(N*N1)*P(:,N1)**2)
  
end subroutine lglnodes
