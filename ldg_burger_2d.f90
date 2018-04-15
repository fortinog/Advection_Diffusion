program burger2d
  use type_defs
  use problem_burger_2d
  use legfuns
  implicit none
  real(dp) :: u(0:q,0:q,nvar,nelemx,nelemy),lap(0:q,0:q,nvar,nelemx,nelemy),&
    f_x(0:q,0:q,nvar,nelemx,nelemy),g_y(0:q,0:q,nvar,nelemx,nelemy)
  real(dp) :: up(0:q,0:q,nvar,nelemx,nelemy),kstage(0:q,0:q,nvar,nelemx,nelemy,4)
  real(dp) :: uexact(0:q,0:q,nvar,nelemx,nelemy)
  real(dp) :: xm(nelemx), x(0:nelemx), ym(nelemy), y(0:nelemy)
  real(dp) :: hx,hy,r_x,x_r,r_y,y_r,utmp,dt,lam_max,t
  integer  :: i,j,k,l,it,i1,i2,ivar,row,col,istart,jstart
  real(dp) :: weights(nint), znodes(nint), leg_mat(nint,0:q),&
    leg_der_mat(nint,0:q)
  real(dp) :: mass_mat_x(0:q), mass_mat_y(0:q), xloc(nint), yloc(nint), uloc(nint,nint)
  character(100) :: str
  real(dp) :: diff_mat_x((q+1)**2,(q+1)**2),diff_mat_y((q+1)**2,(q+1)**2)
  real(dp) :: xplot(0:nplot),cof_2_plot(0:nplot,0:q)
  !
  !
  ! Grid and metric
  hx = (xr-xl)/real(nelemx,dp)
  hy = (yt-yb)/real(nelemy,dp)
  r_x = 2.0_dp/hx
  x_r = hx/2.0_dp
  r_y = 2.0_dp/hy
  y_r = hy/2.0_dp
  do i = 0,nelemx
   x(i) = xl + real(i,dp)*hx
  end do
  xm = 0.5_dp*(x(1:nelemx)+x(0:nelemx-1))
  do i = 0,nelemy
   y(i) = yb + real(i,dp)*hy
  end do
  ym = 0.5_dp*(y(1:nelemy)+y(0:nelemy-1))
  do i = 0,nplot
   xplot(i) = -1.0_dp + 2.0_dp*real(i,dp)/real(nplot,dp)
  end do
  ! for plotting
  do k = 0,q
   do l = 0,nplot
    cof_2_plot(l,k) = legendre(xplot(l),k)
   end do
  end do
  ! For quadrature
  i = nint-1
  call lglnodes(znodes,weights,i)
  do k = 0,q
   do l = 1,nint
    leg_mat(l,k) = legendre(znodes(l),k)
    leg_der_mat(l,k) = legendre_der(znodes(l),k)
   end do
  end do
  do k = 0,q
   mass_mat_x(k) = sum(weights*leg_mat(1:nint,k)**2)
   mass_mat_y(k) = sum(weights*leg_mat(1:nint,k)**2)
  end do
  ! Project initial data
  do ivar = 1,nvar
   do j = 1,nelemy
    yloc = ym(j) + znodes*y_r
    do i = 1,nelemx
     xloc = xm(i) + znodes*x_r
     do k = 1,nint
      do l = 1,nint
       uloc(k,l) = init_data(xloc(k),yloc(l),ivar)
      end do
     end do
     do l = 0,q
      do k = 0,q
       u(k,l,ivar,i,j) = 0.0_dp
       do i2 = 1,nint
        do i1 = 1,nint
         u(k,l,ivar,i,j) = u(k,l,ivar,i,j) &
           +weights(i1)*weights(i2)&
           *leg_mat(i1,k)*leg_mat(i2,l)&
           *uloc(i1,i2)
        end do
       end do
       u(k,l,ivar,i,j) = u(k,l,ivar,i,j)/(mass_mat_x(k)*mass_mat_y(l))
      end do
     end do
    end do
   end do
  end do
  uexact = u
  ! rescale
  mass_mat_x = mass_mat_x*x_r
  mass_mat_y = mass_mat_y*y_r
  !
  diff_mat_x = 0.0_dp
  diff_mat_y = 0.0_dp
  do l =0,q
   do k = 0,q
    row = 1 + k + l*(q+1)
    do j =0,q
     do i = 0,q
      col = 1 + i + j*(q+1)
      diff_mat_x(row,col) = 0.0_dp
      diff_mat_y(row,col) = 0.0_dp
      do i2 = 1,nint
       do i1 = 1,nint
        diff_mat_x(row,col) = diff_mat_x(row,col) &
          +x_r*r_x*y_r*weights(i1)*weights(i2)*(&
          leg_der_mat(i1,k)*leg_mat(i2,l)*&
          leg_mat(i1,i)*leg_mat(i2,j))
        diff_mat_y(row,col) = diff_mat_y(row,col) &
          +y_r*r_y*x_r*weights(i1)*weights(i2)*(&
          leg_mat(i1,k)*leg_der_mat(i2,l)*&
          leg_mat(i1,i)*leg_mat(i2,j))
       end do
      end do
     end do
    end do
   end do
  end do
  !
  !
  it = 0
  do ivar = 1,nvar
   WRITE(str,'("sol",I2.2,"_",I8.8,".txt")') ivar,it
   OPEN(unit=29,file=trim(str),status="REPLACE")
   ! Sweep across each x-line one
   do j = 1,nelemy
    jstart = 1
    if (j.eq.1) jstart = 0
    do l = jstart,nplot
     do i = 1,nelemx
      istart = 1
      if (i.eq.1) istart = 0
      do k = istart,nplot
       utmp = 0.0_dp
       do i2 = 0,q
        do i1 = 0,q
         utmp = utmp + u(i1,i2,ivar,i,j)&
           *cof_2_plot(k,i1)*cof_2_plot(l,i2)
        end do
       end do
       write(29,*) xm(i)+xplot(k)*x_r, ym(j)+xplot(l)*y_r, utmp
      end do
     end do
    end do
   end do
   close(29)
  end do
  lap = 0.0_dp
  t = 0.0_dp
  it = 0
  do while  (t < tend)
   up = u
   call compute_advection(f_x,g_y,u,hx,hy,nint,mass_mat_x,mass_mat_y,lam_max)
   dt = CFL*min(hx,hy)/real(q,dp)**2/lam_max
   dt = min(dt,tend-dt)
   if(nu.gt.0.0_dp) then
    call compute_laplacian(lap,u,x_r,y_r,weights,diff_mat_x,diff_mat_y,&
      mass_mat_x,mass_mat_y,leg_mat)
   end if
   kstage(:,:,:,:,:,1) = f_x+g_y+nu*lap
   u = up + 0.5_dp*dt*kstage(:,:,:,:,:,1)
   call compute_advection(f_x,g_y,u,hx,hy,nint,mass_mat_x,mass_mat_y,lam_max)
   if(nu.gt.0.0_dp) then
    call compute_laplacian(lap,u,x_r,y_r,weights,diff_mat_x,diff_mat_y,&
      mass_mat_x,mass_mat_y,leg_mat)
   end if
   kstage(:,:,:,:,:,2) = f_x+g_y+nu*lap
   u = up + 0.5_dp*dt*kstage(:,:,:,:,:,2)
   call compute_advection(f_x,g_y,u,hx,hy,nint,mass_mat_x,mass_mat_y,lam_max)
   if(nu.gt.0.0_dp) then
    call compute_laplacian(lap,u,x_r,y_r,weights,diff_mat_x,diff_mat_y,&
      mass_mat_x,mass_mat_y,leg_mat)
   end if
   kstage(:,:,:,:,:,3) = f_x+g_y+nu*lap
   u = up + dt*kstage(:,:,:,:,:,3)
   call compute_advection(f_x,g_y,u,hx,hy,nint,mass_mat_x,mass_mat_y,lam_max)
   if(nu.gt.0.0_dp) then
    call compute_laplacian(lap,u,x_r,y_r,weights,diff_mat_x,diff_mat_y,&
      mass_mat_x,mass_mat_y,leg_mat)
   end if
   kstage(:,:,:,:,:,4) = f_x+g_y+nu*lap
   u = up + (dt/6.0_dp)&
     *(kstage(:,:,:,:,:,1)+2.0_dp*kstage(:,:,:,:,:,2)&
     +2.0_dp*kstage(:,:,:,:,:,3)+kstage(:,:,:,:,:,4))
   !
   t = t + dt
   it = it + 1
   if (mod(it,npl_mod) .eq. 0) then
    write(*,*) "Saving at time: ", t
    do ivar = 1,nvar
     WRITE(str,'("sol",I2.2,"_",I8.8,".txt")') ivar,it
     OPEN(unit=29,file=trim(str),status="REPLACE")
     ! Sweep across each x-line one
     do j = 1,nelemy
      jstart = 1
      if (j.eq.1) jstart = 0
      do l = jstart,nplot
       do i = 1,nelemx
        istart = 1
        if (i.eq.1) istart = 0
        do k = istart,nplot
         utmp = 0.0_dp
         do i2 = 0,q
          do i1 = 0,q
           utmp = utmp + u(i1,i2,ivar,i,j)&
             *cof_2_plot(k,i1)*cof_2_plot(l,i2)
          end do
         end do
         write(29,*) xm(i)+xplot(k)*x_r, ym(j)+xplot(l)*y_r, utmp
        end do
       end do
      end do
     end do
     close(29)
    end do
   end if
  end do
  

end program burger2d

subroutine compute_advection(f_x,g_y,u,hx,hy,nint_ps,mass_mat_x,mass_mat_y,alpha)
  !
  use type_defs
  use problem_burger_2d
  use legfuns
  implicit none
  integer :: nint_ps
  real(dp) :: mass_mat_x(0:q), mass_mat_y(0:q)
  real(dp) :: u(0:q,0:q,nvar,nelemx,nelemy),f_x(0:q,0:q,nvar,nelemx,nelemy),&
    g_y(0:q,0:q,nvar,nelemx,nelemy)
  real(dp) :: hx,hy
  ! Naive implementation, convert to ps-nodal values then integrate
  real(dp) :: u_ps(nint_ps,nint_ps,nvar,nelemx,nelemy)
  real(dp) :: f_ps(nint_ps,nint_ps,nvar)
  real(dp) :: ub(nint_ps,4,nvar,nelemx,nelemy),fb_star(nint_ps,4,nvar,nelemx,nelemy),&
    fb(nint_ps,4,nvar,nelemx,nelemy)
  real(dp) :: x_r,y_r,r_x,r_y,mode_c,alpha
  integer  :: i,j,k,l,i1,i2,ivar
  real(dp) :: weights(nint_ps),znodes(nint_ps),leg_mat(nint_ps,0:q),&
    leg_der_mat(nint_ps,0:q)

  real(dp) :: ub_plus(nint_ps,nvar), fb_plus(nint_ps,nvar)
  i = nint_ps-1

  call lglnodes(znodes,weights,i)
  do k = 0,q
   do l = 1,nint_ps
    leg_mat(l,k) = legendre(znodes(l),k)
    leg_der_mat(l,k) = legendre_der(znodes(l),k)
   end do
  end do

  r_x = 2.0_dp/hx
  x_r = hx/2.0_dp
  r_y = 2.0_dp/hy
  y_r = hy/2.0_dp

  !
  ! Build u on quadrature points
  !
  do j = 1,nelemy
   do i = 1,nelemx
    do ivar = 1,nvar
     u_ps(:,:,ivar,i,j) = 0.0_dp
     do i2 = 0,q
      do i1 = 0,q
       mode_c = u(i1,i2,ivar,i,j)
       do l = 1,nint_ps
        do k = 1,nint_ps
         u_ps(k,l,ivar,i,j) = u_ps(k,l,ivar,i,j) &
           + mode_c*leg_mat(k,i1)*leg_mat(l,i2)
        end do
       end do
      end do
     end do
    end do
   end do
  end do

  alpha = maxval(abs(u_ps(:,:,1,:,:)))

  do j = 1,nelemy
   do i = 1,nelemx
    do ivar = 1,nvar
     do l = 1,nint_ps
      ub(l,1,ivar,i,j) = u_ps(1,l,ivar,i,j)
      ub(l,2,ivar,i,j) = u_ps(nint_ps,l,ivar,i,j)
      ub(l,3,ivar,i,j) = u_ps(l,1,ivar,i,j)
      ub(l,4,ivar,i,j) = u_ps(l,nint_ps,ivar,i,j)
     end do
    end do
   end do
  end do
  ! 
  fb = 0.5_dp*ub**2
  !
  f_x = 0.0_dp
  g_y = 0.0_dp
  !
  ! periodic bc ...
  do ivar = 1,nvar
   do j = 1,nelemy
    i = 1
    fb_star(:,1,ivar,i,j) = 0.5_dp*(fb(:,1,ivar,i,j) + fb(:,2,ivar,nelemx,j)) &
      - 0.5_dp*alpha*(ub(:,1,ivar,i,j)-ub(:,2,ivar,nelemx,j))
    fb_star(:,2,ivar,i,j) = 0.5_dp*(fb(:,2,ivar,i,j) + fb(:,1,ivar,i+1,j)) &
      + 0.5_dp*alpha*(ub(:,2,ivar,i,j)-ub(:,1,ivar,i+1,j))
    do i = 2,nelemx-1
     fb_star(:,1,ivar,i,j) = 0.5_dp*(fb(:,1,ivar,i,j) + fb(:,2,ivar,i-1,j)) &
       - 0.5_dp*alpha*(ub(:,1,ivar,i,j)-ub(:,2,ivar,i-1,j))
     fb_star(:,2,ivar,i,j) = 0.5_dp*(fb(:,2,ivar,i,j) + fb(:,1,ivar,i+1,j)) &
       + 0.5_dp*alpha*(ub(:,2,ivar,i,j)-ub(:,1,ivar,i+1,j))
    end do
    i = nelemx
    fb_star(:,1,ivar,i,j) = 0.5_dp*(fb(:,1,ivar,i,j) + fb(:,2,ivar,i-1,j)) &
      - 0.5_dp*alpha*(ub(:,1,ivar,i,j)-ub(:,2,ivar,i-1,j))
    fb_star(:,2,ivar,i,j) = 0.5_dp*(fb(:,2,ivar,i,j) + fb(:,1,ivar,1,j)) &
      + 0.5_dp*alpha*(ub(:,2,ivar,i,j)-ub(:,1,ivar,1,j))
   end do
  end do

  do i = 1,nelemx
   do ivar = 1,nvar
    ! Inner faces
    j = 1
    fb_star(:,3,ivar,i,j) = 0.5_dp*(fb(:,3,ivar,i,j) + fb(:,4,ivar,i,nelemy)) &
      - 0.5_dp*alpha*(ub(:,3,ivar,i,j)-ub(:,4,ivar,i,nelemy))
    fb_star(:,4,ivar,i,j) = 0.5_dp*(fb(:,4,ivar,i,j) + fb(:,3,ivar,i,j+1)) &
      + 0.5_dp*alpha*(ub(:,4,ivar,i,j)-ub(:,3,ivar,i,j+1))
    do j = 2,nelemy-1
     fb_star(:,3,ivar,i,j) = 0.5_dp*(fb(:,3,ivar,i,j) + fb(:,4,ivar,i,j-1)) &
       - 0.5_dp*alpha*(ub(:,3,ivar,i,j)-ub(:,4,ivar,i,j-1))
     fb_star(:,4,ivar,i,j) = 0.5_dp*(fb(:,4,ivar,i,j) + fb(:,3,ivar,i,j+1)) &
       + 0.5_dp*alpha*(ub(:,4,ivar,i,j)-ub(:,3,ivar,i,j+1))
    end do
    j = nelemy
    fb_star(:,3,ivar,i,j) = 0.5_dp*(fb(:,3,ivar,i,j) + fb(:,4,ivar,i,j-1)) &
      - 0.5_dp*alpha*(ub(:,3,ivar,i,j)-ub(:,4,ivar,i,j-1))
   end do
   ! Bottom
   j = 1
   ub_plus(:,1) = 0.0_dp
   fb_plus(:,1) = 0.0_dp
   do ivar = 1,nvar
    fb_star(:,3,ivar,i,j) = 0.5_dp*(fb(:,3,ivar,i,j) + fb_plus(:,ivar)) &
      - 0.5_dp*alpha*(ub(:,3,ivar,i,j)-ub_plus(:,ivar))
   end do

   j = nelemy
   ub_plus(:,1) = 0.0_dp
   fb_plus(:,1) = 0.0_dp
   do ivar = 1,nvar
    fb_star(:,4,ivar,i,j) = 0.5_dp*(fb(:,4,ivar,i,j) + fb_plus(:,ivar)) &
      + 0.5_dp*alpha*(ub(:,4,ivar,i,j)-ub_plus(:,ivar))
   end do
  end do
  !
  do j = 1,nelemy
   do i = 1,nelemx
    do ivar = 1,nvar
     ! First add flux integrals in y along r = +/- 1
     do l =0,q
      do k = 0,q
       f_x(k,l,ivar,i,j) = &
         - y_r*sum(weights*fb_star(1:nint_ps,2,ivar,i,j)&
         *leg_mat(nint_ps,k)*leg_mat(1:nint_ps,l)) &
         + y_r*sum(weights*fb_star(1:nint_ps,1,ivar,i,j)&
         *leg_mat(1,k)*leg_mat(1:nint_ps,l))
      end do
     end do
    end do
   end do
  end do
  !
  ! Integrate against all test functions
  do j = 1,nelemy
   do i = 1,nelemx
    f_ps(:,:,1) = 0.5_dp*u_ps(:,:,1,i,j)**2
    !
    do ivar = 1,nvar
     do l =0,q
      do k = 0,q
       do i2 = 1,nint_ps
        do i1 = 1,nint_ps
         f_x(k,l,ivar,i,j) = f_x(k,l,ivar,i,j) &
           + x_r*r_x*y_r*weights(i1)*weights(i2)&
           *f_ps(i1,i2,ivar)*leg_der_mat(i1,k)*leg_mat(i2,l)
        end do
       end do
       f_x(k,l,ivar,i,j) = f_x(k,l,ivar,i,j)/(mass_mat_x(k)*mass_mat_y(l))
      end do
     end do
    end do
   end do
  end do
  !
  !
  !
  do j = 1,nelemy
   do i = 1,nelemx
    do ivar = 1,nvar
     ! First add flux integrals in y along r = +/- 1
     do l =0,q
      do k = 0,q
       g_y(k,l,ivar,i,j) = &
         - x_r*sum(weights*fb_star(1:nint_ps,4,ivar,i,j)&
         *leg_mat(1:nint_ps,k)*leg_mat(nint_ps,l)) &
         + x_r*sum(weights*fb_star(1:nint_ps,3,ivar,i,j)&
         *leg_mat(1:nint_ps,k)*leg_mat(1,l))
      end do
     end do
    end do
   end do
  end do
  ! Integrate against all test functions
  do j = 1,nelemy
   do i = 1,nelemx
    f_ps(:,:,1) = 0.5_dp*u_ps(:,:,1,i,j)**2
    do ivar = 1,nvar
     do l =0,q
      do k = 0,q
       do i2 = 1,nint_ps
        do i1 = 1,nint_ps
         g_y(k,l,ivar,i,j) = g_y(k,l,ivar,i,j) &
           + y_r*r_y*x_r*weights(i1)*weights(i2)&
           *f_ps(i1,i2,ivar)*leg_mat(i1,k)*leg_der_mat(i2,l)
        end do
       end do
       g_y(k,l,ivar,i,j) = g_y(k,l,ivar,i,j)/(mass_mat_x(k)*mass_mat_y(l))
      end do
     end do
    end do
   end do
  end do

end subroutine compute_advection

subroutine compute_laplacian(lap,u,x_r,y_r,weights,diff_mat_x,diff_mat_y,&
  mass_mat_x,mass_mat_y,leg_mat)
  !
  use type_defs
  use problem_burger_2d
  implicit none
  real(dp) :: u(0:q,0:q,nvar,nelemx,nelemy),gradx(0:q,0:q,nvar,nelemx,nelemy),&
    grady(0:q,0:q,nvar,nelemx,nelemy),lap(0:q,0:q,nvar,nelemx,nelemy)
  real(dp) :: ub(nint,4,nvar,nelemx,nelemy),gradb(nint,4,nvar,nelemx,nelemy)
  real(dp) :: x_r,y_r
  integer  :: i,j,k,l,i1,i2,sz1,sz2,ivar
  real(dp) :: weights(nint), leg_mat(nint,0:q)
  real(dp) :: mass_mat_x(0:q), mass_mat_y(0:q)
  real(dp) :: diff_mat_x((q+1)**2,(q+1)**2),diff_mat_y((q+1)**2,(q+1)**2)
  !
  ! Build u on the boundary 1,2,3,4 are x = -1, x = 1, y = -1, y = 1.
  !
  do ivar = 1,nvar
   do j = 1,nelemy
    do i = 1,nelemx
     ub(:,:,ivar,i,j) = 0.0_dp
     do i2 = 0,q
      do i1 = 0,q
       ub(:,1,ivar,i,j) = ub(:,1,ivar,i,j) + u(i1,i2,ivar,i,j)*leg_mat(1,i1)*leg_mat(:,i2)
       ub(:,2,ivar,i,j) = ub(:,2,ivar,i,j) + u(i1,i2,ivar,i,j)*leg_mat(nint,i1)*leg_mat(:,i2)
       ub(:,3,ivar,i,j) = ub(:,3,ivar,i,j) + u(i1,i2,ivar,i,j)*leg_mat(:,i1)*leg_mat(1,i2)
       ub(:,4,ivar,i,j) = ub(:,4,ivar,i,j) + u(i1,i2,ivar,i,j)*leg_mat(:,i1)*leg_mat(nint,i2)
      end do
     end do
    end do
   end do
  end do
  ! Periodic boundary conditions in the x-direction
  do j = 1,nelemy
   ! Use the left side flux
   do i = 1,nelemx-1
    ub(:,2,:,i,j) = ub(:,1,:,i+1,j)
   end do
   ub(:,2,:,nelemx,j) = ub(:,1,:,1,j)
  end do
  ! Dirichlet boundary conditions in the y-direction
  do i = 1,nelemx
   do j = 1,nelemy-1
    ub(:,4,:,i,j) = ub(:,3,:,i,j+1)
   end do
   ub(:,4,:,i,nelemy) =  0.0_dp
   ub(:,3,:,i,1)      =  0.0_dp
  end do
  !
  sz1 = 1
  sz2 = (q+1)**2
  do j = 1,nelemy
   do i = 1,nelemx
    do ivar = 1,nvar
     ! First add flux integrals
     gradx(0:q,0:q,ivar,i,j) = 0.0_dp
     do l =0,q
      do k = 0,q
       gradx(k,l,ivar,i,j) = &
         +y_r*sum(weights*ub(1:nint,2,ivar,i,j)*leg_mat(nint,k)*leg_mat(1:nint,l)) &
         -y_r*sum(weights*ub(1:nint,1,ivar,i,j)*leg_mat(1,k)*leg_mat(1:nint,l))
      end do
     end do
     ! Then add the derivative in x
     call DGEMV('N',sz2,sz2,-1.d0,diff_mat_x,sz2,&
       u(0:q,0:q,ivar,i,j),sz1,1.d0,gradx(0:q,0:q,ivar,i,j),sz1)
     ! Solve
     do l =0,q
      do k = 0,q
       gradx(k,l,ivar,i,j) = gradx(k,l,ivar,i,j)/(mass_mat_x(k)*mass_mat_y(l))
      end do
     end do
     ! First add flux integrals
     grady(0:q,0:q,ivar,i,j) = 0.0_dp
     do l =0,q
      do k = 0,q
       grady(k,l,ivar,i,j) =&
         +x_r*sum(weights*ub(1:nint,4,ivar,i,j)*leg_mat(1:nint,k)*leg_mat(nint,l)) &
         -x_r*sum(weights*ub(1:nint,3,ivar,i,j)*leg_mat(1:nint,k)*leg_mat(1,l))
      end do
     end do
     ! Then add the derivative in x
     call DGEMV('N',sz2,sz2,-1.d0,diff_mat_y,sz2,&
       u(0:q,0:q,ivar,i,j),sz1,1.d0,grady(0:q,0:q,ivar,i,j),sz1)
     ! Solve
     do l =0,q
      do k = 0,q
       grady(k,l,ivar,i,j) = grady(k,l,ivar,i,j)/(mass_mat_x(k)*mass_mat_y(l))
      end do
     end do
    end do
   end do
  end do
  ! Now compute u_xx+u_yy
  ! Build u on the boundary 1,2,3,4 are x = -1, x = 1, y = -1, y = 1.
  do j = 1,nelemy
   do i = 1,nelemx
    do ivar = 1,nvar
     gradb(:,:,ivar,i,j) = 0.0_dp
     do i2 = 0,q
      do i1 = 0,q
       gradb(:,1,ivar,i,j) = gradb(:,1,ivar,i,j) &
         + gradx(i1,i2,ivar,i,j)*leg_mat(1,i1)*leg_mat(:,i2)
       gradb(:,2,ivar,i,j) = gradb(:,2,ivar,i,j) &
         + gradx(i1,i2,ivar,i,j)*leg_mat(nint,i1)*leg_mat(:,i2)
       gradb(:,3,ivar,i,j) = gradb(:,3,ivar,i,j) &
         + grady(i1,i2,ivar,i,j)*leg_mat(:,i1)*leg_mat(1,i2)
       gradb(:,4,ivar,i,j) = gradb(:,4,ivar,i,j) &
         + grady(i1,i2,ivar,i,j)*leg_mat(:,i1)*leg_mat(nint,i2)
      end do
     end do
    end do
   end do
  end do
  ! Periodic in x
  do j = 1,nelemy
   do i = 2,nelemx
    gradb(:,1,:,i,j) = gradb(:,2,:,i-1,j)
   end do
  end do
  gradb(:,1,:,1,:) = gradb(:,2,:,nelemx,:)
  !
  do i = 1,nelemx
   do j = 2,nelemy
    gradb(:,3,:,i,j) = gradb(:,4,:,i,j-1)
   end do
  end do
  !gradb(:,3,:,:,1) = gradb(:,4,:,:,nelemy)
  !gradb(:,3,:,:,1) = gradb(:,4,:,:,nelemy)
  !
  sz1 = 1
  sz2 = (q+1)**2
  do j = 1,nelemy
   do i = 1,nelemx
    do ivar = 1,nvar
     ! First add flux integrals
     lap(0:q,0:q,ivar,i,j) = 0.0_dp
     do l =0,q
      do k = 0,q
       lap(k,l,ivar,i,j) = &
         + y_r*sum(weights*gradb(1:nint,2,ivar,i,j)*leg_mat(nint,k)*leg_mat(1:nint,l)) &
         - y_r*sum(weights*gradb(1:nint,1,ivar,i,j)*leg_mat(1,k)*leg_mat(1:nint,l)) &
         + x_r*sum(weights*gradb(1:nint,4,ivar,i,j)*leg_mat(1:nint,k)*leg_mat(nint,l)) &
         - x_r*sum(weights*gradb(1:nint,3,ivar,i,j)*leg_mat(1:nint,k)*leg_mat(1,l))
      end do
     end do
     ! Then add the derivative in x
     call DGEMV('N',sz2,sz2,-1.d0,diff_mat_x,sz2,&
       gradx(0:q,0:q,ivar,i,j),sz1,1.d0,lap(0:q,0:q,ivar,i,j),sz1)
     ! Then add the derivative in x
     call DGEMV('N',sz2,sz2,-1.d0,diff_mat_y,sz2,&
       grady(0:q,0:q,ivar,i,j),sz1,1.d0,lap(0:q,0:q,ivar,i,j),sz1)
     ! Solve
     do l =0,q
      do k = 0,q
       lap(k,l,ivar,i,j) = lap(k,l,ivar,i,j)/(mass_mat_x(k)*mass_mat_y(l))
      end do
     end do
    end do
   end do
  end do

end subroutine compute_laplacian

