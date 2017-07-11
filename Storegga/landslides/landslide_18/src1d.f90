! This routine should be a simplified version of src2
! which applies source terms for a 1-d slice of data along the
! edge of a grid.  This is called only from qad where the conservative
! fix-up is applied and is used to apply source terms over partial
! time steps to the coarse grid cell values used in solving Riemann
! problems at the interface between coarse and fine grids.
subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)
      
    use geoclaw_module, only: grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: omega, coordinate_system, manning_coefficient
    use geoclaw_module, only: manning_break, num_manning
    use amr_module, only: xupper, yupper, xlower, ylower
    use bing_module                           

    implicit none

    ! Input
    integer, intent(in) :: meqn, mbc, mx1d, maux
    real(kind=8), intent(in) :: t, dt
    real(kind=8), intent(inout) :: q1d(meqn, mx1d), aux1d(maux, mx1d)

    ! Local storage
    integer :: i, nman
    logical :: found
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    ! Bingham
    real(kind=8), parameter:: hs_min = 0.0d0
    real(kind=8) :: u,v,up,vp,dp,ds,upbar,gxmod,gymod,g,cDx,cDy,ubar,usum
    real(kind=8) :: bx_fric(2),by_fric(2),nu_m,drag_force,cF,cP,tau_y
    real(kind=8) :: hx,hy,x

    g = grav*(1.d0-rho_a/rho_s)
	
    call checkval_bing(mx1d,1,meqn,mbc,q1d,alpha_1)

    hx=(xupper-xlower)/(mx1d+1)

    if (hydrodrag) then

       cF = cF_hyd
       dP = cP_hyd

       do i=1,mx1d

             h=q1d(1,i)

             if (h>depth_tolerance) then

                up=q1d(2,i)
                vp=q1d(3,i)
                u =q1d(4,i)/h
                v =q1d(5,i)/h
                gxmod = aux1d(4,i)
                gymod = aux1d(5,i)

                ubar=dsqrt(u**2+v**2)

                upbar = dsqrt(up**2+vp**2)

                if (upbar>1d-4) then
                   dp=(dsqrt(u**2+v**2)/upbar-alpha_1)*h/(1.d0-alpha_1)
                else
                   dp=h
                endif

                dp = min(h,max(hs_min,dp))
                ds = min(h,max(h-dp,hs_min))

                cDx=cF

                if (q1d(1,i-1)>depth_tolerance.and. &
                q1d(1,i+1)>depth_tolerance) then
                   cDx = cF - cP*min(0.d0, &
                         sign(1.d0,u)*(q1d(1,i+1)+aux1d(1,i+1) & 
                              -q1d(1,i-1) - aux1d(1,i-1) )/hx/2.d0)
                elseif (q1d(1,i+1)<depth_tolerance.and.u>0.d0) then
                   cDx = 1.732d0*cP
                elseif (q1d(1,i-1)<depth_tolerance.and.u<0.d0) then
                   cDx = 1.732d0*cP
                endif

                drag_force=-dt*.5d0*rho_a/rho_s*cDx*u*ubar

                if (abs(q1d(4,i))>=abs(drag_force)) then
                    q1d(4,i)= drag_force
                    q1d(2,i)= q1d(4,i)/(dp+alpha_1*ds)
                else 
                    q1d(4,i)=-q1d(4,i)
                    q1d(2,i)=-q1d(2,i)
                endif

                drag_force=-dt*.5d0*rho_a/rho_s*cDy*v*ubar

                if (abs(q1d(5,i))>=abs(drag_force)) then
                    q1d(5,i)= q1d(5,i)+drag_force
                    q1d(3,i)= q1d(5,i)/(dp+alpha_1*ds)
                else 
                    q1d(5,i)=-q1d(5,i)
                    q1d(3,i)=-q1d(3,i)
                endif

             endif

       end do

    endif

    call checkval_bing(mx1d,1,meqn,mbc,q1d,alpha_1)

    ! Herschel-Bulkley Friction terms
    do i=1,mx1d

          h=q1d(1,i)

          if (h<depth_tolerance) then
             q1d(:5,i) = 0.d0
          else

             up=q1d(2,i)
             vp=q1d(3,i)
             u =q1d(4,i)/h
             v =q1d(5,i)/h

             upbar = dsqrt(up**2+vp**2)

             if (upbar>1d-8) then
                dp=(dsqrt(u**2+v**2)/upbar-alpha_1)*h/(1.d0-alpha_1)
             else
                dp=h
             endif

             dp = min(h,max(hs_min,dp))
             ds = min(h,max(h-dp,hs_min))

             ! Bingham Friction
             if (dp<depth_tolerance) then
                bx_fric(1:2)=0.d0
                by_fric(1:2)=0.d0
                q1d(2:3,i)=0.d0
             else
                bx_fric(1)=1.d0/dp
                by_fric(1)=1.d0/dp
             endif

             tau_y=tauy_r+(tauy_i-tauy_r)*exp(-q1d(6,i)*remold_coeff)

             if (upbar>1d-20) then
                bx_fric(1)=-bx_fric(1)*dt*tau_y*up/upbar/rho_s
                by_fric(1)=-by_fric(1)*dt*tau_y*vp/upbar/rho_s
             else
                bx_fric(1)=0.d0
                by_fric(1)=0.d0
             endif

             bx_fric(1)=bx_fric(1)/cm_coeff
             by_fric(1)=by_fric(1)/cm_coeff

             if (abs(q1d(2,i))>abs(bx_fric(1))) then
                q1d(2,i)= bx_fric(1)
             else
                bx_fric(1)=-q1d(2,i)
                q1d(2,i)= bx_fric(1)
             endif

             if (abs(q1d(3,i))>abs(by_fric(1))) then
                q1d(3,i)= by_fric(1)
             else
                by_fric(1)=-q1d(3,i)
                q1d(3,i)=by_fric(1)
             endif

             !gamma_r = tau_y/visc_dyn
             visc_dyn = tau_y/gamma_r

             if (ds>depth_tolerance) then
                bx_fric(2)= beta*abs(up/gamma_r/ds)**n_param
                by_fric(2)= beta*abs(vp/gamma_r/ds)**n_param
             else
                bx_fric(2)= 0.d0
                by_fric(2)= 0.d0
             endif

             if (upbar>1d-20) then
                bx_fric(2)=-bx_fric(2)*dt*tau_y*up/upbar/rho_s/cm_coeff
                bx_fric(2)= bx_fric(2)+bx_fric(1)*dp
                by_fric(2)=-by_fric(2)*dt*tau_y*vp/upbar/rho_s/cm_coeff
                by_fric(2)= by_fric(2)+by_fric(1)*dp
             else
                bx_fric(:)=0.d0
                by_fric(:)=0.d0
             endif

             if (abs(q1d(4,i))>abs(bx_fric(2))) then
                q1d(4,i)= bx_fric(2)
             else
                q1d(4,i)= -q1d(4,i)
             endif

             if (abs(q1d(5,i))>abs(by_fric(2))) then
                q1d(5,i)= q1d(5,i) + by_fric(2)
             else
                q1d(5,i)=q1d(5,i)
             endif

          endif

    end do

    call checkval_bing(mx1d,1,meqn,mbc,q1d,alpha_1)

    ! Only lat-long coordinate system supported here right now
    if (coriolis_forcing .and. coordinate_system == 2) then

        do i=1,mx1d
            ! aux(3,:,:) stores the y coordinates multiplied by deg2rad
            fdt = 2.d0 * omega * sin(aux1d(3,i)) * dt

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)
    
            q1d(2,i) = q1d(2,i) * a(1,1) + q1d(3,i) * a(1,2)
            q1d(3,i) = q1d(2,i) * a(2,1) + q1d(3,i) * a(2,2)
        enddo
    endif

end subroutine src1d
