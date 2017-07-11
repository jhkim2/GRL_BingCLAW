
subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use geoclaw_module, only: grav,spherical_distance
    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use bing_module

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j
    real(kind=8) :: h, hu, hv

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.d-5

    ! Bingham
    real(kind=8), parameter:: hs_min = 0.0d0
    real(kind=8) :: u,v,up,vp,dp,ds,upbar,gxmod,gymod,g,cDx,cDy,ubar,usum
    real(kind=8) :: bx_fric(2),by_fric(2),nu_m,drag_force,cF,cP,tau_y
    real(kind=8) :: dx_meters,dy_meters,x,y,dgamma,gamma
    real(kind=8) :: xm,xc,xp,ym,yc,yp
    double precision :: qold(5,1-mbc:mx+mbc,1-mbc:my+mbc)

    g = grav*(1.d0-rho_a/rho_s)

    !print*,'Max. speed =',dsqrt(maxval(q(2,:,:)**2+q(3,:,:)**2))
	
    call checkval_bing(mx,my,meqn,mbc,q,alpha_1)

    if (hydrodrag) then

       cF = cF_hyd
       cP = cP_hyd

       do i=1,mx

          xm = xlower + (i - 1.d0) * dx
          xc = xlower + (i - 0.5d0) * dx
          xp = xlower + i * dx

          do j=1,my
             y=ylower+(j-.5d0)*dy

             ym = ylower + (j - 1.d0) * dy
             yc = ylower + (j - 0.5d0) * dy
             yp = ylower + j * dy

             h=q(1,i,j)

             if (h>depth_tolerance) then

                if (coordinate_system == 2) then
                    ! Convert distance in lat-long to meters
                    dx_meters = spherical_distance(xp,yc,xm,yc)
                    dy_meters = spherical_distance(xc,yp,xc,ym)
                else
                    dx_meters = dx
                    dy_meters = dy
                endif

                up=q(2,i,j)
                vp=q(3,i,j)
                u =q(4,i,j)/h
                v =q(5,i,j)/h
                gxmod = aux(4,i,j)
                gymod = aux(5,i,j)

                ubar=dsqrt(u**2+v**2)

                upbar = dsqrt(up**2+vp**2)

                if (upbar>1d-4) then
                   dp=(dsqrt(u**2+v**2)/upbar-alpha_1)*h/(1.d0-alpha_1)
                else
                   dp=h
                endif

                dp = min(h,max(hs_min,dp))
                ds = min(h,max(h-dp,hs_min))

                cDx = cF
                cDy = cF

                if (q(1,i+1,j)<depth_tolerance.and.u>0.d0) then
                   cDx = 1.732d0*cP
                elseif (q(1,i-1,j)<depth_tolerance.and.u<0.d0) then
                   cDx = 1.732d0*cP
                elseif (q(1,i,j+1)<depth_tolerance.and.v>0.d0) then
                   cDy = 1.732d0*cP
                elseif (q(1,i,j-1)<depth_tolerance.and.v<0.d0) then
                   cDy = 1.732d0*cP
                !elseif (q(1,i+1,j+1)<depth_tolerance.and.u>0.d0.and.v>0.d0) then
                !   cDx = 1.732d0*cP
                !   cDy = 1.732d0*cP
                !elseif (q(1,i+1,j-1)<depth_tolerance.and.u>0.d0.and.v<0.d0) then
                !   cDx = 1.732d0*cP
                !   cDy = 1.732d0*cP
                !elseif (q(1,i-1,j+1)<depth_tolerance.and.u<0.d0.and.v>0.d0) then
                !   cDx = 1.732d0*cP
                !   cDy = 1.732d0*cP
                !elseif (q(1,i-1,j-1)<depth_tolerance.and.u<0.d0.and.v<0.d0) then
                !   cDx = 1.732d0*cP
                !   cDy = 1.732d0*cP
                else
                   cDx = cF - cP*min(0.d0, &
                         sign(1.d0,u)*(q(1,i+1,j)-q(1,i-1,j) )/dx_meters/2.d0)
                   cDy = cF - cP*min(0.d0, &
                         sign(1.d0,v)*(q(1,i,j+1)-q(1,i,j-1) )/dy_meters/2.d0)

                endif

                drag_force=-dt*.5d0*rho_a/rho_s*cDx*u*ubar

                if (abs(q(4,i,j))>=abs(drag_force)) then
                    q(4,i,j)= q(4,i,j)+drag_force
                    q(2,i,j)= q(4,i,j)/(dp+alpha_1*ds)
                else 
                    q(4,i,j)= 0.d0
                    q(2,i,j)= 0.d0
                endif

                drag_force=-dt*.5d0*rho_a/rho_s*cDx*up*upbar/dp

                !if (abs(q(2,i,j))>=abs(drag_force)) then
                !    q(2,i,j)= q(2,i,j)+drag_force
                !else 
                !    q(2,i,j)= 0.d0
                !endif

                drag_force=-dt*.5d0*rho_a/rho_s*cDy*v*ubar

                if (abs(q(5,i,j))>=abs(drag_force)) then
                    q(5,i,j)= q(5,i,j)+drag_force
                    q(3,i,j)= q(5,i,j)/(dp+alpha_1*ds)
                else 
                    q(5,i,j)= 0.d0
                    q(3,i,j)= 0.d0
                endif

                drag_force=-dt*.5d0*rho_a/rho_s*cDy*vp*upbar/dp

             endif

          end do
       end do

    endif

    call checkval_bing(mx,my,meqn,mbc,q,alpha_1)

    ! Herschel-Bulkley Friction terms
    do i=1,mx
       do j=1,my
          h=q(1,i,j)

          if (h<depth_tolerance) then
             q(1:5,i,j) = 0.d0
          else

             up=q(2,i,j)
             vp=q(3,i,j)
             u =q(4,i,j)/h
             v =q(5,i,j)/h

             upbar = dsqrt(up**2+vp**2)
             ubar  = dsqrt(u**2+v**2)

             if (upbar>1d-8) then
                dp=(ubar/upbar-alpha_1)*h/(1.d0-alpha_1)
             else
                dp=h
             endif

             dp = min(h,max(hs_min,dp))
             ds = min(h,max(h-dp,hs_min))

             ! Bingham Friction
             if (dp<depth_tolerance) then
                bx_fric(1:2)=0.d0
                by_fric(1:2)=0.d0
             else
                bx_fric(1)=1.d0/dp
                by_fric(1)=1.d0/dp
             endif

             if (remolding) then
                tau_y=tauy_r+(tauy_i-tauy_r)*exp(-q(6,i,j)*remold_coeff)
             else
                tau_y = tauy_i
             endif

             if (upbar>1d-20) then
                bx_fric(1)=-bx_fric(1)*dt*tau_y*up/upbar/rho_s
                by_fric(1)=-by_fric(1)*dt*tau_y*vp/upbar/rho_s
             else
                bx_fric(1)=0.d0
                by_fric(1)=0.d0
             endif

             bx_fric(1)=bx_fric(1)/cm_coeff
             by_fric(1)=by_fric(1)/cm_coeff

             if (abs(q(2,i,j))>abs(bx_fric(1))) then
                q(2,i,j)= q(2,i,j)/(1- bx_fric(1)/up)
             else
                !print*,i,dt,q(1,i,j),q(2,i,j),bx_fric(1)
                !print*,dp
                bx_fric(1)=-q(2,i,j)
                q(2,i,j)=0.d0
             endif

             if (abs(q(3,i,j))>abs(by_fric(1))) then
                q(3,i,j)= q(3,i,j)/(1-by_fric(1)/vp)
             else
                by_fric(1)=-q(3,i,j)
                q(3,i,j)=0.d0
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

             if (abs(q(4,i,j))>abs(bx_fric(2))) then
                q(4,i,j)= q(4,i,j) + bx_fric(2)
             else
                q(4,i,j)=0.d0
             endif

             if (abs(q(5,i,j))>abs(by_fric(2))) then
                q(5,i,j)= q(5,i,j) + by_fric(2)
             else
                q(5,i,j)=0.d0
             endif

             if (remolding) then

                if (ds>0.0001d0) then
                   q(6,i,j)= q(6,i,j)+(n_param+1.d0)/n_param &
                            *dsqrt(up**2+vp**2)/ds*dt
                endif

             endif

          endif

       end do
    end do

    call checkval_bing(mx,my,meqn,mbc,q,alpha_1)

    ! viscous friction term
    if (visc_fric) then    
      do i=1,mx
          do j=1,my
          h = q(1,i,j)
          if (h.gt.1d-3) then
            up = q(2,i,j)
            vp = q(3,i,j)
            hu = q(4,i,j)
            hv = q(5,i,j)
            u = hu/h
            v = hv/h

            upbar = dsqrt(up**2+vp**2)

            if (upbar>1d-4) then
               dp=(dsqrt(u**2+v**2)/upbar-alpha_1)*h/(1.d0-alpha_1)
            else
               dp=h
            endif

            dp = min(h,max(hs_min,dp))
            ds = min(h,max(h-dp,hs_min))
                                
            if (u*q(2,i,j).gt.0.d0) then
               gamma=3.d0*visc_kin*u/h
               dgamma=dt*gamma
                   
               if (abs(dgamma).gt.abs(q(4,i,j))) then
                  q(4,i,j) = 0.d0
               else 
                  q(4,i,j) = q(4,i,j) - dgamma
               endif
               q(2,i,j)= q(4,i,j)/(Dp+alpha_1*Ds)
            endif

            if (v*q(5,i,j).gt.0.d0) then
               gamma=3.d0*visc_kin*v/h
               dgamma=dt*gamma
                   
               if (abs(dgamma).gt.abs(q(5,i,j))) then
                  q(5,i,j) = 0.d0
               else 
                  q(5,i,j) = q(5,i,j) - dgamma
               endif
               q(3,i,j)= q(5,i,j)/(Dp+alpha_1*Ds)
            endif

           endif

         enddo
      enddo
      call checkval_bing(mx,my,meqn,mbc,q,alpha_1)
    endif    
    ! End of viscous friction term

    ! Artificial viscosity
      if (a_vis>0.d0) then
         qold=q(1:5,:,:)
         do i=1,mx
            do j=1,my
               h = q(1,i,j)
               if (h.gt.1d-3.and.q(1,i+1,j)>0.d0.and.q(1,i-1,j)>0.d0.and. &
                  q(1,i,j+1)>0.d0.and.q(1,i,j-1)>0.d0) then
                  nu_m=a_vis * &
                       abs(q(1,i+1,j)+q(1,i-1,j)+q(1,i,j+1)+q(1,i,j-1)-4*h) &
                      /abs(q(1,i+1,j)+q(1,i-1,j)+q(1,i,j+1)+q(1,i,j-1)+4*h)

                  u = qold(4,i,j)/h + nu_m* &
                      (qold(4,i+1,j)/q(1,i+1,j)+qold(4,i-1,j)/q(1,i-1,j) &
                      +qold(4,i,j+1)/q(1,i,j+1)+qold(4,i,j-1)/q(1,i,j-1) &
                      -4*qold(4,i,j)/q(1,i,j))
                  q(4,i,j)=h*u
                  up= qold(2,i,j) + nu_m* &
                     (qold(2,i+1,j)+qold(2,i-1,j)+qold(2,i,j+1) &
                     +qold(2,i,j-1)-4*qold(2,i,j))
                  q(2,i,j)=up
               endif
            end do
         end do
         call checkval_bing(mx,my,meqn,mbc,q,alpha_1)
      endif
    ! End of artificial viscosity

end subroutine src2
