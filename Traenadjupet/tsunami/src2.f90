subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: manning_coefficient
    use geoclaw_module, only: manning_break, num_manning
    use geoclaw_module, only: rad2deg
    use bous_module

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, nman
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    logical :: verbose
    real(kind=8) :: wave_max_slope,eta(1:mx),eta_max

    if (t<=2.d0*dt) then
        aux(4,:,:) = 0.d0
    else
        do j = 1,my
            do i=1,mx
                if (q(1,i,j)>0.d0) then
                    aux(4,i,j) = max(aux(4,i,j),q(1,i,j)+aux(1,i,j))
                endif
            enddo
        enddo
    endif

    ! Friction source term
    if (friction_forcing) then
        do j=1,my
            do i=1,mx
                ! Extract appropriate momentum
                if (q(1,i,j) < depth_tolerance) then
                    q(2:3,i,j) = 0.d0
                else
                    ! Apply friction source term only if in shallower water
                    if (q(1,i,j) <= friction_depth) then
                        do nman = num_manning, 1, -1
                            if (aux(1,i,j) .lt. manning_break(nman)) then
                                coeff = manning_coefficient(nman)
                            endif
                        enddo
                        ! Calculate source term
                        gamma = sqrt(q(2,i,j)**2 + q(3,i,j)**2) * g     &   
                              * coeff**2 / (q(1,i,j)**(7.d0/3.d0))
                        dgamma = 1.d0 + dt * gamma
                        q(2, i, j) = q(2, i, j) / dgamma
                        q(3, i, j) = q(3, i, j) / dgamma
                    endif
                endif
            enddo
        enddo
    endif
    ! End of friction source term

    ! Coriolis source term
    ! TODO: May want to remove the internal calls to coriolis as this could 
    !       lead to slow downs.
    if (coriolis_forcing) then
        do j=1,my
            y = ylower + (j - 0.5d0) * dy
            fdt = coriolis(y) * dt ! Calculate f dependent on coordinate system

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)

            do i=1,mx
                q(2,i,j) = q(2, i, j) * a(1,1) + q(3, i, j) * a(1,2)
                q(3,i,j) = q(2, i, j) * a(2,1) + q(3, i, j) * a(2,2)
            enddo
        enddo
    endif
    ! End of coriolis source term
          
    ! Computation of Dispersive Terms
     
      if (bouss) then

         call src2_bous(meqn,mbc,mx,my, &
                  xlower,ylower,dx,dy,q,maux,aux,t,dt)
     
      endif

      verbose = .True.

      if (verbose) then
         j = 1
         wave_max_slope = 0.d0
         eta = 0.d0
         do i=1,mx
            if (q(1,i,j)>0.d0) then
               wave_max_slope=max(wave_max_slope, &
                  abs( (q(1,i,j)+aux(1,i,j)-(q(1,i-1,j)+aux(1,i-1,j)) )/dx/1.d0))
               eta(i) = q(1,i,j)+aux(1,i,j)
            endif
         end do
         eta_max=maxval(eta)/dx
         !wave_max_slope = wave_max_slope*rad2deg
         !print*,'maximum wave slope =', wave_max_slope
         !print*,'max possible slope =', eta_max/dx
         !print*,'ratio              =', wave_max_slope/eta_max,eta_max*dx,dx
      endif

      return
      end
      
