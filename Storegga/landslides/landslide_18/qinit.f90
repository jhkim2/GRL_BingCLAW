subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type,add_perturbation
    use geoclaw_module, only: sea_level
    use geoclaw_module, only: spherical_distance
    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use bing_module
    
    implicit double precision (a-h,o-z)
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
 
    ! Locals
    integer :: i,j,n_cell
    real(kind=8) :: corner(4,2),s1,s2,s3,s4,s5,b5
    real(kind=8) :: b1,b2,b3,b4,thick_avg,thick_max
    real(kind=8) :: total_vol,td1,td2,dist1,dist2
    real(kind=8) :: h(1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8) :: h_init(1-mbc:mx+mbc,1-mbc:my+mbc)

    real(kind=8) :: xm,xc,xp,ym,yc,yp,dx_meters,dy_meters

    q = 0.d0
       
       corner(1,:)=[2.2 ,64.3]
       corner(2,:)=[3.3 ,64.8]
       corner(3,:)=[6.  ,64.5]
       corner(4,:)=[3.4  ,62.6]

       !corner(1,:)=[3.2 ,64. ]
       !corner(2,:)=[4.2 ,64.6]
       !corner(3,:)=[5.7 ,64.7]
       !corner(4,:)=[3.8 ,62.8]
       
       s1=(corner(2,2)-corner(1,2))/(corner(2,1)-corner(1,1))
       s2=(corner(3,2)-corner(2,2))/(corner(3,1)-corner(2,1))
       s3=(corner(4,2)-corner(3,2))/(corner(4,1)-corner(3,1))
       s4=(corner(1,2)-corner(4,2))/(corner(1,1)-corner(4,1))
       s5 = 0.1d0
      
       b1=corner(1,2)-s1*corner(1,1)
       b2=corner(2,2)-s2*corner(2,1)
       b3=corner(3,2)-s3*corner(3,1)
       b4=corner(4,2)-s4*corner(4,1)
       b5 = 64.2d0-s5*4.0d0 
       !b5 = 64.0d0 -s5*4.2d0 
                     
       thick_avg=0.d0
       thick_max=0.d0
       n_cell=0
       total_vol=0.d0

       do j=1,my
           y = ylower + (j-0.5d0)*dy
           do i=1,mx
               x = xlower + (i-0.5d0)*dx

               if (y < s2*x+b2 .and. y > s4*x+b4 .and. &
                   y < s1*x+b1 .and. x<6. .and. y > 62.6d0 )  then

                   if (aux(1,i,j)>-710.d0.and.aux(1,i,j)<=-210.d0) then
                       q(1,i,j)= -(aux(1,i,j)+210.d0)
                   elseif (aux(1,i,j)>-1100.d0.and.aux(1,i,j)<=-710.d0) then
                       q(1,i,j)= 450.d0
                   elseif (aux(1,i,j)>-1500.d0.and.aux(1,i,j)<=-1100.d0) then
                       q(1,i,j)= 450.d0*(1.d0-abs((aux(1,i,j)+1100.d0)/400.d0)**.5)
                   endif

               endif


               if ((x-3.9d0)**2+((y-64.d0)*1.5d0)**2<0.7**2) then
                   q(1,i,j) = max(q(1,i,j),250.d0)
               endif

           end do
       end do

       !q(1,:,:) = 0.8d0*q(1,:,:)

       ! surface smoothing

       h_init(:,:)=q(1,:,:)
       do k=1,100
           h(:,:)=q(1,:,:)
           do j=1,my
               do i=1,mx
                   if (maxval(h_init(i-1:i+1,j-1:j+1))>0.d0) then
                      q(1,i,j)=(4*h(i,j)+h(i+1,j)+h(i-1,j)+h(i,j+1)+h(i,j-1))/8.d0  &
                         +(4*aux(1,i,j)+aux(1,i+1,j)+aux(1,i-1,j)+aux(1,i,j+1) &
                         +aux(1,i,j-1))/8.d0 -aux(1,i,j)
                   endif
               end do
           end do
       end do

       q(1,:,:) = 0.9d0*q(1,:,:)

       do j=1,my
           y = ylower + (j-0.5d0)*dy
           do i=1,mx
               x = xlower + (i-0.5d0)*dx

               if ((x-4.15d0)**2+((y-63.8d0)*1.5d0)**2<0.7**2) then
                   q(1,i,j) = max(q(1,i,j),250.d0)
               endif

               dist_s1 = abs(0.7d0*x +y-5.3d0*0.7d0-63.3d0)/sqrt(0.7d0**2+1.d0)
               if (q(1,i,j)>0.) then
                   q(1,i,j) = q(1,i,j)*exp(-1.d0*dist_s1**(2.d0)/.6d0) 
               !else
               !    q(1,i,j) = q(1,i,j)*exp(-1.d0*dist_s1**(2.d0)/.6d0) 
               endif

               !if (aux(1,i,j)<-2000.d0) q(1,i,j) = 0.1d0*q(1,i,j)

               !if (aux(1,i,j)+q(1,i,j)<-1300.) q(1,i,j)=0.

               if (y<62.) q(1,i,j)=0.

           end do
       end do

       ! surface smoothing
       h_init(:,:)=q(1,:,:)
       do k=1,10
           h(:,:)=q(1,:,:)
           do j=1,my
               do i=1,mx
                   if (maxval(h_init(i-1:i+1,j-1:j+1))>0.d0) then
                      q(1,i,j)=(4*h(i,j)+h(i+1,j)+h(i-1,j)+h(i,j+1)+h(i,j-1))/8.d0  &
                         +(4*aux(1,i,j)+aux(1,i+1,j)+aux(1,i-1,j)+aux(1,i,j+1) &
                         +aux(1,i,j-1))/8.d0 -aux(1,i,j)
                   endif
               end do
           end do
       end do

       do j=1,my
           y = ylower + (j-0.5d0)*dy
           do i=1,mx
               x = xlower + (i-0.5d0)*dx

               if (aux(1,i,j)+q(1,i,j)<-1350.) q(1,i,j)=0.
               if (aux(1,i,j)>-280.) q(1,i,j)=0.

           end do
       end do

       !q = 0.d0

    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
       call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    endif

       total_vol=0.
       n_cell   =0
       do j=1,my

          ym = ylower + (j - 1.d0) * dy
          yc = ylower + (j - 0.5d0) * dy
          yp = ylower + j * dy

          do i=1,mx
             xm = xlower + (i - 1.d0) * dx
             xc = xlower + (i - 0.5d0) * dx
             xp = xlower + i * dx

             if (coordinate_system == 2) then
                 ! Convert distance in lat-long to meters
                 dx_meters = spherical_distance(xp,yc,xm,yc)
                 dy_meters = spherical_distance(xc,yp,xc,ym)
             else
                 dx_meters = dx
                 dy_meters = dy
             endif

             total_vol = total_vol + q(1,i,j)*dx_meters*dy_meters/1.d9
             if (q(1,i,j)>0.d0) then
                 total_area = total_area + dx_meters*dy_meters
             endif
          end do
       end do
    
       print*,'Maximum thickness =',maxval(q(1,:,:)),'meter'
       print*,'landslide area    =',total_area/1.d6,'km^2'
       print*,'landslide volume  =',total_vol,'km^3'

    call checkval_bing(mx,my,meqn,mbc,q,alpha_1)
    
       ! Test front shear thinning
       do j=1,my
           y = ylower + (j-0.5d0)*dy
           do i=1,mx
               x = xlower + (i-0.5d0)*dx

               dist1 = sqrt((x-3.5d0 )**2 + ((y-64.2d0)*2.d0)**2)

               if (q(1,i,j)>0.) then
                   !q(6,i,j) = 1.d5*(exp(-(dist1*1.d0)*4))
               endif

           end do
       end do

end subroutine qinit
