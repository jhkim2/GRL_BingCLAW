subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type,add_perturbation
    use geoclaw_module, only: sea_level
    use bing_module
    use geoclaw_module, only: spherical_distance
    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    
    implicit double precision (a-h,o-z)
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
 
    ! Locals
    integer :: i,j,n_cell
    real(kind=8) :: corner(6,2),s1,s2,s3,s4,s5,b5,s6,b6
    real(kind=8) :: b1,b2,b3,b4,thick_avg,thick_max
    real(kind=8) :: total_vol,td1,td2,dist1,dist2
    real(kind=8) :: h(1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8) :: h_init(1-mbc:mx+mbc,1-mbc:my+mbc)

    real(kind=8) :: xm,xc,xp,ym,yc,yp,dx_meters,dy_meters

    q = 0.d0
       
       corner(1,:)=[8.9 , 67.3]
       corner(2,:)=[9.2 , 67.45]
       corner(3,:)=[9.3, 67.75]
       corner(4,:)=[8.95 , 68.1]
       corner(5,:)=[7.95 , 68.3]
       corner(6,:)=[7.45 , 67.7]

       s1=(corner(2,2)-corner(1,2))/(corner(2,1)-corner(1,1))
       s2=(corner(3,2)-corner(2,2))/(corner(3,1)-corner(2,1))
       s3=(corner(4,2)-corner(3,2))/(corner(4,1)-corner(3,1))
       s4=(corner(5,2)-corner(4,2))/(corner(5,1)-corner(4,1))
       s5=(corner(6,2)-corner(5,2))/(corner(6,1)-corner(5,1))
       s6=(corner(1,2)-corner(6,2))/(corner(1,1)-corner(6,1))
      
       b1=corner(1,2)-s1*corner(1,1)
       b2=corner(2,2)-s2*corner(2,1)
       b3=corner(3,2)-s3*corner(3,1)
       b4=corner(4,2)-s4*corner(4,1)
       b5=corner(5,2)-s5*corner(5,1)
       b6=corner(6,2)-s6*corner(6,1)
                     
       thick_avg=0.d0
       thick_max=0.d0
       n_cell=0
       total_vol=0.d0

       do j=1,my
           y = ylower + (j-0.5d0)*dy
           do i=1,mx
               x = xlower + (i-0.5d0)*dx

               if (y > s1*x+b1 .and. y > s2*x+b2 .and. y < s3*x+b3 .and. &
                   y < s4*x+b4 .and. y < s5*x+b5 .and. y > s6*x+b6)  then
                   if (aux(1,i,j)>-600.d0.and.aux(1,i,j)<=-400.d0) then
                       !q(1,i,j)= -(aux(1,i,j)+400.d0)
                       q(1,i,j)= 250.d0*(1.d0-((aux(1,i,j)+600.d0)/200.d0)**2)
                   elseif (aux(1,i,j)>-2000.d0.and.aux(1,i,j)<=-600.d0) then
                       q(1,i,j)= 250.d0 !-35.d0/50.d0*(aux(1,i,j)+350.d0)+100.d0
                   elseif (aux(1,i,j)>-3000.d0.and.aux(1,i,j)<=-2000.d0) then
                       !q(1,i,j)= 250.d0*(1.d0-((aux(1,i,j)+2000.d0)/500.d0)**(0.5))
                       q(1,i,j)= 250.d0*(1.d0-abs((aux(1,i,j)+2000.d0)/1000.d0)**.5)
                   endif
               endif
               if ((x-8.1d0)**2+((y-68.2d0)*3.d0)**2<0.1d0) then
                   !q(1,i,j) = max(q(1,i,j),100.d0)
               endif

           end do
       end do

       ! surface smoothing
       h_init(:,:)=q(1,:,:)
       do k=1,20
           h(:,:)=q(1,:,:)
           do j=1,my
               do i=1,mx
                   if (maxval(h_init(i-1:i+1,j-1:j+1))>0.d0) then
                      q(1,i,j)=(2.d0*h(i,j)+h(i+1,j)+h(i-1,j))/4.d0  &
                         +(2.d0*aux(1,i,j)+aux(1,i+1,j)+aux(1,i-1,j))/4.d0 -aux(1,i,j)
                      q(1,i,j) = max(0.d0,q(1,i,j))
                   endif
               end do
           end do
           h(:,:)=q(1,:,:)
           do j=1,my
               do i=1,mx
                   if (maxval(h_init(i-1:i+1,j-1:j+1))>0.d0) then
                      q(1,i,j)=(2.d0*h(i,j)+h(i,j+1)+h(i,j-1))/4.d0  &
                         +(2.d0*aux(1,i,j)+aux(1,i,j+1)+aux(1,i,j-1))/4.d0 -aux(1,i,j)
                      q(1,i,j) = max(0.d0,q(1,i,j))
                   endif
               end do
           end do
       end do

       q(1,:,:) = .7d0*q(1,:,:)

       do k=1,2
           h(:,:)=q(1,:,:)
           do j=1,my
               do i=1,mx
                   if (maxval(h(i-1:i+1,j-1:j+1))>0.d0) then
                      q(1,i,j)=(2.d0*h(i,j)+h(i+1,j)+h(i-1,j))/4.d0  &
                         +(2.d0*aux(1,i,j)+aux(1,i+1,j)+aux(1,i-1,j))/4.d0 -aux(1,i,j)
                      q(1,i,j) = max(0.d0,q(1,i,j))
                   endif
               end do
           end do
           h(:,:)=q(1,:,:)
           do j=1,my
               do i=1,mx
                   if (maxval(h(i-1:i+1,j-1:j+1))>0.d0) then
                      q(1,i,j)=(2.d0*h(i,j)+h(i,j+1)+h(i,j-1))/4.d0  &
                         +(2.d0*aux(1,i,j)+aux(1,i,j+1)+aux(1,i,j-1))/4.d0 -aux(1,i,j)
                      q(1,i,j) = max(0.d0,q(1,i,j))
                   endif
               end do
           end do
       end do

       do j=1,my
           y = ylower + (j-0.5d0)*dy
           do i=1,mx
               x = xlower + (i-0.5d0)*dx

               if (y > s1*x+b1 .and. y > s2*x+b2 .and. y < s3*x+b3 .and. &
                   y < s4*x+b4 .and. y < s5*x+b5 .and. y > s6*x+b6)  then
                   if (aux(1,i,j)>-2400.d0.and.aux(1,i,j)<=-2000.d0) then
                       !q(1,i,j)= 200.d0
                   endif
               endif

               !if (aux(1,i,j)<-2350.d0.and.(x-7.9d0)+y-68.2d0>0.d0) q(1,i,j)=0.d0
               if ( y < 67.3d0) q(1,i,j)=0.d0

           end do
       end do

       ! surface smoothing
       h_init(:,:)=q(1,:,:)
       do k=1,2
           h(:,:)=q(1,:,:)
           do j=1,my
               do i=1,mx
                   if (h_init(i,j)>0.d0) then
                      q(1,i,j)=(2.d0*h(i,j)+h(i+1,j)+h(i-1,j))/4.d0  &
                         +(2.d0*aux(1,i,j)+aux(1,i+1,j)+aux(1,i-1,j))/4.d0 -aux(1,i,j)
                      q(1,i,j) = max(0.d0,q(1,i,j))
                   endif
               end do
           end do
           h(:,:)=q(1,:,:)
           do j=1,my
               do i=1,mx
                   if (h_init(i,j)>0.d0) then
                      q(1,i,j)=(2.d0*h(i,j)+h(i,j+1)+h(i,j-1))/4.d0  &
                         +(2.d0*aux(1,i,j)+aux(1,i,j+1)+aux(1,i,j-1))/4.d0 -aux(1,i,j)
                      q(1,i,j) = max(0.d0,q(1,i,j))
                   endif
               end do
           end do
       end do

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
             if (q(1,i,j)>1.d0) then
                 total_area = total_area + dx_meters*dy_meters/1.d6
             endif
          end do
       end do
    
       print*,'Maximum thickness =',maxval(q(1,:,:)),'meter'
       print*,'landslide area    =',total_area,'km^2'
       print*,'landslide volume  =',total_vol,'km^3'

    call checkval_bing(mx,my,meqn,mbc,q,alpha_1)
    
end subroutine qinit
