!======================================================================
       subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,  &
                      ql,qr,auxl,auxr,fwave,s,amdq,apdq)
!======================================================================
!                                                                     !
!      # Bingham Fluid Solver with Herschel-Bulkley rheology          !
!        Modified Geoclaw rpn2 solver                                 !
!                                                                     !
!           Jihwan Kim, NGI, Oslo Norway, 2014-2015                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use geoclaw_module, only: grav, drytol => dry_tolerance
      use geoclaw_module, only: earth_radius, deg2rad
      use amr_module, only: mcapa
      use bing_module

      implicit none

      !input
      integer maxm,meqn,maux,mwaves,mbc,mx,ixy

      double precision  fwave(meqn, mwaves, 1-mbc:maxm+mbc)
      double precision  s(mwaves, 1-mbc:maxm+mbc)
      double precision  ql(meqn, 1-mbc:maxm+mbc)
      double precision  qr(meqn, 1-mbc:maxm+mbc)
      double precision  apdq(meqn,1-mbc:maxm+mbc)
      double precision  amdq(meqn,1-mbc:maxm+mbc)
      double precision  auxl(maux,1-mbc:maxm+mbc)
      double precision  auxr(maux,1-mbc:maxm+mbc)

      !local only
      integer m,i,mw,maxiter,mu,nv
      double precision wall(3)
      double precision fw(5,5)
      double precision sw(5)

      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL
      double precision bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision s1m,s2m,phiL,phiR
      double precision hstar,hstartest,hstarHLL,sLtest,sRtest
      double precision tw,dxdc,hbar
	  
      !Bing only
      double precision upL,upR,vpL,vpR,r1,dpL,dpR,gmod
      double precision WR(5),EVR(5,5),EVR1(5,5)
      double precision del(5),beta1,beta2
      double precision sfract(5,1-mbc:maxm+mbc),sfract2(5,1-mbc:maxm+mbc)
      double precision lambdaL(5),lambdaR(5)
      integer INFO, IPIV(5)

      logical rare1,rare2,efix

      !Initialize Riemann problem for grid interface
      s(:,:)=0.d0
      fwave(:,:,:)=0.d0

      efix = .True.
      sfract = 0.d0
      sfract2= 0.d0
	  
      beta1= (alpha_2-alpha_1)/(1.d0-alpha_1)
      beta2= (1.d0-alpha_2)/(1.d0-alpha_1)

      !loop through Riemann problems at each grid cell
      do i=2-mbc,mx+mbc

!-----------------------Initializing-----------------------------------
         !inform of a bad riemann problem from the start
         if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
            write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
         endif

!        !set normal direction
         if (ixy.eq.1) then
            mu=2
            nv=3
         else
            mu=3
            nv=2
         endif

         !zero (small) negative values if they exist
         if (qr(1,i-1).lt.0.d0) then
             qr(1:6,i-1)=0.d0
         endif

         if (ql(1,i).lt.0.d0) then
             ql(1:6,i) = 0.d0
         endif

         !skip problem if in a completely dry area
         if (qr(1,i-1) <= drytol .and. ql(1,i) <= drytol) then
            go to 30
         endif

         if (auxR(mu+4,i-1)>0.d0.and.auxL(mu+4,i)>0.d0) go to 30

         !Riemann problem variables
         hL = qr(1,i-1) 
         hR = ql(1,i)
         upL = qr(mu,i-1)
         upR = ql(mu,i)
         vpL = qr(nv,i-1)
         vpR = ql(nv,i)
         huL = qr(mu+2,i-1)
         huR = ql(mu+2,i)
         hvL = qr(nv+2,i-1)
         hvR = ql(nv+2,i)
         bL = auxR(1,i-1)
         bR = auxL(1,i)

         gmod = .5d0*(auxR(mu+2,i-1)+auxL(mu+2,i))
		 		 
         !check for wet/dry boundary
         if (hR.gt.drytol) then
            uR=huR/hR
            vR=hvR/hR
            phiR = 0.5d0*gmod*hR**2 + huR**2/hR
         else
            hR = 0.d0
            huR = 0.d0
            hvR = 0.d0
            uR = 0.d0
            vR = 0.d0
            upR= 0.d0
            vpR= 0.d0
            phiR=0.d0
         endif

         if (hL.gt.drytol) then
            uL=huL/hL
            vL=hvL/hL
            phiL = 0.5d0*gmod*hL**2 + huL**2/hL
         else
            hL=0.d0
            huL=0.d0
            hvL=0.d0
            uL=0.d0
            vL=0.d0
            upL=0.d0
            vpL=0.d0
            phiL=0.d0
         endif

         !determine plug layer depth
         if (sqrt(upL**2+vpL**2)>1d-2) then
             r1= sqrt(uL**2+vL**2)/sqrt(upL**2+vpL**2)
            dpL= (r1-alpha_1)*hL/(1.d0-alpha_1)
         else
            dpL=hL
         endif
         if (upL**2+vpL**2==uL**2+vL**2) dpL=hL
         
         dpL = max(0.d0,min(hL,dpL))

         if (sqrt(upR**2+vpR**2)>1d-2) then
             r1= sqrt(uR**2+vR**2)/sqrt(upR**2+vpR**2)
            dpR= (r1-alpha_1)*hR/(1.d0-alpha_1)
         else
            dpR=hR
         endif
         if (upR**2+vpR**2==uR**2+vR**2) dpR=hR

         dpR = max(0.d0,min(hR,dpR))

         wall(1) = 1.d0
         wall(2) = 1.d0
         wall(3) = 1.d0
		 
         if (hR.le.drytol) then
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,       &
                                       rare1,rare2,1,drytol,gmod)
            hstartest=max(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
               !bR=hstartest+bL
               wall(2)=0.d0
               wall(3)=0.d0
               hR=hL
               dpR=dpL
               huR=-huL
               bR=bL
               uR=-uL
               upR=-upL
               phiR=phiL
               vR=vL
               vpR=vpL
            elseif (hL+bL.lt.bR) then
               bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,       &
                                       rare1,rare2,1,drytol,gmod)
            hstartest=max(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
               !bL=hstartest+bR
               wall(1)=0.d0
               wall(2)=0.d0
               hL=hR
               dpL=dpR
               huL=-huR
               bL=bR
               uL=-uR
               phiL=phiR
               upL=-upR
               vL=vR
               vpL=vpR
            elseif (hR+bR.lt.bL) then
               bL=hR+bR
            endif
         endif

!------------------------------------------------------------------------
         !go to 20
         ! if shear layer exists, solve the full system
         if (hL-dpL>1d-3.or.hR-dpR>1d-3) then
             go to 20
         endif
!------------------------------------------------------------------------
         ! Under development
         ! if there is no plug layer, then solve SWE instead

         !determine wave speeds
         sL=uL-sqrt(gmod*hL) ! 1 wave speed of left state
         sR=uR+sqrt(gmod*hR) ! 2 wave speed of right state

         uhat=(sqrt(gmod*hL)*uL + sqrt(gmod*hR)*uR) &
              /(sqrt(gmod*hR)+sqrt(gmod*hL)) ! Roe average
         chat=sqrt(gmod*0.5d0*(hR+hL)) ! Roe average
         sRoe1=uhat-chat ! Roe wave speed 1 wave
         sRoe2=uhat+chat ! Roe wave speed 2 wave

         sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
         sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

         !--------------------end initializing...finally----------
         !solve Riemann problem.

         maxiter = 1

         call riemann_aug_JCP(maxiter,3,3,hL,hR,huL, &
              huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2, &
                                         drytol,gmod,sw(1:3),fw(1:3,1:3))

        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)
            fw(1,mw)=fw(1,mw)*wall(mw) 
            fw(2,mw)=fw(2,mw)*wall(mw)
            fw(3,mw)=fw(3,mw)*wall(mw)
         enddo

         do mw=1,3
            s(mw,i)=sw(mw)
            fwave(1,mw,i)=fw(1,mw)
            if (abs(fw(1,mw))>1.d-300) then
            ! Not sure on this choice of fwave, but does not make a difference
               fwave(mu,mw,i)=fw(2,mw)/fw(1,mw)
               fwave(nv,mw,i)=fw(3,mw)/fw(1,mw)
            endif
            fwave(mu+2,mw,i)=fw(2,mw)
            fwave(nv+2,mw,i)=fw(3,mw)
         enddo

         go to 25

 20      continue

!-------------------------------------------------------------------------------------
         ! if there is a shear layer, solve the full system

         ! Find the eigenvalues and eigenvectors
         call find_eigval(hL,hR,dpL,dpR,upL,upR,vpL,vpR,   &
                 uL,uR,vL,vR,gmod,WR(1:5),beta1,beta2)

         if (efix) then

            lambdaL = 0.d0
            lambdaR = 0.d0

            if (hL>drytol) then
               call find_eigval(hL,hL,dpL,dpL,upL,upL,vpL,vpL,   &
                 uL,uL,vL,vL,gmod,lambdaL,beta1,beta2)
            endif
            if (hR>drytol) then
               call find_eigval(hR,hR,dpR,dpR,upR,upR,vpR,vpR,   &
                 uR,uR,vR,vR,gmod,lambdaR,beta1,beta2)
            endif

            do mw=1,3
               sR = lambdaR(mw) 
               sL = lambdaL(mw)
               if (sR>1.d-4.and.sL<-1.d-4.and.abs(WR(mw))>1d-6 &
               .and.sL<WR(mw).and.WR(mw)<sR) then
                  WR(mw) = .5d0*(sR+sL)
                  sfract(mw,i) = (sR-WR(mw))/(sR-sL)*sL/WR(mw)
                  sfract2(mw,i) =(1.d0-(sR-WR(mw))/(sR-sL))*sR/WR(mw)
                  sfract(mw,i) = abs(sL)/(sR-sL)
                  sfract2(mw,i) = 1.d0 - sfract(mw,i)
               endif
            enddo

         endif

         call find_eigvec(hL,hR,dpL,dpR,upL,upR,vpL,vpR,   &
                 uL,uR,vL,vR,gmod,WR(1:5),EVR(1:5,1:5),beta1,beta2)
	 
         hbar = (hR+hL)/2.d0

	 sw(1:5) = WR(1:5)
		 
         EVR1(1:5,1:5)= EVR(1:5,1:5)
		 
         del(1) = hR*uR-hL*uL

         del(2) = (upR**2 - upL**2)/2.d0
         del(2) = del(2) + gmod*(hR+bR-hL-bL)/cm_coeff 

         del(3) = .5d0*(upR+upL)*(vpR-vpL)

         del(4) = beta2*(hR*uR*upR-hL*uL*upL)
         del(4) = del(4)+beta1*(hR*upR**2-hL*upL**2)  &
                 + gmod*hbar*(hR+bR-hL-bL)/cm_coeff 

         del(5) = beta2*(hR*vR*upR-hL*vL*upL)
         del(5) = del(5)+beta1*(hR*upR*vpR-hL*upL*vpL)

         !Compute the inverse of the eigenvector matrix	 
         call DGESV( 5, 1, EVR1(1:5,1:5), 5, IPIV, del(1:5), 5, INFO )
		 		 
         if (INFO.ne.0) then
            print*,'error with 5x5 eigenvalues solver'
            print*,i,ixy,INFO
         endif
		 
         do mw=1,5
            fw(1:5,mw) = del(mw)*EVR(1:5,mw)
         enddo 

         !eliminate ghost fluxes for wall
         if (wall(1).eq.0.d0) then
            do mw=1,mwaves
               if (sw(mw).le.0.d0) then
                  sw(mw) = 0.d0
                  fw(1:5,mw) = 0.d0
               endif
            end do
         endif
		 
         if (wall(3).eq.0.d0) then
            do mw=1,mwaves
               if (sw(mw).ge.0.d0) then
                  sw(mw) = 0.d0
                  fw(1:5,mw) = 0.d0
               endif
            end do
         endif

         do mw=1,5
            s(mw,i)=sw(mw)
            fwave(1,mw,i)=fw(1,mw)
            fwave(mu,mw,i)=fw(2,mw)
            fwave(nv,mw,i)=fw(3,mw)
            fwave(mu+2,mw,i)=fw(4,mw)
            fwave(nv+2,mw,i)=fw(5,mw) !+hR*uR*vR - hL*uL*vL 
         enddo

 25      continue

         !solve advection equation for the remolding variable
         s(6,i) = 0.d0
         fwave(:,6,i) = 0.d0
         fwave(6,:,i) = 0.d0

         if (remolding) then
            if (hL<drytol) then
               uhat = uR
            elseif (hR<drytol) then
               uhat = uL
            else
               uhat  =(sqrt(hL)*uL + sqrt(hR)*uR) &
                  /(sqrt(hR)+sqrt(hL)) ! Roe average
               uhat = .5d0*(uR+uL)
            endif
            s(6,i) = alpha_1*uhat
            fwave(6,6,i) = alpha_1*(uR*ql(6,i)-uL*qr(6,i-1))

            if (wall(1).eq.0.d0) then
               if (s(6,i).le.0.d0) then
                  s(6,i) = 0.d0
                  fwave(:,6,i) = 0.d0
               endif
            endif
		 
            if (wall(3).eq.0.d0) then
               if (s(6,i).ge.0.d0) then
                  s(6,i) = 0.d0
                  fwave(:,6,i) = 0.d0
               endif
            endif
         endif

 30      continue

      enddo


!==========Capacity for mapping from latitude longitude to physical space====
      if (mcapa.gt.0) then
         do i=2-mbc,mx+mbc
            if (ixy.eq.1) then
               dxdc=(earth_radius*deg2rad)
            else
               dxdc=earth_radius*cos(auxl(3,i))*deg2rad
            endif

            do mw=1,mwaves
!              if (s(mw,i) .gt. 316.d0) then
!               # shouldn't happen unless h > 10 km!
!                write(6,*) 'speed > 316: i,mw,s(mw,i): ',i,mw,s(mw,i)
!                endif
	           s(mw,i)=dxdc*s(mw,i)
               fwave(1,mw,i)=dxdc*fwave(1,mw,i)
               fwave(2,mw,i)=dxdc*fwave(2,mw,i)
               fwave(3,mw,i)=dxdc*fwave(3,mw,i)
               fwave(4,mw,i)=dxdc*fwave(4,mw,i)
               fwave(5,mw,i)=dxdc*fwave(5,mw,i)
               fwave(6,mw,i)=dxdc*fwave(6,mw,i)
            enddo
         enddo
      endif

!===============================================================================


!============= compute fluctuations=============================================

         do i=2-mbc,mx+mbc
            amdq(1:meqn,i) = 0.d0
            apdq(1:meqn,i) = 0.d0
            do  mw=1,5
               if (efix.and.abs(sfract(mw,i))>0.d0) then
                  amdq(1:5,i) = amdq(1:5,i) + sfract(mw,i)*fwave(1:5,mw,i)
                  apdq(1:5,i) = apdq(1:5,i) + sfract2(mw,i)*fwave(1:5,mw,i)
               else
                  if (s(mw,i) < 0.d0) then
                     amdq(1:5,i) = amdq(1:5,i) + fwave(1:5,mw,i)
                  else if (s(mw,i) > 0.d0) then
                     apdq(1:5,i)  = apdq(1:5,i) + fwave(1:5,mw,i)
                  else
                     amdq(1:5,i) = amdq(1:5,i) + 0.5d0 * fwave(1:5,mw,i)
                     apdq(1:5,i) = apdq(1:5,i) + 0.5d0 * fwave(1:5,mw,i)
                  endif
               endif
            end do
         enddo

         do i=2-mbc,mx+mbc
            if (s(6,i) < 0.d0) then
               amdq(6,i) = amdq(6,i) + fwave(6,6,i)
            else if (s(6,i) > 0.d0) then
               apdq(6,i)  = apdq(6,i) + fwave(6,6,i)
            endif
         enddo

      return
      end subroutine
	  
! =====================================================================
      subroutine find_eigval(hL,hR,dpL,dpR,upL,upR,vpL,vpR,    &
                         uL,uR,vL,vR,g,EVAL,beta1,beta2)
! ======================================================================     
!     Determine eigenvalues and eigenvectors
!-----------------------------------------------------------------------
      use bing_module
      use geoclaw_module, only: drytol => dry_tolerance
      implicit none

      integer i,j,lsup,bubble
      double precision hL,hR,upL,upR,uL,uR,g
      double precision vpL,vpR,vL,vR
      double precision hbar,uhat,temp,uphat,dpL,dpR,vhat,vphat
      double precision A(3,3),EVAL(5)
      double precision beta1,beta2

      INTEGER           INFO
      DOUBLE PRECISION  EVL( 3, 3 ), EVR( 3, 3 ),     &
                        WI( 3 ), WORK( 50 ), WR( 3 )

      g = g/cm_coeff

      uhat  =(sqrt(hL)*uL + sqrt(hR)*uR) &
            /(sqrt(hR)+sqrt(hL)) ! Roe average
      vhat  =(sqrt(hL)*vL + sqrt(hR)*vR) &
            /(sqrt(hR)+sqrt(hL)) ! Roe average

      if (hL<drytol) then
         hbar = hR
         uhat = uR
      elseif (hR<drytol) then
         hbar = hL
         uhat = uL
      else
         hbar =(hR+hL)/2.d0
      endif

      uphat=(upR+upL)/2.d0
      vphat=(vpR+vpL)/2.d0

      ! When uhat and uphat is close to 0, 
      ! numerical error can lead to instability
      if (abs(uphat)<abs(uhat)-0.001d0) then
         uhat=uphat
         !stop
      endif

      if (abs(vphat)<abs(vhat)-0.001d0) then
         vhat=vphat
         !stop
      endif

      ! Make sure uhat is between alpha*uphat and uhat
      uhat=max(abs(uphat)*alpha_1,min(abs(uhat),abs(uphat)))*sign(1.d0,uphat)
      vhat=max(abs(vphat)*alpha_1,min(abs(vhat),abs(vphat)))*sign(1.d0,vphat)

      A(:,:) = 0.d0

      A(1,:) = [ 0.d0,  0.d0, 1.d0 ]
      A(2,:) = [    g, uphat, 0.d0 ]
      A(3,1) = g*hbar + beta1*uphat**2
      A(3,2) = 2.d0*beta1*hbar*uphat + beta2*hbar*uhat
      A(3,3) = beta2*uphat

      do i=1,3
         do j=1,3
         if (abs(A(i,j))<1.d-15) then
            A(i,j)=0.d0
         endif
         enddo
      enddo
 
      ! find eigenvalues and eigenvectors
      call DGEEV('N', 'N', 3, A(1:3,1:3), 3, WR(1:3), WI(1:3), &
                 EVL(1:3,1:3), 1, EVR(1:3,1:3), &
                 3, WORK, 50, INFO )

      if (INFO.ne.0) print*,'error with DGEEV',INFO

      do i=1,3
          if (abs(WI(i))>0.d0) then
              print*,'Warning!!! Non-hyperbolic!!'
              print*,'Real Part of EV=',WR
              print*,'Imaginar Part of EV=',WI
              !print*,g,hbar
              !print*,uL,upL,uR,upR
              !print*,vL,vpL,vR,vpR
              !print*,hL,hR,dpL,dpR
              !print*,uhat,uphat,uhat/uphat
              !print*,vhat,vphat,vhat/vphat
              if (abs(WI(i))>1.d-6) stop
          endif
      enddo

      EVAL=0.d0

      do i=1,3
         EVAL(i) = WR(i)
      end do

      EVAL(4) = uphat

      EVAL(5) = beta2*uphat

      lsup = 3 !lsup is the size of the array to be used

      do while (lsup > 1)
         bubble = 0 !bubble in the greatest element out of order
         do j = 1, (lsup-1)
            if (eval(j) > eval(j+1)) then
               temp = eval(j)
               eval(j) = eval(j+1)
               eval(j+1) = temp
               bubble = j
            endif 
         enddo
         lsup = bubble   
      enddo   

      return
      end subroutine

! =====================================================================
      subroutine find_eigvec(hL,hR,dpL,dpR,upL,upR,vpL,vpR,    &
                         uL,uR,vL,vR,g,EVAL,EVEC,beta1,beta2)
! ======================================================================     
!     Determine eigenvectors
!-----------------------------------------------------------------------
      use bing_module
      use geoclaw_module, only: drytol => dry_tolerance
      implicit none

      integer i,j,lsup,bubble
      double precision hL,hR,upL,upR,uL,uR,g
      double precision vpL,vpR,vL,vR
      double precision hbar,uhat,temp,uphat,dpL,dpR,vhat,vphat
      double precision EVEC(5,5),EVAL(5)
      double precision beta1,beta2

      g = g/cm_coeff

      uhat  =(sqrt(hL)*uL + sqrt(hR)*uR) &
            /(sqrt(hR)+sqrt(hL)) ! Roe average
      vhat  =(sqrt(hL)*vL + sqrt(hR)*vR) &
            /(sqrt(hR)+sqrt(hL)) ! Roe average

      if (hL<drytol.or.hR<drytol) then
         hbar = max(hR,hL)
      else
         hbar =(hR+hL)/2.d0
      endif

      uphat=(upR+upL)/2.d0
      vphat=(vpR+vpL)/2.d0

      ! When uhat and uphat is close to 0, 
      ! numerical error can lead to instability
      if (abs(uphat)<abs(uhat)-0.001d0) then
         uhat=uphat
         !stop
      endif

      if (abs(vphat)<abs(vhat)-0.001d0) then
         vhat=vphat
         !stop
      endif

      ! Make sure uhat is between alpha*uphat and uhat
      uhat=max(abs(uphat)*alpha_1,min(abs(uhat),abs(uphat)))*sign(1.d0,uphat)
      vhat=max(abs(vphat)*alpha_1,min(abs(vhat),abs(vphat)))*sign(1.d0,vphat)

      EVEC=0.d0

      do i = 1,3
         EVEC(1,i) = 1.d0
         if (abs(EVAL(i)-uphat)>0.d0) then
            EVEC(2,i) = g/(EVAL(i)-uphat)
         else
            EVEC(2,i) = 1.d0
         endif
         EVEC(4,i) = EVAL(i)
         EVEC(3,i) = 0.d0 !vphat
         EVEC(5,i) = vhat
      end do

      EVEC(3,4)= 1.d0
      EVEC(5,5)= 1.d0

      return
      end subroutine


