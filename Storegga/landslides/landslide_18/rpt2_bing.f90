! =====================================================
      subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx, &
                     ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
      use geoclaw_module, only: g => grav, tol => dry_tolerance
      use geoclaw_module, only: coordinate_system,earth_radius,deg2rad

      implicit none
!
!     # Riemann solver in the transverse direction 
!     not implemented yet

!-----------------------last modified 1/10/05----------------------

      integer ixy,maxm,meqn,maux,mwaves,mbc,mx,imp

      double precision  ql(meqn,1-mbc:maxm+mbc)
      double precision  qr(meqn,1-mbc:maxm+mbc)
      double precision  asdq(meqn,1-mbc:maxm+mbc)
      double precision  bmasdq(meqn,1-mbc:maxm+mbc)
      double precision  bpasdq(meqn,1-mbc:maxm+mbc)
      double precision  aux1(maux,1-mbc:maxm+mbc)
      double precision  aux2(maux,1-mbc:maxm+mbc)
      double precision  aux3(maux,1-mbc:maxm+mbc)

      double precision  abs_tol
      double precision  hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br
      double precision  uhat,vhat,hhat,roe1,roe3,s1,s2,s3,s1l,s3r
      double precision  delf1,delf2,delf3,dxdcd,dxdcu
      double precision  dxdcm,dxdcp,topo1,topo3,eta

      integer i,m,mw,mu,mv

      return
      end
