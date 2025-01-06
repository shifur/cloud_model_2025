cc D. Posselt 7/28/2008
cc Modified this to be called as a column model from an external driver
cc Stripped out unneccessary common blocks and statistics
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine radrat (twcz,tcont,cosz,tb,qb,rho,p0,p00,pi,dz0,iflag,
     1                   iopcloud,iradave,jj,airsfc)
      include 'dimensions.h'
      include 'radiation.h'
c nb is a number of blocks for calculation
!       integer nbb,nq,lay,nadd,nxi,nw,nx2,nz2,nz4,nz15,nx3,mb1
      integer nbb,nq,nxi,nw,nx2,nz2,nz4,nz15,nx3,mb1
      integer mb2,iflag,npp1
c     parameter(nx=66,ny=10,nz=34,kl2=nz-2,lay=88,nadd=7)
!       parameter (nxf=4000,nyf=75,nzf=37,nt=38640,itt=244) 
!       parameter (lb=2,kb=1)
! ! define decomposition from rmp_switch.h
!       parameter (npes=1,ncol=1,nrow=1)
! ! define partial dimension for computation by decomposition
!       parameter (nx=(nxf-lb*2-1)/ncol+1+lb*2)
!       parameter (ny=1)
!       parameter (nz=nzf)
! ! define partial dimension for fft by transpose decomposition
!       parameter (nyc=1)
!       parameter (nyr= ny)
!       parameter (nxc= nx)
!       parameter (nxr=(nxf-lb*2-1)/nrow+1+lb*2)
!       parameter (nzc=(nzf-kb*2-1)/ncol+1+kb*2)
!       parameter (nzr=(nzf-kb*2-1)/nrow+1+kb*2)
!
!        parameter (nadd=7,lay=88)
      integer nx1,ny1, kles,kl2,kmax, np
!       parameter(kmax=nzf,kles=kmax-1,kl2=kles-1)
c     parameter(nx1=nx-4,iles=nx-1)
c     parameter(nx1=nxf-4,ny1=nyf-4,iles=nxf-1)   ! cccshie 8/23/04
c     parameter(nx1=nx-4,ny1=ny-4,iles=nxf-1)   ! cccshie 8/23/04
!       parameter(iles=nx-1)   ! cccshie 9/15/04
       parameter (nx1=1,ny1=1)    ! cccshie 9/15/04
      parameter(nq=nz+nadd-1,nw=nq-1, np=nz-2+nadd)
cccshie 8/9/04
      integer jj2
      real rflux(nx1,ny1,8)
      common/radflux/rflux
c     m=nx1
c     nw=np
c     npp1=nq
c
      integer iradave,iopcloud
      real    twcz,tcont
      real    p00(nz),dz0(nz),tb(nz),qb(nz),rho(nz),p0(nz),pi(nz)
      dimension airsfc(nx,ny)
      common/radtemp1/ qc(nx,nz),qr(nx,nz),qi(nx,nz),qs(nx,nz),
     1 qg(nx,nz),pt(nx,nz),qv(nx,nz),qh(nx,nz),qsw(nx,nz),qgw(nx,nz),
     2 qhw(nx,nz)
      common/radtemp2/ rsw(nx,nz),rlw(nx,nz)
      common/radtemp3/ tsfc_a(nx),tsst(nx)
!       common/surface1/ sun_4(nx,ny,4)
c-----------------------------------------------------------------------
c local variables
c-----------------------------------------------------------------------
       real    dfdts(nx1,1,nq),st4(nx1,1),tsq(nx1,1),tb_r(nx1,1)
c     tb_r:  sfc air temperature - dan please include this
c     tsq: sfc temperature
c
       real    rsirbm(nx1,1),rsirdf(nx1,1),rsuvbm(nx1,1),
     $         rsuvdf(nx1,1),emiss(nx1,10)
       common/albedo1/ rsirbm,rsirdf,rsuvbm,rsuvdf,emiss
       common/var6/ rldown(nx,ny),swdown(nx,ny),albedo(nx,ny)
       integer ict,icb
       common/cloudp/ ict,icb
c
       real   coolr1(nx,1,nw),heatr1(nx,1,nw)
       real   tausw1(nx1,1,nw,3),reff(nx1,1,nw,3),cwc(nx1,1,nw,3),
     1        taucld(nx1,1,np,3),tauir1(nx1,1,nw,3)
c
       real    cosz,pl(lay),pa(lay)
       real    plq(nx1,1,nq),taq(nx1,1,nw),waq(nx1,1,nw),oaq(nx1,1,nw)
       real    cosz1(nx1,1)
       real    cmc,n2o,ch4,cfc11,cfc12,cfc22
       real    aco,sc0
c      cmc=co2
       real    ta(lay),wa(lay),oa(lay),plu(nadd)
       real    fcld(nx1,1,nw),flx(nx1,1,nq)
       real    flx1(nx1,1,nq),flc(nx1,1,nq),flc1(nx1,1,nq)
       real    fdirir(nx1,1),fdifir(nx1,1)
       real    fdirpar(nx1,1),fdifpar(nx1,1)
       real    fdiruv(nx1,1),fdifuv(nx1,1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      real cwc(m,np,3),taucld(m,np,3),reff(m,np,3)
c
c  cloud water mixing ratio (cwc)                  gm/gm      m*np*3
c        index 1 for ice particles
c        index 2 for liquid drops
c        index 3 for rain drops
c  cloud optical thickness (taucld)                 n/d       m*np*3
c        index 1 for ice particles
c        index 2 for liquid drops
c        index 3 for rain drops
c  effective cloud-particle size (reff)          micrometer   m*np*3
c        index 1 for ice particles
c        index 2 for liquid drops
c        index 3 for rain drops
c  cloud amount (fcld)                            fraction    m*np
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       logical high,overcast,cldwater,trace
       data high/.true./
c     if high is true: more accurate cooling/heating between 0.01-24 mb
c     however, it will take more computation than if it false
c
       data overcast/.true./
c
       data cldwater/.false./
c
       data trace/.false./
c
c   if trace = .true., absorption due to n2o, ch4, cfcs, and the
c   two minor co2 bands in the window region is included.
c   if trace = .false., absorption in minor bands is neglected.
c
c   overcast=.true.  ! no partial cloudiness
c
c   fcld is 1 or 0
c
c   cldwater=.false.  ! optical thickness is specified
c   cldwater=.true .  ! optical thickness is calculated by md chou's
c                       radiation codes
c
c     choose nadd= 1, 3, 4, 7
c
c     data plu/ 0.01                                       /
c      data plu/ 0.01,             4.04,               22.46/
c     data plu/ 0.01      , 3.04      , 7.18        , 22.46/
      data plu/ 0.01, 0.56, 3.04, 4.04, 7.18,  12.38, 22.46/
       integer i,k,irst,km,ib,i1,k1,ihalfp
       data irst/0/
       save

c D. Posselt Set values of some necessary variables
      kmax = nz
      kles = nz - 1
      kl2 = nz - 2

       do k=1,nq
       do i=1,nx1
         flx(i,1,k)=0.0
         flx1(i,1,k)=0.0
         flc(i,1,k)=0.0
         flc1(i,1,k)=0.0
c
         dfdts(i,1,k)=0.0
c
       enddo
       enddo
       do k=1,nw
          do i=1,nx1
             fcld(i,1,k)=0.
             tausw1(i,1,k,1)=0.
             tausw1(i,1,k,2)=0.
             tausw1(i,1,k,3)=0.
             reff(i,1,k,1)=0.
             reff(i,1,k,2)=0.
             reff(i,1,k,3)=0.
             cwc(i,1,k,1)=0.
             cwc(i,1,k,2)=0.
             cwc(i,1,k,3)=0.
             tauir1(i,1,k,1)=0.
             tauir1(i,1,k,2)=0.
             tauir1(i,1,k,3)=0.
          enddo
       enddo
       do k=1,np
          do i=1,nx1
             taucld(i,1,k,1)=0.
             taucld(i,1,k,2)=0.
             taucld(i,1,k,3)=0.
          enddo
       enddo
       do k=1,nw
       do i=1,nx
          coolr1(i,1,k)=0.0
          heatr1(i,1,k)=0.0
       enddo
       enddo
       do i=1,nx1
         st4(i,1)=0.0
         fdirir(i,1)=0.0
         fdifir(i,1)=0.0
         fdirpar(i,1)=0.0
         fdifpar(i,1)=0.0
         fdiruv(i,1)=0.0
         fdifuv(i,1)=0.0
       enddo
       do i=1,lay
          pa(i)=0.0
       enddo
c
       if (irst.ne.0) go to 500
c
c     np=number of atmospheric layers; npp1=surface level
       npp1=np+1
c     solar constant and cosine of solar zenith angle
c     solar constant and surface albedo (0.07 water; 0.19 land)
         sc0=1378.
         rsfcir=0.07
         rsfcuv=0.07
         do i=1,nx1
           rsirbm(i,1)=rsfcir
           rsirdf(i,1)=rsfcir
           rsuvbm(i,1)=rsfcuv
           rsuvdf(i,1)=rsfcuv
         enddo
c
         ict=24
         icb=32
c
c  assign co2 (cmc) and others. units are parts/part
       cmc=300.e-6
       n2o=0.
       ch4=0.
       cfc11=0.
       cfc12=0.
       cfc22=0.
!        write(6,*) 'npp1(3)=',npp1
c      write(6,*) 'npp1(3)=',npp1,'  nadd = ',nadd
c
c    ocean, emiss=0.97 for all bands. the land emissivity in the lw
c    varies significantly depending upon the suface type. there is no
c    "standard" value for land emissivity. i would "empirically" set it
c    to be 0.9 for all bands over land.
c
       do ia=1,10
       do i=1,nx1
         emiss(i,ia)=0.97
       enddo
       enddo
       do k=1,kl2+1
         km=kl2+3-k
         pl(k+nadd)=1.e-3*p00(km)
       end do
       do k=1,nadd
         pl(k)=plu(k)
       end do
       do k=1,kl2
         km=kl2+2-k
         pa(k+nadd)=1.e-3*p0(km)
       end do
      do 13 k=1,nadd
         pa(k)=0.5*(pl(k)+pl(k+1))
         wa(k)=twcz
   13 continue
c
       call fito3 (npp1,pa,ta,oa)
!        print*
!        print*,' k       pa        pl       pl        ao        ta
!      *   wa'
!        write (6,975) (k,pa(k),pl(k),pl(k),oa(k),ta(k),wa(k),k=1,npp1)
!   975 format(i3,3f10.3,e12.4,f12.5,f12.7)
       irst=1
c
cc
c
  500 continue
         rsfcir=0.07
         rsfcuv=0.07
         do i=1,nx1
           rsirbm(i,1)=rsfcir
           rsirdf(i,1)=rsfcir
           rsuvbm(i,1)=rsfcuv
           rsuvdf(i,1)=rsfcuv
         enddo
c
c    ocean, emiss=0.97 for all bands. the land emissivity in the lw
c    varies significantly depending upon the suface type. there is no
c    "standard" value for land emissivity. i would "empirically" set it
c    to be 0.9 for all bands over land.
c
       do ia=1,10
       do i=1,nx1
         emiss(i,ia)=0.97
       enddo
       enddo
c
c
      do i=1,nx1
        do k=1,npp1
          plq(i,1,k)=pl(k)
        enddo
        do k=1,nw
          oaq(i,1,k)=oa(k)
          fcld(i,1,k)=0.0
        enddo
        do k=1,nadd
          waq(i,1,k)=wa(k)
          taq(i,1,k)=ta(k)
        enddo
        do k=2,kles
           km=kmax-k+nadd
          waq(i,1,km)=qv(i,k)
          taq(i,1,km)=pt(i,k)
        enddo
      enddo
        do i=1,nx1
           tb_r(i,1)=tsfc_a(i)
           tsq(i,1)=tsst(i)
           cosz1(i,1)=cosz
        enddo
! D.Posselt 07/23/07 Include j index in opt4 for RAMS microphysics
         call opt4 (iopcloud,iradave,nx1,twcz,rho,dz0,
     +              tausw1,tauir1,fcld,taq,waq,oaq,plq,cwc,reff,jj)
cccshie 8/19/04
c       do j=1,ny1
        jj2=jj
        do i=1,nx1
          rflux(i,jj2,1)=0.
          rflux(i,jj2,2)=0.
          rflux(i,jj2,3)=0.
          rflux(i,jj2,4)=0.
          rflux(i,jj2,5)=0.
          rflux(i,jj2,6)=0.
          rflux(i,jj2,7)=0.
          rflux(i,jj2,8)=0.
        enddo
c       enddo
c ------------------------------------------------------------
         if (cosz.ge.0.005) then
         if (cldwater) then
           do k=1,np
           do i=1,nx1
             taucld(i,1,k,1)=0.
             taucld(i,1,k,2)=0.
             taucld(i,1,k,3)=0.
           enddo
           enddo
         else
           do k=1,np
           do i=1,nx1
             taucld(i,1,k,1)=tausw1(i,1,k,1)
             taucld(i,1,k,2)=tausw1(i,1,k,2)
             taucld(i,1,k,3)=tausw1(i,1,k,3)
           enddo
           enddo
         endif
          do i=1,nx1
             cosz1(i,1)=cosz
          enddo
cccshie 8/19/04
c           call sorad (nx1,np,plq,taq,waq,oaq,cmc,
            call sorad (nx1,np,jj2,plq,taq,waq,oaq,cmc,
     *                  overcast,cldwater,cwc,taucld,reff,fcld,ict,icb,
     *                  cosz1,rsuvbm,rsuvdf,rsirbm,rsirdf,
     *                  flx1,flc1,fdiruv,fdifuv,fdirpar,fdifpar,
     *                  fdirir,fdifir)
         end if
         if (cldwater) then
           do k=1,np
           do i=1,nx1
c
             taucld(i,1,k,1)=0.
             taucld(i,1,k,2)=0.
             taucld(i,1,k,3)=0.
c
           enddo
           enddo
         else
           do k=1,np
           do i=1,nx1
             taucld(i,1,k,1)=tauir1(i,1,k,1)
             taucld(i,1,k,2)=tauir1(i,1,k,2)
             taucld(i,1,k,3)=tauir1(i,1,k,3)
           enddo
           enddo
         endif
cccshie 8/19/04
c         call irrad (nx1,np,plq,taq,waq,oaq,tb_r,  tsq  ,cmc,
          call irrad (nx1,np,jj2,plq,taq,waq,oaq,tb_r,  tsq  ,cmc,
c         call irrad (nx1,ny1,np,jj2,plq,taq,waq,oaq,tb_r,  tsq  ,cmc,
     *               n2o,ch4,cfc11,cfc12,cfc22,emiss,
     *               overcast,cldwater,cwc,taucld,reff,fcld,ict,icb,
     *               high,trace,flx,flc,dfdts,st4)
c
c since flx1 is normalized, they should be multiplied
c by sc0*cosz
c since upward flux should be positive, heatr is multiplied by
c a minus sign
c save solar and long-wave radiative fluxes
c
c
c
         do k=1,nw
           do i=1,nx1
              heatr1(i,1,k)=(flx1(i,1,k+1)-flx1(i,1,k))*8.441874/
     1                       (plq(i,1,k+1)-plq(i,1,k))
              heatr1(i,1,k)=-heatr1(i,1,k)*sc0*cosz1(i,1)
              coolr1(i,1,k)=(flx(i,1,k+1)-flx(i,1,k))*8.441874/
     1                       (plq(i,1,k+1)-plq(i,1,k))
c              print*,'flx',k,i,flx(i,1,k),flx(i,1,k+1),plq(i,1,k),
c     1                     plq(i,1,k+1)
         enddo
         enddo
c
         do k1=nadd+1,np
           k=np+2-k1
           do i=1,nx1
              aco=1./(86400.*pi(k))
             rsw(i,k)=aco*heatr1(i,1,k1)
             rlw(i,k)=aco*coolr1(i,1,k1)
           enddo
         enddo
cccshie 8/19/04
       do i=1,nx1
       rflux(i,jj2,1)=rflux(i,jj2,1)*sc0*cosz1(i,1)  ! downward SW surface (+:down)
       rflux(i,jj2,3)=rflux(i,jj2,3)*sc0*cosz1(i,1)  ! downward SW TOA     (+:down)
       rflux(i,jj2,4)=flx1(i,1,1)*sc0*cosz1(i,1)-rflux(i,jj2,3) ! upward SW TOA (-:up) : (net SW flux at TOA) - downward SW flux TOA
       rflux(i,jj2,6)=flx1(i,1,npp1)*sc0*cosz1(i,1)-rflux(i,jj2,1) ! upward SW surface (-:up) : (net SW flux at surface) - downward SW flux surface
       enddo

1000   continue
cc
c -----------------------------------------------------------------------
!        if(iflag.eq.2) then
!          ihalfp=(iles-1)/2
!          write(6,*) 'cosz=',cosz
!          write(6,*) 'check rsw(i,k) and rlw(i,k) at central point'
!          do k1=nadd+1,np
!            k=np+2-k1
!            write(6,7821) k,pl(k),rsw(ihalfp,k)*86400.,
!      1                           rlw(ihalfp,k)*86400.,
!      1        heatr1(ihalfp,1,k1),coolr1(ihalfp,1,k1)
!          enddo
!        endif
! 7821  format(2x,i6,2x,f10.3,2x,4e20.10)
      return
      end
