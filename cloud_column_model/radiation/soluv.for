c******************************************************************
cccshie 8/19/04
c     subroutine soluv (m,np,wh,oh,dp,overcast,cldwater,
      subroutine soluv (m,np,jj2,wh,oh,dp,overcast,cldwater,
     *          cwp,taucld,reff,ict,icb,fcld,cosz,
     *          taual,ssaal,asyal,rsuvbm,rsuvdf,
     *          flx,flc,fdiruv,fdifuv,fdirpar,fdifpar)
      implicit none
      include 'radiation.h'
      include 'dimensions.h'
c******************************************************************
c  compute solar fluxes in the uv+par region. the spectrum is
c  grouped into 8 bands:
c  
c              band     micrometer
c
c       uv-c    1.     .175 - .225
c               2.     .225 - .245
c                      .260 - .280
c               3.     .245 - .260
c
c       uv-b    4.     .280 - .295
c               5.     .295 - .310
c               6.     .310 - .320
c      
c       uv-a    7.     .320 - .400
c      
c       par     8.     .400 - .700
c
c----- input parameters:                            units      size
c
c  number of soundings (m)                          n/d         1
c  number of atmospheric layers (np)                n/d         1
c  layer scaled-water vapor content (wh)          gm/cm^2      m*np
c  layer ozone content (oh)                      (cm-atm)stp   m*np
c  layer pressure thickness (dp)                    mb         m*np
c  option for scaling cloud optical thickness       n/d         1
c        overcast="true" if scaling is not required
c        overcast="fasle" if scaling is required
c  input option for cloud optical thickness         n/d         1
c        cldwater="true" if taucld is provided
c        cldwater="false" if cwp is provided
c  cloud water amount (cwp)                        gm/m**2     m*np*3
c        index 1 for ice particles
c        index 2 for liquid drops
c        index 3 for rain drops
c  cloud optical thickness (taucld)                 n/d        m*np*3
c       index 1 for ice particles
c       index 2 for liquid drops
c       index 3 for rain drops
c  effective cloud-particle size (reff)          micrometer    m*np*3
c       index 1 for ice paticles
c       index 2 for liquid drops
c       index 3 for rain drops
c  level index separating high and                  n/d        m
c       middle clouds (ict)
c  level indiex separating middle and               n/d        m
c       low clouds (icb)
c  cloud amount (fcld)                            fraction     m*np
c  cosine of solar zenith angle (cosz)              n/d        m
c  aerosol optical thickness (taual)                n/d        m*np*11
c  aerosol single-scattering albedo (ssaal)         n/d        m*np*11
c  aerosol asymmetry factor (asyal)                 n/d        m*np*11
c  uv+par surface albedo for beam                 fraction     m
c       radiation (rsuvbm)
c  uv+par surface albedo for diffuse              fraction     m
c       radiation (rsuvdf)
c
c---- temporary array
c
c  scaled cloud optical thickness                   n/d        m*np
c       for beam radiation (tauclb)
c  scaled cloud optical thickness                   n/d        m*np
c       for diffuse radiation  (tauclf)     
c
c----- output (updated) parameters:
c
c  all-sky net downward flux (flx)               fraction      m*(np+1)
c  clear-sky net downward flux (flc)             fraction      m*(np+1)
c  all-sky direct downward uv flux at
c       the surface (fdiruv)                     fraction      m
c  all-sky diffuse downward uv flux at
c       the surface (fdifuv)                     fraction      m
c  all-sky direct downward par flux at
c       the surface (fdirpar)                    fraction      m
c  all-sky diffuse downward par flux at
c       the surface (fdifpar)                    fraction      m
c
c***********************************************************************
cccshie 8/19/04
c     implicit none
!       parameter (nxf=4000,nyf=75,nzf=37,nt=38640,itt=244) 
!       parameter (lb=2,kb=1)
! define decomposition from rmp_switch.h
!       parameter (npes=1,ncol=1,nrow=1)
! define partial dimension for computation by decomposition
!       parameter (nx=(nxf-lb*2-1)/ncol+1+lb*2)
!       parameter (ny=1)
!       parameter (nz=nzf)
! define partial dimension for fft by transpose decomposition
!       parameter (nyc=1)
!       parameter (nyr= ny)
!       parameter (nxc= nx)
!       parameter (nxr=(nxf-lb*2-1)/nrow+1+lb*2)
!       parameter (nzc=(nzf-kb*2-1)/ncol+1+kb*2)
!       parameter (nzr=(nzf-kb*2-1)/nrow+1+kb*2)
!
!        parameter (nadd=7,lay=88)
      integer nx1, ny1
      parameter (nx1=1,ny1=1)    ! cccshie 9/15/04
c-----input parameters
c     integer m,np,ict,icb
      integer m,np,jj2,ict,icb
      real taucld(m,np,3),reff(m,np,3),fcld(m,np)
      real cwp(m,np,3),wh(m,np),oh(m,np),dp(m,np)
      real taual(m,np,11),ssaal(m,np,11),asyal(m,np,11)
      real rsuvbm(m),rsuvdf(m),cosz(m)
      logical overcast,cldwater
c-----output (updated) parameter
      real flx(m,np+1),flc(m,np+1)
      real fdiruv (m),fdifuv (m)
      real fdirpar(m),fdifpar(m)
c-----static parameters
      integer nband
      parameter (nband=8)
      real hk(nband),wk(nband),zk(nband),ry(nband)
      real aig(3),awg(3),arg(3)
      real aib(2),awb(2),arb(2)
c-----temporary array
      integer i,k,ib
      integer ih1,ih2,im1,im2,is1,is2
!     real dsm(m)
!     real tauclb(m,np),tauclf(m,np),asycl(m,np)
      real taurs,tauoz,tauwv
!     real tausto(m,np),ssatau(m,np),asysto(m,np)
!     real tautob(m,np),ssatob(m,np),asytob(m,np)
!     real tautof(m,np),ssatof(m,np),asytof(m,np)
      real taux,reff1,reff2,g1,g2,g3
!     real rr(m,np+1,2),tt(m,np+1,2),td(m,np+1,2),
!    *     rs(m,np+1,2),ts(m,np+1,2)
!     real fall(m,np+1),fclr(m,np+1),fsdir(m),fsdif(m)
!     real asyclt(m),cc(m,3)
!     real rrt(m,np),ttt(m,np),tdt(m,np),rst(m,np),tst(m,np)
!     real dum1(m,np+1),dum2(m),dum3(m),dum(m,np)
      real   ,  allocatable :: dsm(:)
      real   ,  allocatable :: tauclb(:,:)
      real   ,  allocatable :: tauclf(:,:)
      real   ,  allocatable :: asycl(:,:)
      real   ,  allocatable :: tausto(:,:)
      real   ,  allocatable :: ssatau(:,:)
      real   ,  allocatable :: asysto(:,:)
      real   ,  allocatable :: tautob(:,:)
      real   ,  allocatable :: ssatob(:,:)
      real   ,  allocatable :: asytob(:,:)
      real   ,  allocatable :: tautof(:,:)
      real   ,  allocatable :: ssatof(:,:)
      real   ,  allocatable :: asytof(:,:)
      real   ,  allocatable :: rr(:,:,:)
      real   ,  allocatable :: tt(:,:,:)
      real   ,  allocatable :: td(:,:,:)
      real   ,  allocatable :: rs(:,:,:)
      real   ,  allocatable :: ts(:,:,:)
      real   ,  allocatable :: fall(:,:)
      real   ,  allocatable :: fclr(:,:)
      real   ,  allocatable :: fsdir(:)
      real   ,  allocatable :: fsdif(:)
      real   ,  allocatable :: asyclt(:)
      real   ,  allocatable :: cc(:,:)
      real   ,  allocatable :: rrt(:,:)
      real   ,  allocatable :: ttt(:,:)
      real   ,  allocatable :: tdt(:,:)
      real   ,  allocatable :: rst(:,:)
      real   ,  allocatable :: tst(:,:)
      real   ,  allocatable :: dum1(:,:)
      real   ,  allocatable :: dum2(:)
      real   ,  allocatable :: dum3(:)
      real   ,  allocatable :: dum(:,:)
ccshie 8/19/04
      real rflux(nx1,ny1,8)
      common/radflux/rflux
      real    dwflux(nx1,np+1)
      real    upflux(nx1,np+1)
c-----hk is the fractional extra-terrestrial solar flux in each
c     of the 8 bands. the sum of hk is 0.47074. (table 3)
      data hk/.00057, .00367, .00083, .00417,
     *        .00600, .00556, .05913, .39081/
c-----zk is the ozone absorption coefficient. unit: /(cm-atm)stp
c     (table 3)
      data zk /30.47, 187.2,  301.9,   42.83,
     *         7.09,  1.25,   0.0345,  0.0572/
c-----wk is the water vapor absorption coefficient. unit: cm**2/g
c     (table 3)
      data wk /7*0.0, 0.00075/
c-----ry is the extinction coefficient for rayleigh scattering.
c     unit: /mb. (table 3)
      data ry /.00604, .00170, .00222, .00132,
     *         .00107, .00091, .00055, .00012/ 
c-----coefficients for computing the extinction coefficients of ice, 
c     water, and rain particles, independent of spectral band. (table 4)
      data aib/ 3.33e-4,2.52/
      data awb/-6.59e-3,1.65/
      data arb/ 3.07e-3,0.00/
c-----coefficients for computing the asymmetry factor of ice, water,
c     and rain particles, independent of spectral band. (table 6)
      data aig/.74625,.0010541,-.00000264/
      data awg/.82562,.0052900,-.00014866/
      data arg/.883,0.0,0.0/
      allocate(dsm(m))
      allocate(tauclb(m,np))
      allocate(tauclf(m,np))
      allocate(asycl(m,np))
      allocate(tausto(m,np))
      allocate(ssatau(m,np))
      allocate(asysto(m,np))
      allocate(tautob(m,np))
      allocate(ssatob(m,np))
      allocate(asytob(m,np))
      allocate(tautof(m,np))
      allocate(ssatof(m,np))
      allocate(asytof(m,np))
      allocate(rr(m,np+1,2))
      allocate(tt(m,np+1,2))
      allocate(td(m,np+1,2)) 
      allocate(rs(m,np+1,2))
      allocate(ts(m,np+1,2))
      allocate(fall(m,np+1))
      allocate(fclr(m,np+1))
      allocate(fsdir(m))
      allocate(fsdif(m))
      allocate(asyclt(m))
      allocate(cc(m,3))
      allocate(rrt(m,np))
      allocate(ttt(m,np))
      allocate(tdt(m,np))
      allocate(rst(m,np))
      allocate(tst(m,np))
      allocate(dum1(m,np+1))
      allocate(dum2(m))
      allocate(dum3(m))
      allocate(dum(m,np))
c-----initialize fdiruv, fdifuv, surface reflectances and transmittances.
c     the reflectance and transmittance of the clear and cloudy portions
c     of a layer are denoted by 1 and 2, respectively.
c     cc is the maximum cloud cover in each of the high, middle, and low
c     cloud groups.
c     1/dsm=1/cos(53) = 1.66
      do i=1,m                    
         dsm(i)=0.602
         fdiruv(i)=0.0
         fdifuv(i)=0.0
         rr(i,np+1,1)=rsuvbm(i)
         rr(i,np+1,2)=rsuvbm(i)
         rs(i,np+1,1)=rsuvdf(i)
         rs(i,np+1,2)=rsuvdf(i)
         td(i,np+1,1)=0.0
         td(i,np+1,2)=0.0
         tt(i,np+1,1)=0.0
         tt(i,np+1,2)=0.0
         ts(i,np+1,1)=0.0
         ts(i,np+1,2)=0.0
         cc(i,1)=0.0
         cc(i,2)=0.0
         cc(i,3)=0.0
      enddo
c-----compute cloud optical thickness.  eqs. (4.6) and (4.11)
       do k=1,np
        do i=1,m
cc          reff(i,k,1)=max(reff(i,k,1),1.0)
          reff(i,k,1)=min(reff(i,k,1),150.0)
cc          reff(i,k,2)=max(reff(i,k,2),4.0)
          reff(i,k,2)=min(reff(i,k,2),20.0)
        enddo
       enddo
      if (cldwater) then
       do k=1,np
        do i=1,m
          taucld(i,k,1)=cwp(i,k,1)*(aib(1)+aib(2)/reff(i,k,1))
          taucld(i,k,2)=cwp(i,k,2)*(awb(1)+awb(2)/reff(i,k,2))
          taucld(i,k,3)=cwp(i,k,3)* arb(1)
        enddo
       enddo
      endif
c-----options for scaling cloud optical thickness
      if (overcast) then
       do k=1,np
        do i=1,m
          tauclb(i,k)=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)
          tauclf(i,k)=tauclb(i,k)
        enddo
       enddo
       do k=1,3
        do i=1,m
           cc(i,k)=1.0
        enddo
       enddo
      else
c-----scale cloud optical thickness in each layer from taucld (with
c     cloud amount fcld) to tauclb and tauclf (with cloud amount cc).
c     tauclb is the scaled optical thickness for beam radiation and
c     tauclf is for diffuse radiation (see section 7).
         call cldscale(m,np,cosz,fcld,taucld,ict,icb,
     *                 cc,tauclb,tauclf)
      endif
c-----cloud asymmetry factor for a mixture of liquid and ice particles.
c     unit of reff is micrometers. eqs. (4.8) and (6.4)
      do k=1,np
       do i=1,m
           asyclt(i)=1.0
           taux=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)
         if (taux.gt.0.02 .and. fcld(i,k).gt.0.01) then
c           reff1=reff(i,k,1)
c           reff2=reff(i,k,2)
            reff1=min(reff(i,k,1),150.)
            reff2=min(reff(i,k,2),20.0)
           g1=(aig(1)+(aig(2)+aig(3)*reff1)*reff1)*taucld(i,k,1)
           g2=(awg(1)+(awg(2)+awg(3)*reff2)*reff2)*taucld(i,k,2)
           g3= arg(1)*taucld(i,k,3)
           asyclt(i)=(g1+g2+g3)/taux
         endif
       enddo
         do i=1,m
           asycl(i,k)=asyclt(i)
         enddo
      enddo
c-----integration over spectral bands
      do 100 ib=1,nband
       do k=1,np
        do i=1,m
c-----compute rayleigh, ozone and water vapor optical thicknesses
          taurs=ry(ib)*dp(i,k)
          tauoz=zk(ib)*oh(i,k)
          tauwv=wk(ib)*wh(i,k)
c-----compute clear-sky optical thickness, single scattering albedo,
c     and asymmetry factor (eqs. 6.2-6.4)
          tausto(i,k)=taurs+tauoz+tauwv+taual(i,k,ib)+1.0e-8
          ssatau(i,k)=ssaal(i,k,ib)*taual(i,k,ib)+taurs
          asysto(i,k)=asyal(i,k,ib)*ssaal(i,k,ib)*taual(i,k,ib)
c-----compute reflectance and transmittance of the clear portion of a layer
          tautob(i,k)=tausto(i,k)
          ssatob(i,k)=ssatau(i,k)/tautob(i,k)+1.0e-8
          ssatob(i,k)=min(ssatob(i,k),0.999999)
          asytob(i,k)=asysto(i,k)/(ssatob(i,k)*tautob(i,k))
        enddo
       enddo
c-----for direct incident radiation
         call deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)
c-----diffuse incident radiation is approximated by beam radiation with
c     an incident angle of 53 degrees, eqs. (6.5) and (6.6)
         call deledd (m,np,tautob,ssatob,asytob,dsm,rst,tst,dum)
       do k=1,np
        do i=1,m
           rr(i,k,1)=rrt(i,k)
           tt(i,k,1)=ttt(i,k)
           td(i,k,1)=tdt(i,k)
           rs(i,k,1)=rst(i,k)
           ts(i,k,1)=tst(i,k)
        enddo
       enddo
c-----compute reflectance and transmittance of the cloudy portion of a layer
       do k=1,np
        do i=1,m
c-----for direct incident radiation
c     the effective layer optical properties. eqs. (6.2)-(6.4)
           tautob(i,k)=tausto(i,k)+tauclb(i,k)
           ssatob(i,k)=(ssatau(i,k)+tauclb(i,k))/tautob(i,k)+1.0e-8
           ssatob(i,k)=min(ssatob(i,k),0.999999)
           asytob(i,k)=(asysto(i,k)+asycl(i,k)*tauclb(i,k))
     *                /(ssatob(i,k)*tautob(i,k))
c-----for diffuse incident radiation
           tautof(i,k)=tausto(i,k)+tauclf(i,k)
           ssatof(i,k)=(ssatau(i,k)+tauclf(i,k))/tautof(i,k)+1.0e-8
           ssatof(i,k)=min(ssatof(i,k),0.999999)
           asytof(i,k)=(asysto(i,k)+asycl(i,k)*tauclf(i,k))
     *                /(ssatof(i,k)*tautof(i,k))
        enddo
       enddo
c-----for direct incident radiation
c     note that the cloud optical thickness is scaled differently for direct
c     and diffuse insolation, eqs. (7.3) and (7.4).
         call deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)
c-----diffuse incident radiation is approximated by beam radiation with
c     an incident angle of 53 degrees, eqs. (6.5) and (6.6)
         call deledd (m,np,tautof,ssatof,asytof,dsm,rst,tst,dum)
       do k=1,np
        do i=1,m
           rr(i,k,2)=rrt(i,k)
           tt(i,k,2)=ttt(i,k)
           td(i,k,2)=tdt(i,k)
           rs(i,k,2)=rst(i,k)
           ts(i,k,2)=tst(i,k)
        enddo
       enddo
c-----flux calculations
c     initialize clear-sky flux (fclr), all-sky flux (fall), 
c     and surface downward fluxes (fsdir and fsdif)
        do k=1,np+1
         do i=1,m
           fclr(i,k)=0.0
           fall(i,k)=0.0
         enddo
        enddo
        do i=1,m
           fsdir(i)=0.0
           fsdif(i)=0.0
        enddo
        if (overcast) then
              ih1=1
              ih2=1
              im1=1
              im2=1
              is1=1
              is2=1
c-----for clear-sky fluxes only
           do k=1,np+1
            do i=1,m
              dum1(i,k)=0.0
            enddo
           enddo
           do i=1,m
              dum2(i)=0.
              dum3(i)=0.
           enddo
cccshie 8/19/04 dwflux,upflux, for cloudy-sky or all-sky only,
ccc              so it is dummy here as for clear-sky.
           do k=1,np+1
            do i=1,m
              dwflux(i,k)=0.0
              upflux(i,k)=0.0
            enddo
           enddo
         call cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *                cc,rr,tt,td,rs,ts,fclr,dum1,dum2,dum3,
     *                dwflux,upflux)
c    *                cc,rr,tt,td,rs,ts,fclr,dum1,dum2,dum3)
cccshie 8/19/04 dwflux,upflux, for cloudy-sky when overcast=true
           do k=1,np+1
            do i=1,m
              dwflux(i,k)=0.0
              upflux(i,k)=0.0
            enddo
           enddo
           do k=1,np+1
            do i=1,m
              fall(i,k)=0.0
            enddo
           enddo
           do i=1,m
              fsdir(i)=0.0
              fsdif(i)=0.0
           enddo
              ih1=2
              ih2=2
              im1=2
              im2=2
              is1=2
              is2=2
c-----for cloudy-sky fluxes only
         call cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *                cc,rr,tt,td,rs,ts,dum1,fall,fsdir,fsdif,
     *                dwflux,upflux)
c    *                cc,rr,tt,td,rs,ts,dum1,fall,fsdir,fsdif)
        else
              ih1=1
              ih2=2
              im1=1
              im2=2
              is1=1
              is2=2
c-----for clear- and all-sky fluxes
c     the all-sky flux, fall is the summation inside the brackets
c     of eq. (7.11)
cccshie 8/19/04 dwflux,upflux, for all-sky when overcast=false
           do k=1,np+1
            do i=1,m
              dwflux(i,k)=0.0
              upflux(i,k)=0.0
            enddo
           enddo
         call cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *                cc,rr,tt,td,rs,ts,fclr,fall,fsdir,fsdif,
     *                dwflux,upflux)
c    *                cc,rr,tt,td,rs,ts,fclr,fall,fsdir,fsdif)
        endif
c-----flux integration, eq. (6.1)
       do k=1,np+1
        do i=1,m
          flx(i,k)=flx(i,k)+fall(i,k)*hk(ib)
cccshie 8/19/04
         if(k.eq.1)then
c       rflux(i,jj2,3)=rflux(i,jj2,3)+hk(ib) ! djohnson way in old radiation code, may do it later
        rflux(i,jj2,3)=rflux(i,jj2,3)+dwflux(i,k)*hk(ib)
         endif
         if(k.eq.(np+1))then
        rflux(i,jj2,1)=rflux(i,jj2,1)+dwflux(i,k)*hk(ib)
         endif
        enddo
        do i=1,m
          flc(i,k)=flc(i,k)+fclr(i,k)*hk(ib)
        enddo
       enddo
c-----compute direct and diffuse downward surface fluxes in the uv
c     and par regions
       if(ib.lt.8) then
        do i=1,m
          fdiruv(i)=fdiruv(i)+fsdir(i)*hk(ib)
          fdifuv(i)=fdifuv(i)+fsdif(i)*hk(ib)
        enddo
       else
        do i=1,m
          fdirpar(i)=fsdir(i)*hk(ib)
          fdifpar(i)=fsdif(i)*hk(ib)
        enddo
       endif
 100  continue
      deallocate(dsm)
      deallocate(tauclb)
      deallocate(tauclf)
      deallocate(asycl)
      deallocate(tausto)
      deallocate(ssatau)
      deallocate(asysto)
      deallocate(tautob)
      deallocate(ssatob)
      deallocate(asytob)
      deallocate(tautof)
      deallocate(ssatof)
      deallocate(asytof)
      deallocate(rr)
      deallocate(tt)
      deallocate(td)
      deallocate(rs)
      deallocate(ts)
      deallocate(fall)
      deallocate(fclr)
      deallocate(fsdir)
      deallocate(fsdif)
      deallocate(asyclt)
      deallocate(cc)
      deallocate(rrt)
      deallocate(ttt)
      deallocate(tdt)
      deallocate(rst)
      deallocate(tst)
      deallocate(dum1)
      deallocate(dum2)
      deallocate(dum3)
      deallocate(dum)
      return
      end
