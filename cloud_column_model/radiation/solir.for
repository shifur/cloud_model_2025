c***********************************************************************
cccshie 8/19/04
c     subroutine solir (m,np,wh,dp,overcast,cldwater,
      subroutine solir (m,np,jj2,wh,dp,overcast,cldwater,
     *                  cwp,taucld,reff,ict,icb,fcld,cosz,
     *                  taual,ssaal,asyal,
     *                  rsirbm,rsirdf,flx,flc,fdirir,fdifir)
      implicit none
      include 'dimensions.h'
      include 'radiation.h'
c************************************************************************
c  compute solar flux in the infrared region. the spectrum is divided
c   into three bands:
c
c          band   wavenumber(/cm)  wavelength (micron)
c          1( 9)    14280-8200         0.70-1.22
c          2(10)     8200-4400         1.22-2.27
c          3(11)     4400-1000         2.27-10.0
c
c----- input parameters:                            units      size
c
c  number of soundings (m)                          n/d         1
c  number of atmospheric layers (np)                n/d         1
c  layer scaled-water vapor content (wh)          gm/cm^2      m*np
c  option for scaling cloud optical thickness       n/d         1
c        overcast="true" if scaling is not required
c        overcast="fasle" if scaling is required
c  input option for cloud optical thickness         n/d         1
c        cldwater="true" if taucld is provided
c        cldwater="false" if cwp is provided
c  cloud water concentration (cwp)                gm/m**2      m*np*3
c        index 1 for ice particles
c        index 2 for liquid drops
c        index 3 for rain drops
c  cloud optical thickness (taucld)                 n/d        m*np*3
c        index 1 for ice paticles
c        index 2 for liquid drops
c        index 3 for rain drops
c  effective cloud-particle size (reff)           micrometer   m*np*3
c        index 1 for ice paticles
c        index 2 for liquid drops
c        index 3 for rain drops
c  level index separating high and                  n/d        m
c        middle clouds (ict)
c  level index separating middle and                n/d        m
c        low clouds (icb)
c  cloud amount (fcld)                            fraction     m*np
c  aerosol optical thickness (taual)                n/d        m*np*11
c  aerosol single-scattering albedo (ssaal)         n/d        m*np*11
c  aerosol asymmetry factor (asyal)                 n/d        m*np*11 
c  near ir surface albedo for beam                fraction     m
c        radiation (rsirbm)
c  near ir surface albedo for diffuse             fraction     m
c        radiation (rsirdf)
c
c---- temporary array
c
c  scaled cloud optical thickness                   n/d        m*np
c          for beam radiation (tauclb)
c  scaled cloud optical thickness                   n/d        m*np
c          for diffuse radiation  (tauclf)     
c
c----- output (updated) parameters:
c
c  all-sky flux (downward-upward) (flx)           fraction     m*(np+1)
c  clear-sky flux (downward-upward) (flc)         fraction     m*(np+1)
c  all-sky direct downward ir flux at
c          the surface (fdirir)                   fraction     m
c  all-sky diffuse downward ir flux at
c          the surface (fdifir)                   fraction     m
c
c**********************************************************************
cccshie 8/19/04
c     implicit none
!       parameter (nxf=4000,nyf=75,nzf=37,nt=38640,itt=244) 
!       parameter (nt=38640,itt=244) 
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
      integer ih1,ih2,im1,im2,is1,is2
      real cwp(m,np,3),taucld(m,np,3),reff(m,np,3)
      real fcld(m,np),cosz(m)
      real rsirbm(m),rsirdf(m)
      real taual(m,np,11),ssaal(m,np,11),asyal(m,np,11)
      real dp(m,np),wh(m,np)
      logical overcast,cldwater
c-----output (updated) parameters
      real flx(m,np+1),flc(m,np+1)
      real fdirir(m),fdifir(m)
c-----static parameters
      integer nk,nband
      parameter (nk=10,nband=3)
      real hk(nband,nk),xk(nk),ry(nband)
      real aib(nband,2),awb(nband,2),arb(nband,2)
      real aia(nband,3),awa(nband,3),ara(nband,3)
      real aig(nband,3),awg(nband,3),arg(nband,3)
c-----temporary array
      integer ib,iv,ik,i,k
!     real dsm(m)
!     real tauclb(m,np),tauclf(m,np),cc(m,3)
!     real ssacl(m,np),asycl(m,np)
!     real rr(m,np+1,2),tt(m,np+1,2),td(m,np+1,2),
!    *     rs(m,np+1,2),ts(m,np+1,2)
!     real fall(m,np+1),fclr(m,np+1),fsdir(m),fsdif(m)
      real taurs,tauwv
!     real tausto(m,np),ssatau(m,np),asysto(m,np)
!     real tautob(m,np),ssatob(m,np),asytob(m,np)
!     real tautof(m,np),ssatof(m,np),asytof(m,np)
      real taux,reff1,reff2,w1,w2,w3,g1,g2,g3
!     real ssaclt(m),asyclt(m)
!     real rrt(m,np),ttt(m,np),tdt(m,np),rst(m,np),tst(m,np)
!     real dum1(m,np+1),dum2(m),dum3(m),dum(m,np)
ccshie 8/19/04
      real rflux(nx1,ny1,8)
      common/radflux/rflux
!     real    dwflux(nx1,np+1)
!     real    upflux(nx1,np+1)
      real   ,  allocatable :: dsm(:)
      real   ,  allocatable :: tauclb(:,:)
      real   ,  allocatable :: tauclf(:,:)
      real   ,  allocatable :: cc(:,:)
      real   ,  allocatable :: ssacl(:,:)
      real   ,  allocatable :: asycl(:,:)
      real   ,  allocatable :: rr(:,:,:)
      real   ,  allocatable :: tt(:,:,:)
      real   ,  allocatable :: td(:,:,:)
      real   ,  allocatable :: rs(:,:,:)
      real   ,  allocatable :: ts(:,:,:)
      real   ,  allocatable :: fall(:,:)
      real   ,  allocatable :: fclr(:,:)
      real   ,  allocatable :: fsdir(:)
      real   ,  allocatable :: fsdif(:)
      real   ,  allocatable :: tausto(:,:)
      real   ,  allocatable :: ssatau(:,:)
      real   ,  allocatable :: asysto(:,:)
      real   ,  allocatable :: tautob(:,:)
      real   ,  allocatable :: ssatob(:,:)
      real   ,  allocatable :: asytob(:,:)
      real   ,  allocatable :: tautof(:,:)
      real   ,  allocatable :: ssatof(:,:)
      real   ,  allocatable :: asytof(:,:)
      real   ,  allocatable :: ssaclt(:)
      real   ,  allocatable :: asyclt(:)
      real   ,  allocatable :: rrt(:,:)
      real   ,  allocatable :: ttt(:,:)
      real   ,  allocatable :: tdt(:,:)
      real   ,  allocatable :: rst(:,:)
      real   ,  allocatable :: tst(:,:)
      real   ,  allocatable :: dum1(:,:)
      real   ,  allocatable :: dum2(:)
      real   ,  allocatable :: dum3(:)
      real   ,  allocatable :: dum(:,:)
      real   ,  allocatable :: dwflux(:,:)
      real   ,  allocatable :: upflux(:,:)
c-----water vapor absorption coefficient for 10 k-intervals.
c     unit: cm^2/gm (table 2)
      data xk/            
     1  0.0010, 0.0133, 0.0422, 0.1334, 0.4217,            
     2  1.334,  5.623,  31.62,  177.8,  1000.0/  
c-----water vapor k-distribution function,
c     the sum of hk is 0.52926. unit: fraction (table 2)
      data hk/
     1 .20673,.08236,.01074,  .03497,.01157,.00360,
     2 .03011,.01133,.00411,  .02260,.01143,.00421,
     3 .01336,.01240,.00389,  .00696,.01258,.00326,
     4 .00441,.01381,.00499,  .00115,.00650,.00465,
     5 .00026,.00244,.00245,  .00000,.00094,.00145/
c-----ry is the extinction coefficient for rayleigh scattering.
c     unit: /mb (table 3)
      data ry /.0000156, .0000018, .000000/
c-----coefficients for computing the extinction coefficients of
c     ice, water, and rain particles (table 4)
      data aib/
     1  .000333, .000333, .000333,
     2     2.52,    2.52,    2.52/
      data awb/
     1  -0.0101, -0.0166, -0.0339,
     2     1.72,    1.85,    2.16/
      data arb/
     1   0.00307, 0.00307, 0.00307,
     2   0.0    , 0.0    , 0.0    /
c-----coefficients for computing the single-scattering co-albedo of
c     ice, water, and rain particles (table 5)
      data aia/
     1 -.00000260, .00215346, .08938331,
     2  .00000746, .00073709, .00299387,
     3  .00000000,-.00000134,-.00001038/
      data awa/
     1  .00000007,-.00019934, .01209318,
     2  .00000845, .00088757, .01784739,
     3 -.00000004,-.00000650,-.00036910/
      data ara/
     1  .029,      .342,      .466,
     2  .0000,     .000,      .000,
     3  .0000,     .000,      .000/
c-----coefficients for computing the asymmetry factor of 
c     ice, water, and rain particles (table 6)
      data aig/
     1  .74935228, .76098937, .84090400,
     2  .00119715, .00141864, .00126222,
     3 -.00000367,-.00000396,-.00000385/
      data awg/
     1  .79375035, .74513197, .83530748,
     2  .00832441, .01370071, .00257181,
     3 -.00023263,-.00038203, .00005519/
      data arg/
     1  .891,      .948,      .971,
     2  .0000,     .000,      .000,
     3  .0000,     .000,      .000/
      allocate(dsm(m))
      allocate(tauclb(m,np))
      allocate(tauclf(m,np))
      allocate(cc(m,3))
      allocate(ssacl(m,np))
      allocate(asycl(m,np))
      allocate(rr(m,np+1,2))
      allocate(tt(m,np+1,2))
      allocate(td(m,np+1,2))
      allocate(rs(m,np+1,2))
      allocate(ts(m,np+1,2))
      allocate(fall(m,np+1))
      allocate(fclr(m,np+1))
      allocate(fsdir(m))
      allocate(fsdif(m))
      allocate(tausto(m,np))
      allocate(ssatau(m,np))
      allocate(asysto(m,np))
      allocate(tautob(m,np))
      allocate(ssatob(m,np))
      allocate(asytob(m,np))
      allocate(tautof(m,np))
      allocate(ssatof(m,np))
      allocate(asytof(m,np))
      allocate(ssaclt(m))
      allocate(asyclt(m))
      allocate(rrt(m,np))
      allocate(ttt(m,np))
      allocate(tdt(m,np))
      allocate(rst(m,np))
      allocate(tst(m,np))
      allocate(dum1(m,np+1))
      allocate(dum2(m))
      allocate(dum3(m))
      allocate(dum(m,np))
      allocate(dwflux(nx1,np+1))
      allocate(upflux(nx1,np+1))
c-----initialize surface fluxes, reflectances, and transmittances.
c     the reflectance and transmittance of the clear and cloudy portions
c     of a layer are denoted by 1 and 2, respectively.
c     cc is the maximum cloud cover in each of the high, middle, and low
c     cloud groups.
c     1/dsm=1/cos(53)=1.66
      do i=1,m
         dsm(i)=0.602
         fdirir(i)=0.0
         fdifir(i)=0.0
         rr(i,np+1,1)=rsirbm(i)
         rr(i,np+1,2)=rsirbm(i)
         rs(i,np+1,1)=rsirdf(i)
         rs(i,np+1,2)=rsirdf(i)
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
c-----integration over spectral bands
      do 100 ib=1,nband
       iv=ib+8
c-----compute cloud optical thickness. eqs. (4.6) and (4.11)
       do k=1,np
        do i=1,m
c          reff(i,k,1)=max(reff(i,k,1),1.0)
          reff(i,k,1)=min(reff(i,k,1),150.0)
c          reff(i,k,2)=max(reff(i,k,2),4.0)
          reff(i,k,2)=min(reff(i,k,2),20.0)
        enddo
       enddo
      if (cldwater) then
       do k=1,np
        do i=1,m
          taucld(i,k,1)=cwp(i,k,1)*(aib(ib,1)
     *                    +aib(ib,2)/reff(i,k,1))
          taucld(i,k,2)=cwp(i,k,2)*(awb(ib,1)
     *                    +awb(ib,2)/reff(i,k,2))
          taucld(i,k,3)=cwp(i,k,3)*arb(ib,1)
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
c     tauclf is for diffuse radiation.
       call cldscale(m,np,cosz,fcld,taucld,ict,icb,
     *              cc,tauclb,tauclf)
      endif
c-----compute cloud single scattering albedo and asymmetry factor
c     for a mixture of ice and liquid particles.
c     eqs.(4.6)-(4.8), (6.2)-(6.4)
       do k=1,np
        do i=1,m
           ssaclt(i)=0.99999
           asyclt(i)=1.0
           taux=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)
          if (taux.gt.0.02 .and. fcld(i,k).gt.0.01) then
c           reff1=reff(i,k,1)
c           reff2=reff(i,k,2)
            reff1=min(reff(i,k,1),150.)
            reff2=min(reff(i,k,2),20.0)
           w1=(1.-(aia(ib,1)+(aia(ib,2)+
     *         aia(ib,3)*reff1)*reff1))*taucld(i,k,1)
           w2=(1.-(awa(ib,1)+(awa(ib,2)+
     *         awa(ib,3)*reff2)*reff2))*taucld(i,k,2)
           w3=(1.- ara(ib,1))*taucld(i,k,3)
           ssaclt(i)=(w1+w2+w3)/taux
           g1=(aig(ib,1)+(aig(ib,2)+aig(ib,3)*reff1)*reff1)*w1
           g2=(awg(ib,1)+(awg(ib,2)+awg(ib,3)*reff2)*reff2)*w2
           g3= arg(ib,1)*w3
           asyclt(i)=(g1+g2+g3)/(w1+w2+w3)
          endif
        enddo
         do i=1,m
           ssacl(i,k)=ssaclt(i)
         enddo
         do i=1,m
           asycl(i,k)=asyclt(i)
         enddo
       enddo
c-----integration over the k-distribution function
       do 200 ik=1,nk
        do k=1,np
         do i=1,m
           taurs=ry(ib)*dp(i,k)
           tauwv=xk(ik)*wh(i,k)
c-----compute clear-sky optical thickness, single scattering albedo,
c     and asymmetry factor. eqs.(6.2)-(6.4)
           tausto(i,k)=taurs+tauwv+taual(i,k,iv)+1.0e-8
           ssatau(i,k)=ssaal(i,k,iv)*taual(i,k,iv)+taurs
           asysto(i,k)=asyal(i,k,iv)*ssaal(i,k,iv)*taual(i,k,iv)
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
c-----for direct incident radiation. eqs.(6.2)-(6.4)
           tautob(i,k)=tausto(i,k)+tauclb(i,k)
           ssatob(i,k)=(ssatau(i,k)+ssacl(i,k)*tauclb(i,k))
     *                /tautob(i,k)+1.0e-8
           ssatob(i,k)=min(ssatob(i,k),0.999999)
           asytob(i,k)=(asysto(i,k)+asycl(i,k)*ssacl(i,k)*tauclb(i,k))
     *                /(ssatob(i,k)*tautob(i,k))
c-----for diffuse incident radiation
           tautof(i,k)=tausto(i,k)+tauclf(i,k)
           ssatof(i,k)=(ssatau(i,k)+ssacl(i,k)*tauclf(i,k))
     *                /tautof(i,k)+1.0e-8
           ssatof(i,k)=min(ssatof(i,k),0.999999)
           asytof(i,k)=(asysto(i,k)+asycl(i,k)*ssacl(i,k)*tauclf(i,k))
     *                /(ssatof(i,k)*tautof(i,k))
         enddo
        enddo
c-----for direct incident radiation
          call deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)
c-----diffuse incident radiation is approximated by beam radiation with
c     an incident angle of 53 degrees, eqs.(6.5) and (6.6)
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
c-----initialize clear-sky flux (fclr), all-sky flux (fall), 
c     and surface downward fluxes (fsdir and fsdif)
        do k=1,np+1
         do i=1,m
           fclr(i,k)=0.0
           fall(i,k)=0.0
           dum1(i,k)=0.0
         enddo
        enddo
        do i=1,m
           fsdir(i)=0.0
           fsdif(i)=0.0
           dum2(i)=0.
           dum3(i)=0.
        enddo
        if (overcast) then
              ih1=1
              ih2=1
              im1=1
              im2=1
              is1=1
              is2=1
c-----for clear-sky fluxes only
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
c-----flux integration following eq. (6.1)
       do k=1,np+1
        do i=1,m
          flx(i,k) = flx(i,k)+fall(i,k)*hk(ib,ik)
        enddo
        do i=1,m
          flc(i,k) = flc(i,k)+fclr(i,k)*hk(ib,ik)
        enddo
cccshie 8/19/04
        do i=1,m
         if(k.eq.1)then
c       rflux(i,jj2,3)=rflux(i,jj2,3)+hk(ib,ik) ! downward SW TOA, djohnson's old, may do it later
        rflux(i,jj2,3)=rflux(i,jj2,3)+dwflux(i,k)*hk(ib,ik) ! downward SW TOA
         endif
         if(k.eq.(np+1))then
        rflux(i,jj2,1)=rflux(i,jj2,1)+dwflux(i,k)*hk(ib,ik) ! downward SW surface
         endif
        enddo
       enddo
c-----compute downward surface fluxes in the ir region
       do i=1,m
          fdirir(i) = fdirir(i)+fsdir(i)*hk(ib,ik)
          fdifir(i) = fdifir(i)+fsdif(i)*hk(ib,ik)
       enddo
  200 continue
  100 continue
      deallocate(dsm)
      deallocate(tauclb)
      deallocate(tauclf)
      deallocate(cc)
      deallocate(ssacl)
      deallocate(asycl)
      deallocate(rr)
      deallocate(tt)
      deallocate(td)
      deallocate(rs)
      deallocate(ts)
      deallocate(fall)
      deallocate(fclr)
      deallocate(fsdir)
      deallocate(fsdif)
      deallocate(tausto)
      deallocate(ssatau)
      deallocate(asysto)
      deallocate(tautob)
      deallocate(ssatob)
      deallocate(asytob)
      deallocate(tautof)
      deallocate(ssatof)
      deallocate(asytof)
      deallocate(ssaclt)
      deallocate(asyclt)
      deallocate(rrt)
      deallocate(ttt)
      deallocate(tdt)
      deallocate(rst)
      deallocate(tst)
      deallocate(dum1)
      deallocate(dum2)
      deallocate(dum3)
      deallocate(dum)
      deallocate(dwflux)
      deallocate(upflux)
      return
      end
