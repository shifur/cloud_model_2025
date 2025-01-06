c*******************************************************************
      subroutine cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *           cc,rr,tt,td,rs,ts,fclr,fall,fsdir,fsdif,
     *           dwflux,upflux)   ! cccshie 8/19/04
c    *           cc,rr,tt,td,rs,ts,fclr,fall,fsdir,fsdif)
      include 'radiation.h'
      include 'dimensions.h'
c*******************************************************************
c  compute upward and downward fluxes using a two-stream adding method
c  following equations (6.9)-(6.16).
c
c  clouds are grouped into high, middle, and low clouds which are assumed
c  randomly overlapped. it involves a maximum of 8 sets of calculations.
c  in each set of calculations, each atmospheric layer is homogeneous,
c  either totally filled with clouds or without clouds.
c  input parameters:
c
c   m:   number of soundings
c   np:  number of atmospheric layers
c   ict: the level separating high and middle clouds
c   icb: the level separating middle and low clouds
c   ih1,ih2,im1,im2,is1,is2: indices for three group of clouds
c   cc:  effective cloud covers for high, middle and low clouds
c   rr:  reflection of a layer illuminated by beam radiation
c   tt:  total (direct+diffuse) transmission of a layer illuminated 
c        by beam radiation
c   td:  direct beam transmission
c   rs:  reflection of a layer illuminated by diffuse radiation
c   ts:  transmission of a layer illuminated by diffuse radiation
c
c  output parameters:
c
c     fclr:  clear-sky flux (downward minus upward)
c     fall:  all-sky flux (downward minus upward)
c     fsdir: surface direct downward flux
c     fsdif: surface diffuse downward flux
c
c*********************************************************************c
cccshie 8/19/04
c     implicit none
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
      integer nx1, ny1  ! 9/15/04
       parameter (nx1=nx,ny1=1)    ! cccshie 9/15/04
c-----input parameters
      integer m,np,ict,icb,ih1,ih2,im1,im2,is1,is2
      real rr(m,np+1,2),tt(m,np+1,2),td(m,np+1,2)
      real rs(m,np+1,2),ts(m,np+1,2)
      real cc(m,3)
c-----temporary array
      integer i,k,ih,im,is
!     real rra(m,np+1,2,2),tta(m,np+1,2,2),tda(m,np+1,2,2)
!     real rsa(m,np+1,2,2),rxa(m,np+1,2,2)
!     real ch(m),cm(m),ct(m),flxdn(m,np+1)
!     real fdndir(m),fdndif(m),fupdif
      real denm,xx,yy
      real fupdif
      real   ,  allocatable :: rra(:,:,:,:)
      real   ,  allocatable :: tta(:,:,:,:)
      real   ,  allocatable :: tda(:,:,:,:)
      real   ,  allocatable :: rsa(:,:,:,:)
      real   ,  allocatable :: rxa(:,:,:,:)
      real   ,  allocatable :: ch(:)
      real   ,  allocatable :: cm(:)
      real   ,  allocatable :: ct(:)
      real   ,  allocatable :: flxdn(:,:)
      real   ,  allocatable :: fdndir(:)
      real   ,  allocatable :: fdndif(:)
c-----output parameters
      real fclr(m,np+1),fall(m,np+1)
      real fsdir(m),fsdif(m)
cccshie 8/19/04
      real dwflux(m,np+1)
      real upflux(m,np+1)
cccshie 8/19/04
      real rflux(nx1,ny1,8)
      common/radflux/rflux
      allocate(rra(m,np+1,2,2))
      allocate(tta(m,np+1,2,2))
      allocate(tda(m,np+1,2,2))
      allocate(rsa(m,np+1,2,2))
      allocate(rxa(m,np+1,2,2))
      allocate(ch(m))
      allocate(cm(m))
      allocate(ct(m))
      allocate(flxdn(m,np+1))
      allocate(fdndir(m))
      allocate(fdndif(m))
c-----compute transmittances and reflectances for a composite of
c     layers. layers are added one at a time, going down from the top.
c     tda is the composite direct transmittance illuminated by beam radiation
c     tta is the composite total transmittance illuminated by
c         beam radiation
c     rsa is the composite reflectance illuminated from below
c         by diffuse radiation
c     tta and rsa are computed from eqs. (6.10) and (6.12)
c-----for high clouds
c     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition
      do ih=ih1,ih2
       do i=1,m
          tda(i,1,ih,1)=td(i,1,ih)
          tta(i,1,ih,1)=tt(i,1,ih)
          rsa(i,1,ih,1)=rs(i,1,ih)
          tda(i,1,ih,2)=td(i,1,ih)
          tta(i,1,ih,2)=tt(i,1,ih)
          rsa(i,1,ih,2)=rs(i,1,ih)
       enddo
       do k=2,ict-1
        do i=1,m
          denm = ts(i,k,ih)/( 1.-rsa(i,k-1,ih,1)*rs(i,k,ih))
          tda(i,k,ih,1)= tda(i,k-1,ih,1)*td(i,k,ih)
          tta(i,k,ih,1)= tda(i,k-1,ih,1)*tt(i,k,ih)
     *          +(tda(i,k-1,ih,1)*rsa(i,k-1,ih,1)*rr(i,k,ih)
     *          +tta(i,k-1,ih,1)-tda(i,k-1,ih,1))*denm
          rsa(i,k,ih,1)= rs(i,k,ih)+ts(i,k,ih)
     *                  *rsa(i,k-1,ih,1)*denm
          if(tda(i,k,ih,1).lt.1.e-10) tda(i,k,ih,1)=0.
          if(tta(i,k,ih,1).lt.1.e-10) tta(i,k,ih,1)=0.
          tda(i,k,ih,2)= tda(i,k,ih,1)
          tta(i,k,ih,2)= tta(i,k,ih,1)
          rsa(i,k,ih,2)= rsa(i,k,ih,1)
        enddo
       enddo
c-----for middle clouds
c     im=1 for clear-sky condition, im=2 for cloudy-sky condition
      do im=im1,im2
       do k=ict,icb-1
        do i=1,m
          denm = ts(i,k,im)/( 1.-rsa(i,k-1,ih,im)*rs(i,k,im))
          tda(i,k,ih,im)= tda(i,k-1,ih,im)*td(i,k,im)
          tta(i,k,ih,im)= tda(i,k-1,ih,im)*tt(i,k,im)
     *         +(tda(i,k-1,ih,im)*rsa(i,k-1,ih,im)*rr(i,k,im)
     *         +tta(i,k-1,ih,im)-tda(i,k-1,ih,im))*denm
          rsa(i,k,ih,im)= rs(i,k,im)+ts(i,k,im)
     *                  *rsa(i,k-1,ih,im)*denm
          if(tda(i,k,ih,im).lt.1.e-10) tda(i,k,ih,im)=0.
          if(tta(i,k,ih,im).lt.1.e-10) tta(i,k,ih,im)=0.
        enddo
       enddo
      enddo                 ! end im loop
      enddo                 ! end ih loop
c-----layers are added one at a time, going up from the surface.
c     rra is the composite reflectance illuminated by beam radiation
c     rxa is the composite reflectance illuminated from above
c         by diffuse radiation
c     rra and rxa are computed from eqs. (6.9) and (6.11)
c-----for the low clouds
c     is=1 for clear-sky condition, is=2 for cloudy-sky condition
      do is=is1,is2
       do i=1,m
         rra(i,np+1,1,is)=rr(i,np+1,is)
         rxa(i,np+1,1,is)=rs(i,np+1,is)
         rra(i,np+1,2,is)=rr(i,np+1,is)
         rxa(i,np+1,2,is)=rs(i,np+1,is)
       enddo
       do k=np,icb,-1
        do i=1,m
          denm=ts(i,k,is)/( 1.-rs(i,k,is)*rxa(i,k+1,1,is) )
          rra(i,k,1,is)=rr(i,k,is)+(td(i,k,is)*rra(i,k+1,1,is)
     *        +(tt(i,k,is)-td(i,k,is))*rxa(i,k+1,1,is))*denm
          rxa(i,k,1,is)= rs(i,k,is)+ts(i,k,is)
     *        *rxa(i,k+1,1,is)*denm
          rra(i,k,2,is)=rra(i,k,1,is)
          rxa(i,k,2,is)=rxa(i,k,1,is)
        enddo
       enddo
c-----for middle clouds
      do im=im1,im2
       do k=icb-1,ict,-1
        do i=1,m
          denm=ts(i,k,im)/( 1.-rs(i,k,im)*rxa(i,k+1,im,is) )
          rra(i,k,im,is)= rr(i,k,im)+(td(i,k,im)*rra(i,k+1,im,is)
     *        +(tt(i,k,im)-td(i,k,im))*rxa(i,k+1,im,is))*denm
          rxa(i,k,im,is)= rs(i,k,im)+ts(i,k,im)
     *        *rxa(i,k+1,im,is)*denm
        enddo
       enddo
      enddo                 ! end im loop
      enddo                 ! end is loop
c-----integration over eight sky situations.
c     ih, im, is denotes high, middle and low cloud groups.
      do ih=ih1,ih2
c-----clear portion 
         if(ih.eq.1) then
           do i=1,m
             ch(i)=1.0-cc(i,1)
           enddo
          else
c-----cloudy portion
           do i=1,m
             ch(i)=cc(i,1)
           enddo
          endif
      do im=im1,im2
c-----clear portion
         if(im.eq.1) then
           do i=1,m
              cm(i)=ch(i)*(1.0-cc(i,2))
           enddo
         else
c-----cloudy portion
           do i=1,m
              cm(i)=ch(i)*cc(i,2) 
           enddo
         endif
      do is=is1,is2
c-----clear portion
         if(is.eq.1) then
           do i=1,m
             ct(i)=cm(i)*(1.0-cc(i,3)) 
           enddo
         else
c-----cloudy portion
           do i=1,m
             ct(i)=cm(i)*cc(i,3)
           enddo
         endif
c-----add one layer at a time, going down.
       do k=icb,np
        do i=1,m
          denm = ts(i,k,is)/( 1.-rsa(i,k-1,ih,im)*rs(i,k,is) )
          tda(i,k,ih,im)= tda(i,k-1,ih,im)*td(i,k,is)
          tta(i,k,ih,im)=  tda(i,k-1,ih,im)*tt(i,k,is)
     *         +(tda(i,k-1,ih,im)*rr(i,k,is)
     *         *rsa(i,k-1,ih,im)+tta(i,k-1,ih,im)-tda(i,k-1,ih,im))*denm
          rsa(i,k,ih,im)= rs(i,k,is)+ts(i,k,is)
     *         *rsa(i,k-1,ih,im)*denm
          if(tda(i,k,ih,im).lt.1.e-10) tda(i,k,ih,im)=0.
          if(tta(i,k,ih,im).lt.1.e-10) tta(i,k,ih,im)=0.
        enddo
       enddo
c-----add one layer at a time, going up.
       do k=ict-1,1,-1
        do i=1,m
          denm =ts(i,k,ih)/(1.-rs(i,k,ih)*rxa(i,k+1,im,is))
          rra(i,k,im,is)= rr(i,k,ih)+(td(i,k,ih)*rra(i,k+1,im,is)
     *        +(tt(i,k,ih)-td(i,k,ih))*rxa(i,k+1,im,is))*denm
          rxa(i,k,im,is)= rs(i,k,ih)+ts(i,k,ih)
     *        *rxa(i,k+1,im,is)*denm
        enddo
       enddo
c-----compute fluxes following eq. (6.15) for fupdif and
c     eq. (6.16) for (fdndir+fdndif)
c     fdndir is the direct  downward flux
c     fdndif is the diffuse downward flux
c     fupdif is the diffuse upward flux
      do k=2,np+1
       do i=1,m
         denm= 1./(1.-rsa(i,k-1,ih,im)*rxa(i,k,im,is))
         fdndir(i)= tda(i,k-1,ih,im)
         xx= tda(i,k-1,ih,im)*rra(i,k,im,is)
         yy= tta(i,k-1,ih,im)-tda(i,k-1,ih,im)
         fdndif(i)= (xx*rsa(i,k-1,ih,im)+yy)*denm
         fupdif= (xx+yy*rxa(i,k,im,is))*denm
         flxdn(i,k)= fdndir(i)+fdndif(i)-fupdif
cccshie 8/19/04
         dwflux(i,k) = fdndir(i)+fdndif(i)
         upflux(i,k) = fupdif
       enddo
      enddo
       do i=1,m
         flxdn(i,1)=1.0-rra(i,1,im,is)
c        dwflux(i,1)=1.0-rra(i,1,im,is)  ! cccshie 8/19/04
         dwflux(i,1)=1.0                 ! cccshie 10/19/04 after Ming-Dah's correction
         upflux(i,1)=-rra(i,1,im,is)     ! cccshie 10/19/04 after Ming-Dah's correction
       enddo
c-----summation of fluxes over all sky situations;
c     the term in the brackets of eq. (7.11)
       do k=1,np+1
        do i=1,m
           if(ih.eq.1 .and. im.eq.1 .and. is.eq.1) then
             fclr(i,k)=flxdn(i,k)
           endif
             fall(i,k)=fall(i,k)+flxdn(i,k)*ct(i)
cccshie 8/19/04
         dwflux(i,k) = dwflux(i,k)*ct(i)
         upflux(i,k) = upflux(i,k)*ct(i)
        enddo
       enddo
        do i=1,m
            fsdir(i)=fsdir(i)+fdndir(i)*ct(i)
            fsdif(i)=fsdif(i)+fdndif(i)*ct(i)
        enddo
       enddo                 ! end is loop
       enddo                 ! end im loop
       enddo                 ! end ih loop
      deallocate(rra)
      deallocate(tta)
      deallocate(tda)
      deallocate(rsa)
      deallocate(rxa)
      deallocate(ch)
      deallocate(cm)
      deallocate(ct)
      deallocate(flxdn)
      deallocate(fdndir)
      deallocate(fdndif)
      return
      end
