cc D. Posselt 7/28/2008
cc Modified this to be called as a column model from an external driver
cc Stripped out unneccessary common blocks and statistics
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine pradrat (iflag,iopcloud,iradave,cosz,tb,qb,rho,
     1                    p0,p00,pi,dz0,tland,qcl,qrn,qci,qcs,qcg,
     2                    dpt,dqv,rsw,rlw,rflux1,improve)
c      implicit none
      include 'dimensions.h'

c D.Posselt Variable declarations
      integer kles,improve
      real    qcl(nx,ny,nz),qrn(nx,ny,nz)
      real    qci(nx,ny,nz),qcs(nx,ny,nz),qcg(nx,ny,nz)
      real    dpt(nx,ny,nz),dqv(nx,ny,nz)
      real    rsw(nx,ny,nz),rlw(nx,ny,nz)
      real    cosz
      real    rflux1(nx,ny,8)

c      parameter (nx=66,ny=10,nz=34)
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
       integer nx1, ny1
       parameter (nx1=1,ny1=1)    ! cccshie 9/15/04
       integer nq,nw,np
cccccc
!       common/mpi_parameter/imax,iles1,iles,il2,jmax,jles1,jles,jl2,
!      1    kmax,kles,kl2
!        common/b1cr/ qcl(nx,ny,nz),qrn(nx,ny,nz)
!        common/b1ig/ qci(nx,ny,nz),qcg(nx,ny,nz)
!        common/b1s/ qcs(nx,ny,nz)
!        common/b1tq/ dpt(nx,ny,nz),dqv(nx,ny,nz)
!        common/slwave/ rsw(nx,ny,nz),rlw(nx,ny,nz)
       real tb(nz),qb(nz),rho(nz),p0(nz),pi(nz),p00(nz),dz0(nz),
     1      tland(nx,ny)
c
       common/radtemp1/ qc(nx,nz),qr(nx,nz),qi(nx,nz),qs(nx,nz),
     1 qg(nx,nz),pt(nx,nz),qv(nx,nz),qh(nx,nz),qsw(nx,nz),qgw(nx,nz),
     2 qhw(nx,nz)
       common/radtemp2/ rsw1(nx,nz),rlw1(nx,nz)
       common/radtemp3/ tsfc_a(nx),tsst(nx)
cccshie 8/18/04
      real rflux(nx1,ny1,8)
      common/radflux/rflux
       save
       twcz=1.e-6
       tcont=190.
c       iles=nx-1
c       jles=ny-1
       kles=nz-1
       nq=nz+nadd-1
       nw=nq-1
       np=nz-2+nadd

!       print *, 'tland ',tland
!       print '(a)',' k          pi                  p0               p00'
!       do k=1,nz
!         print '(i3,3(2X,e10.5))',k,pi(k),p0(k),p00(k)
!       enddo
!       print '(a)',' k          tb                  qb               rho'
!       do k=1,nz
!         print '(i3,3(2X,e10.5))',k,tb(k),qb(k),rho(k)
!       enddo
!       print '(a)',' k          dz0'
!       do k=1,nz
!         print '(i3,1(2X,e10.5))',k,dz0(k)
!       enddo
!       print *,'iopcloud, iradave, cosz, twcz, tcont ',
!      1         iopcloud, iradave, cosz, twcz, tcont
c
cc
c
!         print*,'In pradrat cosz: ',cosz
         if (iflag .eq. 0) then
           j=1
           do 100 i=1,1
             im=i
             tsst(im)=tland(i,j)
             tsfc_a(im)=.5*(tland(i,j)+(dpt(i,j,2)+tb(2))*pi(2))
c
             do k=2,kles
               qc(im,k)=qcl(i,j,k)
               qr(im,k)=qrn(i,j,k)
               qi(im,k)=qci(i,j,k)
               qs(im,k)=qcs(i,j,k)
               qg(im,k)=qcg(i,j,k)
               pt(im,k)=(tb(k)+dpt(i,j,k))*pi(k)
                 if (pt(im,k) .le. tcont) pt(im,k)=tcont
               qv(im,k)=qb(k)+dqv(i,j,k)
                 if (qv(im,k) .le. twcz) qv(im,k)=twcz
             enddo
  100     continue
            call radrat (twcz,tcont,cosz,tb,qb,rho,p0,p00,pi,dz0,iflag,
     1                  iopcloud,iradave,j,tland)
            do k=1,nz
              do i=1,nx
                 rsw1(i,k)=0.
                 rlw1(i,k)=0.
              enddo
            enddo
         endif
c
cc
c
         do 2000 j=1,1
c
c
           do 200 i=1,1
             im=i
c
            tsst(im)=tland(i,j)
            tsfc_a(im)=.5*(tland(i,j)+(dpt(i,j,2)+tb(2))*pi(2))
c
             do k=2,kles
               qc(im,k)=qcl(i,j,k)
               qr(im,k)=qrn(i,j,k)
               qi(im,k)=qci(i,j,k)
               qs(im,k)=qcs(i,j,k)
               qg(im,k)=qcg(i,j,k)
               pt(im,k)=(tb(k)+dpt(i,j,k))*pi(k)
                 if (pt(im,k) .le. tcont) pt(im,k)=tcont
               qv(im,k)=qb(k)+dqv(i,j,k)
                 if (qv(im,k) .le. twcz) qv(im,k)=twcz
             enddo
  200     continue
           call radrat (twcz,tcont,cosz,tb,qb,rho,p0,p00,pi,dz0,iflag,
     1	iopcloud,iradave,j,tland)
         do k=2,kles
           do i=1,1
               rsw(i,j,k)=rsw1(i,k)
               rlw(i,j,k)=rlw1(i,k)
             enddo
         enddo
c
 2000   continue

       do i=1,8
         rflux1(1,1,i) = rflux(1,1,i)
       enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       return
       end
