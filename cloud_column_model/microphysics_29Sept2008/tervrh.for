      subroutine tervrh (dpt,qcl,qrn,qci,qcs,qcg,ta1,pi,
     1                   ww1,irsg,rho,fv,improve)
      include 'dimensions.h'

c D.Posselt input variable declarations
      real dpt(nx,ny,nz)
      real qcl(nx,ny,nz), qrn(nx,ny,nz)
      real qci(nx,ny,nz), qcs(nx,ny,nz), qcg(nx,ny,nz)
      real ww1(nx,ny,nz)
      real rho(nz)
      real ta1(nz)
      real pi(nz), fv(nz)
      integer irsg, improve

c D.Posselt variable declarations
      integer kles

c     parameter (nx=66,ny=10,nz=34,nm=nx)
!     parameter (nxf=260,nyf=260,nzf=43,nt=38640,itt=244) 
!       parameter (nxf=36,nyf=6,nzf=43,nt=28640,itt=182) 
!       parameter (lb=2,kb=1)
! ! define decomposition from rmp_switch.h
!       parameter (npes=1,ncol=1,nrow=1)
! ! define partial dimension for computation by decomposition
!       parameter (nx=(nxf-lb*2-1)/ncol+1+lb*2)
!       parameter (ny=(nyf-lb*2-1)/nrow+1+lb*2)
!       parameter (nz=nzf)
! ! define partial dimension for fft by transpose decomposition
!       parameter (nyc=(nyf-lb*2-1)/ncol+1+lb*2)
!       parameter (nyr= ny)
!       parameter (nxc= nx)
!       parameter (nxr=(nxf-lb*2-1)/nrow+1+lb*2)
!       parameter (nzc=(nzf-kb*2-1)/ncol+1+kb*2)
!       parameter (nzr=(nzf-kb*2-1)/nrow+1+kb*2)
c hmhj modified (nm no need here)
cx    parameter (nm=nx)
!       parameter (nz4=nz*4,nz3=nz*3,nz2=nz*2,nb3=nz+nx)
!       common/mpi_parameter/imax,iles1,iles,il2,jmax,jles1,jles,jl2,
!      1    kmax,kles,kl2
!       common/bxyz/ n,isec,nran,kt1,kt2
!       common/b5/ tb(nz),qb(nz),rho1_tmp(nz2),ta(nz),qa(nz),ta1(nz),
!      1  qa1(nz),zx(nz4),am(nz),zq(nz3),wb(nz),zw(nz2),rrho(nz),wbx(nb3)
!       common/b6/ fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),st(nz),sv(nz),
!      1   sq(nz),sc(nz),se(nz),sqa(nz)
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/rterv/ zrc,zgc,zsc,vr0,vr1,vr2,vr3,vgc,vsc
      common/rsnw/ alv,alf,als,t0,t00,avc,afc,asc,rn1,rn2,bnd2,rn3,rn4,
     1  rn5,rn50,rn51,rn52,rn53,rn6,rn60,rn61,rn62,rn63,rn7,rn8,rn9,
     2  rn10,rn101,rn102,rn10a,rn10b,rn10c,rn11,rn12,rn12a(31),
     3  rn12b(31),rn13(31),rn14,rn15,rn15a,rn16,rn171,rn172,rn17a,rn17b,
     4  rn17c,rn18,rn18a,rn19,rn191,rn192,rn19a,rn20,rn20a,rn20b,rn30,
     5  rn30a,rn21,bnd21,rn22,rn23,rn231,rn232,rn25,rn25a(31),rn31,beta,
     6  rn32,rn33,rn331,rn332,rn34,rn35,bnd1
!       common/b1tq/ dpt(nx,ny,nz),dqv(nx,ny,nz)
!       common/b1cr/ qcl(nx,ny,nz),qrn(nx,ny,nz)
!       common/b1ig/ qci(nx,ny,nz),qcg(nx,ny,nz)
!       common/b1s/ qcs(nx,ny,nz)

! Local variables

c hmhj modified (nm no need)
cx    dimension y1(nx,ny),vr(nx,ny),vs(nx,ny),vg(nx,ny),z1(nm)
      dimension y1(nx,ny),vr(nx,ny),vs(nx,ny),vg(nx,ny)
      dimension tns_funT(nx,ny),tair(nx,ny),tairc(nx,ny)
cxp ---
      real    vgcr,vscf
      real aice(7),vice(7)
      data aice/1.e-6, 1.e-5, 1.e-4, 1.e-3, 0.01, 0.1, 1./
c     data vice/10,30,60,70,80,90,100/
      data vice/5,15,30,35,40,45,50/
cxp ---
      save
c     ******************************************************************
c       ijles=nx*(ny-1)-1
c       istart=nx+2
c
c       ijles=max(nx*(ny-1), ny*(nx-1))-1
c       istart=min(nx,ny)+2
c

c D.Posselt: Set values of parameters for column model
      kles = nz-1

       cmin=1.e-15
c$doacross local(k,ij)
      do k=1,nz
      do j=1,ny
      do i=1,nx
	   ww1(i,j,k)=0.
      enddo
      enddo
      enddo
cxp ---
c     ******************************************************************
      IF (IRSG .EQ. 4) THEN           !cloud ice fall speed
       if(improve.ge.3) const_vt=1.49e4      
       if(improve.ge.3) const_d=11.9
       if(improve.ge.3) const_m=1./5.38e7
       do k=2, kles
       do j=1, ny
       do i=1, nx
        ww1(i,j,k)=0.
        y1(i,j)=.5*1.e6*rho(k)*(qci(i,j,k)+qci(i,j,k-1))        ! to g/m**3
        if (y1(i,j) .lt. 1.e-6) then
         ww1(i,j,k)=0.
        else
         if(improve.ge.3)then            !from Hong et al. (2004)
          y1(i,j)=y1(i,j)*1.e-3
          bb1=const_m*y1(i,j)**0.25
          bb2=const_d*bb1**0.5
          ww1(i,j,k)=max(const_vt*bb2**1.31, 0.0)
          ww1(i,j,k)=ww1(i,j,k)*100. !cm/s
          if (ww1(i,j,k) .gt. 20.) ww1(i,j,k)=20.    !lang et al
         else                          !from Starr & Cox
          do ic=1,6
           if (y1(i,j).gt.aice(ic) .and. y1(i,j).le.aice(ic+1)) then
            ww1(i,j,k)=vice(ic)+(vice(ic+1)-vice(ic))*
     &               (y1(i,j)-aice(ic))/(aice(ic+1)-aice(ic))
            if (ww1(i,j,k) .le. 0.0) ww1(i,j,k)=0.
           endif
          enddo
         endif
        endif
       enddo
       enddo
       enddo
      ENDIF
cxp ---
      if(irsg.ne.0) go to 1
c$doacross local(k,a1,km,j,i)
      do 100 k=2,kles
	a1=.5*rho(k)
	km=k-1
       do 10 j=1,ny
       do 10 i=1,nx
	    y1(i,j)=a1*(qrn(i,j,k)+qrn(i,j,km))
	    if (y1(i,j) .gt. cmin) then
	       vs(i,j)=sqrt( y1(i,j) )
	       vg(i,j)=sqrt( vs(i,j) )
	       vr(i,j)=vr0+vr1*vg(i,j)+vr2*vs(i,j)+vr3*vg(i,j)*vs(i,j)
	      ww1(i,j,k)=max(fv(k)*vr(i,j), 0.e0)
	    endif
   10  continue
  100 continue
      return
    1 if(irsg.ne.1) go to 2
      do 200 k=2,kles
	 km=k-1
	 a1=.5*rho(k)
	vscf=vsc*fv(k)
c$doacross local(j,i)
       do 20 j=1,ny
       do 20 i=1,nx
	    y1(i,j)=a1*(qcs(i,j,k)+qcs(i,j,km))
	    if (y1(i,j) .gt. cmin) then
               tns_funT(i,j)=1.
               if(improve.ge.3)then                            !Lang et al. 2007b
                tair(I,J)=(dpt(i,j,k)+ta1(K))*pi(k)
                tairc(i,j)=min(0.,max(-50.,tair(i,j)-t0))
                tns_funT(i,j)=exp(0.060000*tairc(i,j)*0.25*bs)               
               endif
	       ww1(i,j,k)=max(vscf*y1(i,j)**bsq *tns_funT(i,j), 0.e0)
	    endif
   20  continue
  200 continue
      return
    2 do 300 k=2,kles
	km=k-1
	 a1=.5*rho(k)
	vgcr=vgc*fv(k)
c$doacross local(j,i)
       do 30 j=1,ny
       do 30 i=1,nx
	    y1(i,j)=a1*(qcg(i,j,k)+qcg(i,j,km))
	    if (y1(i,j) .gt. cmin) then
	       ww1(i,j,k)=max(vgcr*y1(i,j)**bgq, 0.e0)
	    endif
   30  continue
  300 continue
      return
      end
