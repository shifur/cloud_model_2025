c     4-6 2005   TAO
c
      subroutine akcoef (ak,ak1,u1,dpt1,ta1,dqv1,qa1,pi,qcl1,qci1,
     1                   rho1,rrho,am,am1,ww1,dz,dt)
      implicit none
      include 'dimensions.h'

c Input variables
      real    ak(nx,ny,nz), ak1(nx,ny,nz), u1(nx,ny,nz)
      real    dpt1(nx,ny,nz), ta1(nz), dqv1(nx,ny,nz), qa1(nz)
      real    pi(nz), qcl1(nx,ny,nz), qci1(nx,ny,nz)
      real    rho1(nz), rrho(nz), am(nz), am1(nz), ww1(nx,ny,nz)
      real    dz, dt


c Parameters
      integer iadvh
      parameter (iadvh = 4)
      real    al, cp, ra, ck, ce, a2k
      parameter (al=2.5e10)
      parameter (cp=1.004e7)
      parameter (ra=2.87e6)
      parameter (ck=0.2)
      parameter (ce=0.7)
c      parameter (a2k=0.25e6)
      parameter (a2k=0.15e6)
c
      real  xy1(nx,ny),xy2(nx,ny),xy3(nx,ny),xy4(nx,ny),xy5(nx,ny),
     1   xy6(nx,ny),xy7(nx,ny)
      real  tp(nx,ny),tm(nx,ny),qp(nx,ny),qm(nx,ny),y1(nz),y2(nz),y3(nz)

c Local variables
      integer i,j,k,km,kp
      real y1k,pik,ta1kp,qa1kp,y2kp,y2km,y3k,ta1k
      real amk,am1k,am1kp,a22z,ammkp,ammk,a33k,a44k,a55k,coefk,z2k
      real alcp, cpra, alsq, cpal, rdz, rdz2, rd4z, rd2z, r2dz2
      real coef(nz), z1(nz), z2(nz), bskm(nz)

      save
c
      alcp=al/cp
      cpra=cp*ra
      alsq=.622*al*al
      cpal=.622*1.61*cp*al
      rdz=1./dz
      rd2z=.5*rdz
      rdz2=rdz*rdz
      rd4z=.25*rdz
      r2dz2=.5*rdz2

      do k=1,nz
        coef(k)=ck*ck*.5*dz*dz/(am(k)*am(k))
        z1(k)=3.*980.*coef(k)
        z2(k)=.5*ce/(ck*dz*dz/(am(k)*am(k)))
        bskm(k)=a2k/2.0
!         ckh=2.
!         bskm(k)=a2k/ckh
      enddo

      do k=1,nz
       y1(k)=1./ta1(k)
       y2(k)=1./pi(k)
       y3(k)=-am(k)*z1(k)*rd2z
      enddo
c
c$doacross local(k,j,i)
	 do k=1,nz
	 do j=1,ny
	 do i=1,nx
	    u1(i,j,k)=0.
	 enddo
	 enddo
	 enddo
c    ****************************************
c$doacross local(j,i)
      do j=1,ny
      do i=1,nx
        tm(i,j)=ta1(1)+dpt1(i,j,1)
        qm(i,j)=qa1(1)+dqv1(i,j,1)
        xy1(i,j)=ta1(2)+dpt1(i,j,2)
	  xy2(i,j)=qa1(2)+dqv1(i,j,2)
      enddo
      enddo

      do 2001 k=2,nz-1
	kp=k+1
	km=k-1
	ta1kp=ta1(kp)
	qa1kp=qa1(kp)
	y1k=y1(k)
	pik=pi(k)
	y2kp=y2(kp)
	y2km=y2(km)
	y3k=y3(k)
	ta1k=ta1(k)

c$doacross local(j,i)
      do 2000 j=1,ny
       do i=1,nx
	 xy4(i,j)=qcl1(i,j,k)+qci1(i,j,k)
	 tp(i,j)=ta1kp+dpt1(i,j,kp)
	 qp(i,j)=qa1kp+dqv1(i,j,kp)
	if (xy4(i,j) .lt. 1.e-5) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 u1(i,j,k)=y1k*( tp(i,j)*(1.+.61*qp(i,j))
     1                  -tm(i,j)*(1.+.61*qm(i,j)) )
	else
	  xy3(i,j)=xy1(i,j)*pik
	  xy4(i,j)=cpra*xy3(i,j)*xy3(i,j)
	  xy5(i,j)=(xy4(i,j)+cpal*xy3(i,j)*xy2(i,j))
     1                /((xy4(i,j)+alsq*xy2(i,j))*ta1k)
	 xy6(i,j)=tp(i,j)-tm(i,j)+alcp*(qp(i,j)*y2kp-qm(i,j)*y2km)
	  xy7(i,j)=qp(i,j)-qm(i,j)+qcl1(i,j,kp)
     1                  -qcl1(i,j,km)+qci1(i,j,kp)-qci1(i,j,km)
	 u1(i,j,k)=xy5(i,j)*xy6(i,j)-xy7(i,j)
	endif
        u1(i,j,k)=y3k*u1(i,j,k)
	 tm(i,j)=xy1(i,j)
	 qm(i,j)=xy2(i,j)
	 xy1(i,j)=tp(i,j)
	 xy2(i,j)=qp(i,j)
       enddo
 2000  continue

 2001  continue
c
c$doacross local(j,i)
      do j=1,ny
      do i=1,nx
         qm(i,j)=ak1(i,j,1)*ak1(i,j,1)
         tm(i,j)=ak1(i,j,2)*ak1(i,j,2)
      enddo
      enddo


      do 3000 k=2,nz-1
	kp=k+1
	km=k-1
	amk=am(k)
	am1k=am1(k)
	am1kp=am1(kp)
	a22z=amk*rrho(k)*rd2z
c        amm2=amk*amk*r4dz2
	ammkp=am1kp*rho1(kp)*bskm(kp)
	ammk=am1k*rho1(k)*bskm(k)
	a33k=rd4z*amk
	a44k=amk*amk*rdz2
	a55k=amk*r2dz2
       coefk=coef(k)
       z2k=z2(k)

c$doacross local(j,i)
       do j=1,ny
       do i=1,nx
	 tp(i,j)=ak1(i,j,kp)*ak1(i,j,kp)
       enddo
       enddo

c$doacross local(j,jp,jm,i,ip,im)
       do j=1,ny
!        jp=j
!        jm=j
       do i=1,nx
!        ip=i
!        im=i
       xy1(i,j)=0.0
!      1  (ww1(ip,j,kp)+ww1(ip,j,k)-ww1(im,j,kp)-ww1(im,j,k))*rd4x
!      2  +(uu1(ip,j,kp)+uu1(i,j,kp)-uu1(ip,j,km)-uu1(i,j,km))*a33k
	xy1(i,j)=xy1(i,j)*xy1(i,j)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       xy2(i,j)=0.0 ! u and v wind = 0.0
!      1  (uu1(ip,jp,k)+uu1(i,jp,k)-uu1(ip,jm,k)-uu1(i,jm,k))*rd4y
!      2  +(vv1(ip,j,k)+vv1(ip,jp,k)-vv1(im,j,k)-vv1(im,jp,k))*rd4x
	xy2(i,j)=xy2(i,j)*xy2(i,j)
       xy3(i,j)=0.0
!      1  (ww1(i,jp,kp)+ww1(i,jp,k)-ww1(i,jm,k)-ww1(i,jm,kp))*rd4y
!      2  +(vv1(i,j,kp)+vv1(i,jp,kp)-vv1(i,j,km)-vv1(i,jp,km))
!      3  *a33k
	xy3(i,j)=xy3(i,j)*xy3(i,j)
       qp(i,j)=coefk*(xy1(i,j)+xy2(i,j)+xy3(i,j)
!      1  +((uu1(ip,j,k)-uu1(i,j,k))**2*rdx2
!      2  +(vv1(i,jp,k)-vv1(i,j,k))**2*rdy2
     3  +a44k*(ww1(i,j,kp)-ww1(i,j,k))**2 )*2.!)
       u1(i,j,k)=u1(i,j,k)+qp(i,j)-z2k*tm(i,j)
!      1  +((tm(ip,j)-tm(i,j))-(tm(i,j)-tm(im,j)))*r2dx2
!      2  +((tm(i,jp)-tm(i,j))-(tm(i,j)-tm(i,jm)))*r2dy2
     3  +a55k*(am1kp*(tp(i,j)-tm(i,j))-am1k*(tm(i,j)-qm(i,j)))
       enddo
       enddo

c$doacross local(j,i)
       do j=1,ny
       do i=1,nx
	u1(i,j,k)=u1(i,j,k)-a22z*(
     1      -2.*(ammkp*(ak1(i,j,kp)-ak1(i,j,k))
     2      -ammk*(ak1(i,j,k)-ak1(i,j,km)))*rdz)
       enddo
       enddo
	
ccc   4th order numerical smoother   ccccccccccccccccccccccccccccccccccc
c$doacross local(j,i)
         do j=1,ny
           do i=1,nx
               xy2(i,j)=0.0
            enddo
         enddo
c		

c$doacross local(i,j)
         do i=1,nx
           do j=1,ny
             xy1(i,j)=0.0
           enddo
         enddo
c
c$doacross local(j,i)
	do j=1,ny
	do i=1,nx
	 qm(i,j)=tm(i,j)
	 tm(i,j)=tp(i,j)
	enddo
        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 3000  continue

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     ******************************************************************
c
       call advectak (ak,ak1,u1,ww1,am,am1,rho1,dz,dt)
cc
ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
      end
