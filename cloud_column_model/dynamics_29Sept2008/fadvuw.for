c
cc    4-6 2005  TAO
c
      subroutine fadvuw (x,x1,w,am,am1,dz,dt)
      implicit none
cc    ****   compute anti-diffusive velocities   *******************
      include 'dimensions.h'

cc    input variables
      real    x(nx,ny,nz),x1(nx,ny,nz)
      real    w(nx,ny,nz)
      real    am(nz), am1(nz)
      real    dz, dt

cc    local variables
      integer nm
      parameter (nm=nx+nz)

      real    dttk
      real    dzzt(nz)
      real    dzzt1(nz)
      integer i, j, k, im, ip, jm, jp, km, kp

      real wmod(nx,ny,nz)
      real eps

      real xy1(nx,ny),xy2(nx,ny),xy3(nx,ny),xy4(nx,ny),xy5(nx,ny),
     1  y1(nm),y2(nm),y3(nm),y4(nm),y5(nm)

      save
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       eps=1.e-10

      do k=1,nz
        dzzt1(k) = am(k) * dt / dz
        dzzt(k) = am1(k) * dt / dz  ! Should be am1
      enddo

c$doacross local(k,j,i)
       do k=1,nz
	  do j=1,ny
	     do i=1,nx
		 wmod(i,j,k)=0.
	     enddo
	  enddo
       enddo

c   ***   w-component wind   *******************************************
      do 3000 k=3,nz-1
         kp=k+1
         km=k-1
      do 3000 j=1,ny
        jp=j
        jm=j
       do 300 i=1,nx
        ip=i
        im=i
        y5(i)=x(i,j,k)+x(i,j,km)
c
       if (y5(i) .ge. eps) then
        y4(i)=x(ip,j,k)+x(ip,j,km)
         wmod(i,j,k)=(abs(w(i,j,k))-w(i,j,k)*w(i,j,k)*dzzt(k))
     1                *(x(i,j,k)-x(i,j,km))/(x(i,j,k)+x(i,j,km)+eps)
       endif
  300 continue
 3000 continue
c   *****************************************************************

ccccc   non-oscillatory option (smolarkiewicz and grabowski, 1990)
      do 9000 k=1,nz
      do 9000 j=1,ny
      do 9000 i=1,nx
         w(i,j,k)=0.
 9000 continue
      do 4000 k=2,nz-1
         kp=k+1
         km=k-1
         dttk=dzzt1(k)
       do 4000 j=1,ny
         jp=j
         jm=j
       do 400 i=1,nx
         im=i
         ip=i
         y1(i)=0.
         y2(i)=0.
         y3(i)=0.
         y4(i)=0.
         y5(i)=x(i,j,k)+x(im,j,k)
        if (y5(i) .ge. eps) then
         y1(i)=max(x(im,j,k),x(i,j,k),x(ip,j,k),x1(im,j,k),
     1             x1(i,j,k),x1(ip,j,k))
         y2(i)=min(x(im,j,k),x(i,j,k),x(ip,j,k),x1(im,j,k),
     1             x1(i,j,k),x1(ip,j,k))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        y3(i)=(y1(i)-x(i,j,k))/(dttk*(max(wmod(i,j,k),0.)*x(i,j,km)
     5                       -min(wmod(i,j,kp),0.)*x(i,j,kp))+eps)
        y4(i)=(x(i,j,k)-y2(i))/(dttk*(max(wmod(i,j,kp),0.)*x(i,j,k)
     5                       -min(wmod(i,j,k),0.)*x(i,j,k))+eps)
        endif
  400  continue
 4000 continue
      do 5000 k=2,nz-1
         kp=k+1
         km=k-1
         dttk=dzzt1(k)
c
       do 500 j=1,ny
         jm=j
         jp=j
       do i=1,nx
         im=i
         ip=i

        xy1(i,j)=0.
        xy2(i,j)=0.
        xy3(i,j)=0.
        xy4(i,j)=0.

        xy5(i,j)=x(i,j,k)+x(i,jm,k)
        if (xy5(i,j) .ge. eps) then
         xy1(i,j)=max(x(i,jm,k),x(i,j,k),x(i,jp,k),x1(i,jm,k),
     1                x1(i,j,k),x1(i,jp,k))
         xy2(i,j)=min(x(i,jm,k),x(i,j,k),x(i,jp,k),x1(i,jm,k),
     1                x1(i,j,k),x1(i,jp,k))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       xy3(i,j)=(xy1(i,j)-x(i,j,k))/(dttk*(max(wmod(i,j,k),0.)*x(i,j,km)
     5                       -min(wmod(i,j,kp),0.)*x(i,j,kp))+eps)
       xy4(i,j)=(x(i,j,k)-xy2(i,j))/(dttk*(max(wmod(i,j,kp),0.)*x(i,j,k)
     5                       -min(wmod(i,j,k),0.)*x(i,j,k))+eps)
        endif
       enddo
  500  continue
 5000 continue

       do 6000 j=1,ny
         jp=j
         jm=j
       do 6000 i=1,nx
         im=i
         ip=i
        do 600 k=2,nz-1
         kp=k+1
         km=k-1
          dttk=dzzt1(k)
         y1(k)=0.
         y2(k)=0.
         y3(k)=0.
         y4(k)=0.
         y5(k)=x(i,j,k)+x(i,j,km)
        if (y5(k) .ge. eps) then
         y1(k)=max(x(i,j,km),x(i,j,k),x(i,j,kp),x1(i,j,km),
     1               x1(i,j,k),x1(i,j,kp))
         y2(k)=min(x(i,j,km),x(i,j,k),x(i,j,kp),x1(i,j,km),
     1               x1(i,j,k),x1(i,j,kp))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        y3(k)=(y1(k)-x(i,j,k))/(dttk*(max(wmod(i,j,k),0.)*x(i,j,km)
     5                       -min(wmod(i,j,kp),0.)*x(i,j,kp))+eps)
        y4(k)=(x(i,j,k)-y2(k))/(dttk*(max(wmod(i,j,kp),0.)*x(i,j,k)
     3                        -min(wmod(i,j,k),0.)*x(i,j,k))+eps)
        endif
  600  continue
       do 650 k=3,nz-1
         km=k-1
        if (y5(k) .ge. eps) then
          w(i,j,k)=min(1.,y4(km),y3(k))*max(0.,wmod(i,j,k))
     1            +min(1.,y4(k),y3(km))*min(0.,wmod(i,j,k))
        endif
  650  continue
 6000 continue

c   *****************************************************************
c$doacross local(j,i)
      do j=1,ny
	 do i=1,nx
	    w(i,j,nz)=0.
	 enddo
      enddo

      return
      end
