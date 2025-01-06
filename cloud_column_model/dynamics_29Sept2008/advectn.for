c    4-6 2005  TAO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine advectn (x,x1,y,ww1,am,am1,rho1,dz,dt)
      implicit none
      include 'dimensions.h'
c     ****   compute advection for various scalar variables
      real    ww1(nx,ny,nz)
      real    x(nx,ny,nz),x1(nx,ny,nz),y(nx,ny,nz)
      real    am(nz), am1(nz), rho1(nz)
      real    dz, dt

ccc   local variables
!       common /tmp1/ u1(nx,ny,nz),v1(nx,ny,nz),w1(nx,ny,nz)
      real    w1(nx,ny,nz)
      integer i,j,k
      save
cc    **************************************************************
c$doacross local(k,j,i)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            w1(i,j,k)=ww1(i,j,k)
          enddo
        enddo
      enddo
c
c     positive definite scheme

c$doacross local(k,j,i)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            x1(i,j,k)=x(i,j,k)
          enddo
        enddo
      enddo

c     set x_temp=x  and x1_temp=x1  and store forcing y

ccc      print '(a,2F20.10)','1Min/max qx ',minval(x), maxval(x)
      call fadv (x,y,rho1,am,dz,w1,dt)
ccc      print '(a,2F20.10)','2Min/max qx ',minval(x), maxval(x)
      call fadvuw (x,x1,w1,am,am1,dz,dt)
ccc      print '(a,2F20.10)','3Min/max qx ',minval(x), maxval(x)
      call fadv (x,y,rho1,am,dz,w1,dt)
ccc      print '(a,2F20.10)','4Min/max qx ',minval(x), maxval(x)
c
      return
      end

