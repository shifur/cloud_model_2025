c     4-6 2005  TAO
      subroutine advectak (x,x1,ak_f,ww1,am,am1,rho1,dz,dt)
      implicit none
      include 'dimensions.h'

      real    x(nx,ny,nz),x1(nx,ny,nz),ak_f(nx,ny,nz),ww1(nx,ny,nz)
      real    am(nz), am1(nz), rho1(nz)
      real    dz, dt

cc    local variables
      real    w1(nx,ny,nz)
      integer i,j,k
      save
cc    **************************************************************
c$doacross local(k,j,i)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         if(x(i,j,k) .lt. 1.e-5) x(i,j,k)=0.
         x1(i,j,k)=x(i,j,k)
         w1(i,j,k)=ww1(i,j,k)
      enddo
      enddo
      enddo
c
cc    need to save the forcing
c
      call fadv (x,ak_f,rho1,am,dz,w1,dt)
      call fadvuw (x,x1,w1,am,am1,dz,dt)
      call fadv (x,ak_f,rho1,am,dz,w1,dt)
c
c$doacross local(k,j,i)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         if(x(i,j,k) .lt. 1.e-5) x(i,j,k)=0.
      enddo
      enddo
      enddo

      return
      end

