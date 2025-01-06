c
cc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fadv (x,y,rho1,am,dz,w,dt)
c     ****   compute advection of different tracers
c Note: rho1 should be the density on layer boundaries...
      implicit none
      include 'dimensions.h'

      real    w(nx,ny,nz),x(nx,ny,nz)
      real    y(nx,ny,nz)
      real    rho1(nz), am(nz)
      real    dzr(nz)
      real    dt, dz

c     local variables
      integer nm
      parameter (nm=nx+nz)
      real    xy1(nx,ny),xy2(nx,ny),xy3(nx,ny),
     1   y3(nm),y4(nm)

      real xy7(nx,ny), xy8(nx,ny), xy9(nx,ny)
      integer i, j, k, kp
      save

c     ******   horizontal and vertical advection terms   ***************

      do k=1,nz
        dzr(k) = am(k) / ( dz * rho1(k) )
      enddo

       do 10 j=1,ny
       do 10 i=1,nx
c   10	    xy1(i,j)=x(i,j,1)
   10	    xy1(i,j)=0.0

      do 1000 k=2,nz-1
        kp=k+1
       do 100 j=1,ny
       do 100 i=1,nx
  100    xy2(i,j)=x(i,j,k)

       do 150 j=1,ny
       do 150 i=1,nx
       y3(i)=x(i,j,kp)
       y4(i)=x(i,j,k)
        if (w(i,j,kp) .ge. 0.0) y3(i)=x(i,j,k)
        if (w(i,j,k) .gt. 0.0) y4(i)=xy1(i,j)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      xy7(i,j)=-dzr(k)*(rho1(kp)*w(i,j,kp)*y3(i)-rho1(k)*w(i,j,k)*y4(i))
      xy8(i,j)=0.0
      xy9(i,j)=0.0
  150 continue
       do 200 j=1,ny
       do 200 i=1,nx
        x(i,j,k)=x(i,j,k)+(xy7(i,j)+xy8(i,j)+xy9(i,j))*dt
  200   xy1(i,j)=xy2(i,j)

 1000 continue

      return
      end
