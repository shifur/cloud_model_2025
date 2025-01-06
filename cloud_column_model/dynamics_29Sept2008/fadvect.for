

      subroutine fadvect (x,y,w,rho,rho1,am,dz)
c     ****   compute advection of terminal velocity
      implicit none

      include 'dimensions.h'
      real    x(nx,ny,nz),y(nx,ny,nz),w(nx,ny,nz)
      real    rho(nz),rho1(nz),am(nz)
      real    dz

cc    local variables
      real    dzr(nz)
      integer i,j,k,kp
      real    dzrk
      save
c     ******   vertical advection for terminal velocity   *************

      do k=1,nz
        dzr(k) = am(k) / ( dz * rho(k) ) ! Should be rho, not rho1...
      enddo

c$doacross local(k,kp,dzrk,j,i)
      do k=2,nz-1
	 kp=k+1
	 dzrk=dzr(k)
	 do j=1,ny
	    do i=1,nx
	       y(i,j,k)=y(i,j,k)+dzrk*(rho1(kp)*w(i,j,kp)*x(i,j,kp)
     1                 -rho1(k)*w(i,j,k)*x(i,j,k))
	    enddo
	 enddo
      enddo
      return
      end
