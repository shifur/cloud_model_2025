c     4-6 2005   TAO
c
      subroutine satdt_ak (ak,ak1,ak_f,dt)
      implicit none
      include 'dimensions.h'

      real dt
      real ak_f(nx,ny,nz), ak(nx,ny,nz), ak1(nx,ny,nz)
      real a1
      parameter (a1=2.0e6)
      integer i,j,k

c     ****************************************************************
c
       do k=2,nz-1
c$doacross local(j,i)
	    do j=1,ny
	    do i=1,nx
	       ak(i,j,k)=ak(i,j,k)+ak_f(i,j,k)*dt
              ak1(i,j,k)=ak(i,j,k)
	    enddo
	    enddo

c$doacross local(j,i)
	  do j=1,ny
	  do i=1,nx
	     ak(i,j,k)=min(a1, max(ak(i,j,k), 0.))
	     ak1(i,j,k)=min(a1, max(ak1(i,j,k), 0.))
	  enddo
	  enddo
       enddo

c   *****************************************************************
       return
       end

