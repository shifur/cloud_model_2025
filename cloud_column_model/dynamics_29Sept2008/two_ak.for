      subroutine two_ak (ak,ak2)
      implicit none
      include 'dimensions.h'

c     this needs to be called after finishing all u,v and w
c     ***  set two time ak for scalar variables *********************

      real ak(nx,ny,nz)
      real ak2(nx,ny,nz)
      integer i,j,k

c$doacross local(k,j,i)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         ak2(i,j,k)=2.*ak(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
