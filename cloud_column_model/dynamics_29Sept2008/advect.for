c
cc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine advect (x,x1,w,xxw,ak,dfc,rho,rho1,rrho,qc,
     1                   am,am1,vtp,ismg,dz,dt,dlt1)
c     ****   compute advection of qc, qr, qi, qs and qg
c
      implicit none
      include 'dimensions.h'
c
c   input variables
      integer ismg
      real    dt, dz
      real    vtp(nx,ny,nz)
      real    x(nx,ny,nz)
      real    x1(nx,ny,nz)
      real    rho(nz)
      real    rho1(nz)
      real    rrho(nz)
      real    w(nx,ny,nz)
      real    am(nz)
      real    am1(nz)
      real    xxw(nx,ny,nz)
      real    dfc(nx,ny,nz)
      real    qc(nx,ny,nz)
      real    dlt1
      real    ak(nx,ny,nz)


c   local variables
      real aa(nz)	! hmhj
      real y(nx,ny,nz), y4d(nx,ny,nz), xy3(nx,ny)
      real a0, a1, a2, y3, y4
      real a2k, a0k
      parameter (a2k=0.15e6)

c hmhj
      real yz1(nz),yz2(nz),yz3(nz)
      real tmpd(nx,ny,nz)
      integer i,j,k,km

      save
c
c Zero dummy arrays
      do j=1,ny
         do i=1,nx
            xy3(i,j)=0.0
         enddo
      enddo
      do k=1,nz
        do j=1,ny
          do i=1,nx
            y(i,j,k)=0.
            y4d(i,j,k)=0.
          enddo
        enddo
      enddo

      do k=3,nz
         km=k-1
         a0k=2.0 * a2k
         do j=1,ny
            do i=1,nx
              y4d(i,j,k)=rho1(k)*(-(am1(k)*(x1(i,j,k)-x1(i,j,km))
     1          *(ak(i,j,k)*(1.+xxw(i,j,k))+ak(i,j,km)*(1.+xxw(i,j,km)))
     2          +a0k*(x1(i,j,k)-x1(i,j,km)))/dz)
            enddo
         enddo
      enddo

      do k=3,nz
       km=k-1
       a0k=2.0 * a2k
       do j=1,ny
          do i=1,nx
            y3=qc(i,j,k) * dlt1
            y4=qc(i,j,km) * dlt1
            if(y3 .ge. 1.e-5 .and. y4 .ge. 1.e-5) then
              y4d(i,j,k)=rho1(k)*(dfc(i,j,k)-(a0k/dz)
     1                   *(x1(i,j,k)-x1(i,j,km)))
            endif
          enddo
        enddo
      enddo
c
      do j=1,ny
        do i=1,nx
          y4d(i,j,nz)=0.
        enddo
      enddo
      do k=3,nz
         km=k-1

         do j=1,ny
            do i=1,nx
               y(i,j,km)=(am(km)*rrho(km)/(2.*dz))*(xy3(i,j)-y4d(i,j,k))
               xy3(i,j)=y4d(i,j,k)
            enddo
         enddo

      enddo  ! k
c
c     ******   Terminal Velocity   *********************************
c
      if (ismg.eq.4.or.ismg.eq.5) then
        call fadvect (x,y,vtp,rho,rho1,am,dz)
      endif

ccc      print '(a,2F20.10)','Min/max qx  ',minval(x), maxval(x)

c   Call routine to do positive definite advection
      call advectn (x,x1,y,w,am,am1,rho,dz,dt)

ccc      print '(a,2F20.10)','Min/max qx  ',minval(x), maxval(x)

c      need to save y (forcing) and change x/x1 values

c$doacross local(k,j,i)
      do k=2,nz-1
         do j=1,ny
            do i=1,nx
               x(i,j,k)=x(i,j,k)+y(i,j,k)*dt
            enddo
         enddo
      enddo
c
c
c     ***********************************************************
cc    adjust if there is small negative cloud fields

      if(ismg .gt. 2) then ! if advecting condensate

      do k=2,nz-1
         yz1(k)=0.
         yz2(k)=0.
         yz3(k)=rho(k)*dz/am(k)
      enddo

c$doacross local(k,j,i)
      do k=2,nz-1
         do j=1,ny
            do i=1,nx
               yz1(k)=yz1(k)+max(x(i,j,k),0.0)
               yz2(k)=yz2(k)+max(-x(i,j,k),0.0)
            enddo
         enddo
      enddo

      a1=0.
      a2=0.

      do k=2,nz-1
         a1=a1+yz1(k)*yz3(k)
         a2=a2+yz2(k)*yz3(k)
      enddo

      a0=0.
      if(a1.ne.0.) a0=(a1-a2)/a1

c$doacross local(k,j,i)
      do k=2,nz-1
         do j=1,ny
            do i=1,nx
               x(i,j,k)=max(x(i,j,k),0.0)*a0
            enddo
         enddo
      enddo

cc Zero condensate in bottom and top layer
      do j=1,ny
        do i=1,nx
          x(i,j,1)  = 0.0
          x(i,j,nz) = 0.0
        enddo
      enddo

cc If, in the unlikely event that all values are <= 0, zero the entire array
      if ( maxval(x) .le. 0.0 ) then
        do k=1,nz
          do j=1,ny
            do i=1,nx
              x(i,j,k) = 0.0
            enddo
          enddo
        enddo
      endif

      endif  ! if(ismg .gt. 2) then

ccc      print '(a,2F20.10)','Min/max qx  ',minval(x), maxval(x)


      return
      end

