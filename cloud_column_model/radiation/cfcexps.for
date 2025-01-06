c**********************************************************************
      subroutine cfcexps(ib,m,np,a1,b1,fk1,a2,b2,fk2,dcfc,dt,cfcexp)
c**********************************************************************
c   compute cfc(-11, -12, -22) exponentials for individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  parameters for computing the scaled cfc amounts
c             for temperature scaling (a1,b1,a2,b2)
c  the absorption coefficients for the
c     first k-distribution function due to cfcs (fk1,fk2)
c  layer cfc amounts (dcfc)
c  layer temperature minus 250k (dt)
c
c---- output parameters
c  1 exponential for each layer (cfcexp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
c---- input parameters -----
      real dcfc(m,np),dt(m,np)
c---- output parameters -----
      real cfcexp(m,np)
c---- static data -----
      real a1,b1,fk1,a2,b2,fk2
c---- temporary arrays -----
      real xf
c**********************************************************************
       do k=1,np
        do i=1,m
c-----compute the scaled cfc amount (xf) and exponential (cfcexp)
          if (ib.eq.4) then
           xf=dcfc(i,k)*(1.+(a1+b1*dt(i,k))*dt(i,k))
           cfcexp(i,k)=exp(-xf*fk1)
          else
           xf=dcfc(i,k)*(1.+(a2+b2*dt(i,k))*dt(i,k))
           cfcexp(i,k)=exp(-xf*fk2)
          endif
        enddo
       enddo
      return
      end
