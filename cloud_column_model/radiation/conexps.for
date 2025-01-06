c**********************************************************************
      subroutine conexps(ib,m,np,dcont,xke,conexp)
c**********************************************************************
c   compute exponentials for continuum absorption in individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  layer scaled water vapor amount for continuum absorption (dcont) 
c  absorption coefficients for the first k-distribution function
c     due to water vapor continuum absorption (xke)
c
c---- output parameters
c  1 or 3 exponentials for each layer (conexp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
c---- input parameters ------
      real dcont(m,np)
c---- updated parameters -----
      real conexp(m,np,3)
c---- static data -----
      real xke(9)
c****************************************************************
        do k=1,np
         do i=1,m
           conexp(i,k,1) = exp(-dcont(i,k)*xke(ib))
         enddo
        enddo
       if (ib .eq. 3) then
c-----the absorption coefficients for sub-bands 3b and 3a are, respectively,
c     two and four times the absorption coefficient for sub-band 3c (table 9).
c     note that conexp(i,k,3) is for sub-band 3a. 
         do k=1,np
          do i=1,m
            conexp(i,k,2) = conexp(i,k,1) *conexp(i,k,1)
            conexp(i,k,3) = conexp(i,k,2) *conexp(i,k,2)
          enddo
         enddo
       endif
      return
      end
