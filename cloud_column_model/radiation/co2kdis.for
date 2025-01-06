c**********************************************************************
      subroutine co2kdis(m,np,k,co2exp,tco2,tran)
c**********************************************************************
c   compute co2 transmittances between levels k1 and k2 for
c    m soundings, using the k-distribution method with linear
c    pressure scaling.
c
c---- input parameters
c   number of grid intervals (m)
c   number of levels (np)
c   current level (k)
c   exponentials for co2 absorption (co2exp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to co2 absorption
c     for the various values of the absorption coefficient (tco2)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer m,np,i,k
c---- input parameters -----
      real co2exp(m,np,6,2)
c---- updated parameters -----
      real tco2(m,6,2),tran(m)
c---- temporary arrays -----
      real xc
c-----tco2 is the 6 exp factors between levels k1 and k2. 
c     xc is the total co2 transmittance given by eq. (4.30).
c     the k-distribution functions are given in table 10.
         do i=1,m
c-----band-wings
           tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
           xc=   0.1395 *tco2(i,1,1)
           tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
           xc=xc+0.1407 *tco2(i,2,1)
           tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
           xc=xc+0.1549 *tco2(i,3,1)
           tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
           xc=xc+0.1357 *tco2(i,4,1)
           tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
           xc=xc+0.0182 *tco2(i,5,1)
           tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
           xc=xc+0.0220 *tco2(i,6,1)
c-----band-center region
           tco2(i,1,2)=tco2(i,1,2)*co2exp(i,k,1,2)
           xc=xc+0.0766 *tco2(i,1,2)
           tco2(i,2,2)=tco2(i,2,2)*co2exp(i,k,2,2)
           xc=xc+0.1372 *tco2(i,2,2)
           tco2(i,3,2)=tco2(i,3,2)*co2exp(i,k,3,2)
           xc=xc+0.1189 *tco2(i,3,2)
           tco2(i,4,2)=tco2(i,4,2)*co2exp(i,k,4,2)
           xc=xc+0.0335 *tco2(i,4,2)
           tco2(i,5,2)=tco2(i,5,2)*co2exp(i,k,5,2)
           xc=xc+0.0169 *tco2(i,5,2)
           tco2(i,6,2)=tco2(i,6,2)*co2exp(i,k,6,2)
           xc=xc+0.0059 *tco2(i,6,2)
           tran(i)=tran(i)*xc
         enddo
      return
      end
