c**********************************************************************
      subroutine comkdis(ib,m,np,k,comexp,tcom,tran)
c**********************************************************************
c  compute co2-minor transmittances between levels k1 and k2
c   for m soundings, using the k-distribution method
c   with linear pressure scaling.
c
c---- input parameters
c   spectral band (ib)
c   number of grid intervals (m)
c   number of levels (np)
c   current level (k)
c   exponentials for co2-minor absorption (comexp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to co2-minor absorption
c     for the various values of the absorption coefficient (tcom)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
c---- input parameters -----
      real comexp(m,np,6)
c---- updated parameters -----
      real tcom(m,6),tran(m)
c---- temporary arrays -----
      real xc
c-----tcom is computed from eq. (8.20). 
c     xc is the total co2 transmittance computed from (8.22)
c     the k-distribution functions are given in table 6.
         do i=1,m
c-----band 4
           if (ib.eq.4) then
            tcom(i,1)=tcom(i,1)*comexp(i,k,1)
            xc=   0.12159*tcom(i,1)
            tcom(i,2)=tcom(i,2)*comexp(i,k,2)
            xc=xc+0.24359*tcom(i,2)
            tcom(i,3)=tcom(i,3)*comexp(i,k,3)
            xc=xc+0.24981*tcom(i,3)
            tcom(i,4)=tcom(i,4)*comexp(i,k,4)
            xc=xc+0.26427*tcom(i,4)
            tcom(i,5)=tcom(i,5)*comexp(i,k,5)
            xc=xc+0.07807*tcom(i,5)
            tcom(i,6)=tcom(i,6)*comexp(i,k,6)
            xc=xc+0.04267*tcom(i,6)
c-----band 5
           else
            tcom(i,1)=tcom(i,1)*comexp(i,k,1)
            xc=   0.06869*tcom(i,1)
            tcom(i,2)=tcom(i,2)*comexp(i,k,2)
            xc=xc+0.14795*tcom(i,2)
            tcom(i,3)=tcom(i,3)*comexp(i,k,3)
            xc=xc+   0.19512*tcom(i,3)
            tcom(i,4)=tcom(i,4)*comexp(i,k,4)
            xc=xc+   0.33446*tcom(i,4)
            tcom(i,5)=tcom(i,5)*comexp(i,k,5)
            xc=xc+   0.17199*tcom(i,5)
            tcom(i,6)=tcom(i,6)*comexp(i,k,6)
            xc=xc+   0.08179*tcom(i,6)
           endif
            tran(i)=tran(i)*xc
         enddo
      return
      end
