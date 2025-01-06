c**********************************************************************
      subroutine h2okdis(ib,m,np,k,fkw,gkw,ne,h2oexp,conexp,
     *                   th2o,tcon,tran)
c**********************************************************************
c   compute water vapor transmittance between levels k1 and k2 for
c   m soundings, using the k-distribution method.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of levels (np)
c  current level (k)
c  planck-weighted k-distribution function due to
c    h2o line absorption (fkw)
c  planck-weighted k-distribution function due to
c    h2o continuum absorption (gkw)
c  number of terms used in each band to compute water vapor
c     continuum transmittance (ne)
c  exponentials for line absorption (h2oexp) 
c  exponentials for continuum absorption (conexp) 
c
c---- updated parameters
c  transmittance between levels k1 and k2 due to
c    water vapor line absorption (th2o)
c  transmittance between levels k1 and k2 due to
c    water vapor continuum absorption (tcon)
c  total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
c---- input parameters ------
      real conexp(m,np,3),h2oexp(m,np,6)
      integer ne(9)
      real  fkw(6,9),gkw(6,3)
c---- updated parameters -----
      real th2o(m,6),tcon(m,3),tran(m)
c---- temporary arrays -----
      real trnth2o
c-----tco2 are the six exp factors between levels k1 and k2 
c     tran is the updated total transmittance between levels k1 and k2
c-----th2o is the 6 exp factors between levels k1 and k2 due to
c     h2o line absorption. 
c-----tcon is the 3 exp factors between levels k1 and k2 due to
c     h2o continuum absorption.
c-----trnth2o is the total transmittance between levels k1 and k2 due
c     to both line and continuum absorption.
c-----comoute th2o following eq. (8.20).
         do i=1,m
           th2o(i,1) = th2o(i,1)*h2oexp(i,k,1)
           th2o(i,2) = th2o(i,2)*h2oexp(i,k,2)
           th2o(i,3) = th2o(i,3)*h2oexp(i,k,3)
           th2o(i,4) = th2o(i,4)*h2oexp(i,k,4)
           th2o(i,5) = th2o(i,5)*h2oexp(i,k,5)
           th2o(i,6) = th2o(i,6)*h2oexp(i,k,6)
         enddo
      if (ne(ib).eq.0) then
c-----comoute trnh2o following eq. (8.22). fkw is given in table 4.
         do i=1,m
           trnth2o      =(fkw(1,ib)*th2o(i,1)
     *                  + fkw(2,ib)*th2o(i,2)
     *                  + fkw(3,ib)*th2o(i,3)
     *                  + fkw(4,ib)*th2o(i,4)
     *                  + fkw(5,ib)*th2o(i,5)
     *                  + fkw(6,ib)*th2o(i,6))
          tran(i)=tran(i)*trnth2o
         enddo
      elseif (ne(ib).eq.1) then
c-----comoute trnh2o following eq. (8.22) and (4.27).
         do i=1,m
           tcon(i,1)= tcon(i,1)*conexp(i,k,1)
           trnth2o      =(fkw(1,ib)*th2o(i,1)
     *                  + fkw(2,ib)*th2o(i,2)
     *                  + fkw(3,ib)*th2o(i,3)
     *                  + fkw(4,ib)*th2o(i,4)
     *                  + fkw(5,ib)*th2o(i,5)
     *                  + fkw(6,ib)*th2o(i,6))*tcon(i,1)
          tran(i)=tran(i)*trnth2o
         enddo
      else
c-----for band 3. this band is divided into 3 subbands.
         do i=1,m
           tcon(i,1)= tcon(i,1)*conexp(i,k,1)
           tcon(i,2)= tcon(i,2)*conexp(i,k,2)
           tcon(i,3)= tcon(i,3)*conexp(i,k,3)
c-----comoute trnh2o following eq. (4.29).
           trnth2o      = (  gkw(1,1)*th2o(i,1)
     *                     + gkw(2,1)*th2o(i,2)
     *                     + gkw(3,1)*th2o(i,3)
     *                     + gkw(4,1)*th2o(i,4)
     *                     + gkw(5,1)*th2o(i,5)
     *                     + gkw(6,1)*th2o(i,6) ) * tcon(i,1)
     *                  + (  gkw(1,2)*th2o(i,1)
     *                     + gkw(2,2)*th2o(i,2)
     *                     + gkw(3,2)*th2o(i,3)
     *                     + gkw(4,2)*th2o(i,4)
     *                     + gkw(5,2)*th2o(i,5)
     *                     + gkw(6,2)*th2o(i,6) ) * tcon(i,2)
     *                  + (  gkw(1,3)*th2o(i,1)
     *                     + gkw(2,3)*th2o(i,2)
     *                     + gkw(3,3)*th2o(i,3)
     *                     + gkw(4,3)*th2o(i,4)
     *                     + gkw(5,3)*th2o(i,5)
     *                     + gkw(6,3)*th2o(i,6) ) * tcon(i,3)
          tran(i)=tran(i)*trnth2o
         enddo
      endif
      return
      end
