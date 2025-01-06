c**********************************************************************
      subroutine b10kdis(m,np,k,h2oexp,conexp,co2exp,n2oexp
     *          ,th2o,tcon,tco2,tn2o,tran)
c**********************************************************************
c
c   compute h2o (line and continuum),co2,n2o transmittances between
c   levels k1 and k2 for m soundings, using the k-distribution
c   method with linear pressure scaling.
c
c---- input parameters
c   number of grid intervals (m)
c   number of levels (np)
c   current level (k)
c   exponentials for h2o line absorption (h2oexp)
c   exponentials for h2o continuum absorption (conexp)
c   exponentials for co2 absorption (co2exp)
c   exponentials for n2o absorption (n2oexp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to h2o line absorption
c     for the various values of the absorption coefficient (th2o)
c   transmittance between levels k1 and k2 due to h2o continuum
c     absorption for the various values of the absorption
c     coefficient (tcon)
c   transmittance between levels k1 and k2 due to co2 absorption
c     for the various values of the absorption coefficient (tco2)
c   transmittance between levels k1 and k2 due to n2o absorption
c     for the various values of the absorption coefficient (tn2o)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer m,np,i,k
c---- input parameters -----
      real h2oexp(m,np,6),conexp(m,np,3),co2exp(m,np,6,2)
     *    ,n2oexp(m,np,4)
c---- updated parameters -----
      real th2o(m,6),tcon(m,3),tco2(m,6,2),tn2o(m,4)
     *    ,tran(m)
c---- temporary arrays -----
      real xx
c-----for h2o line. the k-distribution functions are given in table 4.
        do i=1,m
           th2o(i,1)=th2o(i,1)*h2oexp(i,k,1)
           xx=   0.3153*th2o(i,1)
           th2o(i,2)=th2o(i,2)*h2oexp(i,k,2)
           xx=xx+0.4604*th2o(i,2)
           th2o(i,3)=th2o(i,3)*h2oexp(i,k,3)
           xx=xx+0.1326*th2o(i,3)
           th2o(i,4)=th2o(i,4)*h2oexp(i,k,4)
           xx=xx+0.0798*th2o(i,4)
           th2o(i,5)=th2o(i,5)*h2oexp(i,k,5)
           xx=xx+0.0119*th2o(i,5)
           tran(i)=xx
        enddo
c-----for h2o continuum. note that conexp(i,k,3) is for subband 3a.
        do i=1,m
           tcon(i,1)=tcon(i,1)*conexp(i,k,1)
           tran(i)=tran(i)*tcon(i,1)
        enddo
c-----for co2 (table 6)
        do i=1,m
           tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
           xx=    0.2673*tco2(i,1,1)
           tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
           xx=xx+ 0.2201*tco2(i,2,1)
           tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
           xx=xx+ 0.2106*tco2(i,3,1)
           tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
           xx=xx+ 0.2409*tco2(i,4,1)
           tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
           xx=xx+ 0.0196*tco2(i,5,1)
           tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
           xx=xx+ 0.0415*tco2(i,6,1)
           tran(i)=tran(i)*xx
        enddo
c-----for n2o (table 5)
        do i=1,m
           tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
           xx=   0.970831*tn2o(i,1)
           tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
           xx=xx+0.029169*tn2o(i,2)
           tran(i)=tran(i)*(xx-1.0)
        enddo
      return
      end
