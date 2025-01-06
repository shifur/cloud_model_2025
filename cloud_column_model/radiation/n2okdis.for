c**********************************************************************
      subroutine n2okdis(ib,m,np,k,n2oexp,tn2o,tran)
c**********************************************************************
c   compute n2o transmittances between levels k1 and k2 for
c    m soundings, using the k-distribution method with linear
c    pressure scaling.
c
c---- input parameters
c   spectral band (ib)
c   number of grid intervals (m)
c   number of levels (np)
c   current level (k)
c   exponentials for n2o absorption (n2oexp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to n2o absorption
c     for the various values of the absorption coefficient (tn2o)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
c---- input parameters -----
      real n2oexp(m,np,4)
c---- updated parameters -----
      real tn2o(m,4),tran(m)
c---- temporary arrays -----
      real xc
c-----tn2o is computed from eq. (8.20). 
c     xc is the total n2o transmittance computed from (8.22)
c     the k-distribution functions are given in table 5.
        do i=1,m
c-----band 6
          if (ib.eq.6) then
           tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
           xc=   0.940414*tn2o(i,1)
           tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
           xc=xc+0.059586*tn2o(i,2)
c-----band 7
          else
           tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
           xc=   0.561961*tn2o(i,1)
           tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
           xc=xc+0.138707*tn2o(i,2)
           tn2o(i,3)=tn2o(i,3)*n2oexp(i,k,3)
           xc=xc+0.240670*tn2o(i,3)
           tn2o(i,4)=tn2o(i,4)*n2oexp(i,k,4)
           xc=xc+0.058662*tn2o(i,4)
          endif
           tran(i)=tran(i)*xc
        enddo
      return
      end
