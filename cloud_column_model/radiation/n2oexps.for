c**********************************************************************
      subroutine n2oexps(ib,m,np,dn2o,pa,dt,n2oexp)
c**********************************************************************
c   compute n2o exponentials for individual layers 
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  layer n2o amount (dn2o)
c  layer pressure (pa)
c  layer temperature minus 250k (dt)
c
c---- output parameters
c  2 or 4 exponentials for each layer (n2oexp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
c---- input parameters -----
      real dn2o(m,np),pa(m,np),dt(m,np)
c---- output parameters -----
      real n2oexp(m,np,4)
c---- temporary arrays -----
      real xc,xc1,xc2
c-----scaling and absorpton data are given in table 5.
c     transmittances are computed using eqs. (8.18) and (8.19).
       do k=1,np
        do i=1,m
c-----four exponential by powers of 21 for band 6.
          if (ib.eq.6) then
           xc=dn2o(i,k)*(1.+(1.9297e-3+4.3750e-6*dt(i,k))*dt(i,k))
           n2oexp(i,k,1)=exp(-xc*6.31582e-2)
           xc=n2oexp(i,k,1)*n2oexp(i,k,1)*n2oexp(i,k,1)
           xc1=xc*xc
           xc2=xc1*xc1
           n2oexp(i,k,2)=xc*xc1*xc2
c-----four exponential by powers of 8 for band 7
          else
           xc=dn2o(i,k)*(pa(i,k)/500.0)**0.48
     *        *(1.+(1.3804e-3+7.4838e-6*dt(i,k))*dt(i,k))
           n2oexp(i,k,1)=exp(-xc*5.35779e-2)
           xc=n2oexp(i,k,1)*n2oexp(i,k,1)
           xc=xc*xc
           n2oexp(i,k,2)=xc*xc
           xc=n2oexp(i,k,2)*n2oexp(i,k,2)
           xc=xc*xc
           n2oexp(i,k,3)=xc*xc
           xc=n2oexp(i,k,3)*n2oexp(i,k,3)
           xc=xc*xc
           n2oexp(i,k,4)=xc*xc
          endif
        enddo
       enddo
      return
      end
