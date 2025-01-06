c**********************************************************************
      subroutine co2exps(m,np,dco2,pa,dt,co2exp)
c**********************************************************************
c   compute co2 exponentials for individual layers.
c
c---- input parameters
c  number of grid intervals (m)
c  number of layers (np)
c  layer co2 amount (dco2)
c  layer pressure (pa)
c  layer temperature minus 250k (dt)
c
c---- output parameters
c  6 exponentials for each layer (co2exp)
c**********************************************************************
      implicit none
      integer m,np,i,k
c---- input parameters -----
      real dco2(m,np),pa(m,np),dt(m,np)
c---- output parameters -----
      real co2exp(m,np,6,2)
c---- temporary arrays -----
      real xc
c**********************************************************************
        do k=1,np
         do i=1,m
c-----the scakubg oaraneters are given in table 3, and values of
c     the absorption coefficient are given in table 10.
c     scaled co2 amount for band-wings (sub-bands 3a and 3c)
           xc = dco2(i,k)*(pa(i,k)/300.0)**0.5
     1             *(1.+(0.0182+1.07e-4*dt(i,k))*dt(i,k))
c-----six exponentials by powers of 8 (see eqs. 8.18, 8.19 and table 10).
           co2exp(i,k,1,1)=exp(-xc*2.656e-5)
           xc=co2exp(i,k,1,1)*co2exp(i,k,1,1)
           xc=xc*xc
           co2exp(i,k,2,1)=xc*xc
           xc=co2exp(i,k,2,1)*co2exp(i,k,2,1)
           xc=xc*xc
           co2exp(i,k,3,1)=xc*xc
           xc=co2exp(i,k,3,1)*co2exp(i,k,3,1)
           xc=xc*xc
           co2exp(i,k,4,1)=xc*xc
           xc=co2exp(i,k,4,1)*co2exp(i,k,4,1)
           xc=xc*xc
           co2exp(i,k,5,1)=xc*xc
           xc=co2exp(i,k,5,1)*co2exp(i,k,5,1)
           xc=xc*xc
           co2exp(i,k,6,1)=xc*xc
c-----for band-center region (sub-band 3b)
           xc = dco2(i,k)*(pa(i,k)/30.0)**0.85
     1             *(1.+(0.0042+2.00e-5*dt(i,k))*dt(i,k))
           co2exp(i,k,1,2)=exp(-xc*2.656e-3)
           xc=co2exp(i,k,1,2)*co2exp(i,k,1,2)
           xc=xc*xc
           co2exp(i,k,2,2)=xc*xc
           xc=co2exp(i,k,2,2)*co2exp(i,k,2,2)
           xc=xc*xc
           co2exp(i,k,3,2)=xc*xc
           xc=co2exp(i,k,3,2)*co2exp(i,k,3,2)
           xc=xc*xc
           co2exp(i,k,4,2)=xc*xc
           xc=co2exp(i,k,4,2)*co2exp(i,k,4,2)
           xc=xc*xc
           co2exp(i,k,5,2)=xc*xc
           xc=co2exp(i,k,5,2)*co2exp(i,k,5,2)
           xc=xc*xc
           co2exp(i,k,6,2)=xc*xc
         enddo
        enddo
      return
      end
