c**********************************************************************
      subroutine ch4exps(ib,m,np,dch4,pa,dt,ch4exp)
c**********************************************************************
c   compute ch4 exponentials for individual layers
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  layer ch4 amount (dch4)
c  layer pressure (pa)
c  layer temperature minus 250k (dt)
c
c---- output parameters
c  1 or 4 exponentials for each layer (ch4exp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
c---- input parameters -----
      real dch4(m,np),pa(m,np),dt(m,np)
c---- output parameters -----
      real ch4exp(m,np,4)
c---- temporary arrays -----
      real xc
c*****  scaling and absorpton data are given in table 5  *****
       do k=1,np
        do i=1,m
c-----four exponentials for band 6
          if (ib.eq.6) then
           xc=dch4(i,k)*(1.+(1.7007e-2+1.5826e-4*dt(i,k))*dt(i,k))
           ch4exp(i,k,1)=exp(-xc*5.80708e-3)
c-----four exponentials by powers of 12 for band 7
          else
           xc=dch4(i,k)*(pa(i,k)/500.0)**0.65
     *       *(1.+(5.9590e-4-2.2931e-6*dt(i,k))*dt(i,k))
           ch4exp(i,k,1)=exp(-xc*6.29247e-2)
           xc=ch4exp(i,k,1)*ch4exp(i,k,1)*ch4exp(i,k,1)
           xc=xc*xc
           ch4exp(i,k,2)=xc*xc
           xc=ch4exp(i,k,2)*ch4exp(i,k,2)*ch4exp(i,k,2)
           xc=xc*xc
           ch4exp(i,k,3)=xc*xc
           xc=ch4exp(i,k,3)*ch4exp(i,k,3)*ch4exp(i,k,3)
           xc=xc*xc
           ch4exp(i,k,4)=xc*xc
          endif
        enddo
       enddo
      return
      end
