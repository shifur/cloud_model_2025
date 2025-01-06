c**********************************************************************
      subroutine comexps(ib,m,np,dcom,dt,comexp)
c**********************************************************************
c   compute co2-minor exponentials for individual layers using 
c   eqs. (8.18) and (8.19).
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  layer co2 amount (dcom)
c  layer temperature minus 250k (dt)
c
c---- output parameters
c  6 exponentials for each layer (comexp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k,ik
c---- input parameters -----
      real dcom(m,np),dt(m,np)
c---- output parameters -----
      real comexp(m,np,6)
c---- temporary arrays -----
      real xc
c*****  scaling and absorpton data are given in table 6  *****
       do k=1,np
        do i=1,m
          if (ib.eq.4) then
           xc=dcom(i,k)*(1.+(3.5775e-2+4.0447e-4*dt(i,k))*dt(i,k))
          endif
          if (ib.eq.5) then
           xc=dcom(i,k)*(1.+(3.4268e-2+3.7401e-4*dt(i,k))*dt(i,k))
          endif
           comexp(i,k,1)=exp(-xc*1.922e-7)
          do ik=2,6
           xc=comexp(i,k,ik-1)*comexp(i,k,ik-1)
           xc=xc*xc
           comexp(i,k,ik)=xc*comexp(i,k,ik-1)
          enddo
        enddo
       enddo
      return
      end
