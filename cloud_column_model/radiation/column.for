c***********************************************************************
      subroutine column (m,np,pa,dt,sabs0,sabs,spre,stem)
c***********************************************************************
c-----compute column-integrated (from top of the model atmosphere)
c     absorber amount (sabs), absorber-weighted pressure (spre) and
c     temperature (stem).
c     computations follow eqs. (8.24) - (8.26).
c
c--- input parameters
c   number of soundings (m)
c   number of atmospheric layers (np)
c   layer pressure (pa)
c   layer temperature minus 250k (dt)
c   layer absorber amount (sabs0)
c
c--- output parameters
c   column-integrated absorber amount (sabs)
c   column absorber-weighted pressure (spre)
c   column absorber-weighted temperature (stem)
c
c--- units of pa and dt are mb and k, respectively.
c    units of sabs are g/cm**2 for water vapor and (cm-atm)stp
c    for co2 and o3
c***********************************************************************
      implicit none
      integer m,np,i,k
c---- input parameters -----
      real pa(m,np),dt(m,np),sabs0(m,np)
c---- output parameters -----
      real sabs(m,np+1),spre(m,np+1),stem(m,np+1)
c*********************************************************************
        do i=1,m
          sabs(i,1)=0.0
          spre(i,1)=0.0
          stem(i,1)=0.0
        enddo
        do k=1,np
         do i=1,m
           sabs(i,k+1)=sabs(i,k)+sabs0(i,k)
           spre(i,k+1)=spre(i,k)+pa(i,k)*sabs0(i,k)
           stem(i,k+1)=stem(i,k)+dt(i,k)*sabs0(i,k)
         enddo
        enddo
       return
       end
