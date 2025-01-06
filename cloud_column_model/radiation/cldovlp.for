c***********************************************************************
      subroutine cldovlp (m,np,k2,ict,icb,it,im,ib,
     *               cldhi,cldmd,cldlw,fcld,tcldlyr,fclr)
c***********************************************************************
c     compute the fractional clear line-of-sight between levels k1
c     and k2
c
c input parameters
c
c  m:       number of soundings
c  np:      number of layers
c  k2:      index for the level
c  ict:     the level separating high and middle clouds
c  icb:     the level separating middle and low clouds
c  it:      number of cloudy layers in the high-cloud group
c  im:      number of cloudy layers in the middle-cloud group
c  ib:      number of cloudy layers in the low-cloud group
c  fcld:    fractional cloud cover of a layer
c  tcldlyr: transmittance of a cloud layer
c  
c output parameter
c
c  fclr:    clear line-of-sight between levels k1 and k2
c***********************************************************************
      implicit none
      integer m,np,k2,ict,icb
      integer i,j,k,ii,it(m),im(m),ib(m),itx(m,np),imx(m,np),ibx(m,np)
      real cldhi(m),cldmd(m),cldlw(m)
      real fcld(m,np),tcldlyr(m,np),fclr(m)
c***********************************************************************
       do i=1,m
c-----for high clouds
c     "it" is the number of high-cloud layers
        if (k2.le.ict) then
         if(fcld(i,k2-1).gt.0.001) then
          it(i)=it(i)+1
          ii=it(i)
          itx(i,ii)=k2-1
         if (ii .eq. 1) go to 11
c-----rearrange the order of cloud layers with increasing cloud amount
         do k=1,ii-1
           j=itx(i,k)
          if(fcld(i,j).gt.fcld(i,k2-1)) then
           do j=ii-1,k,-1
            itx(i,j+1)=itx(i,j)
           enddo
            itx(i,k)=k2-1
            go to 11
          endif
         enddo
   11   continue
c-----compute equivalent black-body high cloud amount
           cldhi(i)=0.0
          do k=1,ii
           j=itx(i,k)
           cldhi(i)=fcld(i,j)-tcldlyr(i,j)*(fcld(i,j)-cldhi(i))
          enddo
        endif
       endif
c-----for middle clouds
c     "im" is the number of middle-cloud layers
       if (k2.gt.ict .and. k2.le.icb) then
        if(fcld(i,k2-1).gt.0.001) then
         im(i)=im(i)+1
         ii=im(i)
         imx(i,ii)=k2-1
        if (ii .eq. 1) go to 21
c-----rearrange the order of cloud layers with increasing cloud amount
         do k=1,ii-1
            j=imx(i,k)
           if(fcld(i,j).gt.fcld(i,k2-1)) then
            do j=ii-1,k,-1
             imx(i,j+1)=imx(i,j)
            enddo
             imx(i,k)=k2-1
             go to 21
           endif
          enddo
   21   continue
c-----compute equivalent black-body middle cloud amount
           cldmd(i)=0.0
          do k=1,ii
           j=imx(i,k)
           cldmd(i)=fcld(i,j)-tcldlyr(i,j)*(fcld(i,j)-cldmd(i))
          enddo
        endif
       endif
c-----for low clouds
c     "ib" is the number of low-cloud layers
       if (k2.gt.icb) then
        if(fcld(i,k2-1).gt.0.001) then
         ib(i)=ib(i)+1
         ii=ib(i)
         ibx(i,ii)=k2-1
        if (ii .eq. 1) go to 31
c-----rearrange the order of cloud layers with increasing cloud amount
         do k=1,ii-1
          j=ibx(i,k)
           if(fcld(i,j).gt.fcld(i,k2-1)) then
            do j=ii-1,k,-1
             ibx(i,j+1)=ibx(i,j)
            enddo
             ibx(i,k)=k2-1
             go to 31
           endif
          enddo
   31    continue
c-----compute equivalent black-body low cloud amount
           cldlw(i)=0.0
          do k=1,ii
           j=ibx(i,k)
           cldlw(i)=fcld(i,j)-tcldlyr(i,j)*(fcld(i,j)-cldlw(i))
          enddo
        endif
       endif
c-----fclr is the equivalent clear fraction between levels k1 and k2
c     assuming the three cloud groups are randomly overlapped.
c     it follows eqs. (10) and (12).
        fclr(i)=(1.0-cldhi(i))*(1.0-cldmd(i))*(1.0-cldlw(i))   
      enddo
      return
      end
