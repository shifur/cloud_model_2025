c**********************************************************************
      subroutine tablup(k1,k2,m,np,nx,nh,sabs,spre,stem,w1,p1,
     *                  dwe,dpe,coef1,coef2,coef3,tran)
c**********************************************************************
c   compute water vapor, co2 and o3 transmittances between level
c   k1 and and level k2 for m soundings, using table look-up.
c
c   calculations follow eq. (4.16).
c
c---- input ---------------------
c  indices for layer (k1) and level (k2)
c  number of grid intervals (m)
c  number of atmospheric layers (np)
c  number of pressure intervals in the table (nx)
c  number of absorber amount intervals in the table (nh)
c  column-integrated absorber amount (sabs)
c  column absorber amount-weighted pressure (spre)
c  column absorber amount-weighted temperature (stem)
c  first value of absorber amount (log10) in the table (w1) 
c  first value of pressure (log10) in the table (p1) 
c  size of the interval of absorber amount (log10) in the table (dwe)
c  size of the interval of pressure (log10) in the table (dpe)
c  pre-computed coefficients (coef1, coef2, and coef3)
c
c---- updated ---------------------
c  transmittance (tran)
c
c  note:
c   (1) units of sabs are g/cm**2 for water vapor and
c       (cm-atm)stp for co2 and o3.
c   (2) units of spre and stem are, respectively, mb and k.
c   
c**********************************************************************
      implicit none
      integer k1,k2,m,np,nx,nh,i
c---- input parameters -----
      real w1,p1,dwe,dpe
      real sabs(m,np+1),spre(m,np+1),stem(m,np+1)
      real coef1(nx,nh),coef2(nx,nh),coef3(nx,nh)
c---- update parameter -----
      real tran(m)
c---- temporary variables -----
      real x1,x2,x3,we,pe,fw,fp,pa,pb,pc,ax,ba,bb,t1,ca,cb,t2
      integer iw,ip
c**********************************************************************
      do i=1,m
        x1=sabs(i,k2)-sabs(i,k1)
        we=(log10(x1)-w1)/dwe
       if (we .ge. (w1-2.)) then
        x2=(spre(i,k2)-spre(i,k1))/x1
        x3=(stem(i,k2)-stem(i,k1))/x1
c-----normalize we and pe
        pe=(log10(x2)-p1)/dpe
c-----restrict the magnitudes of the normalized we and pe.
        we=min(we,real(nh-1))
        pe=max(pe,0.0)
        pe=min(pe,real(nx-1))
c-----assign iw and ip and compute the distance of we and pe 
c     from iw and ip.
        iw=int(we+1.0)
        iw=min(iw,nh-1)
        iw=max(iw, 2)
        fw=we-float(iw-1)
        ip=int(pe+1.0)
        ip=min(ip,nx-1)
        ip=max(ip, 1)
        fp=pe-float(ip-1)
c-----linear interpolation in pressure
        pa = coef1(ip,iw-1)*(1.-fp)+coef1(ip+1,iw-1)*fp
        pb = coef1(ip,  iw)*(1.-fp)+coef1(ip+1,  iw)*fp
        pc = coef1(ip,iw+1)*(1.-fp)+coef1(ip+1,iw+1)*fp
c-----quadratic interpolation in absorber amount for coef1
        ax = (-pa*(1.-fw)+pc*(1.+fw)) *fw*0.5 + pb*(1.-fw*fw)
c-----linear interpolation in absorber amount for coef2 and coef3
        ba = coef2(ip,  iw)*(1.-fp)+coef2(ip+1,  iw)*fp
        bb = coef2(ip,iw+1)*(1.-fp)+coef2(ip+1,iw+1)*fp
        t1 = ba*(1.-fw) + bb*fw
        ca = coef3(ip,  iw)*(1.-fp)+coef3(ip+1,  iw)*fp
        cb = coef3(ip,iw+1)*(1.-fp)+coef3(ip+1,iw+1)*fp
        t2 = ca*(1.-fw) + cb*fw
c-----update the total transmittance between levels k1 and k2
        tran(i)= (ax + (t1+t2*x3) * x3)*tran(i)
        tran(i)=min(tran(i),0.9999999)
        tran(i)=max(tran(i),0.0000001)
       else
        tran(i)=0.9999999
       endif
      enddo
      return
      end
