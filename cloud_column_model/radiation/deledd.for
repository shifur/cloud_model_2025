c*********************************************************************
      subroutine deledd(m,np,tau,ssc,g0,cza,rr,tt,td)
c*********************************************************************
c
c-----uses the delta-eddington approximation to compute the
c     bulk scattering properties of a single layer
c     coded following king and harshvardhan (jas, 1986)
c
c  inputs:
c       m:  number of soundings
c      np:  number of atmospheric layers
c     tau:  optical thickness
c     ssc:  single scattering albedo
c     g0:   asymmetry factor
c     cza:  cosine of the zenith angle
c
c  outputs:
c
c     rr:  reflection of the direct beam
c     tt:  total (direct+diffuse) transmission of the direct beam
c     td:  direct transmission of the direct beam
c
c*********************************************************************
      implicit none
      real zero,one,two,three,four,fourth,seven,thresh
      parameter (one =1., three=3.)
      parameter (two =2., seven=7.)
      parameter (four=4., fourth=.25)
      parameter (zero=0., thresh=1.e-8)
c-----input parameters
      integer m,np
      real tau(m,np),ssc(m,np),g0(m,np),cza(m)
c-----output parameters
      real rr(m,np),tt(m,np),td(m,np)
c-----temporary parameters
      integer i,k
      real zth,ff,xx,taup,sscp,gp,gm1,gm2,gm3,akk,alf1,alf2,
     *     all,bll,st7,st8,cll,dll,fll,ell,st1,st2,st3,st4
       real taupdzth,akkdtaup
c---------------------------------------------------------------------
      do k=1,np
       do i=1,m
          zth = cza(i)
c  delta-eddington scaling of single scattering albedo,
c  optical thickness, and asymmetry factor,
c  k & h eqs(27-29)
           ff  = g0(i,k)*g0(i,k)
           xx  = one-ff *ssc(i,k)
           taup= tau(i,k)*xx
           sscp= ssc(i,k)*(one-ff)/xx
           gp  = g0(i,k) /(one+g0(i,k))
c  gamma1, gamma2, and gamma3. see table 2 and eq(26) k & h
c  ssc and gp are the d-s single scattering
c  albedo and asymmetry factor.
           xx  =  three*gp 
           gm1 =  (seven - sscp*(four+xx))*fourth
           gm2 = -(one   - sscp*(four-xx))*fourth
c  akk is k as defined in eq(25) of k & h
           akk = sqrt((gm1+gm2)*(gm1-gm2))
           xx  = akk * zth
           st7 = one - xx
           st8 = one + xx
           st3 = st7 * st8
           if (abs(st3) .lt. thresh) then
               zth = zth + 0.001
               xx  = akk * zth
               st7 = one - xx
               st8 = one + xx
               st3 = st7 * st8
           endif
c  extinction of the direct beam transmission
           td(i,k)=0.
           taupdzth=taup/zth
           if (taupdzth .lt. 40. ) td(i,k)  = exp(-taup/zth)
c  alf1 and alf2 are alpha1 and alpha2 from eqs (23) & (24) of k & h
           gm3  = (two - zth*three*gp)*fourth
           xx   = gm1 - gm2
           alf1 = gm1 - gm3 * xx
           alf2 = gm2 + gm3 * xx
c  all is last term in eq(21) of k & h
c  bll is last term in eq(22) of k & h
           xx  = akk * two
           all = (gm3 - alf2 * zth    )*xx*td(i,k)
           bll = (one - gm3 + alf1*zth)*xx
           xx  = akk * gm3
           cll = (alf2 + xx) * st7
           dll = (alf2 - xx) * st8
           xx  = akk * (one-gm3)
           fll = (alf1 + xx) * st8
           ell = (alf1 - xx) * st7
           st2=0.
           akkdtaup=akk*taup
           if (akkdtaup.lt.40.) st2 = exp(-akkdtaup)
           st4 = st2 * st2
           st1 =  sscp / ((akk+gm1 + (akk-gm1)*st4) * st3)
c  rr is r-hat of eq(21) of k & h
c  tt is diffuse part of t-hat of eq(22) of k & h
           rr(i,k) =   ( cll-dll*st4         -all*st2)*st1
           tt(i,k) = - ((fll-ell*st4)*td(i,k)-bll*st2)*st1
           rr(i,k) = max(rr(i,k),zero)
           tt(i,k) = max(tt(i,k),zero)
           tt(i,k) = tt(i,k)+td(i,k)
        enddo
       enddo
      return
      end
