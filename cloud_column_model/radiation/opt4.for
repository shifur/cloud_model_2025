cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine opt4 (iopcloud,iradave,nx1,twcz,rho,dz0,
     +                 tausw,tauir,
c D.Posselt add J index for RAMS variables
     +                fcld,taq,waq,oaq,plq,cwc,reff,jj)
c     define variables and calculate the optical thickness             c
c
      include 'dimensions.h'
      include 'radiation.h'
c     parameter(nx=66,ny=10,nz=34,kl2=nz-2,lay=88)
c     parameter(nadd=7)
!       parameter (nxf=4000,nyf=75,nzf=37,nt=38640,itt=244) 
!       parameter (lb=2,kb=1)
! ! define decomposition from rmp_switch.h
!       parameter (npes=1,ncol=1,nrow=1)
! ! define partial dimension for computation by decomposition
!       parameter (nx=(nxf-lb*2-1)/ncol+1+lb*2)
!       parameter (ny=1)
!       parameter (nz=nzf)
! ! define partial dimension for fft by transpose decomposition
!       parameter (nyc=1)
!       parameter (nyr= ny)
!       parameter (nxc= nx)
!       parameter (nxr=(nxf-lb*2-1)/nrow+1+lb*2)
!       parameter (nzc=(nzf-kb*2-1)/ncol+1+kb*2)
!       parameter (nzr=(nzf-kb*2-1)/nrow+1+kb*2)
!
!        parameter (nadd=7,lay=88)
!       parameter(kmax=nzf,kles=kmax-1,kl2=kles-1)
      parameter(npp1=nz+nadd-1,np=nz-2+nadd ,nq=npp1,nw=nq-1)
c
      integer iopcloud,nx1
      real    rho(nz),dz0(nz)

      real    fcld(nx1,1,nw)
      real    plq(nx1,1,nq),taq(nx1,1,nw),waq(nx1,1,nw),oaq(nx1,1,nw)
      real    tausw(nx1,1,nw,3),tauir(nx1,1,nw,3)
      real    cwc(nx1,1,nw,3),reff(nx1,1,nw,3)
c      real    taucld(nx1,1,np,3)
      integer iradave
      real    tnw,tns,tng,roqs,roqg,roqr
      common/size/ tnw,tns,tng,roqs,roqg,roqr
      common/radtemp1/ qc(nx,nz),qr(nx,nz),qi(nx,nz),qs(nx,nz),
     1 qg(nx,nz),pt(nx,nz),qv(nx,nz),qh(nx,nz),qsw(nx,nz),qgw(nx,nz),
     2 qhw(nx,nz)
c-----------------------------------------------------------------------
c Add RAMS variables
c-----------------------------------------------------------------------
! Integer index in y-direction
      integer jj
! Mixing ratios
      real rcp(nx,ny,nz),rrp(nx,ny,nz),rpp(nx,ny,nz),rsp(nx,ny,nz), 
     1     rap(nx,ny,nz),rgp(nx,ny,nz),rhp(nx,ny,nz),rc2p(nx,ny,nz)
      common/rams_mixr1/ rcp,rrp,rpp,rsp,rap,rgp,rhp,rc2p
! Variables from RAMS micro used in radiation
      real gce_reff(nx,ny,nz,3)
      real rams_reff(nx,ny,nz,8)
      common/rams_reff/ gce_reff, rams_reff
! Flags
      integer i_rams
      integer jnmb(8)
      common/rams_parms/ i_rams, jnmb
! Holders for optical depth and mixing ratios
      real tau_rams(8)
      real q4,q6,q7,q8
c-----------------------------------------------------------------------
c End RAMS variables
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c local variables
c-----------------------------------------------------------------------
      real    aa1(nx,npp1),bb(nx,npp1),aa2(nx,npp1)
      real    y1(nw),y2(nw),y3(nw),y4(nw),y5(nw),y6(nw),y7(nw),y8(nw),
     1        y9(nw),y10(nw),y11(nw)
c     integer i,k,km,kmax,kles
      integer i,k,km               !!! cccshie 7/2/02
      real    b0,b1,b2,b4,b5,b6,b7,b9,cpi,effrad
      real    q1,q2,q3,q5,rnx1,tauqc,tauqg,tauqii,tauqis,tauqr,twco
c D.Posselt Add print flag for diagnostic output
      logical print_more
      save

c D.Posselt 7/28/2008 Define variables
      kmax = nz
      kles = nz - 1
      kl2  = nz - 2
      i_rams = 0
c
       cpi =4.*atan(1.)
c
cc
c
      print_more = .false.
c D.Posselt Diagnostic print prior to computation
      if ( i_rams .eq. 1 .and. print_more) then
        print*,'In opt4, RAMS variables: '
c        call tmax (rcp,'RCP ',0.0,wmax,wmin,3)
c        call tmax (rc2p,'RC2P',0.0,wmax,wmin,3)
c        call tmax (rrp,'RRP ',0.0,wmax,wmin,3)
c        call tmax (rpp,'RPP ',0.0,wmax,wmin,3)
c        call tmax (rsp,'RSP ',0.0,wmax,wmin,3)
c        call tmax (rap,'RAP ',0.0,wmax,wmin,3)
c        call tmax (rgp,'RGP ',0.0,wmax,wmin,3)
c        call tmax (rhp,'RHP ',0.0,wmax,wmin,3)
c        call tmax (rams_reff(:,:,:,1),'CRe ',0.0,wmax,wmin,2)
c        call tmax (rams_reff(:,:,:,8),'C2Re',0.0,wmax,wmin,2)
c        call tmax (rams_reff(:,:,:,2),'RRe ',0.0,wmax,wmin,2)
c        call tmax (rams_reff(:,:,:,3),'PRe ',0.0,wmax,wmin,2)
c        call tmax (rams_reff(:,:,:,4),'SRe ',0.0,wmax,wmin,2)
c        call tmax (rams_reff(:,:,:,5),'ARe ',0.0,wmax,wmin,2)
c        call tmax (rams_reff(:,:,:,6),'GRe ',0.0,wmax,wmin,2)
c        call tmax (rams_reff(:,:,:,7),'HRe ',0.0,wmax,wmin,2)
c        call tmax (gce_reff(:,:,:,1),'Re1 ',0.0,wmax,wmin,2)
c        call tmax (gce_reff(:,:,:,2),'Re2 ',0.0,wmax,wmin,2)
c        call tmax (gce_reff(:,:,:,3),'Re3 ',0.0,wmax,wmin,2)
      endif
       do 200 i=1,nx1
         do k=1,nw
           tauir(i,1,k,1)=0.
           tauir(i,1,k,2)=0.
           tauir(i,1,k,3)=0.
           tausw(i,1,k,1)=0.
           tausw(i,1,k,2)=0.
           tausw(i,1,k,3)=0.
           reff(i,1,k,1)=0.
           reff(i,1,k,2)=0.
           reff(i,1,k,3)=0.
           cwc(i,1,k,1)=0.
           cwc(i,1,k,2)=0.
           cwc(i,1,k,3)=0.
c           taucld(i,1,k,1)=0.   
c           taucld(i,1,k,2)=0.
c           taucld(i,1,k,3)=0.
         enddo
         b5=0.
         b6=0.
         b7=0.
         b9=0.
         do 100 k=2,kles
           km=kmax-k+nadd
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc D.Posselt 8/22/2006 If-test for RAMS microphysics
ccc In RAMS micro, effective radius is computed for each species. These
ccc effective radii are stored in rams_reff, and are used to compute 
ccc optical depth for each species.
ccc
ccc For consistency with the Chou radiation code, input effective radius
ccc is defined for only three categories -- ice, cloud, and rain. Hence,
ccc the effective radii for light ice (rpp,rsp,rap), cloud (rcp,rc2p),
ccc and rain (rrp) are computed as a mass-weighted average of the 
ccc individual constituents and stored in gce_reff...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           if (i_rams .eq. 1 .and. iflag .ne. 0) then
c Fill mixing ratio holders
             q1=rcp(i+2,1,k)
             q8=rc2p(i+2,1,k)
             q2=rrp(i+2,1,k)
             q3=rpp(i+2,1,k)
             q4=rsp(i+2,1,k)
             q5=rap(i+2,1,k)
             q6=rgp(i+2,1,k)
             q7=rhp(i+2,1,k)
c Zero optical depths
             tau_rams(1) = 0.
             tau_rams(2) = 0.
             tau_rams(3) = 0.
             tau_rams(4) = 0.
             tau_rams(5) = 0.
             tau_rams(6) = 0.
             tau_rams(7) = 0.
             tau_rams(8) = 0.
c Test for cloud vs. no-cloud in liquid and ice 
             if (q1+q8+q2       .ge. twcz) fcld(i,1,km) = 1.0
             if (q3+q4+q5+q6+q7 .ge. twcz) fcld(i,1,km) = 1.0
c Small cloud optical depth
             if (q1 .ge. twcz) then
c    print*,'rams_reff cloud: ',i,k,rams_reff(i+2,1,k,1)
               effrad = rams_reff(i+2,1,k,1) / 1.e4
               if (effrad.gt.0.) 
     1           tau_rams(1) = ( rho(k)*dz0(k)*q1 ) / effrad
               CWC(I,1,KM,2)  = q1
             endif
c Large cloud optical depth
             if (q8 .ge. twcz) then
c              print*,'rams_reff cloud2: ',i,k,rams_reff(i+2,1,k,8)
               effrad = rams_reff(i+2,1,k,8) / 1.e4
               if (effrad.gt.0.)
     1           tau_rams(8) = ( rho(k)*dz0(k)*q8 ) / effrad
                 CWC(I,1,KM,2)  = CWC(I,1,KM,2) + q8
             endif
c Compute cloud effective radius for use in radiation--stored in gce_reff
             if (q1 .ge. twcz .or. q8 .ge. twcz) then
               REFF(I,1,KM,2) = gce_reff(i+2,1,k,2)
c              print*,'gce_reff cloud: ',i,k,reff(i,1,km,2)
             endif
c Rain optical depth
             if (q2 .ge. twcz) then
c              print*,'rams_reff rain: ',i,k,rams_reff(i+2,1,k,2)
               effrad = rams_reff(i+2,1,k,2) / 1.e4
               if (effrad.gt.0.) 
     1           tau_rams(2) = ( rho(k)*dz0(k)*q2 ) / effrad
               CWC(I,1,KM,3)  = q2
               REFF(I,1,KM,3) = gce_reff(i+2,1,k,3)
c              print*,'gce_reff rain: ',i,k,reff(i,1,km,3)
             endif
c Pristine ice optical depth
             if (q3 .ge. twcz) then
c              print*,'rams_reff ice: ',i,k,rams_reff(i+2,1,k,3)
c              effrad = rams_reff(i+2,1,k,3) / 1.e4 ! in cgs
c              if (effrad.gt.0.)
c     1          tau_rams(3) = ( rho(k)*dz0(k)*q3 ) / effrad
               CWC(I,1,KM,1)  = CWC(I,1,KM,1) + q3
ccc D.Posselt 8/24/2006 Try Fu and Liou method for ice tau and reff
               effrad=0.0125+(taq(i,1,km)-243.16)*0.00050
               if (taq(i,1,km) .gt. 243.16) effrad=0.0125
               if (taq(i,1,km) .lt. 223.16) effrad=0.0025
               tauqis=1.e4*rho(k)*dz0(k)*q3*
     +                    (-0.006656+ 3.686e-4/effrad)
               tauqii=1.e4*rho(k)*dz0(k)*q3*
     +                    (-0.011500+ 4.110e-4/effrad
     +                        +17.300e-8/(effrad*effrad))
               reff(i,1,km,1) = effrad*1.e4
c               print*,'gce_reff ice: ',i,k,reff(i,1,km,1)
             endif
c Snow optical depth
             if (q4 .ge. twcz) then
c              print*,'rams_reff snow: ',i,k,rams_reff(i+2,1,k,4)
               effrad = rams_reff(i+2,1,k,4) / 1.e4
               if (effrad.gt.0.) 
     1           tau_rams(4) = ( rho(k)*dz0(k)*q4 ) / effrad
               CWC(I,1,KM,1)  = CWC(I,1,KM,1) + q4
             endif
c Aggregate optical depth
             if (q5 .ge. twcz) then
c              print*,'rams_reff aggregates: ',i,k,rams_reff(i+2,1,k,5)
               effrad = rams_reff(i+2,1,k,5) / 1.e4
               if (effrad.gt.0.) 
     1           tauqs = ( rho(k)*dz0(k)*q5 ) / effrad
               CWC(I,1,KM,1)  = CWC(I,1,KM,1) + q5
             endif
c Graupel optical depth
             if (q6 .ge. twcz) then
c              print*,'rams_reff graupel: ',i,k,rams_reff(i+2,1,k,6)
               effrad = rams_reff(i+2,1,k,6) / 1.e4
               if (effrad.gt.0.) 
     1           tau_rams(6) = ( rho(k)*dz0(k)*q6 ) / effrad
               CWC(I,1,KM,1)  = CWC(I,1,KM,1) + q6
             endif
c Hail optical depth
             if (q7 .ge. twcz) then
c              print*,'rams_reff hail: ',i,k,rams_reff(i+2,1,k,7)
               effrad = rams_reff(i+2,1,k,7) / 1.e4
               if (effrad.gt.0.) 
     1           tau_rams(7) = ( rho(k)*dz0(k)*q7 ) / effrad
               CWC(I,1,KM,1)  = CWC(I,1,KM,1) + q7
             endif
c            if (tau_rams(1).gt.0.) 
c     1         print*,'i,k,tau rams cloud:  ',i,k,tau_rams(1)
c            if (tau_rams(8).gt.0.) 
c     1         print*,'i,k,tau rams cloud2: ',i,k,tau_rams(8)
c            if (tau_rams(2).gt.0.) 
c     1         print*,'i,k,tau rams rain:   ',i,k,tau_rams(2)
c            if (tau_rams(3).gt.0.) 
c     1         print*,'i,k,tau rams ice:    ',i,k,tau_rams(3)
c            if (tau_rams(4).gt.0.) 
c     1         print*,'i,k,tau rams snow:   ',i,k,tau_rams(4)
c            if (tau_rams(5).gt.0.) 
c     1         print*,'i,k,tau rams agg:    ',i,k,tau_rams(5)
c            if (tau_rams(6).gt.0.) 
c     1         print*,'i,k,tau rams grau:   ',i,k,tau_rams(6)
c            if (tau_rams(7).gt.0.) 
c     1         print*,'i,k,tau rams hail:   ',i,k,tau_rams(7)
c Fill optical depth arrays. IR optical depth is assumed to be ~ 0.5 SW
c Small cloud droplets
             tausw(i,1,km,2)=1.5*(tau_rams(1)+tau_rams(8))
             tauir(i,1,km,2)=0.5*tausw(i,1,km,2)
c Rain
             tausw(i,1,km,3)=1.5*tau_rams(2)
             tauir(i,1,km,3)=0.5*tausw(i,1,km,3)
c Ice--note here that there is no differentiation between LW and SW
c             tausw(i,1,km,1)=tau_rams(3)+tau_rams(4) + 
c     1                       tau_rams(5)+tau_rams(6)+tau_rams(7)
c             tauir(i,1,km,1)=tau_rams(3)+tau_rams(4) + 
c     1                       tau_rams(5)+tau_rams(6)+tau_rams(7)
ccc D.Posselt 8/24/2006 Try Fu and Liou method for ice tau and reff
             tausw(i,1,km,1)=tauqis+tau_rams(4) + 
     1                       tau_rams(5)+tau_rams(6)+tau_rams(7)
             tauir(i,1,km,1)=tauqii+tau_rams(4) + 
     1                       tau_rams(5)+tau_rams(6)+tau_rams(7)
c            if (tauir(i,1,km,2).gt.0.) 
c     1         print*,'i,k,tau cloud:       ',i,km,tauir(i,1,km,2)
c            if (tauir(i,1,km,3).gt.0.) 
c     1         print*,'i,k,tau rain:        ',i,km,tauir(i,1,km,3)
c            if (tauir(i,1,km,1).gt.0.) 
c     1         print*,'i,k,tau ice:         ',i,km,tauir(i,1,km,1)
c            if (tausw(i,1,km,2).gt.0.) 
c     1         print*,'i,k,tau cloud:       ',i,km,tausw(i,1,km,2)
c            if (tausw(i,1,km,3).gt.0.) 
c     1         print*,'i,k,tau rain:        ',i,km,tausw(i,1,km,3)
c            if (tausw(i,1,km,1).gt.0.) 
c     1         print*,'i,k,tau ice:         ',i,km,tausw(i,1,km,1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc D.Posselt End RAMS microphysics section
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           else
            b0=0.
            b1=0.
            b2=0.
            b4=0.
            b3=0.
            tauqc=0.0
            tauqr=0.0
            tauqs=0.0
            tauqis=0.0
            tauqii=0.0
            tauqg=0.0
           q1=qc(i,k)
           q2=qr(i,k)
           q5=qg(i,k)
           if (iopcloud .eq. 1) then
              q3=qi(i,k)+qs(i,k)
              q4=0.
           else
              q3=qi(i,k)
              q4=qs(i,k)
           endif
c
          if(q1+q2 .ge. twcz) fcld(i,1,km)=1.0
          if(q3+q4+q5 .ge. twcz) fcld(i,1,km)=1.0
c
          if(q1 .ge. twcz) then
            b0=rho(k)*dz0(k)*q1
           effrad=0.0015
           tauqc=b0/effrad
            b0=1.e4*b0
c
            cwc(i,1,km,2)=q1
            reff(i,1,km,2)=effrad*1.e4
c
          endif
          if(q2 .ge. twcz) then
            b1=rho(k)*dz0(k)*q2
           effrad=3./((cpi*tnw*roqr/(rho(k)*q2))**.25)
           tauqr=b1/effrad
            b1=1.e4*b1
c
            cwc(i,1,km,3)=q2
            reff(i,1,km,3)=effrad*1.e4
c
          endif
c
cc         for ice particles
c
          if(q3 .ge. twcz) then
            b2=1.e4*rho(k)*dz0(k)*q3
           effrad=0.0125+(taq(i,1,km)-243.16)*0.00050
           if (taq(i,1,km) .gt. 243.16) effrad=0.0125
           if (taq(i,1,km) .lt. 223.16) effrad=0.0025
           tauqis=b2*(-0.006656+ 3.686e-4/effrad)
           tauqii=b2*(-0.011500+ 4.110e-4/effrad
     +                         +17.300e-8/(effrad*effrad))
c
            cwc(i,1,km,1)=q3
            reff(i,1,km,1)=effrad*1.e4
c
          endif
c
          if(q4 .ge. twcz) then
            b3=rho(k)*dz0(k)*q4
           effrad=3./((cpi*tns*roqs/(rho(k)*q4))**.25)
           tauqs=b3/effrad
            b3=1.e4*b3
c
            cwc(i,1,km,1)=cwc(i,1,km,1)+q4
c            reff(i,1,km,1)=reff(i,1,km,1)+effrad*1.e4
c
          endif
c
          if(q5 .ge. twcz) then
            b4=rho(k)*dz0(k)*q5
           effrad=3./((cpi*tng*roqg/(rho(k)*q5))**.25)
           tauqg=b4/effrad
            b4=1.e4*b4
            cwc(i,1,km,1)=cwc(i,1,km,1)+q5
c            reff(i,1,km,1)=reff(i,1,km,1)+effrad*1.e4
          endif
c
c            b4=b4+b3
c            tauqg=tauqg
c
          b5=b5+b0
          b6=b6+b1
          b7=b7+b2
          b9=b9+b4
c
           tausw(i,1,km,2)=1.5*tauqc
           tauir(i,1,km,2)=0.5*tausw(i,1,km,2)
c
cccccccccc           tausw(i,1,km,3)=1.5*tauqr
c
           tausw(i,1,km,3)=1.5*tauqr
           tauir(i,1,km,3)=0.5*tausw(i,1,km,3)
c
           tausw(i,1,km,1)=tauqis+tauqg+tauqs
           tauir(i,1,km,1)=tauqii+tauqg+tauqs
ccc D.Posselt 07/24/07 End of if-test for RAMS microphysics
	     endif
  100    continue
  200    continue
         do i=1,nx1
           k=2
           km=kmax-k+nadd
           k=3
           km1=kmax-k+nadd
           aa1(i,km)=0.8*tausw(i,1,km,1)+0.2*tausw(i,1,km1,1)
           aa2(i,km)=0.8*tausw(i,1,km,2)+0.2*tausw(i,1,km1,2)
           bb(i,km)=0.8*tausw(i,1,km,3)+0.2*tausw(i,1,km1,3)
           do 110 k=3,kles-1
             km=kmax-k+nadd
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            aa1(i,km)=0.2*tausw(i,1,km-1,1)+0.60*tausw(i,1,km,1)
     1                +0.2*tausw(i,1,km+1,1)
             aa2(i,km)=0.2*tausw(i,1,km-1,2)+0.60*tausw(i,1,km,2)
     1                +0.2*tausw(i,1,km+1,2)
             bb(i,km)=0.2*tausw(i,1,km-1,3)+0.60*tausw(i,1,km,3)
     1                +0.2*tausw(i,1,km+1,3)
  110      continue
           k=kles
           km=kmax-k+nadd
           k=kl2
           km1=kmax-k+nadd
            aa1(i,km)=0.8*tausw(i,1,km,1)+0.2*tausw(i,1,km1,1)
            aa2(i,km)=0.8*tausw(i,1,km,2)+0.2*tausw(i,1,km1,2)
            bb(i,km)=0.8*tausw(i,1,km,3)+0.2*tausw(i,1,km1,3)
           do 120 k=2,kles
              km=kmax-k+nadd
             tausw(i,1,km,1)=aa1(i,km)
             tausw(i,1,km,2)=aa2(i,km)
             tausw(i,1,km,3)=bb(i,km)
  120      continue
        enddo
         do i=1,nx1
           k=2
           km=kmax-k+nadd
           k=3
           km1=kmax-k+nadd
           aa1(i,km)=0.8*tauir(i,1,km,1)+0.2*tauir(i,1,km1,1)
           aa2(i,km)=0.8*tauir(i,1,km,2)+0.2*tauir(i,1,km1,2)
           bb(i,km)=0.8*tauir(i,1,km,3)+0.2*tauir(i,1,km1,3)
           do 210 k=3,kles-1
             km=kmax-k+nadd
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            aa1(i,km)=0.2*tauir(i,1,km-1,1)+0.60*tauir(i,1,km,1)
     1                +0.2*tauir(i,1,km+1,1)
             aa2(i,km)=0.2*tauir(i,1,km-1,2)+0.60*tauir(i,1,km,2)
     1                +0.2*tauir(i,1,km+1,2)
             bb(i,km)=0.2*tauir(i,1,km-1,3)+0.60*tauir(i,1,km,3)
     1                +0.2*tauir(i,1,km+1,3)
  210      continue
           k=kles
           km=kmax-k+nadd
           k=kl2
           km1=kmax-k+nadd
            aa1(i,km)=0.8*tauir(i,1,km,1)+0.2*tauir(i,1,km1,1)
            aa2(i,km)=0.8*tauir(i,1,km,2)+0.2*tauir(i,1,km1,2)
            bb(i,km)=0.8*tauir(i,1,km,3)+0.2*tauir(i,1,km1,3)
           do 220 k=2,kles
              km=kmax-k+nadd
             tauir(i,1,km,1)=aa1(i,km)
             tauir(i,1,km,2)=aa2(i,km)
             tauir(i,1,km,3)=bb(i,km)
  220      continue
        enddo
        do i=1,nx1
           do k=1,nw
             if (tauir(i,1,k,1) .le. 0.01) tauir(i,1,k,1)=0.0
             if (tauir(i,1,k,2) .le. 0.01) tauir(i,1,k,2)=0.0
             if (tauir(i,1,k,3) .le. 0.01) tauir(i,1,k,3)=0.0
             if (tausw(i,1,k,1) .le. 0.01) tausw(i,1,k,1)=0.0
             if (tausw(i,1,k,2) .le. 0.01) tausw(i,1,k,2)=0.0
             if (tausw(i,1,k,3) .le. 0.01) tausw(i,1,k,3)=0.0
             if (reff(i,1,k,1) .le. 1.e-10) reff(i,1,k,1)=0.0
             if (reff(i,1,k,2) .le. 1.e-10) reff(i,1,k,2)=0.0
             if (reff(i,1,k,3) .le. 1.e-10) reff(i,1,k,3)=0.0
             if (cwc(i,1,k,1) .le. 1.e-6) cwc(i,1,k,1)=0.0
             if (cwc(i,1,k,2) .le. 1.e-6) cwc(i,1,k,2)=0.0
             if (cwc(i,1,k,3) .le. 1.e-6) cwc(i,1,k,3)=0.0
           enddo
        enddo
c
cc
c
      if (iradave .eq. 1) then
         rnx1=1./float(nx1)
        do k=1,nw
          y1(k)=0.
          y2(k)=0.
          y3(k)=0.
          y4(k)=0.
          y5(k)=0.
          y6(k)=0.
          y7(k)=0.
          y8(k)=0.
          y9(k)=0.
          y10(k)=0.
          y11(k)=0.
        enddo
        do 300 k=2,kles
          km=kmax-k+nadd
          do i=1,nx1
             y1(km)=y1(km)+taq(i,1,km)
             y2(km)=y2(km)+waq(i,1,km)
             y3(km)=y3(km)+tausw(i,1,km,1)
             y4(km)=y4(km)+tausw(i,1,km,2)
             y5(km)=y5(km)+tausw(i,1,km,3)
             y6(km)=y6(km)+reff(i,1,km,1)
             y7(km)=y7(km)+reff(i,1,km,2)
             y8(km)=y8(km)+reff(i,1,km,3)
             y9(km)=y9(km)+tauir(i,1,km,1)
             y10(km)=y10(km)+tauir(i,1,km,2)
             y11(km)=y11(km)+tauir(i,1,km,3)
          enddo
  300   continue
         do k=2,kles
          km=kmax-k+nadd
          taq(1,1,km)=y1(km)*rnx1
          waq(1,1,km)=y2(km)*rnx1
          tausw(1,1,km,1)=y3(km)*rnx1
          tausw(1,1,km,2)=y4(km)*rnx1
          tausw(1,1,km,3)=y5(km)*rnx1
          reff(1,1,km,1)=y6(km)*rnx1
          reff(1,1,km,2)=y7(km)*rnx1
          reff(1,1,km,3)=y8(km)*rnx1
          tauir(1,1,km,1)=y9(km)*rnx1          
          tauir(1,1,km,2)=y10(km)*rnx1
          tauir(1,1,km,3)=y11(km)*rnx1
        enddo
      endif
      return
      end
