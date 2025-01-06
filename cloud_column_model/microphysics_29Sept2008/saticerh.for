c
cc D. Posselt 10/01/2008
cc Modified this to be called as a column model from an external driver
cc Stripped out unneccessary common blocks and statistics
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------
c
c D. Posselt 11/24/2010
c Modified the way mixed phase collisions are done. These are controlled
c by dlt2, dlt3, and dlt4, which previously were set to 0.0 or 1.0 
c depending on the air temperature and threshold values of cloud, rain, 
c or snow. I am trying a new technique, which scales these from 0.0 at
c the freezing level linearly back to 1.0 at -10 deg. C
c
c Turning on or off this behavior is controlled via the parameter
c
c imixed_dp = 1 (on for all) = 2 (on for dlt2 only)
c
c For now, this is a locally set variable--defined at compile time
c
c-----------------------------------------------------------------------

      subroutine saticerh (dt,dpt,dqv,qcl,qrn,qci,qcs,qcg,ww1,
     1                     rho,rrho,ta1,qa1,pi,p0,fv,improve,
     2                     tair_out, tairc_out, theta_out, 
     3                     delta2, delta3, delta4)
      include 'dimensions.h'
c     (r&h)  compute ice phase microphysics and saturation processes
c     parameter (nx=66,ny=10,nz=34,nm=nx,nt=2880)

c D.Posselt input variable declarations
      real dt, d2t, d22t
      real dpt(nx,ny,nz), dqv(nx,ny,nz)
      real qcl(nx,ny,nz), qrn(nx,ny,nz)
      real qci(nx,ny,nz), qcs(nx,ny,nz), qcg(nx,ny,nz)
      real rho(nz), rrho(nz)
      real ta1(nz), qa1(nz)
      real p0(nz), pi(nz), fv(nz)
      real ww1(nx,ny,nz)

      real theta_out(nx,ny,nz), tair_out(nx,ny,nz), tairc_out(nx,ny,nz)
      real delta2(nx,ny,nz), delta3(nx,ny,nz), delta4(nx,ny,nz)

c D.Posselt variable declarations
      integer kles
      integer ijkadv, id
      integer new_ice_sat
      real    rijl2

!DP New mixed phase parameterization. 
!DP 0 = original, 1 = scale all dlt by temperature, 2 = scale dlt2 only...
      integer imixed_dp
      parameter (imixed_dp = 0)
!DP Deg C below zero to end scaling (e.g., 10: scaling done between 273 and 263 K)
      real    mixed_scale
      parameter (mixed_scale = 10.)

      parameter (nm=nx,nz2=2*nz,nz3=3*nz,nz4=4*nz,nm2=2*nm)
!     parameter (nxy=nx*ny,nb=nx*ny*(nz-29),nb2=nx*ny*(nz-25),nt2=2*nt)
      parameter (nxy=nx*ny,nb=nx*ny*(nz-29),nb2=nx*ny*(nz-25))
      parameter (nb3=nz+nx,nb5=nx*ny*(nz-29))
!       parameter (nb3=nz+nx,nb4=5*nt,nb5=nx*ny*(nz-29))
!       common/mpi_parameter/imax,iles1,iles,il2,jmax,jles1,jles,jl2,
!      1    kmax,kles,kl2
!       common/bxyz/ n,isec,nran,kt1,kt2
!       common/option/ lipps,ijkadv,istatmin,iwater,itoga,imlifting,lin,
!      1   irf,iadvh,irfg,ismg,id
!       common/timestat/ ndt_stat,idq
!       common/iice/ new_ice_sat
      common/icemass/ ami50,ami40
!       common/bt/ dt,d2t,rijl2,dts,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,
!      1    psfc,fcor,sec,aminut,rdt
      common/cont/ c38,c358,c610,c149,c879,c172,c409,c76,c218,c580,c141
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/size/ tnw,tns,tng,roqs,roqg,roqr
      common/rterv/ zrc,zgc,zsc,vr0,vr1,vr2,vr3,vgc,vsc
      common/rsnw/ alv,alf,als,t0,t00,avc,afc,asc,rn1,rn2,bnd2,rn3,rn4,
     1  rn5,rn50,rn51,rn52,rn53,rn6,rn60,rn61,rn62,rn63,rn7,rn8,rn9,
     2  rn10,rn101,rn102,rn10a,rn10b,rn10c,rn11,rn12,rn12a(31),
     3  rn12b(31),rn13(31),rn14,rn15,rn15a,rn16,rn171,rn172,rn17a,rn17b,
     4  rn17c,rn18,rn18a,rn19,rn191,rn192,rn19a,rn20,rn20a,rn20b,rn30,
     5  rn30a,rn21,bnd21,rn22,rn23,rn231,rn232,rn25,rn25a(31),rn31,beta,
     6  rn32,rn33,rn331,rn332,rn34,rn35,bnd1

      common /BergCon/BergCon1(31),BergCon2(31)
     1               ,BergCon3(31),BergCon4(31)

c
!       common/b1tq/ dpt(nx,ny,nz),dqv(nx,ny,nz)
!       common/b1cr/ qcl(nx,ny,nz),qrn(nx,ny,nz)
!       common/b1ig/ qci(nx,ny,nz),qcg(nx,ny,nz)
!       common/b2tq/ dpt1(nx,ny,nz),dqv1(nx,ny,nz)
!       common/b2cr/ qcl1(nx,ny,nz),qrn1(nx,ny,nz)
!       common/b2ig/ qci1(nx,ny,nz),qcg1(nx,ny,nz)
!       common/b1s/ qcs(nx,ny,nz)
!       common/b2s/ qcs1(nx,ny,nz)
!       common/b4wp/ ww1(nx,ny,nz)
!       common/slwave/ rsw(nx,ny,nz),rlw(nx,ny,nz)
c
c      common/bw/ trahx(nx,ny,nz),trahy(nx,ny,nz),trav(nx,ny,nz)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/bsat1/ pt(nx,ny),qv(nx,ny),qc(nx,ny),qr(nx,ny),qi(nx,ny),
     1 qs(nx,ny),qg(nx,ny),tair(nx,ny),tairc(nx,ny),rtair(nx,ny),
     2 dep(nx,ny),dd(nx,ny),dd1(nx,ny),qvs(nx,ny),dm(nx,ny),
     3 rq(nx,ny),rsub1(nx,ny),col(nx,ny),cnd(nx,ny),ern(nx,ny),
     4 dlt1(nx,ny),dlt2(nx,ny),dlt3(nx,ny),dlt4(nx,ny),zr(nx,ny),
     5 vr(nx,ny),zs(nx,ny),vs(nx,ny),dbz(nx,ny),dda(nb)
      common/badv/ vg(nx,ny),zg(nx,ny),ps(nx,ny),pg(nx,ny),prn(nx,ny),
     1 psn(nx,ny),pwacs(nx,ny),wgacr(nx,ny),pidep(nx,ny),pint(nx,ny),
     2 qsi(nx,ny),ssi(nx,ny),esi(nx,ny),esw(nx,ny),qsw(nx,ny),
     3 pr(nx,ny),ssw(nx,ny),pihom(nx,ny),pidw(nx,ny),pimlt(nx,ny),
     4 psaut(nx,ny),qracs(nx,ny),psaci(nx,ny),psacw(nx,ny),
     5 qsacw(nx,ny),praci(nx,ny),pmlts(nx,ny),pmltg(nx,ny),
     6 asss(nx,ny),dde(nb5)
      common/bsat/ praut(nx,ny),pracw(nx,ny),psfw(nx,ny),psfi(nx,ny),
     1 dgacs(nx,ny),dgacw(nx,ny),dgaci(nx,ny),dgacr(nx,ny),
     2 pgacs(nx,ny),wgacs(nx,ny),qgacw(nx,ny),wgaci(nx,ny),
     3 qgacr(nx,ny),pgwet(nx,ny),pgaut(nx,ny),pracs(nx,ny),
     4 psacr(nx,ny),qsacr(nx,ny),pgfr(nx,ny),psmlt(nx,ny),
     5 pgmlt(nx,ny),psdep(nx,ny),pgdep(nx,ny),piacr(nx,ny),y5(nx,ny),
     6 ddb(nb2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ta(nz),qa(nz),ta1(nz),
!      1  qa1(nz),zx(nz4),am(nz),zq(nz3),wb(nz),zw(nz2),rrho(nz),wbx(nb3)
!       common/b6/ fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),st(nz),sv(nz),
!      1   sq(nz),sc(nz),se(nz),sqa(nz)
      common/b66b/ s_dep(nz),s_sub(nz),s_qrs(nz),s_qrl(nz),s_mel(nz),
     1   s_frz(nz)
      common/brh1/ srro(nz),qrro(nz),sqc(nz),sqr(nz),sqi(nz),sqs(nz),
     1   sqg(nz),stqc(nz),stqr(nz),stqi(nz),stqs(nz),stqg(nz)
!      1   sqg(nz),stqc(nz),stqr(nz),stqi(nz),stqs(nz),stqg(nz),tttd(nb4)
      common/bsts/ thom(nz,4,7),tdw(nz,4,7),tmlt(nz,4,7),saut(nz,4,7),
     1 saci(nz,4,7),sacw(nz,4,7),raci(nz,4,7),tacr(nz,4,7),raut(nz,4,7),
     2 racw(nz,4,7),sfw(nz,4,7),sfi(nz,4,7),gacs(nz,4,7),gacw(nz,4,7),
     3 gaci(nz,4,7),gacr(nz,4,7),gwet(nz,4,7),gaut(nz,4,7),racs(nz,4,7),
     4 sacr(nz,4,7),gfr(nz,4,7),smlt(nz,4,7),gmlt(nz,4,7),sdep(nz,4,7),
     5 ssub(nz,4,7),gsub(nz,4,7),pern(nz,4,7),d3ri(nz,4,7),d3ir(nz,4,7),
     6 d2sr(nz,4,7),d2rs(nz,4,7),gdry(nz,4,7),coc(nz,4,7),coe(nz,4,7),
     7 smf0(nz,4,7),qc0(nz,4,7),qr0(nz,4,7),qi0(nz,4,7),qs0(nz,4,7),
     8 qg0(nz,4,7),sqc0(nz,4,7),sqr0(nz,4,7),sqi0(nz,4,7),sqs0(nz,4,7),
     9 sqg0(nz,4,7),erns(nz,4,7),wgrs(nz,4,7),qsws(nz,4,7),tb00(nz,4),
     1 qb00(nz,4)
      common/bsts1/ tut1(nz,4,7),tut2(nz,4,7),tvt1(nz,4,7),tvt2(nz,4,7),
     1 tstf(nz,4,7),tstf1(nz,4,7),tstf2(nz,4,7),tsqf(nz,4,7),
     2 qqq(nz,4,7),tsqf1(nz,4,7),tsqf2(nz,4,7),tsqq(nz,4,7),
     3 tsqq1(nz,4,7)
      common/bsts3/ qv0(nz,4,7),tt0(nz,4,7),sqv0(nz,4,7),stt0(nz,4,7),
     1 sgpt(nz,4,7),ssgpt(nz,4,7),snqhd(nz,4,7),snqvd(nz,4,7),
     2 q1t(nz,4,7),snhdh(nz,4,7),sqhdt(nz,4,7),sqvdt(nz,4,7)
      common/bsts4/ srsw(nz,4,7),srlw(nz,4,7),sqtdt(nz,4,7),sqhl(nz,4,7)
      common/bsts7/ pim(nz,4,7),cfr(nz,4,7)

      common/bsts40/ fcld(nz,4,7)

      common/bcs/cf1(61,nz),cf2(61,nz),cf3(61,nz),cf4(61,nz),cf5(61,nz),
     1   cf6(61,nz),cf7(61,nz),cf8(61,nz),cf9(61,nz),cf10(61,nz),
     2   cf11(61,nz),cf12(61,nz),cfnum(61,nz),cfs1(61,nz),cfs2(61,nz),
     3   cfs3(61,nz),cfs4(61,nz),cfs5(61,nz),cfs6(61,nz),cfs7(61,nz),
     4   cfs8(61,nz),cfs9(61,nz),cfs10(61,nz),cfs11(61,nz),cfs12(61,nz),
     5   cfsnum(61,nz),cfz(41,nz),cfw(61,nz),scu1(nz),sed1(nz)

!       common/bsts5/ aco5,aco15,aan5,aan15,anv(nt2),cnv(nt2),spn(nt2)
      common/bsts6/ lconv5,lanvl5,lnspt5,ivv(nx,ny),ics5(nx,ny,4),
     1  ibz(nx,ny,4)
      common/bi/ it(nx,ny),ics(nx,ny,4)
      common/rstat/ csttt(nx,ny),cstt(nx,ny)
      common/bls/ y0(nx,ny),ts0new(nx,ny),qss0new(nx,ny)
      common /micro/physc(nx,ny,nz),physe(nx,ny,nz),physd(nx,ny,nz),
     1              physs(nx,ny,nz),physm(nx,ny,nz),physf(nx,ny,nz)

      dimension y1(nx,ny),y2(nx,ny),y3(nx,ny),y4(nx,ny),
     1                b1(nz),b2(nz)
c    1  y7(nm),b0(nm),b1(nm),b2(nm),y6(nm)
c hmhj modified above

      dimension dda0(nx,ny), ddb0(nx,ny) ! cccshie 7/1/02

c    q budget for each point

c_tao
c     common/q_bugt/ q1_g_h(nx,ny,nz),q1_g_v(nx,ny,nz),
c    1               q1a_g_h(nx,ny,nz),q1a_g_v(nx,ny,nz),
c    2               q1_d_h(nx,ny,nz),q1_d_v(nx,ny,nz),
c    3               q1a_d_h(nx,ny,nz),q1a_d_v(nx,ny,nz),
c    4               q2_g_h(nx,ny,nz),q2_g_v(nx,ny,nz),
c    5               q2a_g_h(nx,ny,nz),q2a_g_v(nx,ny,nz),
c    6               q2_d_h(nx,ny,nz),q2_d_v(nx,ny,nz),
c    7               q2a_d_h(nx,ny,nz),q2a_d_v(nx,ny,nz),
c    8               q1_hyd(nxy,nz),q2_hyd(nxy,nz),
c    9               q1a_hyd(nxy,nz),q2a_hyd(nxy,nz),
c    9               q1_rad(nx,ny,nz),q1a_rad(nx,ny,nz),
c    9               ibudsec,rbud

c_tao

      integer itaobraun
      real tairccri, cn0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       dimension fv(1),rby(7),aa1(31),aa2(31)
      dimension rby(7),aa1(31),aa2(31)
      DIMENSION tairN(nx,ny),tairI(nx,ny),tns_funT(nx,ny)
      DIMENSION pihms(nx,ny),pihmg(nx,ny),pssub(nx,ny),pgsub(nx,ny)
      DIMENSION pimm(nx,ny),pcfr(nx,ny)
c     real tttbud(nxy),qqqbud(nxy)
      data aa1/.7939e-7,.7841e-6,.3369e-5,.4336e-5,.5285e-5,.3728e-5,
     1   .1852e-5,.2991e-6,.4248e-6,.7434e-6,.1812e-5,.4394e-5,.9145e-5,
     2   .1725e-4,.3348e-4,.1725e-4,.9175e-5,.4412e-5,.2252e-5,.9115e-6,
     3   .4876e-6,.3473e-6,.4758e-6,.6306e-6,.8573e-6,.7868e-6,.7192e-6,
     4   .6513e-6,.5956e-6,.5333e-6,.4834e-6/
      data aa2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,
     1   .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,
     2   .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,
     3   .4382,.4361/
      data rby/2.,1.,0.,0.,0.,-.5,-1./
      real, allocatable :: tmpd(:,:,:),sumd(:)
      save

      allocate(tmpd(nx,ny,55))
      allocate(sumd(55))

c D.Posselt: Set values of parameters for column model
      kles = nz-1
      ndt_stat = dt
      id = 1
      ijkadv = 1
      d22t = dt + dt
      rijl2 = 1.
      do k=1,nz
!         fv(k)   = sqrt(rho(2) * rrho(k))
        srro(k) = 1. / sqrt(rho(k))
        qrro(k) = sqrt ( srro(k) )
      enddo

c. D. Posselt: Check input data
!       print '(a)','Input to saticerh'
!       print '(a)','dpt, ta, dqv, qa, qcl, qrn, qci, qcs, qcg'
!       do k=1,nz
!         print '(i5,9(2x,e10.5))', 
!      1        k,dpt(1,1,k),ta1(k),dqv(1,1,k),qa1(k), 
!      2        qcl(1,1,k),qrn(1,1,k),qci(1,1,k),qcs(1,1,k),qcg(1,1,k)
!       enddo

cc    ***   three classes of ice-phase   *******************************


!       d22t=d2t

      if(ijkadv .eq. 1) then
	    d2t=dt
      else
	    d2t=d2t
      endif
ctao	  
	  rd2t=1./d2t

!       print*,'dt, d2t, rd2t', dt,d2t,rd2t
ctao
c
c      ijles=nx*(ny-1)-1
c      istart=nx+2
c
c      ijles=max(nx*(ny-1), ny*(nx-1))-1
c      istart=min(nx, ny)+2
c
ctao
      do k=1,nz
      do j=1,1
      do i=1,1
        physc(i,j,k)=0.
        physe(i,j,k)=0.
        physd(i,j,k)=0.
        physs(i,j,k)=0.
        physm(i,j,k)=0.
        physf(i,j,k)=0.
      enddo
      enddo
      enddo
ctao

      cmin=1.e-40
      cmin1=1.e-20
      cmin2=1.e-40
c      cmax2=1.e20
      ucor=3071.29/tnw**.75
      ucos=687.97*roqs**.25/tns**.75
      ucog=687.97*roqg**.25/tng**.75
      uwet=4.464**.95

!       print*,'tnw,tns,tng,roqs,roqg,ucor,ucos,ucog,uwet: ',
!      1        tnw,tns,tng,roqs,roqg,ucor,ucos,ucog,uwet
	  
C  HALLET-MOSSOP RIME SPLINTERING parameters
      if(improve.ge.3) ihalmos=1
      xnsplnt=350.     ! peak # splinters per milligram of rime
      xmsplnt=4.4e-8   ! mass of a splinter (from Ferrier 1994)
      hmtemp1=-2.
      hmtemp2=-4.
      hmtemp3=-6.
      hmtemp4=-8.

	  
      ft=dt/d2t
      rft=rijl2*ft
!       a0=.5*istatmin*rijl2
ctao	  
! 	  a200=ndt_stat*rijl2
ctao	
c$doacross local(i,j)
      do j=1,1
c hmhj modified
c     do i=1,ny
      do i=1,1
	   it(i,j)=1
ctao	   
! 	   asss(i,j)=csttt(i,j)*a200
ctao	   
      enddo
      enddo

      rt0=1./(t0-t00)
      bs3=bs+3.
      bg3=bg+3.
      bsh5=2.5+bsh
      bgh5=2.5+bgh
      bs6=6.+bs
      betah=.5*beta
      rdt=1./d2t
      r10t=rn10*d2t
      r11t=rn11*d2t
      r19t=rn19*d2t
      r19at=rn19a*d2t
      r20t=rn20*d2t
      r23t=rn23*d2t
	  r25a=rn25
      r30t=rn30*d2t
      r33t=rn33*d2t

!       print*,'rt0,bs3,bg3,bsh5,bgh5,bs6,betah,rdt ',
!      1        rt0,bs3,bg3,bsh5,bgh5,bs6,betah,rdt
!       print*,'r10t,r11t,r19t,r19at,r20t,r23t,r25a,r30t,r33t ',
!      1        r10t,r11t,r19t,r19at,r20t,r23t,r25a,r30t,r33t
!       print*
!       print*
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc    rh-type ice scheme

      do 1000 k=2,kles
c         if (ijkadv .eq. 1) then
c           tb0=ta(k)
c           qb0=qa(k)
c         else
	 tb0=ta1(k)
	 qb0=qa1(k)
c         endif
c         p00=p0(k)
	 rp0=3.799052e3/p0(k)
	 pi0=pi(k)
	 pir=1./(pi(k))
crh       pr0=1./p0(k)
	 r00=rho(k)
c         r0s=sqrt(rho(k))
	 rr0=rrho(k)
	 rrs=srro(k)
	 rrq=qrro(k)
	 fv0=fv(k)
	 fvs=sqrt(fv(k))
	 cp409=c409*pi0
	 cv409=c409*avc
	 cp580=c580*pi0
	 cs580=c580*asc
crh       alvr=r00*alv
	 afcp=afc*pir
	 avcp=avc*pir
	 ascp=asc*pir
	 zrr=1.e5*zrc*rrq
	 zsr=1.e5*zsc*rrq
	 zgr=1.e5*zgc*rrq
	 vscf=vsc*fv0
	 vgcf=vgc*fv0
cs        vgcf=vgc*rrs
c	 r1r=rn1*rr0
	 r3f=rn3*fv0
	 r4f=rn4*fv0
	 r5f=rn5*fv0
	 r6f=rn6*fv0
	 r7rf=rn7*rr0*fv0
	 r8rf=rn8*rr0*fv0
	 r9rf=rn9*rr0*fv0
	 r101r=rn101*rr0
	 r102rf=rn102*rrs*fvs
	 r12r=rn12*r00
	 r14f=rn14*fv0
	 r15f=rn15*fv0
	 r16rf=rn16*rr0*fv0
	 r18r=rn18*rr0
	 r191r=rn191*rr0
	 r192rf=rn192*rrs*fvs
	 r22f=rn22*fv0
	 r231r=rn231*rr0
	 r232rf=rn232*rrs*fvs
	 r25rt=rn25*rr0*d2t
	 r31r=rn31*rr0
	 r32rt=rn32*d2t*rrs
	 r331r=rn331*rr0
	 r332rf=rn332*rrs*fvs
	 r34f=rn34*fv0
	 scc=0.
	 see=0.

c$doacross local(j,i)
	 do 150 j=1,1
c hmhj
c      do i=1,1
       do i=1,1

c            tttbud(i,j)=dpt(i,j,k)
c            qqqbud(i,j)=dqv(i,j,k)

	    pt(i,j)=dpt(i,j,k)
	    qv(i,j)=dqv(i,j,k)
	    qc(i,j)=qcl(i,j,k)
	    qr(i,j)=qrn(i,j,k)
	    qi(i,j)=qci(i,j,k)
	    qs(i,j)=qcs(i,j,k)
	    qg(i,j)=qcg(i,j,k)
c            if (qv(i,j)+qb0 .le. 0.) qv(i,j)=-qb0
	    if (qc(i,j) .le. cmin) qc(i,j)=0.0
	    if (qr(i,j) .le. cmin) qr(i,j)=0.0
	    if (qi(i,j) .le. cmin) qi(i,j)=0.0
	    if (qs(i,j) .le. cmin) qs(i,j)=0.0
	    if (qg(i,j) .le. cmin) qg(i,j)=0.0
c            xx0(i,j)=pt(i,j)
c            xx00(i,j)=qv(i,j)
	    tair(i,j)=(pt(i,j)+tb0)*pi0
	    tairc(i,j)=tair(i,j)-t0
	    zr(i,j)=zrr
	    vr(i,j)=0.0
	    zs(i,j)=zsr
	    vs(i,j)=0.0
	    zg(i,j)=zgr
	    vg(i,j)=0.0

	    if (qr(i,j) .gt. cmin) then
	       dd(i,j)=r00*qr(i,j)
	       y1(i,j)=sqrt(dd(i,j))
	       y2(i,j)=sqrt(y1(i,j))
	       zr(i,j)=zrc/y2(i,j)
	       vr(i,j)=fv0*(vr0+vr1*y2(i,j)+vr2*y1(i,j)+
     1                    vr3*y1(i,j)*y2(i,j))
	       vr(i,j)=max(vr(i,j), 0.0)
	    endif

	    if (qs(i,j) .gt. cmin) then
	       dd(i,j)=r00*qs(i,j)
	       y1(i,j)=dd(i,j)**.25
               tns_funT(i,j)=1.
               if(improve.ge.3)then
                y3(i,j)=min(0.,max(-50.,tair(i,j)-t0))                   !2007b
                tns_funT(i,j)=exp(-0.060000*y3(i,j)*0.25)                !2007b
               endif
               ZS(I,J)=ZSC/Y1(I,J)* tns_funT(i,j)
               if(improve.ge.3)then
                y3(i,j)=min(0.,max(-50.,tair(i,j)-t0))                   !2007b
                tns_funT(i,j)=exp(0.060000*y3(i,j)*0.25*bs)              !2007b
               endif
               VS(I,J)=MAX(VSCF*DD(I,J)**BSQ * tns_funT(i,j), 0.)
	    endif

	    if (qg(i,j) .gt. cmin) then
	       dd(i,j)=r00*qg(i,j)
	       y1(i,j)=dd(i,j)**.25
	       zg(i,j)=zgc/y1(i,j)
	       vg(i,j)=max(vgcf*dd(i,j)**bgq, 0.0)
	    endif

	    if (qr(i,j) .le. cmin1) vr(i,j)=0.0
	    if (qs(i,j) .le. cmin1) vs(i,j)=0.0
	    if (qg(i,j) .le. cmin1) vg(i,j)=0.0
c
c*  1 * psaut : autoconversion of qi to qs                        ***1**
c*  3 * psaci : accretion of qi to qs                             ***3**
c*  4 * psacw : accretion of qc by qs (riming) (qsacw for psmlt)  ***4**
c* 34 * pwacs : collection of qs by qc                            **34**
c*  5 * praci : accretion of qi by qr                             ***5**
c*  6 * piacr : accretion of qr or qg by qi                       ***6**

       psaut(i,j)=0.0
       psaci(i,j)=0.0
       praci(i,j)=0.0
       piacr(i,j)=0.0
       psacw(i,j)=0.0
       pwacs(i,j)=0.0
       qsacw(i,j)=0.0
       if(tair(i,j).lt.t0) then
c	  y1(i,j)=rdt*(qi(i,j)-r1r*exp(beta*tairc(i,j)))
c	 psaut(i,j)=max(y1(i,j),0.0)
c             rn1=1.e-3
c             bnd1=6.e-4 ! D.Posselt bnd1 and rn1 now set in consatrh
            esi(i,j)=exp(.025*tairc(i,j))
            if(improve.ge.3) esi(i,j)=0.15
           psaut(i,j)=max(rn1*esi(i,j)*(qi(i,j)-bnd1*fv0*fv0) ,0.0)
       endif

       if(tair(i,j).lt.t0) then
        tns_funT(i,j)=1.
       if(improve.ge.3) tns_funT(i,j)=min(20.,exp(-0.060000*tairc(i,j)))  !Lang et al. 2007b
        esi(i,j)=1.
        if(improve.ge.3) esi(i,j)=0.15
	psaci(i,j)=r3f*qi(i,j)/zs(i,j)**bs3 *tns_funT(i,j) * esi(i,j)    !Lang et al. 2007b
	psacw(i,j)=r4f*qc(i,j)/zs(i,j)**bs3 *tns_funT(i,j)
        if(ihalmos.eq.1)then                                             !Lang et al. 2007b
          y2(i,j)=0.
          if((tairc(i,j).le.hmtemp1).and.(tairc(i,j).ge.hmtemp4))
     1                                              y2(i,j)=0.5
          if((tairc(i,j).ge.hmtemp2).and.(tairc(i,j).le.hmtemp3))
     1                                              y2(i,j)=1.
          pihms(i,j)=psacw(i,j)*y2(i,j)*xnsplnt*1000.*xmsplnt
          psacw(i,j)=psacw(i,j)-pihms(i,j)
        endif
	pwacs(i,j)=r34f*qc(i,j)/zs(i,j)**bs6 *tns_funT(i,j)
	 y1(i,j)=1./zr(i,j)
	 y2(i,j)=y1(i,j)*y1(i,j)
	 y3(i,j)=y1(i,j)*y2(i,j)
	 dd(i,j)=r5f*qi(i,j)*y3(i,j)*(rn50+rn51*y1(i,j)+rn52*y2(i,j)
     1                            +rn53*y3(i,j))
	praci(i,j)=max(dd(i,j),0.0)
	   y4(i,j)=y3(i,j)*y3(i,j)
	 dd1(i,j)=r6f*qi(i,j)*y4(i,j)*(rn60+rn61*y1(i,j)+rn62*y2(i,j)
     1                             +rn63*y3(i,j))
	piacr(i,j)=max(dd1(i,j),0.0)
       else
	qsacw(i,j)=r4f*qc(i,j)/zs(i,j)**bs3 *tns_funT(i,j)
       endif

c DP Turn off snow-to-rain generation processes
c      qsacw(i,j) = 0.0

c* 21 * praut   autoconversion of qc to qr                        **21**
c* 22 * pracw : accretion of qc by qr                             **22**
	 praut(i,j)=max(rn21*(qc(i,j)-bnd21),0.0)
         y1(i,j)=1./zr(i,j)
	 y2(i,j)=y1(i,j)*y1(i,j)
	 y3(i,j)=y1(i,j)*y2(i,j)
	 y4(i,j)=r22f*qc(i,j)*y3(i,j)*(rn50+rn51*y1(i,j)+rn52*y2(i,j)
     1                             +rn53*y3(i,j))
	pracw(i,j)=max(y4(i,j),0.0)
c* 12 * psfw : bergeron processes for qs (koening, 1971)          **12**
c* 13 * psfi : bergeron processes for qs                          **13**

	  psfw(i,j)=0.0
	  psfi(i,j)=0.0
	   if(tair(i,j).lt.t0) then
	     y1(i,j)=max( min(tairc(i,j), -1.), -31.)
	     it(i,j)=int(abs(y1(i,j)))
	     y1(i,j)=rn12a(it(i,j))
	     y2(i,j)=rn12b(it(i,j))
	     y3(i,j)=rn13(it(i,j))
	    psfw(i,j)=max(d2t*y1(i,j)*(y2(i,j)+r12r*qc(i,j))*qi(i,j),0.0)
	    psfi(i,j)=y3(i,j)*qi(i,j)

            if(improve.eq.-1)then                             !xp
              tmp=BergCon1(it(i,j))*qi(i,j)
     1         +BergCon2(it(i,j))*rrho(k)*rn25*exp(beta*(tair(i,j)-t0))
              psfi(i,j)=max(tmp*d2t,0.0)
            endif

            if(improve.ge.2)then                             !Lang et al. 2007b
             y4(i,j)=1./(tair(i,j)-c358)
             y5(i,j)=1./(tair(i,j)-c76)
             qsw(i,j)=rp0*exp(c172-c409*y4(i,j))
             qsi(i,j)=rp0*exp(c218-c580*y5(i,j))
             hfact=(qv(i,j)+qb0-qsi(i,j))/(qsw(i,j)-qsi(i,j))
             if(hfact.gt.1.) hfact=1.
clang        if(qc(i,j).le.cmin) hfact=0.
             sfact=1.
             if(improve.ge.3)then
              SSI(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
              xssi=min(ssi(i,j), 0.20)                            !cap at 20% super wrt ice
              r_nci=min(1.e-3*exp(-.639+12.96*xssi),1.)           !meyers et al. 1992
              dd(i,j)=min((r00*qi(i,j)/r_nci),ami40)              !mean cloud ice mass
              yy1=1.-aa2(it(i,j))
              sfact=(AMI50**YY1-AMI40**YY1)/(AMI50**YY1-dd(i,j)**YY1)
             endif
             if(hfact.gt.0.)then
               psfi(i,j)=psfi(i,j)*hfact*sfact
             else
               psfi(i,j)=0.
             endif
            endif
	   endif
cttt***** qg=qg+min(pgdry,pgwet)
c*  9 * pgacs : accretion of qs by qg (dgacs,wgacs: dry and wet)  ***9**
c* 14 * dgacw : accretion of qc by qg (qgacw for pgmlt)           **14**
c* 15 * dgaci : accretion of qi by qg (wgaci for wet growth)      **15**
c* 16 * dgacr : accretion of qr to qg (qgacr for pgmlt)           **16**
	 y1(i,j)=abs( vg(i,j)-vs(i,j) )
	 y2(i,j)=zs(i,j)*zg(i,j)
	 y3(i,j)=5./y2(i,j)
	 y4(i,j)=.08*y3(i,j)*y3(i,j)
	 y5(i,j)=.05*y3(i,j)*y4(i,j)
	 y2(i,j)=y1(i,j)*(y3(i,j)/zs(i,j)**5+y4(i,j)/zs(i,j)**3+y5(i,j)
     1                  /zs(i,j))

	pgacs(i,j)=r9rf*y2(i,j) *tns_funT(i,j)
	dgacs(i,j)=pgacs(i,j)
	if(improve.ge.1) dgacs(i,j)=0.0
crh     wgacs(i,j)=10.*r9rf*y2(i,j)
	wgacs(i,j)=0.0
	 y1(i,j)=1./zg(i,j)**bg3
	dgacw(i,j)=r14f*qc(i,j)*y1(i,j)
        if(ihalmos.eq.1)then                                    !Lang et al. 2007b
          y2(i,j)=0.
          if((tairc(i,j).le.hmtemp1).and.(tairc(i,j).ge.hmtemp4))
     1                                               y2(i,j)=0.5
          if((tairc(i,j).ge.hmtemp2).and.(tairc(i,j).le.hmtemp3))
     1                                               y2(i,j)=1.
          pihmg(i,j)=dgacw(i,j)*y2(i,j)*xnsplnt*1000.*xmsplnt
          dgacw(i,j)=dgacw(i,j)-pihmg(i,j)
        endif
	qgacw(i,j)=dgacw(i,j)
	dgaci(i,j)=r15f*qi(i,j)*y1(i,j)
	if(improve.ge.1) dgaci(i,j)=0.0
crh     wgaci(i,j)=r15af*qi(i,j)*y1(i,j)
	wgaci(i,j)=0.0

	 y1(i,j)=abs( vg(i,j)-vr(i,j) )
	 y2(i,j)=zr(i,j)*zg(i,j)
	 y3(i,j)=5./y2(i,j)
	 y4(i,j)=.08*y3(i,j)*y3(i,j)
	 y5(i,j)=.05*y3(i,j)*y4(i,j)
	 dd(i,j)=r16rf*y1(i,j)*(y3(i,j)/zr(i,j)**5+y4(i,j)/zr(i,j)**3
     1                     +y5(i,j)/zr(i,j))
	dgacr(i,j)=max(dd(i,j),0.0)
	qgacr(i,j)=dgacr(i,j)

	 if (tair(i,j) .ge. t0) then
	  dgacs(i,j)=0.0
crh       wgacs(i,j)=0.0
	  dgacw(i,j)=0.0
	  dgaci(i,j)=0.0
crh       wgaci(i,j)=0.0
	  dgacr(i,j)=0.0
	 else
	  pgacs(i,j)=0.0
	  qgacw(i,j)=0.0
	  qgacr(i,j)=0.0
	 endif
       enddo   !!! cccshie added by shie 7/1/02 do i=3,jles
  150  continue

c*******pgdry : dgacw+dgaci+dgacr+dgacs                           ******
c* 17 * pgwet : wet growth of qg                                  **17**
crh       pgwet(ij)=0.0
crh       if (tair(ij) .lt. t0) then
crh         y1(ij)=1./(alf+rn17c*tairc(ij))
crh         y2(ij)=rp0-(qv(ij)+qb0)
crh         y3(ij)=.78/zg(ij)**2+r17arf/zg(ij)**bgh5
crh         y4(ij)=rn171*y2(ij)-r172r*tairc(ij)
crh         dd(ij)=y1(ij)*(y4(ij)*y3(ij)+(wgaci(ij)
crh  1                               +wgacs(ij))*(alf+rn17b*tairc(ij)))
crh        pgwet(ij)=max(dd(ij), 0.0)
crh       endif
c******** shed process (wgacr=pgwet-dgacw-wgaci-wgacs)
crh     wgacr(ij)=pgwet(ij)-dgacw(ij)-wgaci(ij)-wgacs(ij)
crh      y2(ij)=dgacw(ij)+dgaci(ij)+dgacr(ij)+dgacs(ij)
crh      if (pgwet(ij) .ge. y2(ij)) then
crh       wgacr(ij)=0.0
crh       wgaci(ij)=0.0
crh       wgacs(ij)=0.0
crh      else
crh       dgacr(ij)=0.0
crh       dgaci(ij)=0.0
crh       dgacs(ij)=0.0
crh      endif
c*******pgdry : dgacw+dgaci+dgacr+dgacs                           ******
c* 15 * dgaci : accretion of qi by qg (wgaci for wet growth)      **15**
c* 17 * pgwet : wet growth of qg                                  **17**
c********   handling the negative cloud water (qc)    ******************
c********   handling the negative cloud ice (qi)      ******************

c$doacross local(j,i)
      do 200 j=1,1
      do 200 i=1,1
	pgwet(i,j)=0.0
	   y1(i,j)=qc(i,j)/d2t
	  psacw(i,j)=min(y1(i,j), psacw(i,j))
	  pihms(i,j)=min(y1(i,j), pihms(i,j))
	  praut(i,j)=min(y1(i,j), praut(i,j))
	  pracw(i,j)=min(y1(i,j), pracw(i,j))
	  psfw(i,j)= min(y1(i,j), psfw(i,j))
	  dgacw(i,j)=min(y1(i,j), dgacw(i,j))
	  pihmg(i,j)=min(y1(i,j), pihmg(i,j))
	  qsacw(i,j)=min(y1(i,j), qsacw(i,j))
	  qgacw(i,j)=min(y1(i,j), qgacw(i,j))

	y1(i,j)=d2t*(psacw(i,j)+praut(i,j)+pracw(i,j)+psfw(i,j)
     1          +dgacw(i,j)+qsacw(i,j)+qgacw(i,j)+pihms(i,j)+pihms(i,j))

	qc(i,j)=qc(i,j)-y1(i,j)
c
	if (qc(i,j) .lt. 0.0) then
	   y2(i,j)=1.
	    if (y1(i,j) .ne. 0.) y2(i,j)=qc(i,j)/y1(i,j)+1.
	   psacw(i,j)=psacw(i,j)*y2(i,j)
	   praut(i,j)=praut(i,j)*y2(i,j)
	   pracw(i,j)=pracw(i,j)*y2(i,j)
	   psfw(i,j)=psfw(i,j)*y2(i,j)
	   dgacw(i,j)=dgacw(i,j)*y2(i,j)
	   qsacw(i,j)=qsacw(i,j)*y2(i,j)
	   qgacw(i,j)=qgacw(i,j)*y2(i,j)
	   qc(i,j)=0.0
	 endif
c
	    y1(i,j)=qi(i,j)/d2t
	   psaut(i,j)=min(y1(i,j), psaut(i,j))
	   psaci(i,j)=min(y1(i,j), psaci(i,j))
	   praci(i,j)=min(y1(i,j), praci(i,j))
	   psfi(i,j)= min(y1(i,j), psfi(i,j))
	   dgaci(i,j)=min(y1(i,j), dgaci(i,j))
	   wgaci(i,j)=min(y1(i,j), wgaci(i,j))

	y1(i,j)=d2t*(psaut(i,j)+psaci(i,j)+praci(i,j)+psfi(i,j)
     1           +dgaci(i,j)+wgaci(i,j))

	qi(i,j)=qi(i,j)-y1(i,j)
c
	 if (qi(i,j) .lt. 0.0) then
	   y2(i,j)=1.
	    if (y1(i,j) .ne. 0.0) y2(i,j)=qi(i,j)/y1(i,j)+1.
	   psaut(i,j)=psaut(i,j)*y2(i,j)
	   psaci(i,j)=psaci(i,j)*y2(i,j)
	   praci(i,j)=praci(i,j)*y2(i,j)
	   psfi(i,j)=psfi(i,j)*y2(i,j)
	   dgaci(i,j)=dgaci(i,j)*y2(i,j)
	   wgaci(i,j)=wgaci(i,j)*y2(i,j)
	   qi(i,j)=0.0
	 endif
c
	 wgacr(i,j)=qgacr(i,j)+qgacw(i,j)
!DP New mixed phase parameterization -- scale by temperature
      if ( imixed_dp .eq. 1 ) then
! At temperatures greater than/equal to freezing
        if (tair(i,j) .ge. t0) then
          dlt3(i,j) = 0.
          dlt4(i,j) = 0.
! At temperatures less than to freezing - mixed_scale
        else if (tair(i,j) .lt. t0-mixed_scale) then
          dlt3(i,j) = 1.
          dlt4(i,j) = 1.
! At temperatures between freezing and freezing - mixed_scale
        else
          dlt3(i,j) = ( t0-tair(i,j) ) / mixed_scale
          dlt4(i,j) = dlt3(i,j)
        endif

!DP Old parameterization
      else
	dlt3(i,j)=0.0
	 if (qr(i,j) .lt. 1.e-4) dlt3(i,j)=1.
	dlt4(i,j)=1.
c DP Implement Steve Lang's change here--raises the threshold for dlt3/4...
        if (qc(i,j) .gt. 5.e-4) dlt4(i,j)=0.0
!improve if (qc(i,j) .gt. 5.e-4) dlt4(i,j)=0.0
!          if (qc(i,j) .gt. 1.e-3) dlt4(i,j)=0.0
	 if (qs(i,j) .le. 1.e-4) dlt4(i,j)=1.
	  if (tair(i,j) .ge. t0) then
	   dlt3(i,j)=0.0
	   dlt4(i,j)=0.0
	  endif
!DP End if statement for old parameterization
      endif
      delta3(1,1,k) = dlt3(i,j)
      delta4(1,1,k) = dlt4(i,j)
	pr(i,j)=d2t*(qsacw(i,j)+praut(i,j)+pracw(i,j)+wgacr(i,j)
     1           -qgacr(i,j))
	ps(i,j)=d2t*(psaut(i,j)+psaci(i,j)+dlt4(i,j)*psacw(i,j)
     1           +psfw(i,j)+psfi(i,j)+dlt3(i,j)*praci(i,j))
	pg(i,j)=d2t*((1.-dlt3(i,j))*praci(i,j)+dgaci(i,j)
     1          +wgaci(i,j)+dgacw(i,j)+(1.-dlt4(i,j))*psacw(i,j)) !!! cccshie 7
c*  7 * pracs : accretion of qs by qr                             ***7**
c*  8 * psacr : accretion of qr by qs (qsacr for psmlt)           ***8**

	  y1(i,j)=abs( vr(i,j)-vs(i,j) )
	  y2(i,j)=zr(i,j)*zs(i,j)
	  y3(i,j)=5./y2(i,j)
	  y4(i,j)=.08*y3(i,j)*y3(i,j)
	  y5(i,j)=.05*y3(i,j)*y4(i,j)

	pracs(i,j)=r7rf*y1(i,j)*(y3(i,j)/zs(i,j)**5+y4(i,j)/zs(i,j)**3
     1                       +y5(i,j)/zs(i,j)) *tns_funT(i,j)
	qracs(i,j)=min(d2t*pracs(i,j), qs(i,j))
	psacr(i,j)=r8rf*y1(i,j)*(y3(i,j)/zr(i,j)**5+y4(i,j)/zr(i,j)**3
     1                       +y5(i,j)/zr(i,j)) *tns_funT(i,j)
	qsacr(i,j)=psacr(i,j)

c DP Turn off rain generation by snow accretion
c      qracs(i,j) = 0.0

	 if (tair(i,j) .ge. t0) then
	  pgaut(i,j)=0.0
	  pracs(i,j)=0.0
	  psacr(i,j)=0.0
	 else
	  qsacr(i,j)=0.0
	  qracs(i,j)=0.0
	 endif
c*  2 * pgaut : autoconversion of qs to qg                        ***2**
c* 18 * pgfr : freezing of qr to qg                               **18**
	pgfr(i,j)=0.0
	pgaut(i,j)=0.0
	 if (tair(i,j) .lt. t0) then
	   y1(i,j)=exp(rn18a*(t0-tair(i,j)))
	  pgfr(i,j)=max(r18r*(y1(i,j)-1.)/zr(i,j)**7., 0.0)
	 endif
  200 continue

c********   handling the negative rain water (qr)    *******************
c********   handling the negative snow (qs)          *******************

c$doacross local(j,i)
      do 250 j=1,1
      do 250 i=1,1
	  y1(i,j)=qr(i,j)/d2t
	 piacr(i,j)=min(y1(i,j), piacr(i,j))
	 dgacr(i,j)=min(y1(i,j), dgacr(i,j))
	 psacr(i,j)=min(y1(i,j), psacr(i,j))
	 pgfr(i,j)= min(y1(i,j), pgfr(i,j))
	 y1(i,j)=(piacr(i,j)+dgacr(i,j)+psacr(i,j)+pgfr(i,j))*d2t
	qr(i,j)=qr(i,j)+pr(i,j)+qracs(i,j)-y1(i,j)
	if (qr(i,j) .lt. 0.0) then
	  y2(i,j)=1.
	   if (y1(i,j) .ne. 0.0) y2(i,j)=qr(i,j)/y1(i,j)+1.
	  piacr(i,j)=piacr(i,j)*y2(i,j)
	  dgacr(i,j)=dgacr(i,j)*y2(i,j)
	  pgfr(i,j)=pgfr(i,j)*y2(i,j)
	  psacr(i,j)=psacr(i,j)*y2(i,j)
	  qr(i,j)=0.0
	endif
!DP New mixed phase parameterization scales by temperature--see above.
      if ( imixed_dp .ne. 0 ) then
! At temperatures greater than/equal to freezing
        if (tair(i,j) .ge. t0) then
          dlt2(i,j) = 0.
! At temperatures less than freezing - mixed_scale: no mixed phase
        else if (tair(i,j) .lt. t0-mixed_scale) then
          dlt2(i,j) = 1.
! At temperatures between freezing and freezing - mixed_scale: mixed phase
        else
          dlt2(i,j) = ( t0-tair(i,j) ) / mixed_scale
        endif

!DP Old parameterization
      else
!DP If-statement for old delta function settings
	dlt2(i,j)=1.
	 if (qr(i,j) .gt. 1.e-4) dlt2(i,j)=0.
	 if (qs(i,j) .le. 1.e-4) dlt2(i,j)=1.
	 if (tair(i,j) .ge. t0) dlt2(i,j)=0.
      endif
!DP Store delta2 for analysis
      delta2(i,j,k) = dlt2(i,j)
	  y1(i,j)=qs(i,j)/d2t
	 pgacs(i,j)=min(y1(i,j), pgacs(i,j))
	 dgacs(i,j)=min(y1(i,j), dgacs(i,j))
	 wgacs(i,j)=min(y1(i,j), wgacs(i,j))
	 pgaut(i,j)=min(y1(i,j), pgaut(i,j))
	 pracs(i,j)=min(y1(i,j), pracs(i,j))
	 pwacs(i,j)=min(y1(i,j), pwacs(i,j))
	prn(i,j)=d2t*((1.-dlt3(i,j))*piacr(i,j)+dgacr(i,j)+pgfr(i,j)
     1               +(1.-dlt2(i,j))*psacr(i,j))
	ps(i,j)=ps(i,j)+d2t*(dlt3(i,j)*piacr(i,j)+dlt2(i,j)*psacr(i,j))
	 pracs(i,j)=(1.-dlt2(i,j))*pracs(i,j)
	 pwacs(i,j)=(1.-dlt4(i,j))*pwacs(i,j)

      psn(i,j)=d2t*(pgacs(i,j)+dgacs(i,j)+wgacs(i,j)+pgaut(i,j)
     1            +pracs(i,j)+pwacs(i,j))

	qs(i,j)=qs(i,j)+ps(i,j)-qracs(i,j)-psn(i,j)
	 if (qs(i,j) .lt. 0.0) then
	   y2(i,j)=1.
	    if (psn(i,j) .ne. 0.) y2(i,j)=qs(i,j)/psn(i,j)+1.
	   pgacs(i,j)=pgacs(i,j)*y2(i,j)
	   dgacs(i,j)=dgacs(i,j)*y2(i,j)
	   wgacs(i,j)=wgacs(i,j)*y2(i,j)
	   pgaut(i,j)=pgaut(i,j)*y2(i,j)
	   pracs(i,j)=pracs(i,j)*y2(i,j)
	   pwacs(i,j)=pwacs(i,j)*y2(i,j)
	   qs(i,j)=0.0
	 endif
      psn(i,j)=d2t*(pgacs(i,j)+dgacs(i,j)+wgacs(i,j)+pgaut(i,j)
     1             +pracs(i,j)+pwacs(i,j))
       qg(i,j)=qg(i,j)+pg(i,j)+prn(i,j)+psn(i,j)
	y1(i,j)=d2t*(psacw(i,j)+psfw(i,j)+dgacw(i,j)+piacr(i,j)
     1           +dgacr(i,j)+psacr(i,j)+pgfr(i,j))-qracs(i,j)
       pt(i,j)=pt(i,j)+afcp*y1(i,j)
c* 11 * psmlt : melting of qs                                     **11**
c* 19 * pgmlt : melting of qg to qr                               **19**

	  psmlt(i,j)=0.0
	  pgmlt(i,j)=0.0
c DP Turn this off as a test...
c      if ( 1 .eq. 2 ) then
c      print*,'We have serious problems'
	tair(i,j)=(pt(i,j)+tb0)*pi0
	if (tair(i,j) .ge. t0) then
	  tairc(i,j)=tair(i,j)-t0
	  dd(i,j)=r11t*tairc(i,j)*(r101r/zs(i,j)**2+r102rf
     1                         /zs(i,j)**bsh5)
	 psmlt(i,j)=min(qs(i,j),max(dd(i,j),0.0))
	   y2(i,j)=r191r/zg(i,j)**2+r192rf/zg(i,j)**bgh5
	  dd1(i,j)=tairc(i,j)*(r19t*y2(i,j)+r19at*(qgacw(i,j)
     1                                       +qgacr(i,j)))
	 pgmlt(i,j)=min(qg(i,j),max(dd1(i,j),0.0))
	 pt(i,j)=pt(i,j)-afcp*(psmlt(i,j)+pgmlt(i,j))
	 qr(i,j)=qr(i,j)+psmlt(i,j)+pgmlt(i,j)
	 qs(i,j)=qs(i,j)-psmlt(i,j)
	 qg(i,j)=qg(i,j)-pgmlt(i,j)
	endif
c      endif

c* 24 * pihom : homogeneous freezing of qc to qi (t < t00)        **24**
c* 25 * pidw  : deposition growth of qc to qi ( t0 < t <= t00)    **25**
c* 26 * pimlt : melting of qi to qc (t >= t0)                     **26**
c****** PIMM  : IMMERSION FREEZING OF QC TO QI (T < T0)           ******
c****** PCFR  : CONTACT NUCLEATION OF QC TO QI (T < T0)           ******

	if (qc(i,j).le.cmin) qc(i,j)=0.0
	if (qi(i,j).le.cmin) qi(i,j)=0.0
         tair(i,j)=(pt(i,j)+tb0)*pi0
	 if(tair(i,j).le.t00) then
	  pihom(i,j)=qc(i,j)
	 else
	  pihom(i,j)=0.0
	 endif
	 if(tair(i,j).ge.t0) then
	  pimlt(i,j)=qi(i,j)
	 else
	  pimlt(i,j)=0.0
	 endif
	 pidw(i,j)=0.0
	 if (tair(i,j).lt.t0 .and. tair(i,j).gt.t00) then
          if(improve.ge.3)then                                    !Lang et al. 2007b
            if(tairc(i,j).le.-5.)then                             !meyers
             TAIRC(I,J)=TAIR(I,J)-T0
             Y1(I,J)=MAX( MIN(TAIRC(I,J), -1.), -31.)
             IT(I,J)=INT(ABS(Y1(I,J)))
             Y3(I,J)=AA2(IT(I,J))
              y4(i,j)=1./(tair(i,j)-c358)
              qsw(i,j)=rp0*exp(c172-c409*y4(i,j))
              rtair(i,j)=1./(tair(i,j)-c76)
              y2(i,j)=exp(c218-c580*rtair(i,j))
              qsi(i,j)=rp0*y2(i,j)
              SSI(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
              xssi=min(ssi(i,j), 0.20)                            !cap at 20% super wrt ice
              r_nci=min(1.e-3*exp(-.639+12.96*xssi),1.)           !meyers et al. 1992
              if(r_nci.gt.15.) r_nci=15.                          !cap at 15000/liter
              dd(i,j)=(r00*qi(i,j)/r_nci)**y3(i,j)                !meyers
              PIDW(i,j)=min(RR0*D2T*y2(i,j)*r_nci*dd(i,j),qc(i,j))    !meyers
            endif
            pimm(i,j)=0.0
            pcfr(i,j)=0.0
            IF(QC(I,J).gt.0.)THEN
             xncld=qc(i,j)/4.e-9                                   !cloud number
             esat=0.6112*exp(17.67*tairc(i,j)/(tairc(i,j)+243.5))*10.
             rv=0.622*esat/(p0(k)/1000.-esat)
             rlapse_m=980.616*(1.+2.5e6*rv/287./tair(i,j))/
     1           (1004.67+2.5e6*2.5e6*rv*0.622/287./tair(i,j)/tair(i,j))
             delT=rlapse_m*0.5*(ww1(i,j,k)+ww1(i,j,k+1))
             if(delT.lt.0.) delT=0.
             Bhi=1.01e-2                                          !pollen (Deihl et al. 2006)
             pimm(i,j)=xncld*Bhi*4.e-9*exp(-tairc(i,j))*delT*d2t*4.e-9

             CPI=4.*ATAN(1.)
             Rc=1.e-3                                             !cloud droplet radius 10 microns
             Ra=1.e-5                                             !aerosol radius 0.1 microns
             Cna=500.                                             !contact nuclei conc per cc
             xccld=xncld*r00                                      !cloud number concentration
             Xknud=7.37*tair(i,j)/288./p0(k)/Ra                   !Knudsen number
             cunnF=1.257+0.400*exp(-1.10/Xknud)                   !Cunningham correction
             DIFFar=1.3804e-16*tair(i,j)/6./cpi/1.718e-4/Ra*cunnF
             pcfr(i,j)=4.e-9*4.*cpi*Rc*DIFFar*xccld*Cna*rr0*d2t
            ENDIF
          else
	   tairc(i,j)=tair(i,j)-t0
	   y1(i,j)=max( min(tairc(i,j), -1.), -31.)
	   it(i,j)=int(abs(y1(i,j)))
	   y2(i,j)=aa1(it(i,j))
	   y3(i,j)=aa2(it(i,j))
	   y4(i,j)=exp(abs(.5*tairc(i,j)))
	   dd(i,j)=(r00*qi(i,j)/(r25a*y4(i,j)))**y3(i,j)
	   pidw(i,j)=min(r25rt*y2(i,j)*y4(i,j)*dd(i,j),qc(i,j))

            if(improve.eq.-1)then                             !xp
              tmp=BergCon3(it(i,j))*qi(i,j)
     1         +BergCon4(it(i,j))*rrho(k)*rn25*exp(beta*(tair(i,j)-t0))
              pidw(i,j)=min(tmp*d2t,qc(i,j))
            endif

	  endif
	 endif

	 y1(i,j)=pihom(i,j)-pimlt(i,j)+pidw(i,j)+pimm(i,j)+pcfr(i,j)

        if(y1(i,j).gt.qc(i,j))then
           y1(i,j)=qc(i,j)
           y2(i,j)=1.
           y3(i,j)=pihom(i,j)+pidw(i,j)+pimm(i,j)+pcfr(i,j)
           if(y3(i,j).ne.0.) y2(i,j)=(qc(i,j)+pimlt(i,j))/y3(i,j)
           pihom(i,j)=pihom(i,j)*y2(i,j)
           pidw(i,j)=pidw(i,j)*y2(i,j)
           pimm(i,j)=pimm(i,j)*y2(i,j)
           pcfr(i,j)=pcfr(i,j)*y2(i,j)
        endif

	pt(i,j)=pt(i,j)+afcp*y1(i,j)
	qc(i,j)=qc(i,j)-y1(i,j)
	qi(i,j)=qi(i,j)+y1(i,j)
c* 31 * pint  : initiation of qi                                  **31**
c* 32 * pidep : deposition of qi                                  **32**
	pint(i,j)=0.0
	pidep(i,j)=0.0
        if (improve.ge.3) then                                    !Lang et al. 2007b
            tair(i,j)=(pt(i,j)+tb0)*pi0
            if (tair(i,j).lt.t0) THEN
              if (qi(i,j).le.cmin) qi(i,j)=0.
              tairc(i,j)=tair(i,j)-t0
              rtair(i,j)=1./(tair(i,j)-c76)
              y2(i,j)=exp(c218-c580*rtair(i,j))
              qsi(i,j)=rp0*y2(i,j)
              esi(i,j)=C610*y2(i,j)
              SSI(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
              y1(i,j)=1./tair(i,j)
              y3(i,j)=SQRT(qi(i,j))
              dd(i,j)=y1(i,j)*(RN10A*y1(i,j)-RN10B)+RN10C*tair(i,j)/
     1                                                         esi(i,j)
              dm(i,j)=max(qv(i,j)+qb0-qsi(i,j),0.0)
              rsub1(i,j)=cs580*qsi(i,j)*rtair(i,j)*rtair(i,j)
              dep(i,j)=dm(i,j)/(1.+rsub1(i,j))
              if(tairc(i,j).le.-5.)then                              !meyers
                y4(i,j)=1./(tair(i,j)-c358)
                qsw(i,j)=rp0*exp(c172-c409*y4(i,j))
                xssi=min(ssi(i,j), 0.20)                            !cap at 20% super wrt ice
                r_nci=min(1.e-3*exp(-.639+12.96*xssi),1.)           !meyers et al. 1992
                if(r_nci.gt.15.) r_nci=15.                          !cap at 15000/liter
                pidep(i,j)=max(R32RT*1.e-4*SSI(i,j)*sqrt(r_nci)*y3(i,j)/     !meyers
     1                                                       dd(i,j),0.)
                ami50=4.8e-7
                dd(i,j)=max(1.e-9*r_nci/r00-qci(i,j,k)*1.e-9/ami50,0.) !test vs old qci
                pint(i,j)=max(min(dd(i,j),dm(i,j)),0.)
                pint(i,j)=min(pint(i,j)+pidep(i,j),dep(i,j))
                if (pint(i,j).le.cmin) pint(i,j)=0.
                pt(i,j)=pt(i,j)+ascp*pint(i,j)
                qv(i,j)=qv(i,j)-pint(i,j)
                qi(i,j)=qi(i,j)+pint(i,j)
              endif
            endif
          endif

c      itaobraun=0 ! using original way for pint and pidep
       itaobraun=1 ! using scott braun's way for pint and pidep
       if(improve.ge.3) itaobraun=0

       if ( itaobraun.eq.0 ) then      ! tao's original
       cn0=1.e-8
cc     beta=-.6
       elseif ( itaobraun.eq.1 ) then  ! scott's
       cn0=1.e-6
c      cn0=1.e-8  ! 4/22/02 special, still use tao's
cc     beta=-.46
       endif

        if ( itaobraun.eq.1 ) then
           tair(i,j)=(pt(i,j)+tb0)*pi0
           if (tair(i,j) .lt. t0) then
c             if (qi(i,j) .le. cmin) qi(i,j)=0.
              if (qi(i,j) .le. cmin2) qi(i,j)=0.
               tairc(i,j)=tair(i,j)-t0
               rtair(i,j)=1./(tair(i,j)-c76)
               y2(i,j)=exp(c218-c580*rtair(i,j))
              qsi(i,j)=rp0*y2(i,j)
               esi(i,j)=c610*y2(i,j)
              ssi(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
                        ami50=3.76e-8

cccshie with scott braun's help, insert "pidep" and change "betah", "c0" in rou
ccc              "consat" (2d), "consatrh" (3d)
ccc        if ( itaobraun.eq.1 ) --> betah=0.5*beta=-.46*0.5=-0.23;   cn0=1.e-6
ccc        if ( itaobraun.eq.0 ) --> betah=0.5*beta=-.6*0.5=-0.30;    cn0=1.e-8

             y1(i,j)=1./tair(i,j)

cccshie with scott braun's help, insert a restriction on ice collection that ic
ccc              should be stopped at -30 c (with cn0=1.e-6, beta=-.46)

             tairccri=tairc(i,j)          ! in degree c
             if(tairccri.le.-30.) tairccri=-30.

c            y2(i,j)=exp(betah*tairc(i,j))
             y2(i,j)=exp(betah*tairccri)
             y3(i,j)=sqrt(qi(i,j))
             dd(i,j)=y1(i,j)*(rn10a*y1(i,j)-rn10b)+rn10c*tair(i,j)
     1                                          /esi(i,j)
          pidep(i,j)=max(r32rt*ssi(i,j)*y2(i,j)*y3(i,j)/dd(i,j), 0.e0)

           r_nci=min(cn0*exp(beta*tairc(i,j)),1.)      ! cccshie 4/18/02
c          r_nci=min(1.e-6*exp(-.46*tairc(i,j)),1.)

           dd(i,j)=max(1.e-9*r_nci/r00-qi(i,j)*1.e-9/ami50, 0.)
                dm(i,j)=max( (qv(i,j)+qb0-qsi(i,j)), 0.0)
                rsub1(i,j)=cs580*qsi(i,j)*rtair(i,j)*rtair(i,j)
              dep(i,j)=dm(i,j)/(1.+rsub1(i,j))
              pint(i,j)=max(min(dd(i,j), dm(i,j)), 0.)

c             pint(i,j)=min(pint(i,j), dep(i,j))
              pint(i,j)=min(pint(i,j)+pidep(i,j), dep(i,j))  ! cccshie 4/15/02

c              if (pint(i,j) .le. cmin) pint(i,j)=0.
               if (pint(i,j) .le. cmin2) pint(i,j)=0.
              pt(i,j)=pt(i,j)+ascp*pint(i,j)
              qv(i,j)=qv(i,j)-pint(i,j)
              qi(i,j)=qi(i,j)+pint(i,j)
           endif
        endif  ! if ( itaobraun.eq.1 ) ! scott's
        if ( itaobraun.eq.0 .and. improve.le.2) then
	 tair(i,j)=(pt(i,j)+tb0)*pi0
	 if (tair(i,j) .lt. t0) then
	 if (qi(i,j) .le. cmin2) qi(i,j)=0.
	   tairc(i,j)=tair(i,j)-t0
	   dd(i,j)=r31r*exp(beta*tairc(i,j))
	    rtair(i,j)=1./(tair(i,j)-c76)
	     y2(i,j)=exp(c218-c580*rtair(i,j))
	    qsi(i,j)=rp0*y2(i,j)
	    esi(i,j)=c610*y2(i,j)
	    ssi(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
	   dm(i,j)=max( (qv(i,j)+qb0-qsi(i,j)), 0.)
	    rsub1(i,j)=cs580*qsi(i,j)*rtair(i,j)*rtair(i,j)
	  dep(i,j)=dm(i,j)/(1.+rsub1(i,j))
	  pint(i,j)=max(min(dd(i,j), dm(i,j)), 0.)
	     y1(i,j)=1./tair(i,j)
	     y2(i,j)=exp(betah*tairc(i,j))
	     y3(i,j)=sqrt(qi(i,j))
	   dd(i,j)=y1(i,j)*(rn10a*y1(i,j)-rn10b)+rn10c*tair(i,j)
     1                                      /esi(i,j)
	  pidep(i,j)=max(r32rt*ssi(i,j)*y2(i,j)*y3(i,j)/dd(i,j), 0.)
	  pint(i,j)=pint(i,j)+pidep(i,j)
	   pint(i,j)=min(pint(i,j),dep(i,j))
cc             if (pint(i,j) .le. cmin2) pint(i,j)=0.
	  pt(i,j)=pt(i,j)+ascp*pint(i,j)
	  qv(i,j)=qv(i,j)-pint(i,j)
	  qi(i,j)=qi(i,j)+pint(i,j)
	 endif
        endif  ! if ( itaobraun.eq.0 ) ! tao's original
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  250  continue

        if (improve.ge.1) new_ice_sat=2
        if (improve.ge.3) new_ice_sat=3

c*****   tao et al (1989) saturation technique  ***********************
        if (new_ice_sat .eq. 0) then ! cccshie by tao 5/3/01, shie 11/16/01 3d
c$doacross local(j,i)
       do 275 j=1,1
       do 275 i=1,1
	 tair(i,j)=(pt(i,j)+tb0)*pi0
	cnd(i,j)=rt0*(tair(i,j)-t00)
	dep(i,j)=rt0*(t0-tair(i,j))
	  y1(i,j)=1./(tair(i,j)-c358)
	  y2(i,j)=1./(tair(i,j)-c76)
	 qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
	 qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
	  dd(i,j)=cp409*y1(i,j)*y1(i,j)
	  dd1(i,j)=cp580*y2(i,j)*y2(i,j)
	 if (qc(i,j).le.cmin) qc(i,j)=cmin
	 if (qi(i,j).le.cmin) qi(i,j)=cmin
	 if (tair(i,j).ge.t0) then
	  dep(i,j)=0.0
	  cnd(i,j)=1.
	  qi(i,j)=0.0
	 endif
	 if (tair(i,j).lt.t00) then
	  cnd(i,j)=0.0
	  dep(i,j)=1.
	  qc(i,j)=0.0
	 endif
	  y5(i,j)=avcp*cnd(i,j)+ascp*dep(i,j)
	   y1(i,j)=qc(i,j)*qsw(i,j)/(qc(i,j)+qi(i,j))
	   y2(i,j)=qi(i,j)*qsi(i,j)/(qc(i,j)+qi(i,j))
	  y4(i,j)=dd(i,j)*y1(i,j)+dd1(i,j)*y2(i,j)
	 qvs(i,j)=y1(i,j)+y2(i,j)
	 rsub1(i,j)=(qv(i,j)+qb0-qvs(i,j))/(1.+y4(i,j)*y5(i,j))
	cnd(i,j)=cnd(i,j)*rsub1(i,j)
	dep(i,j)=dep(i,j)*rsub1(i,j)
	 if (qc(i,j).le.cmin) qc(i,j)=0.
	 if (qi(i,j).le.cmin) qi(i,j)=0.
cc    ******   condensation or evaporation of qc  ******
	 cnd(i,j)=max(-qc(i,j),cnd(i,j))
cc    ******   deposition or sublimation of qi    ******
	 dep(i,j)=max(-qi(i,j),dep(i,j))
	pt(i,j)=pt(i,j)+avcp*cnd(i,j)+ascp*dep(i,j)
	qv(i,j)=qv(i,j)-cnd(i,j)-dep(i,j)
	qc(i,j)=qc(i,j)+cnd(i,j)
	qi(i,j)=qi(i,j)+dep(i,j)
  275  continue
       endif ! cccshie by tao 5/3/01, shie 11/16/01 3d
cccshie 11/16/01 !  call tao et al (1989) saturation technique twice
        if (new_ice_sat .eq. 1) then
c$doacross local(j,j)
           do j=1,1
           do i=1,1
              tair(i,j)=(pt(i,j)+tb0)*pi0
c             cnd1(i,j)=rt0*(tair(i,j)-t00)
c             dep1(i,j)=rt0*(t0-tair(i,j))
              cnd(i,j)=rt0*(tair(i,j)-t00)  ! cccshie by tao 5/3/01
              dep(i,j)=rt0*(t0-tair(i,j))   ! cccshie by tao 5/3/01
                y1(i,j)=1./(tair(i,j)-c358)
                y2(i,j)=1./(tair(i,j)-c76)
              qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
              qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
              dd(i,j)=cp409*y1(i,j)*y1(i,j)
              dd1(i,j)=cp580*y2(i,j)*y2(i,j)
              y5(i,j)=avcp*cnd(i,j)+ascp*dep(i,j) ! cccshie by tao 5/3/01
              y1(i,j)=rt0*(tair(i,j)-t00)*qsw(i,j) ! cccshie by tao 5/3/01
              y2(i,j)=rt0*(t0-tair(i,j))*qsi(i,j)  ! cccshie by tao 5/3/01
c             if (qc(i,j).le.cmin) qc(i,j)=cmin  ! cccshie by tao 5/3/01
c             if (qi(i,j).le.cmin) qi(i,j)=cmin  ! cccshie by tao 5/3/01
              if (tair(i,j).ge.t0) then
c                 dep1(i,j)=0.0
c                 cnd1(i,j)=1.
c                 qi(i,j)=0.0
                  dep(i,j)=0.0 ! cccshie by tao 5/3/01
                  cnd(i,j)=1.  ! cccshie by tao 5/3/01
                  y2(i,j)=0.   ! cccshie by tao 5/3/01
                  y1(i,j)=qsw(i,j) ! cccshie by tao 5/3/01
              endif
              if (tair(i,j).lt.t00) then
                 cnd(i,j)=0.0  ! cccshie by tao 5/3/01
                 dep(i,j)=1.   ! cccshie by tao 5/3/01
                 y2(i,j)=qsi(i,j) ! cccshie by tao 5/3/01
                 y1(i,j)=0.     ! cccshie by tao 5/3/01
c                cnd1(i,j)=0.0
c                dep1(i,j)=1.
c                qc(i,j)=0.0
              endif
c             y5(i,j)=avcp*cnd1(i,j)+ascp*dep1(i,j) ! cccshie by tao 5/3/01
c             y1(i,j)=qc(i,j)*qsw(i,j)/(qc(i,j)+qi(i,j)) ! cccshie by tao 5/3/0
c             y2(i,j)=qi(i,j)*qsi(i,j)/(qc(i,j)+qi(i,j)) ! cccshie by tao 5/3/0
              y4(i,j)=dd(i,j)*y1(i,j)+dd1(i,j)*y2(i,j)
              qvs(i,j)=y1(i,j)+y2(i,j)
              rsub1(i,j)=(qv(i,j)+qb0-qvs(i,j))/(1.+y4(i,j)*y5(i,j))
             cnd(i,j)=cnd(i,j)*rsub1(i,j) ! cccshie by tao 5/3/01
             dep(i,j)=dep(i,j)*rsub1(i,j) ! cccshie by tao 5/3/01
c            cnd1(i,j)=cnd1(i,j)*rsub1(i,j)
c            dep1(i,j)=dep1(i,j)*rsub1(i,j)
c             if (qc(i,j).le.cmin) qc(i,j)=0. ! cccshie by tao 5/3/01
c             if (qi(i,j).le.cmin) qi(i,j)=0. ! cccshie by tao 5/3/01
cc    ******   condensation or evaporation of qc  ******
c            cnd1(i,j)=max(-qc(i,j),cnd1(i,j))
             cnd(i,j)=max(-qc(i,j),cnd(i,j)) ! cccshie by tao 5/3/01
cc    ******   deposition or sublimation of qi    ******
             dep(i,j)=max(-qi(i,j),dep(i,j)) ! cccshie by tao 5/3/01
             pt(i,j)=pt(i,j)+avcp*cnd(i,j)+ascp*dep(i,j) ! cccshie by tao 5/3/0
             qv(i,j)=qv(i,j)-cnd(i,j)-dep(i,j) ! cccshie by tao 5/3/01
             qc(i,j)=qc(i,j)+cnd(i,j) ! cccshie by tao 5/3/01
             qi(i,j)=qi(i,j)+dep(i,j) ! cccshie by tao 5/3/01
c            dep1(i,j)=max(-qi(i,j),dep1(i,j))
c            pt(i,j)=pt(i,j)+avcp*cnd1(i,j)+ascp*dep1(i,j)
c            qv(i,j)=qv(i,j)-cnd1(i,j)-dep1(i,j)
c            qc(i,j)=qc(i,j)+cnd1(i,j)
c            qi(i,j)=qi(i,j)+dep1(i,j)
           enddo
           enddo
        endif ! if (new_ice_sat .eq. 1)
c
cc
c	
	    if (new_ice_sat .eq. 2) then
          do j=1,1
	  do i=1,1
          dep(i,j)=0.0
          cnd(i,j)=0.0
          tair(i,j)=(pt(i,j)+tb0)*pi0
          if (tair(i,j) .ge. 253.16) then
              y1(i,j)=1./(tair(i,j)-c358)
              qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
              dd(i,j)=cp409*y1(i,j)*y1(i,j)
              dm(i,j)=qv(i,j)+qb0-qsw(i,j)
              cnd(i,j)=dm(i,j)/(1.+avcp*dd(i,j)*qsw(i,j))
cc    ******   condensation or evaporation of qc  ******
              cnd(i,j)=max(-qc(i,j), cnd(i,j))
             pt(i,j)=pt(i,j)+avcp*cnd(i,j)
             qv(i,j)=qv(i,j)-cnd(i,j)
             qc(i,j)=qc(i,j)+cnd(i,j)
         endif
          if (tair(i,j) .le. 258.16) then
cc             cnd(i,j)=0.0
           y2(i,j)=1./(tair(i,j)-c76)
           qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
          dd1(i,j)=cp580*y2(i,j)*y2(i,j)
         dep(i,j)=(qv(i,j)+qb0-qsi(i,j))/(1.+ascp*dd1(i,j)*qsi(i,j))
cc    ******   deposition or sublimation of qi    ******
             dep(i,j)=max(-qi(i,j),dep(i,j))
             pt(i,j)=pt(i,j)+ascp*dep(i,j)
             qv(i,j)=qv(i,j)-dep(i,j)
             qi(i,j)=qi(i,j)+dep(i,j)
         endif
       enddo
       enddo
      endif
c
cc
        if (new_ice_sat .eq. 3) then                              !Lang et al. 2007b
          do j=1,1
        do i=1,1
          dep(i,j)=0.0
          cnd(i,j)=0.0
          tair(i,j)=(pt(i,j)+tb0)*pi0
          if (tair(i,j).ge.t00) THEN
            iter=0
            y2(i,j)=999.
            tairI(i,j)=tair(i,j)
            do while(iter.lt.10.and.y2(i,j).gt.0.01)
             y1(i,j)=1./(tair(i,j)-c358)
             qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
             dd(i,j)=cp409*y1(i,j)*y1(i,j)
             dm(i,j)=qv(i,j)+qb0-qsw(i,j)
             cnd(i,j)=dm(i,j)/(1.+avcp*dd(i,j)*qsw(i,j))
             tairN(i,j)=tairI(i,j)+avcp*cnd(i,j)*pi0
             iter=iter+1
             if(iter.eq.10) write(6,*) 'no convergence-w'
             y2(i,j)=abs(tairN(i,j)-tair(i,j))
             tair(i,j)=0.5*(tairN(i,j)+tair(i,j))
           enddo
             y1(i,j)=1./(tairN(i,j)-c358)
             qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
             dd(i,j)=cp409*y1(i,j)*y1(i,j)
             dm(i,j)=qv(i,j)+qb0-qsw(i,j)
             cnd(i,j)=dm(i,j)/(1.+avcp*dd(i,j)*qsw(i,j))
cc    ******   condensation or evaporation of qc  ******
            cnd(i,j)=max(-qc(i,j),cnd(i,j))
            pt(i,j)=pt(i,j)+avcp*cnd(i,j)
            qv(i,j)=qv(i,j)-cnd(i,j)
            qc(i,j)=qc(i,j)+cnd(i,j)
          endif
          if (tair(i,j).le.273.16) THEN
cc    ******   deposition or sublimation of qi    ******
            y1(i,j)=1./(tair(i,j)-c358)
            qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
            y2(i,j)=1./(tair(i,j)-c76)
            qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
            y3(i,j)=1.+min((qsw(i,j)-qsi(i,j))/qsi(i,j), 0.20)             !lang 2007
            if(tair(i,j).le.268.16.and.(qv(i,j)+qb0.gt.y3(i,j)
     1                                                   *qsi(i,j)))then
             iter=0
             y4(i,j)=999.
             tairI(i,j)=tair(i,j)
             do while(iter.lt.10.and.y4(i,j).gt.0.01)
                y2(i,j)=1./(tair(i,j)-c76)
                qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
                dd1(i,j)=cp580*y2(i,j)*y2(i,j)
                dep(i,j)=(qv(i,j)+qb0-y3(i,j)*qsi(i,j))/
     1                               (1.+ascp*dd1(i,j)*y3(i,j)*qsi(i,j))
                tairN(i,j)=tairI(i,j)+ascp*dep(i,j)*pi0
                iter=iter+1
                if(iter.eq.10) write(6,*) 'no convergence-i'
                y4(i,j)=abs(tairN(i,j)-tair(i,j))
                tair(i,j)=0.5*(tairN(i,j)+tair(i,j))
             enddo
             y2(i,j)=1./(tairN(i,j)-c76)
             qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
             dd1(i,j)=cp580*y2(i,j)*y2(i,j)
             dep(i,j)=(qv(i,j)+qb0-y3(i,j)*qsi(i,j))/(1.+ascp*dd1(i,j)
     1                                                        *qsi(i,j))
            elseif(qv(i,j)+qb0.lt.qsi(i,j).and.qi(i,j).gt.cmin)then
             iter=0
             y4(i,j)=999.
             tairI(i,j)=tair(i,j)
             do while(iter.lt.10.and.y4(i,j).gt.0.01)
                y2(i,j)=1./(tair(i,j)-c76)
                qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
                dd1(i,j)=cp580*y2(i,j)*y2(i,j)
                dep(i,j)=(qv(i,j)+qb0-qsi(i,j))/
     1                                       (1.+ascp*dd1(i,j)*qsi(i,j))
                tairN(i,j)=tairI(i,j)+ascp*dep(i,j)*pi0
                iter=iter+1
                if(iter.eq.10) write(6,*) 'no convergence-i'
                y4(i,j)=abs(tairN(i,j)-tair(i,j))
                tair(i,j)=0.5*(tairN(i,j)+tair(i,j))
             enddo
             y2(i,j)=1./(tairN(i,j)-c76)
             qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
             dd1(i,j)=cp580*y2(i,j)*y2(i,j)
             dep(i,j)=(qv(i,j)+qb0-qsi(i,j))/(1.+ascp*dd1(i,j)*qsi(i,j))
             dep(i,j)=max(-qi(i,j),dep(i,j))
            endif
            pt(i,j)=pt(i,j)+ascp*dep(i,j)
            qv(i,j)=qv(i,j)-dep(i,j)
            qi(i,j)=qi(i,j)+dep(i,j)
          endif
         enddo
         enddo
      endif

        if (new_ice_sat.eq.9) then  !non-iterative
          do j=1,1
            do i=1,1
            dep(i,j)=0.0
            cnd(i,j)=0.0
            tair(i,j)=(pt(i,j)+tb0)*pi0
            if (tair(i,j).ge.t00) THEN
              y1(i,j)=1./(tair(i,j)-c358)
              qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
              dd(i,j)=cp409*y1(i,j)*y1(i,j)
              dm(i,j)=qv(i,j)+qb0-qsw(i,j)
              cnd(i,j)=dm(i,j)/(1.+avcp*dd(i,j)*qsw(i,j))
cc    ******   condensation or evaporation of qc  ******
              cnd(i,j)=max(-qc(i,j),cnd(i,j))
              pt(i,j)=pt(i,j)+avcp*cnd(i,j)
              qv(i,j)=qv(i,j)-cnd(i,j)
              qc(i,j)=qc(i,j)+cnd(i,j)
            endif
            if (tair(i,j).le.273.16) THEN
cc    ******   deposition or sublimation of qi    ******
             y1(i,j)=1./(tair(i,j)-c358)
             qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
             y2(i,j)=1./(tair(i,j)-c76)
             qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
             y3(i,j)=1.+min((qsw(i,j)-qsi(i,j))/qsi(i,j), 0.20)             !lang 2007
             y4(i,j)=qsi(i,j)*y3(i,j)
           if(tair(i,j).le.268.16.and.(qv(i,j)+qb0.gt.y4(i,j)))then
             dd1(i,j)=cp580*y2(i,j)*y2(i,j)
             dep(i,j)=(qv(i,j)+qb0-y4(i,j))/(1.+ascp*dd1(i,j)*y4(i,j))
            elseif(qv(i,j)+qb0.lt.qsi(i,j).and.qi(i,j).gt.cmin)then
             dd1(i,j)=cp580*y2(i,j)*y2(i,j)
             dep(i,j)=(qv(i,j)+qb0-qsi(i,j))/(1.+ascp*dd1(i,j)*qsi(i,j))
             dep(i,j)=max(-qi(i,j),dep(i,j))
            endif
              pt(i,j)=pt(i,j)+ascp*dep(i,j)
              qv(i,j)=qv(i,j)-dep(i,j)
              qi(i,j)=qi(i,j)+dep(i,j)
            endif
            enddo
            enddo
        endif

c
cc
c
c* 10 * psdep : deposition of qs                                  **10**
c* 20 * pgdep : deposition of qg                                  **20**
c$doacross local(j,i)
      do 280 j=1,1
      do 280 i=1,1
       psdep(i,j)=0.0
       pgdep(i,j)=0.0
       pssub(i,j)=0.0
       pgsub(i,j)=0.0
       tair(i,j)=(pt(i,j)+tb0)*pi0
       if (tair(i,j) .lt. t0) then
         if(qc(i,j)+qi(i,j).gt.1.e-5) then
	  dlt1(i,j)=1.
	 else
	  dlt1(i,j)=0.
	 endif
	 rtair(i,j)=1./(tair(i,j)-c76)
	 y2(i,j)=exp(c218-c580*rtair(i,j))
         qsi(i,j)=rp0*y2(i,j)
	 esi(i,j)=c610*y2(i,j)
         if(improve.ge.3)then                                    !Lang et al. 2007b
         SSI(I,J)=(QV(I,J)+QB0)/QSI(I,J)-1.
         IF(DLT1(I,J).EQ.1.) SSI(I,J)=max(SSI(I,J),0.)
         IF(DLT1(I,J).EQ.0.) SSI(I,J)=min(SSI(I,J),0.)
          DM(I,J)=QV(I,J)+QB0-QSI(I,J)
          RSUB1(I,J)=CS580*QSI(I,J)*RTAIR(I,J)*RTAIR(I,J)
          DD1(I,J)=DM(I,J)/(1.+RSUB1(I,J))
           Y3(I,J)=1./TAIR(I,J)
          DD(I,J)=Y3(I,J)*(RN10A*Y3(I,J)-RN10B)+RN10C*TAIR(I,J)/ESI(I,J)
           TAIRC(I,J)=TAIR(I,J)-T0
           tns_funT(i,j)=min(20.,exp(-0.060000*tairc(i,j)))
           Y4(I,J)=R10T*SSI(I,J)*(R101R/ZS(I,J)**2+R102RF/ZS(I,J)**BSH5)
     1                         /DD(I,J) * tns_funT(i,j)              ! Lang et al. 2007b
         PSDEP(I,J)=max(-QS(I,J), Y4(I,J))                           !lang 2007
          DD(I,J)=Y3(I,J)*(RN20A*Y3(I,J)-RN20B)+RN10C*TAIR(I,J)/ESI(I,J)
          Y2(I,J)=R191R/ZG(I,J)**2+R192RF/ZG(I,J)**BGH5
         PGDEP(I,J)=MAX(-qg(i,j), R20T*SSI(I,J)*Y2(I,J)/DD(I,J))     !lang 2007
C     ******************************************************************
         IF(DLT1(I,J).EQ.1.)THEN
           Y1(I,J)=MIN(PSDEP(I,J)+PGDEP(I,J), max(0.,DD1(I,J)) )
           IF(PSDEP(I,J).ge.DD1(I,J))THEN
            PSDEP(I,J)=DD1(I,J)
            PGDEP(I,J)=0.
           ENDIF
           IF(DD1(I,J).gt.PSDEP(I,J).and.(PSDEP(I,J)+PGDEP(I,J)).gt.
     1      DD1(I,J)) PGDEP(I,J)=Y1(I,J)-PSDEP(I,J)
         ENDIF
         IF(DLT1(I,J).EQ.0.)THEN
           Y1(I,J)=MAX(PSDEP(I,J)+PGDEP(I,J), min(0.,DD1(I,J)) )
           IF(DD1(I,J).gt.(PSDEP(I,J)+PGDEP(I,J)))THEN
            Y3(I,J)=(PSDEP(I,J)+PGDEP(I,J))
            IF(Y3(I,J).ne.0.0)THEN
             PSDEP(I,J)=PSDEP(I,J)/Y3(I,J)*DD1(I,J)
             PGDEP(I,J)=PGDEP(I,J)/Y3(I,J)*DD1(I,J)
            ENDIF
           ENDIF
         ENDIF
          PSSUB(i,j)=min(PSDEP(i,j), 0.)
          PSDEP(i,j)=max(PSDEP(i,j), 0.)
          PGSUB(i,j)=min(PGDEP(i,j), 0.)
          PGDEP(i,j)=max(PGDEP(i,j), 0.)
         else                      ! original
	  ssi(i,j)=dlt1(i,j)*((qv(i,j)+qb0)/qsi(i,j)-1.)
          dm(i,j)=qv(i,j)+qb0-qsi(i,j)
          rsub1(i,j)=cs580*qsi(i,j)*rtair(i,j)*rtair(i,j)
          dd1(i,j)=max(dm(i,j)/(1.+rsub1(i,j)),0.0)
          y3(i,j)=1./tair(i,j)
          dd(i,j)=y3(i,j)*(rn10a*y3(i,j)-rn10b)+rn10c*tair(i,j)
     1                                       /esi(i,j)
          y4(i,j)=r10t*ssi(i,j)*(r101r/zs(i,j)**2+r102rf/zs(i,j)
     1                                       **bsh5)/dd(i,j)
	  psdep(i,j)=max(y4(i,j), 0.0)
          dd(i,j)=y3(i,j)*(rn20a*y3(i,j)-rn20b)+rn10c*tair(i,j)
     1                                       /esi(i,j)
          y2(i,j)=r191r/zg(i,j)**2+r192rf/zg(i,j)**bgh5
	  pgdep(i,j)=max(r20t*ssi(i,j)*y2(i,j)/dd(i,j),0.0)
c     ******************************************************************
          y1(i,j)=min(psdep(i,j)+pgdep(i,j),dd1(i,j))
          pgdep(i,j)=y1(i,j)-psdep(i,j)
         endif

	 pt(i,j)=pt(i,j)+ascp*y1(i,j)
	 qv(i,j)=qv(i,j)-y1(i,j)
	 qs(i,j)=qs(i,j)+psdep(i,j)+pssub(i,j)
	 qg(i,j)=qg(i,j)+pgdep(i,j)+pgsub(i,j)
       endif


c* 23 * ern : evaporation of qr                                   **23**
	ern(i,j)=0.0
	if (qr(i,j) .gt. 0.0) then
	 tair(i,j)=(pt(i,j)+tb0)*pi0
	  rtair(i,j)=1./(tair(i,j)-c358)
	   y2(i,j)=exp( c172-c409*rtair(i,j) )
	  esw(i,j)=c610*y2(i,j)
	  qsw(i,j)=rp0*y2(i,j)
	  ssw(i,j)=(qv(i,j)+qb0)/qsw(i,j)-1.
	  dm(i,j)=qv(i,j)+qb0-qsw(i,j)
	   rsub1(i,j)=cv409*qsw(i,j)*rtair(i,j)*rtair(i,j)
	  dd1(i,j)=max(-dm(i,j)/(1.+rsub1(i,j)),0.0)
c           y1(i,j)=r00*qrn(i,j,k)
c         ern(i,j)=(((1.6+124.9*y1(i,j)**.2046)*y1(i,j)**.525)
c    1          /(2.55e6/(p00*qsw(i,j))+5.4e5))*(-dm(i,j)
c    2                                        /(r00*qsw(i,j)))*d2t
	    y3(i,j)=1./tair(i,j)
	   dd(i,j)=y3(i,j)*(rn30a*y3(i,j)-rn10b)+rn10c*tair(i,j)
     1                                          /esw(i,j)
	  y1(i,j)=-r23t*ssw(i,j)*(r231r/zr(i,j)**2+r232rf/
     1                                        zr(i,j)**3)/dd(i,j)
	  ern(i,j)=min(dd1(i,j),qr(i,j),max(y1(i,j),0.0))
	  pt(i,j)=pt(i,j)-avcp*ern(i,j)
	  qv(i,j)=qv(i,j)+ern(i,j)
	  qr(i,j)=qr(i,j)-ern(i,j)
	 endif
  280  continue
c* 30 * pmltg : evaporation of melting qg                         **30**
c* 33 * pmlts : evaporation of melting qs                         **33**

c$doacross local(j.i)
      do 300 j=1,1
      do 300 i=1,1
	pmlts(i,j)=0.0
	pmltg(i,j)=0.0
	 tair(i,j)=(pt(i,j)+tb0)*pi0
	if (tair(i,j) .ge. t0) then
c           rtair(i,j)=1./(tair(i,j)-c358)
	   rtair(i,j)=1./(t0-c358)
	    y2(i,j)=exp( c172-c409*rtair(i,j) )
	   esw(i,j)=c610*y2(i,j)
	   qsw(i,j)=rp0*y2(i,j)
	   ssw(i,j)=1.-(qv(i,j)+qb0)/qsw(i,j)
	   dm(i,j)=qsw(i,j)-qv(i,j)-qb0
	   rsub1(i,j)=cv409*qsw(i,j)*rtair(i,j)*rtair(i,j)
	  dd1(i,j)=max(dm(i,j)/(1.+rsub1(i,j)),0.0)
	    y3(i,j)=1./tair(i,j)
	   dd(i,j)=y3(i,j)*(rn30a*y3(i,j)-rn10b)+rn10c*tair(i,j)
     1                                       /esw(i,j)
	   y1(i,j)=r30t*ssw(i,j)*(r191r/zg(i,j)**2+r192rf
     1             /zg(i,j)**bgh5)
     1                          /dd(i,j)
	  pmltg(i,j)=min(qg(i,j),max(y1(i,j),0.0))
	   y1(i,j)=r33t*ssw(i,j)*(r331r/zs(i,j)**2+r332rf
     1                                     /zs(i,j)**bsh5)/dd(i,j)
	  pmlts(i,j)=min(qs(i,j),max(y1(i,j),0.0))
	   y1(i,j)=min(pmltg(i,j)+pmlts(i,j),dd1(i,j))
	  pmltg(i,j)=y1(i,j)-pmlts(i,j)
	  pt(i,j)=pt(i,j)-ascp*y1(i,j)
	  qv(i,j)=qv(i,j)+y1(i,j)
	  qs(i,j)=qs(i,j)-pmlts(i,j)
	  qg(i,j)=qg(i,j)-pmltg(i,j)
	endif
c       if (qv(i,j)+qb0 .le. 0.) qv(i,j)=-qb0
	 if(qc(i,j).lt.cmin) qc(i,j)=0.0
	 if(qr(i,j).lt.cmin) qr(i,j)=0.0
	 if(qi(i,j).lt.cmin) qi(i,j)=0.0
	 if(qs(i,j).lt.cmin) qs(i,j)=0.0
	 if(qg(i,j).lt.cmin) qg(i,j)=0.0

ccc    store the forcing term and updated in main program??

       dpt(i,j,k)=pt(i,j)
       dqv(i,j,k)=qv(i,j)
       qcl(i,j,k)=qc(i,j)
       qrn(i,j,k)=qr(i,j)
       qci(i,j,k)=qi(i,j)
       qcs(i,j,k)=qs(i,j)
       qcg(i,j,k)=qg(i,j)

  300  continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc     henry:  please take a look  (start)
	   sddd=0.
	   ssss=0.
	   shhh=0.
	   sccc=0.
	   smmm=0.
	   sfff=0.

       do 305 j=1,1
       do 305 i=1,1
      	 dd(i,j)=max(-cnd(i,j), 0.)
	     cnd(i,j)=max(cnd(i,j), 0.)
	     dd1(i,j)=max(-dep(i,j), 0.)
	     dep(i,j)=max(dep(i,j), 0.)
  305   continue

c D. Posselt Eliminate following statistics section
!        do j=1,1
!         do i=1,1
! 	 tmpd(i,j,1)=cnd(i,j)
! 	 tmpd(i,j,2)=dd(i,j)+ern(i,j)
! 	 tmpd(i,j,3)=dep(i,j)+pint(i,j)+psdep(i,j)+pgdep(i,j)
! 	 tmpd(i,j,4)=dd1(i,j)+pmlts(i,j)+pmltg(i,j)
! 	 tmpd(i,j,5)=rsw(i,j,k)*dt
! 	 tmpd(i,j,6)=rlw(i,j,k)*dt
! 	 tmpd(i,j,7)=psmlt(i,j)+pgmlt(i,j)+pimlt(i,j)
! 	 tmpd(i,j,8)=dt*(psacw(i,j)+piacr(i,j)+psfw(i,j)
!      1             +pgfr(i,j)+dgacw(i,j)+dgacr(i,j)+psacr(i,j))
!      2             -qracs(i,j)+pihom(i,j)+pidw(i,j)
! 	 
!         enddo
!        enddo
!        do nsum=1,8
!           sumd(nsum)=0.0
!        enddo
! 
!        call RealSum3D(tmpd,sumd,8)
! 
!        scc=sumd(1)
!        see=sumd(2)
!        sddd=sumd(3)
!        ssss=sumd(4)
!        shhh=sumd(5)
!        sccc=sumd(6)
!        smmm=sumd(7)
!        sfff=sumd(8)
! c
! 	s_dep(k)=s_dep(k)+sddd
! 	s_sub(k)=s_sub(k)+ssss
! 	s_qrs(k)=s_qrs(k)+shhh
! 	s_qrl(k)=s_qrl(k)+sccc
! 	s_mel(k)=s_mel(k)+smmm
! 	s_frz(k)=s_frz(k)+sfff
c
! 	sc(k)=scc+sc(k)
! 	se(k)=see+se(k)
cc     henry:  please take a look  (end)

cc    ***   statistics for convective and anvil regimes   ***********

       if (id .eq. 1) then
         
		   rdts=rd2t

c$doacross local(j,i)
	do 310 j=1,1 
      do 310 i=1,1
	   cnd(i,j)=cnd(i,j)*rdts
	   dd(i,j)=dd(i,j)*rdts
	   pint(i,j)=pint(i,j)*rdts
	   pidw(i,j)=pidw(i,j)*rdts
	   pimlt(i,j)=pimlt(i,j)*rdts
	   pihom(i,j)=pihom(i,j)*rdts
	   psmlt(i,j)=psmlt(i,j)*rdts
	   pgmlt(i,j)=pgmlt(i,j)*rdts
	   psdep(i,j)=psdep(i,j)*rdts
	   pgdep(i,j)=pgdep(i,j)*rdts
	   pmltg(i,j)=pmltg(i,j)*rdts
	   pmlts(i,j)=pmlts(i,j)*rdts
	   dd1(i,j)=dd1(i,j)*rdts
	   dep(i,j)=dep(i,j)*rdts
	   ern(i,j)=ern(i,j)*rdts
	   qracs(i,j)=qracs(i,j)*rdts
ctao
		physc(i,j,k)=cnd(i,j)
		physe(i,j,k)=dd(i,j)+ern(i,j)
		physd(i,j,k)=dep(i,j)+pint(i,j)+psdep(i,j)+pgdep(i,j)
		physs(i,j,k)=dd1(i,j)+pmlts(i,j)+pmltg(i,j)
		physm(i,j,k)=psmlt(i,j)+pgmlt(i,j)+pimlt(i,j)+qracs(i,j)
		physf(i,j,k)=psacw(i,j)+piacr(i,j)+psfw(i,j)
     1              +pgfr(i,j)+dgacw(i,j)+dgacr(i,j)+psacr(i,j)
     2              +pihom(i,j)+pidw(i,j)
ctao

! 	  dda0(i,j)=rsw(i,j,k) ! cccshie 7/1/02
! 	  ddb0(i,j)=rlw(i,j,k) ! cccshie 7/1/02
! 	 y1(i,j)=qc(i,j)+qr(i,j)+qi(i,j)+qs(i,j)+qg(i,j)
! c
! 	 dm(i,j)=a0*(rho1(k)*ww1(i,j,k)+rho1(k+1)*ww1(i,j,k+1)+
!      1              y0(i,j)*(rho1(k)*wb(k)+rho1(k+1)*wb(k+1)))
! 	 rq(i,j)=.005*(rho1(k)*(ww1(i,j,k)+wb(k))+
!      1               rho1(k+1)*(ww1(i,j,k+1)+wb(k+1)))/r00
  310    continue

c D. Posselt Eliminate following statistics section

! c    do 1050 kc=1,7
! 	   kc=4
! 	 do mt=1,4
! c$doacross local(j,i)
! 	 do j=1,1
!        do i=1,1
! 	  ibz(i,j,mt)=0
! 	   if(ics5(i,j,mt).eq.1) ibz(i,j,mt)=1
! 	 enddo
! 	 enddo
!        enddo
! 
! c$doacross local(j,i)
! 	do 315 j=1,1 
!       do 315 i=1,1
! 	    ibz(i,j,1)=1
!   315    continue
! 
! 	 do mt=1,4
! 
! 	 do 330 j=1,1
! 	 do 330 i=1,1
! 	   if(kc.eq.4) go to 330
! 	   if(kc.le.3) go to 36
! 	    if (rq(i,j).gt.rby(kc)) ibz(i,j,mt)=0
! 	    go to 330
!    36       if (rq(i,j).lt.rby(kc)) ibz(i,j,mt)=0
!   330    continue
! 	 enddo
! 
! c$doacross local(mt,ij,sww,scc,see,a1,a2,a3,a4,a5,a6,a7,a8,a9,
! c$&    a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,
! c$&    a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,
! c$&    a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50,a51,a52)
! 	 do 350 mt=1,4
! 
! 	  sww=0.0
! 	  scc=0.0
! 	  see=0.0
! 	  a1=0.0
! 	  a2=0.0
! 	  a3=0.0
! 	  a4=0.0
! 	  a5=0.0
! 	  a6=0.0
! 	  a7=0.0
! 	  a8=0.0
! 	  a9=0.0
! 	  a10=0.0
! 	  a11=0.0
! 	  a12=0.0
! 	  a13=0.0
! 	  a14=0.0
! 	  a15=0.0
! 	  a16=0.0
! 	  a17=0.0
! 	  a18=0.0
! 	  a19=0.0
! 	  a20=0.0
! 	  a21=0.0
! 	  a22=0.0
! 	  a23=0.0
! 	  a24=0.0
! 	  a25=0.0
! 	  a26=0.0
! 	  a27=0.0
! 	  a28=0.0
! 	  a29=0.0
! 	  a30=0.0
! 	  a31=0.0
! 	  a32=0.0
! 	  a33=0.0
! 	  a34=0.0
! 	  a35=0.0
! 	  a36=0.0
! 	  a37=0.0
! 	  a38=0.0
! 	  a39=0.0
! 	  a40=0.0
! 	  a41=0.0
! 	  a42=0.0
! 	  a43=0.0
! 	  a44=0.0
! 	  a45=0.
! 	  a46=0.
! 	  a47=0.
! 	  a48=0.
! 	  a49=0.
! 	  a50=0.
! 	  a51=0.
! 	  a52=0.
! 	  a54=0.
! 	  a55=0.
! 
!           do kk=1,55
!             do j=1,1
!               do i=1,1
!                 tmpd(i,j,kk)=0.0
!               enddo
!             enddo
!           enddo
!           do j=1,1
!             do i=1,1
! 
! 	     if(ibz(i,j,mt).eq.1) then
! 	       tmpd(i,j,51)=dm(i,j)
! 	       tmpd(i,j,52)=cnd(i,j)*asss(i,j)
! 	       tmpd(i,j,53)=(dd(i,j)+ern(i,j))*asss(i,j)
! 	       tmpd(i,j, 1)=(pihom(i,j)+pidw(i,j))*asss(i,j)
! 	       tmpd(i,j, 2)=pint(i,j)*asss(i,j)
! 	       tmpd(i,j, 3)=pgfr(i,j)*asss(i,j)
! 	       tmpd(i,j, 4)=psaut(i,j)*asss(i,j)
! 	       tmpd(i,j, 5)=psaci(i,j)*asss(i,j)
! 	       tmpd(i,j, 6)=psacw(i,j)*asss(i,j)
! 	       tmpd(i,j, 7)=praci(i,j)*asss(i,j)
! 	       tmpd(i,j, 8)=piacr(i,j)*asss(i,j)
! 	       tmpd(i,j, 9)=praut(i,j)*asss(i,j)
! 	       tmpd(i,j,10)=pracw(i,j)*asss(i,j)
! 	       tmpd(i,j,11)=psfw(i,j)*asss(i,j)
! 	       tmpd(i,j,12)=psfi(i,j)*asss(i,j)
! 	       tmpd(i,j,13)=(pgacs(i,j)+dgacs(i,j))*asss(i,j)
! 	       tmpd(i,j,14)=dgacw(i,j)*asss(i,j)
! 	       tmpd(i,j,15)=dgaci(i,j)*asss(i,j)
! 	       tmpd(i,j,16)=dgacr(i,j)*asss(i,j)
! 	       tmpd(i,j,17)=pmltg(i,j)*asss(i,j)
! 	       tmpd(i,j,18)=dep(i,j)*asss(i,j)
! 	       tmpd(i,j,19)=pracs(i,j)*asss(i,j)
! 	       tmpd(i,j,20)=psacr(i,j)*asss(i,j)
! 	       tmpd(i,j,21)=pmlts(i,j)*asss(i,j)
! 	       tmpd(i,j,22)=psmlt(i,j)*asss(i,j)
! 	       tmpd(i,j,23)=pgmlt(i,j)*asss(i,j)
! 	       tmpd(i,j,24)=psdep(i,j)*asss(i,j)
! 	       tmpd(i,j,25)=pimlt(i,j)*asss(i,j)
!                tmpd(i,j,26)=pgdep(i,j)*asss(i,j)
!                tmpd(i,j,27)=dd1(i,j)*asss(i,j)
!                tmpd(i,j,28)=praci(i,j)*dlt3(i,j)*asss(i,j)
!                tmpd(i,j,29)=piacr(i,j)*dlt3(i,j)*asss(i,j)
! !xp            tmpd(i,j,30)=psacr(i,j)*dlt2(i,j)*asss(i,j)
!        	       tmpd(i,j,30)=psacr(i,j)*asss(i,j)
!                tmpd(i,j,31)=qracs(i,j)*asss(i,j)
!                tmpd(i,j,32)=psacw(i,j)*dlt4(i,j)*asss(i,j)
!                tmpd(i,j,33)=qcl(i,j,k)*asss(i,j)
!                tmpd(i,j,34)=qrn(i,j,k)*asss(i,j)
!                tmpd(i,j,35)=qci(i,j,k)*asss(i,j)
!                tmpd(i,j,36)=qcs(i,j,k)*asss(i,j)
!                tmpd(i,j,37)=qcg(i,j,k)*asss(i,j)
!                tmpd(i,j,38)=ern(i,j)*asss(i,j)
!                tmpd(i,j,39)=wgacr(i,j)*asss(i,j)
!                tmpd(i,j,40)=qsacw(i,j)*asss(i,j)
!                tmpd(i,j,41)=dda0(i,j)*asss(i,j) 
!                tmpd(i,j,42)=ddb0(i,j)*asss(i,j) 
!                tmpd(i,j,43)=(qv(i,j)+qa1(k)-qa(k))*asss(i,j)
!                tmpd(i,j,44)=(pt(i,j)+ta1(k)-ta(k))*asss(i,j)
!                tmpd(i,j,45)=asss(i,j)
!                tmpd(i,j,46)=y1(i,j)*asss(i,j)
!                tmpd(i,j,47)=(psacw(i,j)+psfw(i,j)+dgacw(i,j)+piacr(i,j)
!      1                     +dgacr(i,j)+psacr(i,j)+pgfr(i,j)-qracs(i,j)
!      2                     +pihom(i,j)-pimlt(i,j)+pidw(i,j))*asss(i,j)
!                tmpd(i,j,48)=(y1(i,j)-qcl1(i,j,k)-qrn1(i,j,k)-qci1(i,j,k)
!      1                     -qcs1(i,j,k)-qcg1(i,j,k))*asss(i,j)
!                tmpd(i,j,49)=(qv(i,j)-dqv1(i,j,k))*asss(i,j)
!                tmpd(i,j,50)=(pt(i,j)-dpt1(i,j,k))*asss(i,j)
!                tmpd(i,j,54)=pimm(i,j)*asss(i,j)
!                tmpd(i,j,55)=pcfr(i,j)*asss(i,j)
!             endif
!            enddo
!           enddo
!           do nsum=1,55
!             sumd(nsum)=0.0
!           enddo
! 
!           call RealSum3D(tmpd,sumd,55)
! 
!           sww=sumd(51)
!           scc=sumd(52)
!           see=sumd(53)
!           a1 =sumd( 1)
!           a2 =sumd( 2)
!           a3 =sumd( 3)
!           a4 =sumd( 4)
!           a5 =sumd( 5)
!           a6 =sumd( 6)
!           a7 =sumd( 7)
!           a8 =sumd( 8)
!           a9 =sumd( 9)
!           a10=sumd(10)
!           a11=sumd(11)
!           a12=sumd(12)
!           a13=sumd(13)
!           a14=sumd(14)
!           a15=sumd(15)
!           a16=sumd(16)
!           a17=sumd(17)
!           a18=sumd(18)
!           a19=sumd(19)
!           a20=sumd(20)
!           a21=sumd(21)
!           a22=sumd(22)
!           a23=sumd(23)
!           a24=sumd(24)
!           a25=sumd(25)
!           a26=sumd(26)
!           a27=sumd(27)
!           a28=sumd(28)
!           a29=sumd(29)
!           a30=sumd(30)
!           a31=sumd(31)
!           a32=sumd(32)
!           a33=sumd(33)
!           a34=sumd(34)
!           a35=sumd(35)
!           a36=sumd(36)
!           a37=sumd(37)
!           a38=sumd(38)
!           a39=sumd(39)
!           a40=sumd(40)
!           a41=sumd(41)
!           a42=sumd(42)
!           a43=sumd(43)
!           a44=sumd(44)
!           a45=sumd(45)
!           a46=sumd(46)
!           a47=sumd(47)
!           a48=sumd(48)
!           a49=sumd(49)
!           a50=sumd(50)
!           a54=sumd(54)
!           a55=sumd(55)
! 
! 
! 	smf0(k,mt,kc)=sww+smf0(k,mt,kc)
! 	coc(k,mt,kc)=scc+coc(k,mt,kc)
! 	coe(k,mt,kc)=see+coe(k,mt,kc)
! 	thom(k,mt,kc)=thom(k,mt,kc)+a1
! 	tdw(k,mt,kc)=tdw(k,mt,kc)+a2
! 	tmlt(k,mt,kc)=tmlt(k,mt,kc)+a3
! 	saut(k,mt,kc)=saut(k,mt,kc)+a4
! 	saci(k,mt,kc)=saci(k,mt,kc)+a5
! 	sacw(k,mt,kc)=sacw(k,mt,kc)+a6
! 	raci(k,mt,kc)=raci(k,mt,kc)+a7
! 	tacr(k,mt,kc)=tacr(k,mt,kc)+a8
! 	raut(k,mt,kc)=raut(k,mt,kc)+a9
! 	racw(k,mt,kc)=racw(k,mt,kc)+a10
! 	sfw(k,mt,kc)=sfw(k,mt,kc)+a11
! 	sfi(k,mt,kc)=sfi(k,mt,kc)+a12
! 	gacs(k,mt,kc)=gacs(k,mt,kc)+a13
! 	gacw(k,mt,kc)=gacw(k,mt,kc)+a14
! 	gaci(k,mt,kc)=gaci(k,mt,kc)+a15
! 	gacr(k,mt,kc)=gacr(k,mt,kc)+a16
! 	gwet(k,mt,kc)=gwet(k,mt,kc)+a17
! 	gaut(k,mt,kc)=gaut(k,mt,kc)+a18
! 	racs(k,mt,kc)=racs(k,mt,kc)+a19
! 	sacr(k,mt,kc)=sacr(k,mt,kc)+a20
! 	gfr(k,mt,kc)=gfr(k,mt,kc)+a21
! 	smlt(k,mt,kc)=smlt(k,mt,kc)+a22
! 	gmlt(k,mt,kc)=gmlt(k,mt,kc)+a23
! 	sdep(k,mt,kc)=sdep(k,mt,kc)+a24
! 	ssub(k,mt,kc)=ssub(k,mt,kc)+a25
! 	gsub(k,mt,kc)=gsub(k,mt,kc)+a26
! 	pern(k,mt,kc)=pern(k,mt,kc)+a27
! 	d3ri(k,mt,kc)=d3ri(k,mt,kc)+a28
! 	d3ir(k,mt,kc)=d3ir(k,mt,kc)+a29
! 	d2sr(k,mt,kc)=d2sr(k,mt,kc)+a30
! 	d2rs(k,mt,kc)=d2rs(k,mt,kc)+a31
! 	gdry(k,mt,kc)=gdry(k,mt,kc)+a32
! 	erns(k,mt,kc)=erns(k,mt,kc)+a38
! 	wgrs(k,mt,kc)=wgrs(k,mt,kc)+a39
! 	qsws(k,mt,kc)=qsws(k,mt,kc)+a40
! 	srsw(k,mt,kc)=srsw(k,mt,kc)+a41
! 	srlw(k,mt,kc)=srlw(k,mt,kc)+a42
! 	q1t(k,mt,kc)=q1t(k,mt,kc)+a47
! 	qc0(k,mt,kc)=a33
! 	qr0(k,mt,kc)=a34
! 	qi0(k,mt,kc)=a35
! 	qs0(k,mt,kc)=a36
! 	qg0(k,mt,kc)=a37
! 	qv0(k,mt,kc)=a43
! 	tt0(k,mt,kc)=a44
! 	sgpt(k,mt,kc)=a45
! 	tsqq(k,mt,kc)=a46
! 	sqhdt(k,mt,kc)=sqhdt(k,mt,kc)+a48
! 	sqvdt(k,mt,kc)=sqvdt(k,mt,kc)+a49
! 	sqtdt(k,mt,kc)=sqtdt(k,mt,kc)+a50
!         pim(K,MT,KC)=pim(K,MT,KC)+a54
!         cfr(K,MT,KC)=cfr(K,MT,KC)+a55
! 	sqc0(k,mt,kc)=sqc0(k,mt,kc)+qc0(k,mt,kc)
! 	sqr0(k,mt,kc)=sqr0(k,mt,kc)+qr0(k,mt,kc)
! 	sqi0(k,mt,kc)=sqi0(k,mt,kc)+qi0(k,mt,kc)
! 	sqs0(k,mt,kc)=sqs0(k,mt,kc)+qs0(k,mt,kc)
! 	sqg0(k,mt,kc)=sqg0(k,mt,kc)+qg0(k,mt,kc)
! 	sqv0(k,mt,kc)=sqv0(k,mt,kc)+qv0(k,mt,kc)
! 	stt0(k,mt,kc)=stt0(k,mt,kc)+tt0(k,mt,kc)
! 	ssgpt(k,mt,kc)=ssgpt(k,mt,kc)+sgpt(k,mt,kc)
! 	tsqq1(k,mt,kc)=tsqq1(k,mt,kc)+tsqq(k,mt,kc)
!   350   continue

       endif

cc    ********************************************
c D. Posselt Eliminate following statistics section
!        if(id.eq.1) then
c--------------------------------------------------------------------
c    condensation:  cnd(ij)
c    evaporation:   dd(ij)+ern(ij)
c    deposition:    dep(ij)+psdep(ij)+pgdep(ij)+pint(ij)
c    sublimation:   dd1(ij)+pmlts(ij)+pmltg(ij)
c    melting:       psmlt(ij)+pgmlt(ij)+pimlt(ij)+qracs(ij)
c    freezing:      pihom(ij)+pidw(ij)+psacw(ij)+psfw(ij)+dgacw(ij)
c                   +piacr(ij)+dgacr(ij)+psacr(ij)+pgfr(ij)
c    mass flux:     dm(ij)
c    cloud water:   qc(ij)
c    rain:
c    cloud ice
c    snow
c    hail/graupel:
c----------------------------------------------------------------------
! c$doacross local(j,i,a1,a2,a3,a11,a22,a33,zdry,a44,zwet)
!        do 42 j=1,1
!        do 42 i=1,1
! 	cnd(i,j)=cnd(i,j)*rft
! 	ern(i,j)=(ern(i,j)+dd(i,j))*rft
! 	y1(i,j)=(dep(i,j)+psdep(i,j)+pgdep(i,j)+pint(i,j))*rft
! 	y2(i,j)=(dd1(i,j)+pmlts(i,j)+pmltg(i,j))*rft
! 	y3(i,j)=(psmlt(i,j)+pgmlt(i,j)+pimlt(i,j)+qracs(i,j))*rft
! 	y4(i,j)=(pihom(i,j)+pidw(i,j)+psacw(i,j)+psfw(i,j)
!      1	+dgacw(i,j)+piacr(i,j)+dgacr(i,j)+psacr(i,j)
!      2	+pgfr(i,j))*rft
! 	y5(i,j)=dm(i,j)
! 	a1=1.e6*r00*qr(i,j)
! 	a2=1.e6*r00*qs(i,j)
! 	a3=1.e6*r00*qg(i,j)
! 	a11=ucor*(max(1.e-5,a1))**1.75
!         tns_funT(i,j)=1.
!         tairc(i,j)=min(0.,max(-50.,tair(i,j)-t0))                   !2007b
!         if(improve.ge.3) tns_funT(i,j)=exp(0.060000*tairc(i,j)*0.75)  !2007b
! 	a22=ucos*(max(1.e-5,a2))**1.75 *tns_funT(i,j)               !2007b
! 	a33=ucog*(max(1.e-5,a3))**1.75
! 	zdry=max(1.,a11+a22+a33)
! 	dbz(i,j)=10.*log10(zdry)
! 	if (tair(i,j).ge.273.16) then
! 	  a44=a11+uwet*(a22+a33)**.95
! 	  zwet=max(1.,a44)
! 	  dbz(i,j)=10.*log10(zwet)
! 	endif
! 	dbz(i,j)=dbz(i,j)*asss(i,j)
! 	qc(i,j)=qc(i,j)*asss(i,j)
! 	qr(i,j)=qr(i,j)*asss(i,j)
! 	qi(i,j)=qi(i,j)*asss(i,j)
! 	qs(i,j)=qs(i,j)*asss(i,j)
! 	qg(i,j)=qg(i,j)*asss(i,j)
! c store microphysical process rates
! c            if (asss(i,j).ne.0.) THEN
!               physc(i,j,k)=cnd(i,j)*asss(i,j)
!               physe(i,j,k)=ern(i,j)*asss(i,j)
!               physd(i,j,k)=y1(i,j)*asss(i,j)
!               physs(i,j,k)=y2(i,j)*asss(i,j)
!               physm(i,j,k)=y3(i,j)*asss(i,j)
!               physf(i,j,k)=y4(i,j)*asss(i,j)
! c            endif
! 
! 
!    42  continue
! 
!        do j=1,1
!        do i=1,1
!         tmpd(i,j,1)=0.0
!         tmpd(i,j,2)=0.0
! 	if(rq(i,j) .ge. 0.) tmpd(i,j,1)=cnd(i,j)
! 	if(rq(i,j) .lt. 0.) tmpd(i,j,2)=ern(i,j)
!        enddo
!        enddo
!        sumd(1)=0.0
!        sumd(2)=0.0
! 
!        call RealSum3D(tmpd,sumd,2)
! 
!        scu1(k)=scu1(k)+sumd(1)
!        sed1(k)=sed1(k)+sumd(2)
! 
!        do 40 j=1,1
!        do 40 i=1,1
! 
!         IF(dbz(i,j).GE.0.) IZZ=int(dbz(i,j)/2.+0.5)
!         IF(dbz(i,j).LT.0.) IZZ=int(dbz(i,j)/2.-0.5)
! 
!         WW=0.01*(ww1(i,j,k)+ww1(i,j,k+1))/2.     ! W at theta levels
!         IF(WW.GE.0.) IWW=int(WW+0.5)
!         IF(WW.LT.0.) IWW=int(WW-0.5)
! 
! 	 if(iww.ge.-20 .and. iww.le.40) then
!             iw1=iww+21
! 	    if (ics5(i,j,2) .eq. 1) then
! 	      cf1(iw1,k)=cf1(iw1,k)+cnd(i,j)
! 	      cf2(iw1,k)=cf2(iw1,k)+ern(i,j)
! 	      cf3(iw1,k)=cf3(iw1,k)+y1(i,j)
! 	      cf4(iw1,k)=cf4(iw1,k)+y2(i,j)
! 	      cf5(iw1,k)=cf5(iw1,k)+y3(i,j)
! 	      cf6(iw1,k)=cf6(iw1,k)+y4(i,j)
! 	      cf7(iw1,k)=cf7(iw1,k)+y5(i,j)
! 	      cf8(iw1,k)=cf8(iw1,k)+qc(i,j)
! 	      cf9(iw1,k)=cf9(iw1,k)+qr(i,j)
! 	      cf10(iw1,k)=cf10(iw1,k)+qi(i,j)
! 	      cf11(iw1,k)=cf11(iw1,k)+qs(i,j)
! 	      cf12(iw1,k)=cf12(iw1,k)+qg(i,j)
! 	      cfnum(iw1,k)=cfnum(iw1,k)+1
! 	    elseif (ics5(i,j,3) .eq. 1) then
! 	      cfs1(iw1,k)=cfs1(iw1,k)+cnd(i,j)
! 	      cfs2(iw1,k)=cfs2(iw1,k)+ern(i,j)
! 	      cfs3(iw1,k)=cfs3(iw1,k)+y1(i,j)
! 	      cfs4(iw1,k)=cfs4(iw1,k)+y2(i,j)
! 	      cfs5(iw1,k)=cfs5(iw1,k)+y3(i,j)
! 	      cfs6(iw1,k)=cfs6(iw1,k)+y4(i,j)
! 	      cfs7(iw1,k)=cfs7(iw1,k)+y5(i,j)
! 	      cfs8(iw1,k)=cfs8(iw1,k)+qc(i,j)
! 	      cfs9(iw1,k)=cfs9(iw1,k)+qr(i,j)
! 	      cfs10(iw1,k)=cfs10(iw1,k)+qi(i,j)
! 	      cfs11(iw1,k)=cfs11(iw1,k)+qs(i,j)
! 	      cfs12(iw1,k)=cfs12(iw1,k)+qg(i,j)
! 	      cfsnum(iw1,k)=cfsnum(iw1,k)+1
! 	    endif
!          endif
! 
! 	 if(izz.ge.-5 .and. izz.le.35) then   ! restrict to -10 dbZ echo region
!             izz=izz+6
! 	      cfz(iw1,k)=cfz(iw1,k)+1.
!        	    if(iww.ge.-20 .and. iww.le.40) then
!             iw1=iww+21
! 	      cfw(iw1,k)=cfw(iw1,k)+1.
!             endif
!          endif
! 
!    40  continue
! 
!        endif

c DP Fill output arrays
        theta_out(1,1,k) = tb0
        tair_out(1,1,k) = tair(1,1)

 1000 continue

c     ****************************************************************
c D. Posselt Eliminate following statistics section
!       if (id .eq. 1) then
!          do j=1,1
!           do i=1,1
! 	    tmpd(i,j,1)=ics(i,j,2)
! 	    tmpd(i,j,2)=ics(i,j,3)
!           enddo
!          enddo
!          sumd(1)=0.0
!          sumd(2)=0.0
! 
!          call RealSum3D(tmpd,sumd,2)
! 
!          lc=sumd(1)
!          ls=sumd(2)
! 
! 	if (lc .eq. 0) lc=1000000
! 	if (ls .eq. 0) ls=1000000
!        do 400 mt=1,3
! 	a1=rijl2
! 	 if (mt .eq. 2) a1=1./float(lc)
! 	 if (mt .eq. 3) a1=1./float(ls)
! 
!        do k=kt1,kt2
!          do j=1,1
!          do i=1,1
! 	      if(ics5(i,j,mt).eq.1) then
! 	       tmpd(i,j,1)=dpt(i,j,k)
! 	       tmpd(i,j,2)=dqv(i,j,k)
!               else
! 	       tmpd(i,j,1)=0.0
! 	       tmpd(i,j,2)=0.0
! 	      endif
!          enddo
!          enddo
!          sumd(1) = 0.0
!          sumd(2) = 0.0
! 
!          call RealSum3D(tmpd,sumd,2)
! 
!          b1(k)=sumd(1)
!          b2(k)=sumd(2)
!        enddo
! 
!        do 430 k=kt1,kt2
! 	    tb00(k,mt)=b1(k)*a1
! 	    qb00(k,mt)=b2(k)*a1
!   430  continue
!   400 continue
!       endif

      d2t=d22t

c      call satdt

      deallocate(tmpd)
      deallocate(sumd)

      return
      end
