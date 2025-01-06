c*****************************************************************
      subroutine rflx(m,np,swc,u1,du,nu,swh,w1,dw,nw,tbl,df)
c*****************************************************************
c-----compute the reduction of clear-sky downward solar flux
c     due to co2 absorption.
      implicit none
c-----input parameters
      integer m,np,nu,nw
      real u1,du,w1,dw
      real swc(m,np+1),swh(m,np+1),tbl(nu,nw)
c-----output (undated) parameter
      real df(m,np+1)
c-----temporary array
      integer i,k,ic,iw 
      real clog,wlog,dc,dd,x0,x1,x2,y0,y1,y2
c-----table look-up for the reduction of clear-sky solar
         x0=u1+float(nu)*du
         y0=w1+float(nw)*dw
         x1=u1-0.5*du
         y1=w1-0.5*dw
      do k= 2, np+1
       do i= 1, m
          clog=min(swc(i,k),x0)
          clog=max(swc(i,k),x1)
          wlog=min(swh(i,k),y0)
          wlog=max(swh(i,k),y1)
          ic=int( (clog-x1)/du+1.)
          iw=int( (wlog-y1)/dw+1.)
          if(ic.lt.2)ic=2
          if(iw.lt.2)iw=2
          if(ic.gt.nu)ic=nu
          if(iw.gt.nw)iw=nw
          dc=clog-float(ic-2)*du-u1
          dd=wlog-float(iw-2)*dw-w1   
          x2=tbl(ic-1,iw-1)+(tbl(ic-1,iw)-tbl(ic-1,iw-1))/dw*dd
          y2=x2+(tbl(ic,iw-1)-tbl(ic-1,iw-1))/du*dc
          df(i,k)=df(i,k)+y2
       enddo      
      enddo  
      return
      end
