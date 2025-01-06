cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine terp1 (n,x,f,w,y,intr,tab,itab)
      dimension x(n),f(n),w(n),tab(3),itab(3)
      save
      call search (y,x,n,i)
      call interp (n,x,f,w,y,i,intr,tab,itab)
      return
      end
