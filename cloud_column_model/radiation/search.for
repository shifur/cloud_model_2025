cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine search (xbar,x,n,i)
      dimension x(n)
      save
      data b/.69314718/
      if(xbar.gt.x(2)) go to 101
      i=1
      return
  101 continue
      if(xbar.lt.x(n-1)) go to 102
      i=n-1
      return
  102 continue
      m=int((log(real(n)))/b)
      i=2**m
      if(i.ge.n) i=i/2
      k=i
      nm1=n-1
  103 continue
      k=k/2
      if(xbar.ge.x(i)) go to 104
      i=i-k
      go to 103
  104 continue
      if(xbar.le.x(i+1)) return
      i=min0(i+k,nm1)
      go to 103
      end
