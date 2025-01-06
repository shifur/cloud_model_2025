cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trip (n,a,b,c,y,z,intr)
      dimension a(n),b(n),c(n),y(n),z(n)
      save
      ii(i)=(i-1)*intr+1
      bn=b(n)
      yn=y(n)
      v=c(n)
      y(1)=y(1)/b(1)
      a(1)=a(1)/b(1)
      b(1)=c(1)/b(1)
      nm2=n-2
      do 101 j=2,nm2
      den=b(j)-a(j)*b(j-1)
      b(j)=c(j)/den
      y(j)=(y(j)-a(j)*y(j-1))/den
      a(j)=-a(j)*a(j-1)/den
      bn=bn-v*a(j-1)
      yn=yn-v*y(j-1)
      v=-v*b(j-1)
  101 continue
      den=b(n-1)-a(n-1)*b(n-2)
      b(n-1)=(c(n-1)-a(n-1)*a(n-2))/den
      y(n-1)=(y(n-1)-a(n-1)*y(n-2))/den
      bn=bn-v*a(n-2)
      yn=yn-v*y(n-2)
      v=a(n)-v*b(n-2)
      nm1=n-1
      in=ii(n)
      inm=ii(nm1)
      z(in)=(yn-v*y(nm1))/(bn-v*b(nm1))
      z(inm)=y(nm1)-b(nm1)*z(in)
      do 102 j=2,nm1
      k=n-j
      ik=ii(k)
      z(ik)=y(k)-b(k)*z(ik+intr)-a(k)*z(in)
  102 continue
      return
      end
