***  For comments, see sbart.r, written by Finbarr O'Sullivan
***  ratfor output of 14 May 1987 version
***  with bug fix in setreg, 27 Sep 1990
      subroutine sbart(xs,ys,ws,n,knot,nk,coef,sz,lev,crit,icrit,spar,
     &ispar,lspar,uspar,tol,isetup,xwy,hs0,hs1,hs2,hs3,sg0,sg1,sg2,sg3,
     &abd,p1ip,p2ip,ld4,ldnk,ier)
      real xs(n),ys(n),ws(n),knot(nk+4),coef(nk),sz(n),lev(n),crit,spar,
     &lspar,uspar,tol,xwy(nk),hs0(nk),hs1(nk),hs2(nk),hs3(nk), sg0(nk),
     &sg1(nk),sg2(nk),sg3(nk),abd(ld4,nk),p1ip(ld4,nk),p2ip(ldnk,nk)
      integer n,nk,isetup,icrit,ispar,ld4,ldnk,ier
      realt1,t2,ratio, a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w, fu,fv,fw,
     &fx,x,ax,bx
      integer i
      if(.not.(isetup.eq.0))goto 23000
      call sgram(sg0,sg1,sg2,sg3,knot,nk)
      call stxwx(xs,ys,ws,n,knot,nk,xwy,hs0,hs1,hs2,hs3)
      t1=0. 
      t2=0.
      do 23002 i=3,nk-3 
      t1 = t1 + hs0(i) 
23002 continue
      do 23004 i=3,nk-3 
      t2 = t2 + sg0(i) 
23004 continue
      ratio = t1/t2
      isetup = 1 
23000 continue
      if(.not.(ispar.eq.1))goto 23006
      call sslvrg(xs,ys,ws,n,knot,nk,coef,sz,lev,crit,icrit,spar,ratio,
     &xwy,hs0,hs1,hs2,hs3,sg0,sg1,sg2,sg3,abd,p1ip,p2ip,ld4,ldnk,ier)
      return
23006 continue
      ax=lspar 
      bx=uspar
      c = 0.5*(3. - sqrt(5.0))
      eps = 1.00
10    eps = eps/2.00
      tol1 = 1.0 + eps
      if(.not.(tol1 .gt. 1.00))goto 23008
      go to 10
23008 continue
      eps = sqrt(eps)
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0
      spar = x
      call sslvrg(xs,ys,ws,n,knot,nk,coef,sz,lev,crit,icrit,spar,ratio,
     &xwy,hs0,hs1,hs2,hs3,sg0,sg1,sg2,sg3,abd,p1ip,p2ip,ld4,ldnk,ier)
      fx = crit
      fv = fx
      fw = fx
20    xm = 0.5*(a + b)
      tol1 = eps*abs(x) + tol/3.0
      tol2 = 2.0*tol1
      if(.not.(abs(x - xm) .le. (tol2 - 0.5*(b - a))))goto 23010
      go to 90
23010 continue
      if(.not.(abs(e) .le. tol1))goto 23012
      go to 40
23012 continue
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.00*(q - r)
      if(.not.(q .gt. 0.0))goto 23014
      p = -p
23014 continue
      q = abs(q)
      r = e
      e = d
30    if(.not.(abs(p) .ge. abs(0.5*q*r)))goto 23016
      go to 40
23016 continue
      if(.not.(p .le. q*(a - x)))goto 23018
      go to 40
23018 continue
      if(.not.(p .ge. q*(b - x)))goto 23020
      go to 40
23020 continue
      d = p/q
      u = x + d
      if(.not.((u - a) .lt. tol2))goto 23022
      d = sign(tol1, xm - x)
23022 continue
      if(.not.((b - u) .lt. tol2))goto 23024
      d = sign(tol1, xm - x)
23024 continue
      go to 50
40    if(.not.(x .ge. xm))goto 23026
      e = a - x
23026 continue
      if(.not.(x .lt. xm))goto 23028
      e = b - x
23028 continue
      d = c*e
50    if(.not.(abs(d) .ge. tol1))goto 23030
      u = x + d
23030 continue
      if(.not.(abs(d) .lt. tol1))goto 23032
      u = x + sign(tol1, d)
23032 continue
      spar = u
      call sslvrg(xs,ys,ws,n,knot,nk,coef,sz,lev,crit,icrit,spar,ratio,
     &xwy,hs0,hs1,hs2,hs3,sg0,sg1,sg2,sg3,abd,p1ip,p2ip,ld4,ldnk,ier)
      fu = crit
      if(.not.(fu .gt. fx))goto 23034
      go to 60
23034 continue
      if(.not.(u .ge. x))goto 23036
      a = x
23036 continue
      if(.not.(u .lt. x))goto 23038
      b = x
23038 continue
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
60    if(.not.(u .lt. x))goto 23040
      a = u
23040 continue
      if(.not.(u .ge. x))goto 23042
      b = u
23042 continue
      if(.not.(fu .le. fw))goto 23044
      go to 70
23044 continue
      if(.not.(w .eq. x))goto 23046
      go to 70
23046 continue
      if(.not.(fu .le. fv))goto 23048
      go to 80
23048 continue
      if(.not.(v .eq. x))goto 23050
      go to 80
23050 continue
      if(.not.(v .eq. w))goto 23052
      go to 80
23052 continue
      go to 20
70    v = w
      fv = fw
      w = u
      fw = fu
      go to 20
80    v = u
      fv = fu
      go to 20
90    continue 
      spar = x 
      crit = fx
      return
23007 continue
      return
      end
      subroutine sgram(sg0,sg1,sg2,sg3,tb,nb)
      real sg0(nb),sg1(nb),sg2(nb),sg3(nb),tb(nb+4),vnikx(4,3),work(16),
     &yw1(4),yw2(4),wpt
      integer nb,ileft,ilo,mflag,i,ii,jj
      do 23054 i=1,nb
      sg0(i)=0.
      sg1(i)=0.
      sg2(i)=0.
      sg3(i)=0.
23054 continue
      ilo = 1
      do 23056 i=1,nb 
      call interv(tb(1),(nb+1),tb(i),ileft,mflag)
      call bsplvd (tb,4,tb(i),ileft,work,vnikx,3)
      do 23058 ii=1,4 
      yw1(ii) = vnikx(ii,3) 
23058 continue
      call bsplvd (tb,4,tb(i+1),ileft,work,vnikx,3)
      do 23060 ii=1,4 
      yw2(ii) = vnikx(ii,3) - yw1(ii) 
23060 continue
      wpt = tb(i+1) - tb(i)
      if(.not.(ileft.ge.4))goto 23062
      do 23064 ii=1,4 
      jj=ii
      sg0(ileft-4+ii) = sg0(ileft-4+ii) +wpt* (yw1(ii)*yw1(jj) + (yw2(
     &ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50 +yw2(ii)*yw2(jj)*.3330 )
      jj=ii+1
      if(.not.(jj.le.4))goto 23066
      sg1(ileft+ii-4) = sg1(ileft+ii-4) +wpt* (yw1(ii)*yw1(jj) + (yw2(
     &ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50 +yw2(ii)*yw2(jj)*.3330 )
23066 continue
      jj=ii+2
      if(.not.(jj.le.4))goto 23068
      sg2(ileft+ii-4) = sg2(ileft+ii-4) +wpt* (yw1(ii)*yw1(jj) + (yw2(
     &ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50 +yw2(ii)*yw2(jj)*.3330 )
23068 continue
      jj=ii+3
      if(.not.(jj.le.4))goto 23070
      sg3(ileft+ii-4) = sg3(ileft+ii-4) +wpt* (yw1(ii)*yw1(jj) + (yw2(
     &ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50 +yw2(ii)*yw2(jj)*.3330 )
23070 continue
23064 continue
      goto 23063
23062 continue
      if(.not.(ileft.eq.3))goto 23072
      do 23074 ii=1,3 
      jj=ii
      sg0(ileft-3+ii) = sg0(ileft-3+ii) +wpt* (yw1(ii)*yw1(jj) + (yw2(
     &ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50 +yw2(ii)*yw2(jj)*.3330 )
      jj=ii+1
      if(.not.(jj.le.3))goto 23076
      sg1(ileft+ii-3) = sg1(ileft+ii-3) +wpt* (yw1(ii)*yw1(jj) + (yw2(
     &ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50 +yw2(ii)*yw2(jj)*.3330 )
23076 continue
      jj=ii+2
      if(.not.(jj.le.3))goto 23078
      sg2(ileft+ii-3) = sg2(ileft+ii-3) +wpt* (yw1(ii)*yw1(jj) + (yw2(
     &ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50 +yw2(ii)*yw2(jj)*.3330 )
23078 continue
23074 continue
      goto 23073
23072 continue
      if(.not.(ileft.eq.2))goto 23080
      do 23082 ii=1,2 
      jj=ii
      sg0(ileft-2+ii) = sg0(ileft-2+ii) +wpt* (yw1(ii)*yw1(jj) + (yw2(
     &ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50 +yw2(ii)*yw2(jj)*.3330 )
      jj=ii+1
      if(.not.(jj.le.2))goto 23084
      sg1(ileft+ii-2) = sg1(ileft+ii-2) +wpt* (yw1(ii)*yw1(jj) + (yw2(
     &ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50 +yw2(ii)*yw2(jj)*.3330 )
23084 continue
23082 continue
      goto 23081
23080 continue
      if(.not.(ileft.eq.1))goto 23086
      do 23088 ii=1,1 
      jj=ii
      sg0(ileft-1+ii) = sg0(ileft-1+ii) +wpt* (yw1(ii)*yw1(jj) + (yw2(
     &ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50 +yw2(ii)*yw2(jj)*.3330 )
23088 continue
23086 continue
23081 continue
23073 continue
23063 continue
23056 continue
      return
      end
      subroutine sinerp(abd,ld4,nk,p1ip,p2ip,ldnk,flag)
      realabd(ld4,nk),p1ip(ld4,nk),p2ip(ldnk,nk),wjm3(3),wjm2(2),wjm1(1)
     &,c0,c1,c2,c3
      integerflag,ld4,nk,ldnk,i,j,k
      wjm3(1)=0. 
      wjm3(2)=0. 
      wjm3(1)=0.
      wjm2(1)=0. 
      wjm2(2)=0.
      wjm1(1)=0.
      do 23090 i=1,nk 
      j=nk-i+1
      c0 = 1./abd(4,j)
      if(.not.(j.le.nk-3))goto 23092
      c1 = abd(1,j+3)*c0
      c2 = abd(2,j+2)*c0
      c3 = abd(3,j+1)*c0 
      goto 23093
23092 continue
      if(.not.(j.eq.nk-2))goto 23094
      c1 = 0.
      c2 = abd(2,j+2)*c0
      c3 = abd(3,j+1)*c0 
      goto 23095
23094 continue
      if(.not.(j.eq.nk-1))goto 23096
      c1 = 0.
      c2 = 0.
      c3 = abd(3,j+1)*c0 
      goto 23097
23096 continue
      if(.not.(j.eq.nk))goto 23098
      c1 = 0.
      c2 = 0.
      c3 = 0.
23098 continue
23097 continue
23095 continue
23093 continue
      p1ip(1,j) = 0. - (c1*wjm3(1)+c2*wjm3(2)+c3*wjm3(3))
      p1ip(2,j) = 0. - (c1*wjm3(2)+c2*wjm2(1)+c3*wjm2(2))
      p1ip(3,j) = 0. - (c1*wjm3(3)+c2*wjm2(2)+c3*wjm1(1))
      p1ip(4,j) = c0**2 +c1**2*wjm3(1)+2.*c1*c2*wjm3(2)+2.*c1*c3*wjm3(3)
     & +c2**2*wjm2(1)+2.*c2*c3*wjm2(2) +c3**2*wjm1(1)
      wjm3(1)=wjm2(1) 
      wjm3(2)=wjm2(2) 
      wjm3(3)=p1ip(2,j)
      wjm2(1)=wjm1(1) 
      wjm2(2)=p1ip(3,j)
      wjm1(1)=p1ip(4,j)
23090 continue
      if(.not.(flag.eq.0))goto 23100
      return
23100 continue
      do 23102 i=1,nk 
      j=nk-i+1
      k=1
23104 if(.not.(k.le.4.and.j+k-1.le.nk))goto 23106
      p2ip(j,j+k-1) = p1ip(5-k,j) 
      k=k+1
      goto 23104
23106 continue
23102 continue
      do 23107 i=1,nk 
      j=nk-i+1
      k=j-4
23109 if(.not.(k.ge.1))goto 23111
      c0 = 1./abd(4,k) 
      c1 = abd(1,k+3)*c0
      c2 = abd(2,k+2)*c0 
      c3 = abd(3,k+1)*c0
      p2ip(k,j) = 0. - ( c1*p2ip(k+3,j) +c2*p2ip(k+2,j) +c3*p2ip(k+1,j) 
     &) 
      k=k-1
      goto 23109
23111 continue
23107 continue
      return
23101 continue
      end
      subroutine sslvrg(x,y,w,n,knot,nk,coef,sz,lev,crit,icrit,spar,
     &ratio,xwy,hs0,hs1,hs2,hs3,sg0,sg1,sg2,sg3,abd,p1ip,p2ip,ld4,ldnk,
     &info)
      real x(n),y(n),w(n),knot(nk+4),coef(nk),sz(n),lev(n),crit,ratio,
     &spar,xwy(nk),hs0(nk),hs1(nk),hs2(nk),hs3(nk), sg0(nk),sg1(nk),sg2(
     &nk),sg3(nk),abd(ld4,nk),p1ip(ld4,nk),p2ip(ldnk,nk),lambda,b0,b1,
     &b2,b3,eps,vnikx(4,1),work(16),xv,bvalu,rss,df
!     &b2,b3,eps,vnikx(4,1),work(16),xv,bvalu2,rss,df
      integer n,nk,icrit,ld4,ldnk,i,icoef,ileft,ilo,info,j,mflag
      ilo = 1 
      eps = .1e-10
      lambda = ratio*16.**(-2. + spar*(6.))
      do 23112 i=1,nk 
      coef(i) = xwy(i) 
23112 continue
      do 23114 i=1,nk 
      abd(4,i) = hs0(i)+lambda*sg0(i) 
23114 continue
      do 23116 i=1,(nk-1) 
      abd(3,i+1) = hs1(i)+lambda*sg1(i) 
23116 continue
      do 23118 i=1,(nk-2) 
      abd(2,i+2) = hs2(i)+lambda*sg2(i) 
23118 continue
      do 23120 i=1,(nk-3) 
      abd(1,i+3) = hs3(i)+lambda*sg3(i) 
23120 continue
      call spbfa(abd,ld4,nk,3,info)
      if(.not.(info.ne.0))goto 23122
      return
23122 continue
      call spbsl(abd,ld4,nk,3,coef)
      icoef = 1
      do 23124 i=1,n 
      xv = x(i)
      sz(i) = bvalu(knot,coef,nk,4,xv,0)
!      sz(i) = bvalu2(knot,coef,nk,4,xv,0)
23124 continue
      if(.not.(icrit.eq.0))goto 23126
      return
23126 continue
      call sinerp(abd,ld4,nk,p1ip,p2ip,ldnk,0)
      do 23128 i=1,n 
      xv = x(i)
      call interv(knot(1),(nk+1),xv,ileft,mflag)
      if(.not.(mflag.eq.-1))goto 23130
      ileft = 4 
      xv = knot(4)+eps 
23130 continue
      if(.not.(mflag.eq.1))goto 23132
      ileft = nk 
      xv = knot(nk+1)-eps 
23132 continue
      j=ileft-3
      call bsplvd(knot,4,xv,ileft,work,vnikx,1)
      b0=vnikx(1,1)
      b1=vnikx(2,1)
      b2=vnikx(3,1)
      b3=vnikx(4,1)
      lev(i) = (p1ip(4,j)*b0**2 + 2.*p1ip(3,j)*b0*b1 +2.*p1ip(2,j)*b0*
     &b2 + 2.*p1ip(1,j)*b0*b3 +p1ip(4,j+1)*b1**2 + 2.*p1ip(3,j+1)*b1*b2 
     &+2.*p1ip(2,j+1)*b1*b3 +p1ip(4,j+2)*b2**2 + 2.*p1ip(3,j+2)*b2*b3 +
     &p1ip(4,j+3)*b3**2 )*w(i)**2 
23128 continue
      if(.not.(icrit.eq.1))goto 23134
      rss = 0. 
      df = 0.
      do 23136 i=1,n 
      rss = rss + ((y(i)-sz(i))*w(i))**2
23136 continue
      do 23138 i=1,n 
      df = df + 1.-lev(i)
23138 continue
      crit = (rss/n)/((df/n)**2) 
      goto 23135
23134 continue
      crit = 0.
      do 23140 i=1,n 
      crit = crit +(((y(i)-sz(i))*w(i))/(1-lev(i)))**2 
23140 continue
23135 continue
      return
23127 continue
      end
      subroutine stxwx(x,z,w,k,xknot,n,y,hs0,hs1,hs2,hs3)
      real z(k),w(k),x(k),xknot(n+4),y(n),hs0(n),hs1(n),hs2(n),hs3(n),
     &eps,vnikx(4,1),work(16)
      integer k,n,j,i,ilo,ileft,mflag
      do 23142 i=1,n 
      y(i)=0. 
      hs0(i)=0. 
      hs1(i)=0.
      hs2(i)=0. 
      hs3(i)=0. 
23142 continue
      ilo=1 
      eps = .1e-9
      do 23144 i=1,k 
      call interv(xknot(1),(n+1),x(i),ileft,mflag)
      if(.not.(mflag.eq.-1))goto 23146
      write(6,'("Error in hess ",i2)')mflag
      stop
23146 continue
      if(.not.(mflag.eq. 1))goto 23148
      if(.not.(x(i).le.(xknot(ileft)+eps)))goto 23150
      ileft=ileft-1
      goto 23151
23150 continue
      write(6,'("Error in hess ",i2)')mflag
      stop
23151 continue
23148 continue
      call bsplvd (xknot,4,x(i),ileft,work,vnikx,1)
      j= ileft-4+1
      y(j) = y(j)+w(i)**2*z(i)*vnikx(1,1)
      hs0(j)=hs0(j)+w(i)**2*vnikx(1,1)**2
      hs1(j)=hs1(j)+w(i)**2*vnikx(1,1)*vnikx(2,1)
      hs2(j)=hs2(j)+w(i)**2*vnikx(1,1)*vnikx(3,1)
      hs3(j)=hs3(j)+w(i)**2*vnikx(1,1)*vnikx(4,1)
      j= ileft-4+2
      y(j) = y(j)+w(i)**2*z(i)*vnikx(2,1)
      hs0(j)=hs0(j)+w(i)**2*vnikx(2,1)**2
      hs1(j)=hs1(j)+w(i)**2*vnikx(2,1)*vnikx(3,1)
      hs2(j)=hs2(j)+w(i)**2*vnikx(2,1)*vnikx(4,1)
      j= ileft-4+3
      y(j) = y(j)+w(i)**2*z(i)*vnikx(3,1)
      hs0(j)=hs0(j)+w(i)**2*vnikx(3,1)**2
      hs1(j)=hs1(j)+w(i)**2*vnikx(3,1)*vnikx(4,1)
      j= ileft-4+4
      y(j) = y(j)+w(i)**2*z(i)*vnikx(4,1)
      hs0(j)=hs0(j)+w(i)**2*vnikx(4,1)**2 
23144 continue
      return
      end
      subroutine sknotl(x,n,knot,k)
      real x(n),knot(n+6),a1,a2,a3,a4
      integer n,k,ndk,j
      a1 = log(50.)/log(2.) 
      a2 = log(100.)/log(2.)
      a3 = log(140.)/log(2.) 
      a4 = log(200.)/log(2.)
      if(.not.(n.lt.50))goto 23152
      ndk = n 
      goto 23153
23152 continue
      if(.not.(n.ge.50 .and. n.lt.200))goto 23154
      ndk = 2.**(a1+(a2-a1)*(n-50.)/150.) 
      goto 23155
23154 continue
      if(.not.(n.ge.200 .and. n.lt.800))goto 23156
      ndk = 2.**(a2+(a3-a2)*(n-200.)/600.) 
      goto 23157
23156 continue
      if(.not.(n.ge.800 .and. n.lt.3200))goto 23158
      ndk = 2.**(a3+(a4-a3)*(n-800.)/2400.) 
      goto 23159
23158 continue
      if(.not.(n.ge.3200))goto 23160
      ndk = 200. + (n-3200)**.2 
23160 continue
23159 continue
23157 continue
23155 continue
23153 continue
      k = ndk + 6
      do 23162 j=1,3 
      knot(j) = x(1) 
23162 continue
      do 23164 j=1,ndk 
      knot(j+3) = x( 1 + (j-1)*(n-1)/(ndk-1) ) 
23164 continue
      do 23166 j=1,3 
      knot(ndk+3+j) = x(n) 
23166 continue
      return
      end
      subroutine setreg(x,y,w,n,xw,nx,min,range,knot,nk)
      real x(n),y(n),w(n),xw(n),min,range,knot(n+6)
      integer n,nk,nx
      real eps
      integer ycnt,i,k
      call scopy(n,x,1,xw,1)
!      call sortpr(x,n,w)
!      call sortpr(xw,n,y)
      range = x(n)-x(1) 
      min = x(1) 
      eps = .1e-9
      do 23168 i=1,n 
      x(i) = (x(i)-min)/range 
23168 continue
      call scopy(n,x,1,xw,1)
      nx = 1 
      x(nx) = x(1) 
      w(nx) = w(1) 
      y(nx) = y(1)*w(1) 
      i=2
23170 if(.not.(i.le.n))goto 23172
      if(.not.(xw(i).gt.x(nx)+eps))goto 23173
      if(.not.(w(nx).gt.0.0))goto 23175
      y(nx) = y(nx)/w(nx)
23175 continue
      nx = nx + 1
      x(nx) = x(i)
      y(nx) = y(i)*w(i)
      w(nx) = w(i) 
      goto 23174
23173 continue
      y(nx) = y(nx)+y(i)*w(i)
      w(nx) = w(nx) + w(i) 
23174 continue
      i=i+1
      goto 23170
23172 continue
      if(.not.(w(nx).gt.0.0))goto 23177
      y(nx) = y(nx)/w(nx)
23177 continue
      i=1
23179 if(.not.(i.le.nx))goto 23181
      if(.not.(w(i).gt.0))goto 23182
      w(i)=sqrt(w(i)) 
23182 continue
      i=i+1
      goto 23179
23181 continue
      call sknotl(x,nx,knot,k) 
      nk = k-4
      return
      end
