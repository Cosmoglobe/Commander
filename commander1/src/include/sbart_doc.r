#     This computes the univariate smoothing spline, automatically choosing
#  the smoothing parameter by minimizing generalized cross validation.  It
#  was written by Finbarr O'Sullivan using the scheme described in
#  "Comments on Dr. Silverman's Paper", J. Royal Statistical Society B
#  (1985) 47, pp.39-40.
#  This version was installed in netlib 14 Mar 1987.
#     Remember to check whether your program needs to square or sqrt the
#  weights before calling this (or any other) least squares routine.
#    The code has been slightly modified by Trevor Hastie and Eric Grosse
#  to use de Boor's spline routines from netlib/pppack.  Also, auxiliary
#  routines "setreg" (for sorting and standardizing the range of the data)
#  and "sknotl" (for filling the knot array) have been appended.


subroutine sbart(xs,ys,ws,n,knot,nk,
		  coef,sz,lev,
		  crit,icrit,spar,ispar,lspar,uspar,tol,
		  isetup,
		  xwy,
		  hs0,hs1,hs2,hs3,
		  sg0,sg1,sg2,sg3,
		  abd,p1ip,p2ip,ld4,ldnk,ier)


       # A Cubic B-spline Smoothing routine.

#
#          The algorithm minimises:
#
#      (1/n) * sum ws(i)**2 * (ys(i)-sz(i))**2 + lambda* int ( sz"(xs) )**2 dxs
#
#        lambda is a function of the spar which is assumed to be between
#        0 and 1


	      # Input

#   n               number of data points
#  ys(n)	    vector of length n containing the observations
#  ws(n)            vector containing the weights given to each data point
#  xs(n)            vector containing the ordinates of the observations
 

#  nk               number of b-spline coefficients to be estimated
#                   nk <= n+2
#  knot(nk+4)       vector of knot points defining the cubic b-spline basis.


#  spar             penalised likelihood smoothing parameter
#  ispar            indicator saying if spar is supplied or to be estimated
#  lspar, uspar     lower and upper values for spar 0.,1. are good values
#  tol              used in Golden Search routine

#  isetup           setup indicator

#  icrit            indicator saying which cross validation score
#		    is to be computed

#  ld4              the leading dimension of abd (ie ld4=4)
#  ldnk             the leading dimension of p2ip (not referenced)


                # Output

#   coef(nk)       vector of spline coefficients
#   sz(n)          vector of smoothed z-values
#   lev(n)         vector of leverages
#   crit           either ordinary of generalized CV score
#   ier            error indicator
#                  ier = 0 ___  everything fine
#                  ier = 1 ___  spar too small or too big
#                               problem in cholesky decomposition



         # Working arrays/matrix
#   xwy		      X'Wy
#   hs0,hs1,hs2,hs3   the diagonals of the X'WX matrix
#   sg0,sg1,sg2,sg3   the diagonals of the Gram matrix
#   abd(ld4,nk)       [ X'WX+lambda*SIGMA] in diagonal form
#   p1ip(ld4,nk)       inner products between columns of L inverse
#   p2ip(ldnk,nk)      all inner products between columns of L inverse
#                          L'L = [X'WX+lambdaSIGMA]  NOT REFERENCED

real     	xs(n),ys(n),ws(n),
	 	knot(nk+4),
	 	coef(nk),sz(n),lev(n),
	 	crit,spar,lspar,uspar,tol,
		xwy(nk),
	 	hs0(nk),hs1(nk),hs2(nk),hs3(nk),
         	sg0(nk),sg1(nk),sg2(nk),sg3(nk),
	 	abd(ld4,nk),p1ip(ld4,nk),p2ip(ldnk,nk)


integer         n,nk,isetup,icrit,ispar,ld4,ldnk,ier




		  # Local variables

real		  t1,t2,ratio,
                  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w,
                  fu,fv,fw,fx,x,
		  ax,bx

integer           i








     #  Compute SIGMA, X' W**2 X, X' W**2 z, trace ratio, s0, s1.

                     # SIGMA     -> sg0,sg1,sg2,sg3
                     # X' W**2 X -> hs0,hs1,hs2,hs3
     	             # X' W**2 Z -> xwy

	   if(isetup==0){
           call sgram(sg0,sg1,sg2,sg3,knot,nk)
           call stxwx(xs,ys,ws,n,knot,nk,xwy,hs0,hs1,hs2,hs3)
	# check 
	# print 999;999 format(" got through check ")


		    t1=0. ; t2=0.
		    do i=3,nk-3 { t1 = t1 + hs0(i) }
		    do i=3,nk-3 { t2 = t2 + sg0(i) }
		    ratio = t1/t2

		    isetup = 1 }

	# check 1
	# print 1999;1999 format(" got through check 1")


     # Compute estimate


     if(ispar==1) { # Value of spar supplied

		call sslvrg(xs,ys,ws,n,knot,nk,
		  		coef,sz,lev,crit,icrit,
		  		spar,ratio,
				xwy,
		  		hs0,hs1,hs2,hs3,
		  		sg0,sg1,sg2,sg3,
		  		abd,p1ip,p2ip,ld4,ldnk,ier)
	# check 2
	# print 2999;2999 format(" got through check 2")

		    return }

     else         {
     
     # Use Forsythe Malcom and Moler routine to minimise criterion
     
      ax=lspar ; bx=uspar  # f denotes the value of the criterion

#
#      an approximation  x  to the point where	f  attains a minimum  on
#  the interval  (ax,bx)  is determined.
#
#
#  input..
#
#  ax	 left endpoint of initial interval
#  bx	 right endpoint of initial interval
#  f	 function subprogram which evaluates  f(x)  for any  x
#	 in the interval  (ax,bx)
#  tol	 desired length of the interval of uncertainty of the final
#	 result ( .ge. 0.0)
#
#
#  output..
#
#  fmin  abcissa approximating the point where	f  attains a minimum
#
#
#      the method used is a combination of  golden  section  search  and
#  successive parabolic interpolation.	convergence is never much slower
#  than  that  for  a  fibonacci search.  if  f  has a continuous second
#  derivative which is positive at the minimum (which is not  at  ax  or
#  bx),  then  convergence  is	superlinear, and usually of the order of
#  about  1.324....
#      the function  f	is never evaluated at two points closer together
#  than  eps*abs(fmin) + (tol/3), where eps is	approximately the square
#  root  of  the  relative  machine  precision.   if   f   is a unimodal
#  function and the computed values of	 f   are  always  unimodal  when
#  separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates
#  the abcissa of the global minimum of  f  on the interval  ax,bx  with
#  an error less than  3*eps*abs(fmin) + tol.  if   f	is not unimodal,
#  then fmin may approximate a local, but perhaps non-global, minimum to
#  the same accuracy.
#      this function subprogram is a slightly modified	version  of  the
#  algol  60 procedure	localmin  given in richard brent, algorithms for
#  minimization without derivatives, prentice - hall, inc. (1973).
#
#
#      real  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
#      real  fu,fv,fw,fx,x
#
#  c is the squared inverse of the golden ratio
#
      c = 0.5*(3. - sqrt(5.0))
#
#  eps is approximately the square root of the relative machine
#  precision.
#
      eps = 1.00
   10 eps = eps/2.00
      tol1 = 1.0 + eps
      if (tol1 .gt. 1.00) go to 10
      eps = sqrt(eps)
#
#  initialization
#
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0

		spar = x
		call sslvrg(xs,ys,ws,n,knot,nk,
		  		coef,sz,lev,crit,icrit,
		  		spar,ratio,
				xwy,
		  		hs0,hs1,hs2,hs3,
		  		sg0,sg1,sg2,sg3,
		  		abd,p1ip,p2ip,ld4,ldnk,ier)

      fx = crit
      fv = fx
      fw = fx
#
#  main loop starts here
#
   20 xm = 0.5*(a + b)
      tol1 = eps*abs(x) + tol/3.0
      tol2 = 2.0*tol1
#
#  check stopping criterion
#
      if (abs(x - xm) .le. (tol2 - 0.5*(b - a))) go to 90
#
# is golden-section necessary
#
      if (abs(e) .le. tol1) go to 40
#
#  fit parabola
#
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.00*(q - r)
      if (q .gt. 0.0) p = -p
      q =  abs(q)
      r = e
      e = d
#
#  is parabola acceptable
#
   30 if (abs(p) .ge. abs(0.5*q*r)) go to 40
      if (p .le. q*(a - x)) go to 40
      if (p .ge. q*(b - x)) go to 40
#
#  a parabolic interpolation step
#
      d = p/q
      u = x + d
#
#  f must not be evaluated too close to ax or bx
#
      if ((u - a) .lt. tol2) d = sign(tol1, xm - x)
      if ((b - u) .lt. tol2) d = sign(tol1, xm - x)
      go to 50
#
#  a golden-section step
#
   40 if (x .ge. xm) e = a - x
      if (x .lt. xm) e = b - x
      d = c*e
#
#  f must not be evaluated too close to x
#
   50 if (abs(d) .ge. tol1) u = x + d
      if (abs(d) .lt. tol1) u = x + sign(tol1, d)

		spar = u
		call sslvrg(xs,ys,ws,n,knot,nk,
		  		coef,sz,lev,crit,icrit,
		  		spar,ratio,
				xwy,
		  		hs0,hs1,hs2,hs3,
		  		sg0,sg1,sg2,sg3,
		  		abd,p1ip,p2ip,ld4,ldnk,ier)

      fu = crit
#
#  update  a, b, v, w, and x
#
      if (fu .gt. fx) go to 60
      if (u .ge. x) a = x
      if (u .lt. x) b = x
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
   60 if (u .lt. x) a = u
      if (u .ge. x) b = u
      if (fu .le. fw) go to 70
      if (w .eq. x) go to 70
      if (fu .le. fv) go to 80
      if (v .eq. x) go to 80
      if (v .eq. w) go to 80
      go to 20
   70 v = w
      fv = fw
      w = u
      fw = fu
      go to 20
   80 v = u
      fv = fu
      go to 20
#
#  end of main loop
#
   90 continue ; spar = x ; crit = fx
      return         }





return
end
subroutine sgram(sg0,sg1,sg2,sg3,tb,nb)

real             sg0(nb),sg1(nb),sg2(nb),sg3(nb),tb(nb+4),
		 vnikx(4,3),work(16),yw1(4),yw2(4),
		 wpt

integer   nb,ileft,ilo,mflag,
	  i,ii,jj		# indices


 		#PURPOSE

# 	Calculation of the cubic B-spline smoothness prior
# 	for "usual" interior knot setup.


#	Uses BSPVD and INTRV in the CMLIB




#	sgm[0-3](nb)    Symmetric matrix
#                       whose (i,j)'th element contains the integral of
#                       B''(i,.) B''(j,.) , i=1,2 ... nb and j=i,...nb.
#                       Only the upper four diagonals are computed.


	#Initialise the sigma vectors
		do i=1,nb{ sg0(i)=0.;sg1(i)=0.;sg2(i)=0.;sg3(i)=0.}
	
		ilo = 1

	do i=1,nb {

                # Calculate a linear approximation to the
		# second derivative of the non-zero B-splines
		# over the interval [tb(i),tb(i+1)].

			#call intrv(tb(1),(nb+1),tb(i),ilo,ileft,mflag)
			call interv(tb(1),(nb+1),tb(i),ileft,mflag)


		        #Left end second derivatives
			#call bspvd (tb,4,3,tb(i),ileft,4,vnikx,work)
			call bsplvd (tb,4,tb(i),ileft,work,vnikx,3)

		        # Put values into yw1
			do ii=1,4 { yw1(ii) = vnikx(ii,3) }

	
		        #Right end second derivatives
			#call bspvd (tb,4,3,tb(i+1),ileft,4,vnikx,work)
			call bsplvd (tb,4,tb(i+1),ileft,work,vnikx,3)


		# Slope*(length of interval) in Linear Approximation to B''
			do ii=1,4 { yw2(ii) = vnikx(ii,3) - yw1(ii)  }
		        wpt = tb(i+1) - tb(i)


# Calculate Contributions to the simga vectors

if(ileft>=4){
do ii=1,4 {

jj=ii
sg0(ileft-4+ii) = sg0(ileft-4+ii) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )
jj=ii+1
if(jj<=4) {sg1(ileft+ii-4) = sg1(ileft+ii-4) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}
jj=ii+2
if(jj<=4) {sg2(ileft+ii-4) = sg2(ileft+ii-4) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}
jj=ii+3
if(jj<=4) {sg3(ileft+ii-4) = sg3(ileft+ii-4) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}

          }
	  }

else if(ileft==3){
do ii=1,3 {

jj=ii
sg0(ileft-3+ii) = sg0(ileft-3+ii) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )
jj=ii+1
if(jj<=3) {sg1(ileft+ii-3) = sg1(ileft+ii-3) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}
jj=ii+2
if(jj<=3) {sg2(ileft+ii-3) = sg2(ileft+ii-3) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}

          }
	  }

else if(ileft==2){
do ii=1,2 {

jj=ii
sg0(ileft-2+ii) = sg0(ileft-2+ii) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )
jj=ii+1
if(jj<=2) {sg1(ileft+ii-2) = sg1(ileft+ii-2) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}

          }
	  }

else if(ileft==1){
do ii=1,1 {

jj=ii
sg0(ileft-1+ii) = sg0(ileft-1+ii) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )

          }
	  }}


return
end
subroutine sinerp(abd,ld4,nk,p1ip,p2ip,ldnk,flag)

real	abd(ld4,nk),p1ip(ld4,nk),p2ip(ldnk,nk),
	wjm3(3),wjm2(2),wjm1(1),c0,c1,c2,c3

integer	flag,ld4,nk,ldnk,i,j,k




	# Purpose :  Computes Inner Products between columns of L(-1)
	#	     where L = abd is a Banded Matrix with 3 subdiagonals

	# The algorithm works in two passes:
	#
	#		Pass 1 computes (cj,ck) k=j,j-1,j-2,j-3 ,j=nk, .. 1
	#		Pass 2 computes (cj,ck) k<=j-4  (If flag == 1 ).
	#
	#		A refinement of Elden's trick is used.
	#




			# Pass 1


			wjm3(1)=0. ; wjm3(2)=0. ; wjm3(1)=0.
			wjm2(1)=0. ; wjm2(2)=0.
			wjm1(1)=0.

		do i=1,nk { j=nk-i+1

			             c0 = 1./abd(4,j)
			if(j<=nk-3) {c1 = abd(1,j+3)*c0
			             c2 = abd(2,j+2)*c0
				     c3 = abd(3,j+1)*c0 }
		   else if(j==nk-2) {c1 = 0.
			             c2 = abd(2,j+2)*c0
				     c3 = abd(3,j+1)*c0 }
		   else if(j==nk-1) {c1 = 0.
			             c2 = 0.
				     c3 = abd(3,j+1)*c0 }
		   else if(j==nk)   {c1 = 0.
			             c2 = 0.
				     c3 = 0.}


		p1ip(1,j) = 0. - (c1*wjm3(1)+c2*wjm3(2)+c3*wjm3(3))
		p1ip(2,j) = 0. - (c1*wjm3(2)+c2*wjm2(1)+c3*wjm2(2))
		p1ip(3,j) = 0. - (c1*wjm3(3)+c2*wjm2(2)+c3*wjm1(1))

		p1ip(4,j) = c0**2 + 
		c1**2*wjm3(1)+2.*c1*c2*wjm3(2)+2.*c1*c3*wjm3(3) +
		c2**2*wjm2(1)+2.*c2*c3*wjm2(2) +
		c3**2*wjm1(1)

			wjm3(1)=wjm2(1) ; wjm3(2)=wjm2(2) ; wjm3(3)=p1ip(2,j)
			wjm2(1)=wjm1(1) ; wjm2(2)=p1ip(3,j)
			wjm1(1)=p1ip(4,j)

				}


	if(flag==0) {return}


		   # Pass 2


	else	    { # Compute p2ip

			do i=1,nk { j=nk-i+1
			for(k=1;k<=4 & j+k-1<=nk;k=k+1)
			{ p2ip(j,j+k-1) = p1ip(5-k,j) }}
		

			do i=1,nk { j=nk-i+1
			for(k=j-4;k>=1;k=k-1){

			c0 = 1./abd(4,k) ; c1 = abd(1,k+3)*c0
			c2 = abd(2,k+2)*c0 ; c3 = abd(3,k+1)*c0
				    
			p2ip(k,j) = 0. - ( c1*p2ip(k+3,j)  + 
					     c2*p2ip(k+2,j)  + 
					     c3*p2ip(k+1,j)  )  }
				  }

			return

			}


end
subroutine sslvrg(x,y,w,n,knot,nk,
		  coef,sz,lev,
		  crit,icrit,
		  spar,ratio,
		  xwy,
		  hs0,hs1,hs2,hs3,
		  sg0,sg1,sg2,sg3,
		  abd,p1ip,p2ip,ld4,ldnk,info)


real     x(n),y(n),w(n),
	 knot(nk+4),
	 coef(nk),sz(n),lev(n),
	 crit,
	 ratio,spar,
	 xwy(nk),
	 hs0(nk),hs1(nk),hs2(nk),hs3(nk),
         sg0(nk),sg1(nk),sg2(nk),sg3(nk),
	 abd(ld4,nk),p1ip(ld4,nk),p2ip(ldnk,nk),
	 lambda,
	 b0,b1,b2,b3,eps,
	 vnikx(4,1),work(16),
	# xv,bvalu,rss,df
	 xv,bvalu2,rss,df



integer  n,nk,icrit,ld4,ldnk,i,icoef,ileft,ilo,info,j,mflag

         ilo = 1 ; eps = .1e-10

    # Purpose : Solves the smoothing problem and computes the
    #           criterion function (OCV or GCV).


	# The coeficients of estimated smooth

	                      lambda = ratio*16.**(-2. + spar*(6.))

          do i=1,nk     { coef(i)    = xwy(i) }
          do i=1,nk     { abd(4,i)   = hs0(i)+lambda*sg0(i) }
          do i=1,(nk-1) { abd(3,i+1) = hs1(i)+lambda*sg1(i) }
          do i=1,(nk-2) { abd(2,i+2) = hs2(i)+lambda*sg2(i) }
          do i=1,(nk-3) { abd(1,i+3) = hs3(i)+lambda*sg3(i) }


	  call spbfa(abd,ld4,nk,3,info)
	  if(info.ne.0) { return } 
	  call spbsl(abd,ld4,nk,3,coef)


	# Value of smooth at the data points

                     icoef = 1
          do i=1,n {    xv = x(i)
			 #sz(i) = bvalu(knot,coef,
			 #               nk,4,0,xv,icoef,work(1)) }
			 sz(i) = bvalu2(knot,coef,
			                nk,4,xv,0)
						 }


	# Compute the criterion function if requested

	  if(icrit==0) { return}

	  else         { # Ordinary or Generalized CV

	                 # Get Leverages First

	      call sinerp(abd,ld4,nk,p1ip,p2ip,ldnk,0)

     do i=1,n { xv = x(i)
	        #call intrv(knot(1),(nk+1),xv,ilo,ileft,mflag)
	        call interv(knot(1),(nk+1),xv,ileft,mflag)

	        if(mflag==-1) { ileft = 4   ; xv = knot(4)+eps }
	        if(mflag==1)  { ileft = nk  ; xv = knot(nk+1)-eps }
 	        j=ileft-3
                #call bspvd(knot,4,1,xv,ileft,4,vnikx,work)
                call bsplvd(knot,4,xv,ileft,work,vnikx,1)

	        b0=vnikx(1,1);b1=vnikx(2,1);b2=vnikx(3,1);b3=vnikx(4,1)

	lev(i) = (p1ip(4,j)*b0**2   + 2.*p1ip(3,j)*b0*b1 +
				   2.*p1ip(2,j)*b0*b2 + 2.*p1ip(1,j)*b0*b3 +
		  p1ip(4,j+1)*b1**2 + 2.*p1ip(3,j+1)*b1*b2 +
				   2.*p1ip(2,j+1)*b1*b3 +
		  p1ip(4,j+2)*b2**2 + 2.*p1ip(3,j+2)*b2*b3 +
		  p1ip(4,j+3)*b3**2 )*w(i)**2      }


			  


	# Evaluate Criterion

	      if(icrit==1) { # Generalized CV
	       
			                rss = 0. ; df = 0.
			     do i=1,n { rss = rss + ((y(i)-sz(i))*w(i))**2}
			     do i=1,n { df  = df  + 1.-lev(i)}
			     crit = (rss/n)/((df/n)**2)       }

	      else         { # Ordinary CV
	       
			                crit = 0.
			     do i=1,n { crit = crit +
			                (((y(i)-sz(i))*w(i))/(1-lev(i)))**2 }}


                        return }

end
subroutine stxwx(x,z,w,k,xknot,n,y,hs0,hs1,hs2,hs3)

real             z(k),w(k),x(k),
		 xknot(n+4),
		 y(n),
		 hs0(n),hs1(n),hs2(n),hs3(n),
                
		 eps,vnikx(4,1),work(16)     # local

integer          k,n,

		 j,i,ilo,ileft,mflag       # local




     # Initialise the output vectors

     do i=1,n { y(i)=0. ; hs0(i)=0. ; hs1(i)=0.
			  hs2(i)=0. ; hs3(i)=0. }


     # Compute X'WX -> hs0,hs1,hs2,hs3  and X'WZ -> y

	ilo=1  ; eps = .1e-9

     do i=1,k {

	      #call intrv(xknot(1),(n+1),x(i),ilo,ileft,mflag)
	      call interv(xknot(1),(n+1),x(i),ileft,mflag)
	# check 
	# print 999;999 format(" got through check stxwx1")


	      if(mflag==-1) {write(6,'("Error in hess ",i2)')mflag;stop}
	      if(mflag== 1) {if(x(i)<=(xknot(ileft)+eps)){ileft=ileft-1}
			   else{write(6,'("Error in hess ",i2)')mflag;stop}}


              #call bspvd (xknot,4,1,x(i),ileft,4,vnikx,work)
              call bsplvd (xknot,4,x(i),ileft,work,vnikx,1)
	# check 2
	# print 2999;2999 format(" got through check2 stxwx ")



	      j= ileft-4+1
	       y(j) =  y(j)+w(i)**2*z(i)*vnikx(1,1)
	      hs0(j)=hs0(j)+w(i)**2*vnikx(1,1)**2
	      hs1(j)=hs1(j)+w(i)**2*vnikx(1,1)*vnikx(2,1)
	      hs2(j)=hs2(j)+w(i)**2*vnikx(1,1)*vnikx(3,1)
	      hs3(j)=hs3(j)+w(i)**2*vnikx(1,1)*vnikx(4,1)

	      j= ileft-4+2
	       y(j) =  y(j)+w(i)**2*z(i)*vnikx(2,1)
	      hs0(j)=hs0(j)+w(i)**2*vnikx(2,1)**2
	      hs1(j)=hs1(j)+w(i)**2*vnikx(2,1)*vnikx(3,1)
	      hs2(j)=hs2(j)+w(i)**2*vnikx(2,1)*vnikx(4,1)

	      j= ileft-4+3
	       y(j) =  y(j)+w(i)**2*z(i)*vnikx(3,1)
	      hs0(j)=hs0(j)+w(i)**2*vnikx(3,1)**2
	      hs1(j)=hs1(j)+w(i)**2*vnikx(3,1)*vnikx(4,1)

	      j= ileft-4+4
	       y(j) =  y(j)+w(i)**2*z(i)*vnikx(4,1)
	      hs0(j)=hs0(j)+w(i)**2*vnikx(4,1)**2  }

return
end
subroutine sknotl(x,n,knot,k)

real        x(n),knot(n+6),a1,a2,a3,a4
integer     n,k,ndk,j


     # Allocate knots acording to the cutoffs given below


	# Cutoff constants

	 a1 = log(50.)/log(2.)  ; a2 = log(100.)/log(2.)
	 a3 = log(140.)/log(2.) ; a4 = log(200.)/log(2.)

	# Cutoff Criteria

        if(n<50)                    { ndk = n }
        else if (n>=50  & n<200)    { ndk = 2.**(a1+(a2-a1)*(n-50.)/150.)  }
        else if (n>=200 & n<800)    { ndk = 2.**(a2+(a3-a2)*(n-200.)/600.)  }
        else if (n>=800 & n<3200)   { ndk = 2.**(a3+(a4-a3)*(n-800.)/2400.)  }
        else if (n>=3200)           { ndk = 200. + (n-3200)**.2 }
		 
		 k = ndk + 6

     
     # Allocate Knots  ( note no account is taken of any weighting vector )

	do j=1,3   { knot(j) = x(1) }
	do j=1,ndk { knot(j+3) = x( 1 + (j-1)*(n-1)/(ndk-1) ) }
	do j=1,3   { knot(ndk+3+j) = x(n) }

return
end
subroutine setreg(x,y,w,n,xw,nx,min,range,knot,nk)

# sbart uses the square of the weights; the following rectifies that.
# Also, the data is sorted (by a routine you may need to change for
# your local machine) and standardized so that spar=1 gives a smooth fit
# with sum(leverage values) about 2-4.

real       x(n),y(n),w(n),xw(n),
	   min,range,knot(n+6)

integer    n,nk,nx


       # Local

real      eps
integer   ycnt,i,k


  		call scopy(n,x,1,xw,1)
  		#call ssort(x,w,n,2) 
  		#call ssort(xw,y,n,2)
		call sortpr(x,n,w)
		call sortpr(xw,n,y)

            	range = x(n)-x(1) ; min = x(1) ; eps = .1e-9
  		do i=1,n { x(i) = (x(i)-min)/range }
  		call scopy(n,x,1,xw,1)


  	nx = 1 ; x(nx) = x(1) ; w(nx) = w(1) ; y(nx) = y(1)*w(1) ;

  	for(i=2;i<=n;i=i+1)

      	{ if(xw(i)>x(nx)+eps)
	   	{ if(w(nx)>0.0) y(nx) = y(nx)/w(nx)
	         nx = nx + 1
	         x(nx) = x(i)
	         y(nx) = y(i)*w(i)
	         w(nx) = w(i) }
        else
	    { y(nx) = y(nx)+y(i)*w(i)
	      w(nx) = w(nx) + w(i) }
        }
	if(w(nx)>0.0) y(nx) = y(nx)/w(nx)	

	for(i=1;i<=nx;i=i+1)
		{ if (w(i)>0) w(i)=sqrt(w(i)) }


      call sknotl(x,nx,knot,k) ; nk = k-4




return
end


