CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C   Univariate Adaptive weights smoothing                        C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine locuniw(n,y,yhat,yhatnew,sk,s2hat,sknew,controls,
     1      newcontr,minsk,ndelta,lambda,eta,gamma,kern)
C
C     performs one step of univariate adaptive weights smoothing
C
C     called from R routine awsuni if in sequential mode (graph=T)
C     and from locunial otherwise
C
C  
C    n     - number of design points
C    y     - observations
C    yhat  - estimates from last step
C    yhatnew - new estimates (return)
C    sk    - old variance estimates (variance of yhat)
C    sknew - new variance estimates (variance of yhatnew) (return)
C    s2hat - variance estimate (variance of y)
C    controls - lower und upper bounds in control step (input and return)
C    newcontr - either 0 or 1 determining whether  controls is to be
C               updated or not (in dyadic control sceme, always 1 if
C                               control is activ in all steps) 
C    minsk  - minimal variances obtained
C    ndelta - radius of neighbourhood in design points
C    lambda - main smoothing parameter (should be around 3)
C    eta    - main control parameter (should be around 4)   
C    gamma  - allow for increase of variances (over minsk) by factor gamma 
C    kern   - vector of values of K(x) in seq(0,6,.3) (K(x)=exp(-x))
C
       implicit logical (a-z)
       integer n,ndelta,i,j,jmin,jmax,iz,newcontr
       real y(n),yhat(n),yhatnew(n),controls(2,n),s2hat(n),sknew(n),
     1         gamma,z,yhati,sk(n),minsk(n),kern(21),wj,
     2         lambda,eta,ss2,lamski,swkp,yh,az,ss2w
C select neighbourhood
       do 1 i=1,n
         swkp=0
         yh=0.0
         ss2=0.0
         lamski=lambda*sk(i)
         yhati=yhat(i)
         jmin=max0(1,i-ndelta)
         jmax=min0(n,i+ndelta)
         do 11 j=jmin,jmax
            z=yhati-yhat(j)
            z=z*z/lamski
            iz=z
            if(iz.ge.20) goto 11
            az=z-iz
            wj=kern(iz+1)*(1-az)+kern(iz+2)*az
            swkp=swkp+wj
            ss2 = ss2 + wj*wj*s2hat(j)
            yh = yh + wj*y(j)
11          continue
C now we have the new wlj's and sk's 
C still to compute yhat
C now check if the new yhat can be accepted 
C and set new controls if (newcontrols != 0)
         yh=yh/swkp
         ss2=ss2/swkp/swkp
         if(controls(1,i).gt.yh.or.controls(2,i).lt.yh) goto 1
         if(ss2.gt.minsk(i)*gamma) goto 1
         if(newcontr.eq.0) goto 111
         ss2w=eta*sqrt(ss2)
         controls(1,i)=amax1(controls(1,i),yh-ss2w)
         controls(2,i)=amin1(controls(2,i),yh+ss2w)
111      yhatnew(i)=yh
         sknew(i)=ss2
1        continue
       return
       end

       subroutine locunial(kstar,n,y,yhat,yhatnew,sk,s2hat,sknew,
     1      controls,newcontr,minsk,ndelta,lambda,eta,gamma,kern,
     2      break)
C
C     performs univariate adaptive weights smoothing
C
C     called from R routine awsuni if not in sequential mode (graph=F)
C
C    kstar - number of iterations
C    n     - number of design points
C    y     - observations
C    yhat  - initial estimates 
C    yhatnew - new estimates (return)
C    sk    - initial variance estimates (variance of yhat)
C    sknew - new variance estimates (variance of yhatnew) (return)
C    s2hat - variance estimate (variance of y)
C    controls - lower und upper bounds in control step (input and return)
C    newcontr - vector of length kstar, with  either 0's or 1's 
C                determining whether  controls is to be updated or not 
C               (in dyadic control sceme, contains only 1's if
C                               control is activ in all steps) 
C    minsk  - minimal variances obtained
C    ndelta - vector of radii of neighbourhoods (in design points)
C    lambda - main smoothing parameter (should be around 3)
C    eta    - main control parameter (should be around 4)   
C    gamma  - allow for increase of variances (over minsk) by factor gamma 
C    kern   - vector of values of K(x) in seq(0,6,.3) (K(x)=exp(-x))
C    break  - terminate iteration if ||yhatnew-yhat||^2 < break 
C
       implicit logical (a-z)
       integer n,kstar,newcontr(kstar),i,ndelta(kstar),k
       real y(n),yhat(n),yhatnew(n),controls(2,n),s2hat(n),sknew(n),
     1    gamma,sk(n),minsk(n),kern(21),lambda,eta,break,z,zd,s
       do 1 k=2,kstar
          call locuniw(n,y,yhat,yhatnew,sk,s2hat,sknew,controls,
     1      newcontr(k),minsk,ndelta(k)-1,lambda,eta,gamma,kern)
          s=0.0
          do 11 i=1,n
             z=yhatnew(i)
             zd=z-yhat(i)
             s=s+zd*zd
             yhat(i)=z
             sk(i)=sknew(i)
             minsk(i)=amin1(minsk(i),sk(i))
11           continue 
          if(s.le.break) return    
1      continue
       return
       end

       subroutine locuini2(n,y,yhat,sk,ndelta,yh,sh,s2hat)
C
C     get initial estimates in case of nontrivial (1 point)
C     initial neighborhoods
C
C     should only be used in extreme cases (small signal/noise)
C
C     
C    n     - number of design points
C    y     - observations
C    yhat  - initial estimates (return)
C    sk    - initial variance estimates (variance of yhat)(return)
C    ndelta - radius of initial neighbourhood
C    yh, sh  - working arrays
C    s2hat  - variance estimates 
C
       implicit logical (a-z)
       integer n,ndelta,jmin,jmax,nj,j0,n1,n2,nbest,i,j,k
       real y(n),yhat(n),testj,yh(2*ndelta+1),s2hat(n),yj,yh1,z,
     1      testmax,sh(n),sj,sk(n)
       do 1 i=1,n
         jmin=max0(1,i-ndelta)
         jmax=min0(n,i+ndelta)
         nj=jmax-jmin+1
         do 11 j=1,nj
            sh(j)=0.0
11          yh(j)=0.0              
         j0=1
         do 12 j=jmin,jmax
            yj=y(j)
            sj=s2hat(j)
            do 121 k=1,j0
               sh(k)=sh(k)+sj
121            yh(k)=yh(k)+yj 
12          j0=j0+1
C        now test for inhomogeneous neighbourhoods
         yh1=yh(1)
         testmax=0.0
         do 13 n1=1,nj-1
            n2=nj-n1
            z=(yh1-yh(n1+1))/n1-yh(n1+1)/n2
            testj=n1*n2*z*z/nj
            if(testj.le.testmax) goto 13
            testmax=testj
            nbest=n1
13          continue
         if(testmax.gt.3*s2hat(i)) goto 14
         yhat(i)=yh(1)/nj
         sk(i)=sh(1)/nj/nj
         goto 1
14       if(jmin+nbest-1.lt.i) goto 15
         yhat(i)=(yh(1)-yh(nbest+1))/nbest
         sk(i)=(sh(1)-sh(nbest+1))/nbest/nbest
         goto 1
15       yhat(i)=yh(nbest+1)/(nj-nbest)
         sk(i)=sh(nbest+1)/(nj-nbest)
1        continue
       return
       end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C    Bivariate Adaptive weights smoothing                        C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine locbiw(nx,ny,y,yhat,yhatnew,sk,sknew,controls,
     1      newcontr,minsk,radiusq,shat,lam,eta,gamma,kern)
C
C     performs one step of bivariate adaptive weights smoothing
C
C     called from R routine awsuni if in sequential mode (graph=T)
C     and from locbiall otherwise
C  
C    nx, ny     - dimension of the image
C    y     - observations (image)
C    yhat  - estimates from last step
C    yhatnew - new estimates (return)
C    sk    - old variance estimates (variance of yhat)
C    sknew - new variance estimates (variance of yhatnew) (return)
C    controls - lower und upper bounds in control step (input and return)
C    newcontr - either 0 or 1 determining whether  controls is to be
C               updated or not (in dyadic control sceme, always 1 if
C                               control is activ in all steps) 
C    minsk  - minimal variances obtained
C    radiussq - radius^2 of neighbourhood 
C    shat  - variance estimate (variance of y)
C    lam   - main smoothing parameter (should be around 3)
C    eta    - main control parameter (should be around 4)   
C    gamma  - allow for increase of variances (over minsk) by factor gamma 
C    kern   - vector of values of K(x) in seq(0,6,.3) (K(x)=exp(-x))
C
       implicit logical (a-z)
       integer nx,ny,newcontr,ndeltaj,iz,
     1         ix,iy,jx,jy,ndelta,jxmin,jxmax,jymin,jymax
       real yhat(nx,ny),y(nx,ny),controls(2,nx,ny),radiusq,
     1      gamma,z,yh,sxy,yhati,yhatnew(nx,ny),
     2      lam,eta,shat(nx,ny),minsk(nx,ny),ss2w,az,
     2      sknew(nx,ny),sk(nx,ny),ski,kern(21),wj,swkp,ssk
C select neighbourhood
       ndelta=sqrt(radiusq)
       do 1 ix=1,nx
          do 1 iy=1,ny
             yhati=yhat(ix,iy)
             ski=sk(ix,iy)
             sxy=lam*ski
             swkp=0
             ssk=0.0
             yh=0.0
             jxmin=max0(1,ix-ndelta)
             jxmax=min0(nx,ix+ndelta)
             do 11 jx=jxmin,jxmax
                ndeltaj=sqrt(radiusq-(jx-ix)*(jx-ix))
                jymin=max0(1,iy-ndeltaj)
                jymax=min0(ny,iy+ndeltaj)
                do 11 jy=jymin,jymax
C   only use observations which contribute to variance reduction
                   if(shat(jx,jy).gt.2*shat(ix,iy)) goto 11
                   z=(yhati-yhat(jx,jy))
                   z=z*z/sxy
                   iz=z
                   if(iz.ge.20) goto 11
                   az=z-iz
                   wj=kern(iz+1)*(1-az)+kern(iz+2)*az
                   swkp=swkp+wj
                   ssk = ssk + wj*wj*shat(jx,jy)
                   yh = yh + wj*y(jx,jy)
11          continue
C now we have the new yh and sk for point (ix,iy)
C still to compute yhat
C now check if the new yhat can be accepted 
C and set new controls if (newcontrols != 0)
         yh=yh/swkp
         ssk=ssk/swkp/swkp
         if(controls(1,ix,iy).gt.yh.or.controls(2,ix,iy).lt.yh) goto 1
         if(ssk.gt.minsk(ix,iy)*gamma) goto 1
         if(newcontr.eq.0) goto 111
         ss2w=eta*sqrt(ssk)
         controls(1,ix,iy)=amax1(controls(1,ix,iy),yh-ss2w)
         controls(2,ix,iy)=amin1(controls(2,ix,iy),yh+ss2w)
111      yhatnew(ix,iy)=yh
         sknew(ix,iy)=ssk
1        continue
       return
       end

       subroutine locbiall(kstar,nx,ny,y,yhat,yhatnew,sk,sknew,
     1  controls,newcontr,minsk,radiusq,shat,lam,eta,gamma,kern)
C
C     performs bivariate adaptive weights smoothing
C
C     called from R routine awsuni if not in sequential mode (graph=F)
C
C    kstar - number of iterations
C    nx, ny   - dimensions of the image
C    y     - observations (image)
C    yhat  - initial estimates 
C    yhatnew - new estimates (return)
C    sk    - initial variance estimates (variance of yhat)
C    sknew - new variance estimates (variance of yhatnew) (return)
C    controls - lower und upper bounds in control step (input and return)
C    newcontr - vector of length kstar, with  either 0's or 1's 
C                determining whether  controls is to be updated or not 
C               (in dyadic control sceme, contains only 1's if
C                               control is activ in all steps) 
C    minsk  - minimal variances obtained
C    radiussq - radius^2 of neighbourhoods 
C    shat   - variance estimate (variance of y)
C    lambda - main smoothing parameter (should be around 3)
C    eta    - main control parameter (should be around 4)   
C    gamma  - allow for increase of variances (over minsk) by factor gamma 
C    kern   - vector of values of K(x) in seq(0,6,.3) (K(x)=exp(-x))
C
       implicit logical (a-z)
       integer nx,ny,kstar,newcontr(kstar),k,ix,iy
       real yhat(nx,ny),y(nx,ny),controls(2,nx,ny),gamma,
     2      yhatnew(nx,ny),lam,eta,shat(nx,ny),minsk(nx,ny),
     2      sknew(nx,ny),sk(nx,ny),kern(21),radiusq(kstar),eps
       eps=.0001
       do 1 k=2,kstar
          call locbiw(nx,ny,y,yhat,yhatnew,sk,sknew,controls,
     1      newcontr(k),minsk,radiusq(k)+eps,shat,lam,eta,gamma,kern)
C      prepare for next iteration
          do 11 ix=1,nx
             do 11 iy=1,ny
                yhat(ix,iy)=yhatnew(ix,iy)
                sk(ix,iy)=sknew(ix,iy)
                minsk(ix,iy)=amin1(minsk(ix,iy),sk(ix,iy))
11              continue
1      continue
       return
       end

       subroutine getnubi(radiusq,weight,nu,n)
C
C    calculate number of points in a neighbourhood with 
C    given radius^2  
C    weight - excentricities of ellipses used as neighbourhoods
C    n - length of radiusq, i.e. number of neighbouhoods in sequence
C
       implicit logical (a-z)
       integer n,nu(n),jx,jy,ndelta,jxmin,jxmax,jymin,jymax,i,ndeltaj
       real radiusq(n),weight(2),w1
C Calculate size of neighborhood
       w1=weight(1)
       do 1 i=1,n
          nu(i)=0
          ndelta=sqrt(radiusq(i)/w1)
             jxmin=-ndelta
             jxmax=ndelta
             do 11 jx=jxmin,jxmax
               ndeltaj=sqrt((radiusq(i)-jx*jx*w1)/weight(2))
               jymin=-ndeltaj
               jymax=ndeltaj
               do 11 jy=jymin,jymax
                   nu(i)=nu(i)+1
11          continue
1      continue
       return
       end

       subroutine locbinis(nx,ny,y,yhat,s2hat,sk,radiusq)
C
C     get initial estimates in case of nontrivial (1 point)
C     initial neighborhoods (bivariate)
C
C     should only be used in extreme cases (small signal/noise)
C
C     
C    nx, ny  - dimensions of the image
C    y     - observations
C    yhat  - initial estimates (return)
C    sk    - initial variance estimates (variance of yhat)(return)
C    radiussq - radius^2 of initial neighbourhood
C    yh, sh  - working arrays
C    s2hat  - variance estimates 
C
       implicit logical (a-z)
       integer nx,ny,ix,iy,jx,jy,ndelta,swkp,jxmin,jxmax,k,
     1         jymin,jymax,nt(132),ii,jj,nind,swkp0,ndeltaj
       logical ind(132)
       real yhat(nx,ny),y(nx,ny),radiusq,yh,s2hat(nx,ny),st2(132),
     1      sk(nx,ny),s2,s20,yht(132),test,test0,yakt,s2hatakt,yh2,yh0
C    y     I: observations
C    yhat  O: initial estimate
C    sk    O: sum_{U_1} w(,)
C select neighbourhood
       ndelta=sqrt(radiusq)
       nind=0
       if(radiusq.lt.1) goto 100
       nind=8
       if(radiusq.lt.2) goto 100
       nind=16
       if(radiusq.lt.4) goto 100
       nind=64
       if(radiusq.lt.5) goto 100
       nind=132
100    continue
       do 1 ix=1,nx
          do 1 iy=1,ny
             do 19 k=1,nind
                nt(k)=0
                yht(k)=0.0
                st2(k)=0.0
19                ind(k)=.TRUE.
             swkp=0
             s2=0.0
             yh=0.0
             yh2=0.0
             jxmin=max0(1,ix-ndelta)
             jxmax=min0(nx,ix+ndelta)
             do 11 jx=jxmin,jxmax
                ii=ix-jx
                ndeltaj=sqrt(radiusq-(jx-ix)*(jx-ix))
                jymin=max0(1,iy-ndeltaj)
                jymax=min0(ny,iy+ndeltaj)
                do 11 jy=jymin,jymax
                   yakt=y(jx,jy)
                   s2hatakt=s2hat(jx,jy)
                   jj=iy-jy
                   swkp=swkp+1
                   s2=s2+s2hatakt
                   yh=yh+yakt
                   if(radiusq.lt.1) goto 11
                   if(ii.gt.-.5) ind(1)=.FALSE.
                   if(jj.gt.-.5) ind(2)=.FALSE.
                   if(ii+jj.gt.-.5) ind(3)=.FALSE.
                   if(ii-jj.gt.-.5) ind(4)=.FALSE.
                   if(ii.lt..5) ind(5)=.FALSE.
                   if(jj.lt..5) ind(6)=.FALSE.
                   if(ii+jj.lt..5) ind(7)=.FALSE.
                   if(ii-jj.lt..5) ind(8)=.FALSE.
                   if(radiusq.lt.2) goto 12
                   if(ii+2*jj.gt.-.5) ind(9)=.FALSE.
                   if(jj+2*ii.gt.-.5) ind(10)=.FALSE.
                   if(ii-2*jj.gt.-.5) ind(11)=.FALSE.
                   if(jj-2*ii.gt.-.5) ind(12)=.FALSE.
                   if(ii+2*jj.lt..5) ind(13)=.FALSE.
                   if(jj+2*ii.lt..5) ind(14)=.FALSE.
                   if(ii-2*jj.lt..5) ind(15)=.FALSE.
                   if(jj-2*ii.lt..5) ind(16)=.FALSE.
                   if(radiusq.lt.4) goto 12
                   if(ii.gt.-1.5) ind(17)=.FALSE.
                   if(jj.gt.-1.5) ind(18)=.FALSE.
                   if(ii+jj.gt.-1.5) ind(19)=.FALSE.
                   if(ii-jj.gt.-1.5) ind(20)=.FALSE.
                   if(ii.lt.1.5) ind(21)=.FALSE.
                   if(jj.lt.1.5) ind(22)=.FALSE.
                   if(ii+jj.lt.1.5) ind(23)=.FALSE.
                   if(ii-jj.lt.1.5) ind(24)=.FALSE.
                   if(ii+2*jj.gt.-1.5) ind(25)=.FALSE.
                   if(jj+2*ii.gt.-1.5) ind(26)=.FALSE.
                   if(ii-2*jj.gt.-1.5) ind(27)=.FALSE.
                   if(jj-2*ii.gt.-1.5) ind(28)=.FALSE.
                   if(ii+2*jj.lt.1.5) ind(29)=.FALSE.
                   if(jj+2*ii.lt.1.5) ind(30)=.FALSE.
                   if(ii-2*jj.lt.1.5) ind(31)=.FALSE.
                   if(jj-2*ii.lt.1.5) ind(32)=.FALSE.
                   if(ii+2*jj.gt.-2.5) ind(33)=.FALSE.
                   if(jj+2*ii.gt.-2.5) ind(34)=.FALSE.
                   if(ii-2*jj.gt.-2.5) ind(35)=.FALSE.
                   if(jj-2*ii.gt.-2.5) ind(36)=.FALSE.
                   if(ii+2*jj.lt.2.5) ind(37)=.FALSE.
                   if(jj+2*ii.lt.2.5) ind(38)=.FALSE.
                   if(ii-2*jj.lt.2.5) ind(39)=.FALSE.
                   if(jj-2*ii.lt.2.5) ind(40)=.FALSE.
                   if(ii+3*jj.gt.-.5) ind(41)=.FALSE.
                   if(jj+3*ii.gt.-.5) ind(42)=.FALSE.
                   if(ii-3*jj.gt.-.5) ind(43)=.FALSE.
                   if(jj-3*ii.gt.-.5) ind(44)=.FALSE.
                   if(ii+3*jj.lt..5) ind(45)=.FALSE.
                   if(jj+3*ii.lt..5) ind(46)=.FALSE.
                   if(ii-3*jj.lt..5) ind(47)=.FALSE.
                   if(jj-3*ii.lt..5) ind(48)=.FALSE.
                   if(ii+3*jj.gt.-1.5) ind(49)=.FALSE.
                   if(jj+3*ii.gt.-1.5) ind(50)=.FALSE.
                   if(ii-3*jj.gt.-1.5) ind(51)=.FALSE.
                   if(jj-3*ii.gt.-1.5) ind(52)=.FALSE.
                   if(ii+3*jj.lt.1.5) ind(53)=.FALSE.
                   if(jj+3*ii.lt.1.5) ind(54)=.FALSE.
                   if(ii-3*jj.lt.1.5) ind(55)=.FALSE.
                   if(jj-3*ii.lt.1.5) ind(56)=.FALSE.
                   if(ii+3*jj.gt.-2.5) ind(57)=.FALSE.
                   if(jj+3*ii.gt.-2.5) ind(58)=.FALSE.
                   if(ii-3*jj.gt.-2.5) ind(59)=.FALSE.
                   if(jj-3*ii.gt.-2.5) ind(60)=.FALSE.
                   if(ii+3*jj.lt.2.5) ind(61)=.FALSE.
                   if(jj+3*ii.lt.2.5) ind(62)=.FALSE.
                   if(ii-3*jj.lt.2.5) ind(63)=.FALSE.
                   if(jj-3*ii.lt.2.5) ind(64)=.FALSE.
                   if(radiusq.lt.5) goto 12
                   if(ii+jj.gt.-2.5) ind(65)=.FALSE.
                   if(ii-jj.gt.-2.5) ind(66)=.FALSE.
                   if(ii+jj.lt.2.5) ind(67)=.FALSE.
                   if(ii-jj.lt.2.5) ind(68)=.FALSE.
                   if(ii+2*jj.gt.-3.5) ind(69)=.FALSE.
                   if(jj+2*ii.gt.-3.5) ind(70)=.FALSE.
                   if(ii-2*jj.gt.-3.5) ind(71)=.FALSE.
                   if(jj-2*ii.gt.-3.5) ind(72)=.FALSE.
                   if(ii+2*jj.lt.3.5) ind(73)=.FALSE.
                   if(jj+2*ii.lt.3.5) ind(74)=.FALSE.
                   if(ii-2*jj.lt.3.5) ind(75)=.FALSE.
                   if(jj-2*ii.lt.3.5) ind(76)=.FALSE.
                   if(ii+3*jj.gt.-3.5) ind(77)=.FALSE.
                   if(jj+3*ii.gt.-3.5) ind(78)=.FALSE.
                   if(ii-3*jj.gt.-3.5) ind(79)=.FALSE.
                   if(jj-3*ii.gt.-3.5) ind(80)=.FALSE.
                   if(ii+3*jj.lt.3.5) ind(81)=.FALSE.
                   if(jj+3*ii.lt.3.5) ind(82)=.FALSE.
                   if(ii-3*jj.lt.3.5) ind(83)=.FALSE.
                   if(jj-3*ii.lt.3.5) ind(84)=.FALSE.
                   if(2*ii+3*jj.gt.-1) ind(85)=.FALSE.
                   if(2*jj+3*ii.gt.-1) ind(86)=.FALSE.
                   if(2*ii-3*jj.gt.-1) ind(87)=.FALSE.
                   if(2*jj-3*ii.gt.-1) ind(88)=.FALSE.
                   if(2*ii+3*jj.lt.1) ind(89)=.FALSE.
                   if(2*jj+3*ii.lt.1) ind(90)=.FALSE.
                   if(2*ii-3*jj.lt.1) ind(91)=.FALSE.
                   if(2*jj-3*ii.lt.1) ind(92)=.FALSE.
                   if(2*ii+3*jj.gt.-2.5) ind(93)=.FALSE.
                   if(2*jj+3*ii.gt.-2.5) ind(94)=.FALSE.
                   if(2*ii-3*jj.gt.-2.5) ind(95)=.FALSE.
                   if(2*jj-3*ii.gt.-2.5) ind(96)=.FALSE.
                   if(2*ii+3*jj.lt.2.5) ind(97)=.FALSE.
                   if(2*jj+3*ii.lt.2.5) ind(98)=.FALSE.
                   if(2*ii-3*jj.lt.2.5) ind(99)=.FALSE.
                   if(2*jj-3*ii.lt.2.5) ind(100)=.FALSE.
                   if(2*ii+3*jj.gt.-3.5) ind(101)=.FALSE.
                   if(2*jj+3*ii.gt.-3.5) ind(102)=.FALSE.
                   if(2*ii-3*jj.gt.-3.5) ind(103)=.FALSE.
                   if(2*jj-3*ii.gt.-3.5) ind(104)=.FALSE.
                   if(2*ii+3*jj.lt.3.5) ind(105)=.FALSE.
                   if(2*jj+3*ii.lt.3.5) ind(106)=.FALSE.
                   if(2*ii-3*jj.lt.3.5) ind(107)=.FALSE.
                   if(2*jj-3*ii.lt.3.5) ind(108)=.FALSE.
                   if(2*ii+3*jj.gt.-4.5) ind(109)=.FALSE.
                   if(2*jj+3*ii.gt.-4.5) ind(110)=.FALSE.
                   if(2*ii-3*jj.gt.-4.5) ind(111)=.FALSE.
                   if(2*jj-3*ii.gt.-4.5) ind(112)=.FALSE.
                   if(2*ii+3*jj.lt.4.5) ind(113)=.FALSE.
                   if(2*jj+3*ii.lt.4.5) ind(114)=.FALSE.
                   if(2*ii-3*jj.lt.4.5) ind(115)=.FALSE.
                   if(2*jj-3*ii.lt.4.5) ind(116)=.FALSE.
                   if(2*ii+3*jj.gt.-5.5) ind(117)=.FALSE.
                   if(2*jj+3*ii.gt.-5.5) ind(118)=.FALSE.
                   if(2*ii-3*jj.gt.-5.5) ind(119)=.FALSE.
                   if(2*jj-3*ii.gt.-5.5) ind(120)=.FALSE.
                   if(2*ii+3*jj.lt.5.5) ind(121)=.FALSE.
                   if(2*jj+3*ii.lt.5.5) ind(122)=.FALSE.
                   if(2*ii-3*jj.lt.5.5) ind(123)=.FALSE.
                   if(2*jj-3*ii.lt.5.5) ind(124)=.FALSE.
                   if(2*ii+3*jj.gt.-6.5) ind(125)=.FALSE.
                   if(2*jj+3*ii.gt.-6.5) ind(126)=.FALSE.
                   if(2*ii-3*jj.gt.-6.5) ind(127)=.FALSE.
                   if(2*jj-3*ii.gt.-6.5) ind(128)=.FALSE.
                   if(2*ii+3*jj.lt.6.5) ind(129)=.FALSE.
                   if(2*jj+3*ii.lt.6.5) ind(130)=.FALSE.
                   if(2*ii-3*jj.lt.6.5) ind(131)=.FALSE.
                   if(2*jj-3*ii.lt.6.5) ind(132)=.FALSE.
12              do 13 k=1,nind
                   if(ind(k)) goto 13
                   yht(k)=yht(k)+yakt
                   nt(k)=nt(k)+1
                   st2(k)=st2(k)+s2hatakt
13              continue
11          continue
C   now compute statistic for selection of U
         if(nind.eq.0) goto 15 
            test0=-(yh*yh)/swkp+2*s2/swkp
            swkp0=swkp
            s20=s2
            yh0=yh
            do 14 k=1,nind
               if(nt(k).eq.swkp0) goto 14
               test=-yht(k)*yht(k)/nt(k)+st2(k)/nt(k)-
     1          (yh0-yht(k))*(yh0-yht(k))/(swkp0-nt(k))+
     2          (s20-st2(k))/(swkp0-nt(k))   
               if(test.ge.test0) goto 14
               test0=test
               yh=yht(k)
               swkp=nt(k)
               s2=st2(k)
14          continue              
15       sk(ix,iy)=s2/swkp/swkp
         yhat(ix,iy)=yh/swkp
1        continue
       return
       end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C    Trivariate Adaptive weights smoothing                       C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine loctriw(nx,ny,nz,y,yhat,yhatnew,sk,sknew,controls,
     1       newcontr,minsk,radiusq,shat,lam,eta,gamma,weight,kern)
C
C     performs one step of trivariate adaptive weights smoothing
C
C     called  from loctrial
C  
C    nx, ny, nz    - dimension of the image
C    y     - observations (image)
C    yhat  - estimates from last step
C    yhatnew - new estimates (return)
C    sk    - old variance estimates (variance of yhat)
C    sknew - new variance estimates (variance of yhatnew) (return)
C    controls - lower und upper bounds in control step (input and return)
C    newcontr - either 0 or 1 determining whether  controls is to be
C               updated or not (in dyadic control sceme, always 1 if
C                               control is activ in all steps) 
C    minsk  - minimal variances obtained
C    radiussq - radius^2 of neighbourhood 
C    shat  - variance estimate (variance of y)
C    lam   - main smoothing parameter (should be around 3)
C    eta    - main control parameter (should be around 4)   
C    gamma  - allow for increase of variances (over minsk) by factor gamma 
C    weight - excentricities of ellipsoids used as neighbourhoods
C         (this allows to deal with different spacings in vertical direction)
C    kern   - vector of values of K(x) in seq(0,6,.3) (K(x)=exp(-x))
C
       implicit logical (a-z)
       integer nx,ny,nz,jzmin,jzmax,iz,jz,izz,newcontr,ndeltaj,
     1         ix,iy,jx,jy,ndelta,jxmin,jxmax,jymin,jymax
       real yhat(nx,ny,nz),y(nx,ny,nz),controls(2,nx,ny,nz),radiusq,
     1      gamma,z,az,yh,sxy,yhati,skig,weight(3),w1,w2,kern(21),
     2      yhatnew(nx,ny,nz),lam,eta,shat(nx,ny,nz),sk(nx,ny,nz),ski,
     3      minsk(nx,ny,nz),sknew(nx,ny,nz),ssk,wj,swkp,sskw
C    y     I: observations
C    lshat I: lambda * sigmahat^2
C    gamma I: gamma
C    yhat  I: estimates of from step k-1
C          O:              after step k
C    sk      I: sum_{U_k-1} w(l,)
C            O: sum_{U_k} w(l,)
C    yhatold: I: old estimates for comparisons
C    skold:      variances of old estimates
C    maxsk:      maximal sk's from previous steps
C select neighbourhood
       w1=weight(1)
       ndelta=sqrt(radiusq/w1)
       do 1 ix=1,nx
         do 1 iy=1,ny
           do 1 iz=1,nz
             yhati=yhat(ix,iy,iz)
             ski=sk(ix,iy,iz)
             sxy=lam*ski
             skig=ski*gamma
             swkp=0
             ssk=0.0
             yh=0.0
             jxmin=max0(1,ix-ndelta)
             jxmax=min0(nx,ix+ndelta)
             do 11 jx=jxmin,jxmax
               w2=weight(2)
               ndeltaj=sqrt((radiusq-(jx-ix)*(jx-ix)*w1)/w2)
               jymin=max0(1,iy-ndeltaj)
               jymax=min0(ny,iy+ndeltaj)
               do 11 jy=jymin,jymax
                 ndeltaj=sqrt((radiusq-(jx-ix)*(jx-ix)*w1-
     1                            (jy-iy)*(jy-iy)*w2)/weight(3))
                 jzmin=max0(1,iz-ndeltaj)
                 jzmax=min0(nz,iz+ndeltaj)
                 do 11 jz=jzmin,jzmax
                   z=(yhati-yhat(jx,jy,jz))
                   z=z*z/sxy
                   izz=z
                   if(izz.ge.20) goto 11
                   az=z-izz
                   wj=kern(izz+1)*(1-az)+kern(izz+2)*az
                   swkp=swkp+wj
                   ssk = ssk + wj*wj*shat(jx,jy,jz)
                   yh = yh + wj*y(jx,jy,jz)
11          continue
C now we have the new yh and sk for point (ix,iy,iz)
C still to compute yhat
C now check if the new yhat can be accepted 
C and set new controls if (newcontrols != 0)
         yh=yh/swkp
         ssk=ssk/swkp/swkp
         if(controls(1,ix,iy,iz).gt.yh.or.controls(2,ix,iy,iz).lt.yh) 
     1   goto 1
         if(ssk.gt.minsk(ix,iy,iz)*gamma) goto 1
         if(newcontr.eq.0) goto 111
         sskw=eta*sqrt(ssk)
         controls(1,ix,iy,iz)=amax1(controls(1,ix,iy,iz),yh-sskw)
         controls(2,ix,iy,iz)=amin1(controls(2,ix,iy,iz),yh+sskw)
111      yhatnew(ix,iy,iz)=yh
         sknew(ix,iy,iz)=ssk
1        continue
       return
       end

       subroutine loctrial(kstar,nx,ny,nz,y,yhat,yhatnew,sk,sknew,
     1 controls,newcontr,minsk,radiusq,shat,lam,eta,gamma,weight,kern)
C
C     performs trivariate adaptive weights smoothing
C
C     called from R routine awsuni if not in sequential mode (graph=F)
C
C    kstar - number of iterations
C    nx, ny, nz   - dimensions of the image
C    y     - observations (image)
C    yhat  - initial estimates 
C    yhatnew - new estimates (return)
C    sk    - initial variance estimates (variance of yhat)
C    sknew - new variance estimates (variance of yhatnew) (return)
C    controls - lower und upper bounds in control step (input and return)
C    newcontr - vector of length kstar, with  either 0's or 1's 
C                determining whether  controls is to be updated or not 
C               (in dyadic control sceme, contains only 1's if
C                               control is activ in all steps) 
C    minsk  - minimal variances obtained
C    radiussq - radius^2 of neighbourhoods 
C    shat   - variance estimate (variance of y)
C    lambda - main smoothing parameter (should be around 3)
C    eta    - main control parameter (should be around 4)   
C    gamma  - allow for increase of variances (over minsk) by factor gamma 
C    weight - excentricities of ellipsoids used as neighbourhoods
C         (this allows to deal with different spacings in vertical direction)
C    kern   - vector of values of K(x) in seq(0,6,.3) (K(x)=exp(-x))
C
       implicit logical (a-z)
       integer nx,ny,nz,kstar,newcontr(kstar),k,ix,iy,iz
       real yhat(nx,ny,nz),y(nx,ny,nz),controls(2,nx,ny,nz),gamma,
     2   yhatnew(nx,ny,nz),lam,eta,shat(nx,ny,nz),minsk(nx,ny,nz),
     3   sknew(nx,ny,nz),sk(nx,ny,nz),kern(21),radiusq(kstar),
     4   weight(3),eps
       eps=.0001
       do 1 k=2,kstar
         call loctriw(nx,ny,nz,y,yhat,yhatnew,sk,sknew,controls,
     1     newcontr(k),minsk,radiusq(k)+eps,shat,lam,eta,gamma,weight,
     1     kern)
C      prepare for next iteration
          do 11 ix=1,nx
             do 11 iy=1,ny
                do 11 iz=1,nz
                yhat(ix,iy,iz)=yhatnew(ix,iy,iz)
                sk(ix,iy,iz)=sknew(ix,iy,iz)
                minsk(ix,iy,iz)=amin1(minsk(ix,iy,iz),sk(ix,iy,iz))
11              continue
1      continue
       return
       end


       subroutine loctriw0(nx,ny,nz,y,yhat,yhatnew,sk,sknew,controls,
     1       newcontr,radiusq,shat,lam,eta,weight,kern)
       implicit logical (a-z)
       integer nx,ny,nz,jzmin,jzmax,iz,jz,izz,newcontr,
     1         ix,iy,jx,jy,ndelta,jxmin,jxmax,jymin,jymax,ndeltaj
       real yhat(nx,ny,nz),y(nx,ny,nz),controls(2,nx,ny,nz),radiusq,
     1      z,yh,sxy,yhati,weight(3),w1,w2,kern(21),az,
     2      yhatnew(nx,ny,nz),lam,eta,shat,sk(nx,ny,nz),ski,
     3      sknew(nx,ny,nz),ssk,wj,swkp,sskw
C
C     performs one step of trivariate adaptive weights smoothing
C     in case of homogeneous variance shat 
C     just do save memory compared to loctriw
C
C     called  from loctria0   
C  
C    nx, ny, nz    - dimension of the image
C    y     - observations (image)
C    yhat  - estimates from last step
C    yhatnew - new estimates (return)
C    sk    - old variance estimates (variance of yhat)
C    sknew - new variance estimates (variance of yhatnew) (return)
C    controls - lower und upper bounds in control step (input and return)
C    newcontr - either 0 or 1 determining whether  controls is to be
C               updated or not (in dyadic control sceme, always 1 if
C                               control is activ in all steps) 
C    minsk  - minimal variances obtained
C    radiussq - radius^2 of neighbourhood 
C    shat  - homogeneous variance estimate (variance of y)
C    lam   - main smoothing parameter (should be around 3)
C    eta    - main control parameter (should be around 4)   
C    gamma  - allow for increase of variances (over minsk) by factor gamma 
C    weight - excentricities of ellipsoids used as neighbourhoods
C         (this allows to deal with different spacings in vertical direction)
C    kern   - vector of values of K(x) in seq(0,6,.3) (K(x)=exp(-x))
C
C select neighbourhood
       w1=weight(1)
       ndelta=sqrt(radiusq/w1)
       do 1 ix=1,nx
         do 1 iy=1,ny
           do 1 iz=1,nz
             yhati=yhat(ix,iy,iz)
             ski=sk(ix,iy,iz)
             sxy=lam*ski
             swkp=0
             ssk=0.0
             yh=0.0
             jxmin=max0(1,ix-ndelta)
             jxmax=min0(nx,ix+ndelta)
             do 11 jx=jxmin,jxmax
               w2=weight(2)
               ndeltaj=sqrt((radiusq-(jx-ix)*(jx-ix)*w1)/w2)
               jymin=max0(1,iy-ndeltaj)
               jymax=min0(ny,iy+ndeltaj)
               do 11 jy=jymin,jymax
                 ndeltaj=sqrt((radiusq-(jx-ix)*(jx-ix)*w1-
     1                        (jy-iy)*(jy-iy)*w2)/weight(3))
                 jzmin=max0(1,iz-ndeltaj)
                 jzmax=min0(nz,iz+ndeltaj)
                 do 11 jz=jzmin,jzmax
                   z=(yhati-yhat(jx,jy,jz))
                   z=z*z/sxy
                   izz=z
                   if(izz.ge.20) goto 11
                   az=z-izz
                   wj=kern(izz+1)*(1-az)+kern(izz+2)*az
                   swkp=swkp+wj
                   ssk = ssk + wj*wj*shat
                   yh = yh + wj*y(jx,jy,jz)
11          continue
C now we have the new yh and sk for point (ix,iy,iz)
C still to compute yhat
C now check if the new yhat can be accepted 
C and set new controls if (newcontrols != 0)
         yh=yh/swkp
         ssk=ssk/swkp/swkp
         if(controls(1,ix,iy,iz).gt.yh.or.controls(2,ix,iy,iz).lt.yh) 
     1   goto 1
         if(newcontr.eq.0) goto 111
         sskw=eta*sqrt(ssk)
         controls(1,ix,iy,iz)=amax1(controls(1,ix,iy,iz),yh-sskw)
         controls(2,ix,iy,iz)=amin1(controls(2,ix,iy,iz),yh+sskw)
111      yhatnew(ix,iy,iz)=yh
         sknew(ix,iy,iz)=ssk
1        continue
       return
       end

       subroutine loctria0(kstar,nx,ny,nz,y,yhat,yhatnew,sk,sknew,
     1 controls,newcontr,radiusq,shat,lam,eta,weight,kern)
C
C     performs trivariate adaptive weights smoothing
C     in case of homogeneous variance shat 
C     just do save memory compared to loctrial
C
C     called from R routine awsuni if not in sequential mode (graph=F)
C
C    kstar - number of iterations
C    nx, ny, nz   - dimensions of the image
C    y     - observations (image)
C    yhat  - initial estimates 
C    yhatnew - new estimates (return)
C    sk    - initial variance estimates (variance of yhat)
C    sknew - new variance estimates (variance of yhatnew) (return)
C    controls - lower und upper bounds in control step (input and return)
C    newcontr - vector of length kstar, with  either 0's or 1's 
C                determining whether  controls is to be updated or not 
C               (in dyadic control sceme, contains only 1's if
C                               control is activ in all steps) 
C    minsk  - minimal variances obtained
C    radiussq - radius^2 of neighbourhoods 
C    shat   - homogeneous variance estimate (variance of y)
C    lambda - main smoothing parameter (should be around 3)
C    eta    - main control parameter (should be around 4)   
C    gamma  - allow for increase of variances (over minsk) by factor gamma 
C    weight - excentricities of ellipsoids used as neighbourhoods
C         (this allows to deal with different spacings in vertical direction)
C    kern   - vector of values of K(x) in seq(0,6,.3) (K(x)=exp(-x))
C
       implicit logical (a-z)
       integer nx,ny,nz,kstar,newcontr(kstar),k,ix,iy,iz
       real yhat(nx,ny,nz),y(nx,ny,nz),controls(2,nx,ny,nz),
     2   yhatnew(nx,ny,nz),lam,eta,shat,
     3   sknew(nx,ny,nz),sk(nx,ny,nz),kern(21),radiusq(kstar),
     4   weight(3),eps
       eps=.0001
       do 1 k=2,kstar
         call loctriw0(nx,ny,nz,y,yhat,yhatnew,sk,sknew,controls,
     1     newcontr(k),radiusq(k)+eps,shat,lam,eta,weight,kern)
C      prepare for next iteration
          do 11 ix=1,nx
             do 11 iy=1,ny
                do 11 iz=1,nz
                yhat(ix,iy,iz)=yhatnew(ix,iy,iz)
                sk(ix,iy,iz)=sknew(ix,iy,iz)
11              continue
1      continue
       return
       end
