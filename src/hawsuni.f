CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize estimates in local constant univariate aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ihawsuni(y,n,hinit,bi,bi2,ai,kern,sigma2)
C   
C   y        observed values of regression function
C   n        number of observations
C   hinit    initial bandwidth
C   bi       \sum \Psi^T Wi \Psi  (output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kern     discretized localization kernel 
C
      implicit logical (a-z)
      integer n
      real*8 y(n),hinit,bi(n),bi2(n),ai(n),kern(102),sigma2(n)
      integer i,j,ja,je,ih,iz
      real*8 z,wj,az,swj,swjy,swj2
      ih=hinit
      do 1 i=1,n
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         swj=0.d0
         swj2=0.d0
         swjy=0.d0
         do 2 j=ja,je
            z=(i-j)/hinit
            z=z*z*1.d2
            if(z.ge.1.d2) goto 2
            iz=z
            az=z-iz
            wj=kern(iz+1)*(1-az)+kern(iz+2)*az
            wj=wj/sigma2(j)
            swj=swj+wj
            swj2=swj2+wj*wj
            swjy=swjy+wj*y(j)
2        continue
         ai(i)=swjy
         bi(i)=swj
         bi2(i)=swj2
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant univariate aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lhawsuni(y,n,hakt,lamakt,theta,bi,bi2,ai,kernl,kerns,
     1                    sym,sigma2)
C   
C   y        observed values of regression function
C   n        number of observations
C   hakt     actual bandwidth
C   lamakt   lambda*sigma2 
C   theta    estimates    (output)
C   bi       \sum \Psi^T Wi \Psi  (output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n
      logical sym
      real*8 y(n),hakt,theta(n),bi(n),bi2(n),ai(n),kernl(102),
     1      kerns(102),lamakt,lam,sigma2(n)
      integer i,j,ja,je,ih,iz
      real*8 z,wj,az,swj,swj2,swjy,bii,thetai
      ih=hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i=1,n
         thetai=theta(i)
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         swj=0.d0
         swj2=0.d0
         swjy=0.d0
         do 2 j=ja,je
C  first stochastic term
            bii=bi(i)
            if(sym) bii=bii+bi(j)
            z=(thetai-theta(j))
            z=z*z*bii/lam
            if(z.ge.1.d2) goto 2
            iz=z
            az=z-iz
            wj=kerns(iz+1)*(1-az)+kerns(iz+2)*az
            z=(i-j)/hakt
            z=z*z*1.d2
            if(z.ge.1.d2) goto 2
            iz=z
            az=z-iz
            wj=wj*(kernl(iz+1)*(1-az)+kernl(iz+2)*az)
            wj=wj/sigma2(j)
            swj=swj+wj
            swj2=swj2+wj*wj
            swjy=swjy+wj*y(j)
2        continue
         ai(i)=swjy
         bi(i)=swj
         bi2(i)=swj2
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform all iterations in local constant univariate aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ghawsuni(y,n,hinit,hincr,hmax,lamakt,eta,theta, 
     1                   bi,bi2,ai,kernl,kerns,biold,sym,sigma2)
C   
C   y        observed values of regression function
C   n        number of observations
C   hinit    initial bandwidth
C   hincr    factor used to increase bandwiths
C   hmax     maximal bandwidth
C   lamakt   lambda*sigma2 
C   eta      memory parameter
C   theta    estimates    (output)
C   bi       \sum \Psi^T Wi \Psi  (output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   biold    working array to store old values of bi
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n
      logical sym
      real*8 y(n),hinit,hincr,hmax,theta(n),bi(n),bi2(n),ai(n),
     1       kernl(102),kerns(102),lamakt,eta,biold(n),sigma2(n)
      integer i
      real*8 hakt,onemeta,z
      hakt=hinit*hincr
      onemeta=1.d0-eta
1     call lhawsuni(y,n,hakt,lamakt,theta,bi,bi2,ai,kernl,kerns,sym,
     1              sigma2)
      do 11 i=1,n
         z=onemeta*ai(i)+eta*biold(i)*theta(i)
         bi(i)=onemeta*bi(i)+eta*biold(i)
         theta(i)=z/bi(i)
         biold(i)=bi(i)
11    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Univariate local polynomial aws 
C
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize estimates in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine iphawsun(n,dp1,dp2,x,y,hinit,bi,ai,theta,
     1                    kernl,dmat,sigma2)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     x          design points (ordered)
C     y          observed values at design points
C     hinit      initial bandwidth
C     bi         \sum \Psi^T Wi \Psi  (output)
C     ai         \sum \Psi^T Wi Y     (output)
C     theta      initial estimates    (output)
C     kernl      discretized localization kernel 
C     dmat       working array
C
      implicit logical (a-z)
      integer n,dp1,dp2,i,j,k,info,je,ja,iz
      real*8 bi(dp2,n),ai(dp1,n),theta(dp1,n),dmat(dp1,dp1),ha2,d
      real*8 x(n),xij,xi,z,epij,y(n),hinit,ha,kernl(102),az,sigma2(n)
C     loop over i=1,n
C     use points within (xi-hinit,xi+hinit) but at least dp1 points
      do 1 i=1,n
         xi=x(i)
         ha=hinit
C     first get ja and je 
1099     do 1000 j=i,1,-1
            if(x(j).le.xi-ha) goto 1001
            ja=j
1000     continue
1001     do 1002 j=i,n
            if(x(j).gt.xi+ha) goto 1003
            je=j
1002     continue
1003     if(je-ja.gt.dp1) goto 1004
             ha=ha*1.25
C            not enough points in neighborhood to estimate parameters
C            increase ha
             goto 1099
1004     continue
C        first fill bi(,i) and ai(,i)  
         ha2=ha*ha
         do 11 j=ja,je
            xij=(x(j)-xi)
            z=1.d2*xij*xij/ha2
            iz=z
            az=z-iz
            epij=(1.-az)*kernl(iz+1)+az*kernl(iz+2)
            epij=epij/sigma2(j)
            z=1.d0
            do 12 k=1,dp2
               bi(k,i)=bi(k,i)+z*epij
               if(k.le.dp1) ai(k,i)=ai(k,i)+z*y(j)*epij
               z=z*xij
12         continue
11      continue
C       expand bi as p times p matrix
        do 13 j=1,dp1
        do 14 k=1,dp1
           if(j.gt.k) then 
              dmat(j,k)=0.0d0
           else
              dmat(j,k)=bi(j+k-1,i)
           end if
14          continue
13       continue
C     compute inverse of bi by choleski decomposition
         if(dp1.eq.1) goto 15
         call invers(dmat,dp1,info)
         goto 16
15       dmat(1,1)=1.d0/dmat(1,1)
C     now dmat contains inverse of B_i 
C     now calculate theta as B_i^{-1} Z_i
16       do 20 j=1,dp1
           d=0.0d0
           do 21 k=1,dp1
              d=d+dmat(j,k)*ai(k,i)
21         continue
           theta(j,i)=d
20      continue
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lphawsun(n,dp1,dp2,x,y,theta,bi,bin,si0,bi0,
     1    ain,lam,tau,h,kernl,kerns,cb,dmat,dmati,dmat0,
     2    thij,thji,psix,psiy,sym,sigma2)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     x          design points (ordered)
C     y          observed values at design points
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     si0        \sum \Psi^T Wi0 \Psi [1,]    (input)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ain        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     tau        tau                      (extension penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic and extension kernel (exponential)
C     cb         binomial coefficients     
C     dmat, dmati,  dmat0,       working arrays  dp1 times dp1
C     thij,thji                     vector for parameter differences
C     psix,psiy                working memory 
C     sym      asymmetric or symmetric test (logical)
C     
C      implicit logical (a-z)
      integer n,dp1,dp2,i,j,k,l,info,iz,je,ja,nwij 
      logical sym     
      real*8 x(n),y(n),psix(dp2),psiy(dp1),theta(dp1,n),kerns(102),
     1 bi(dp2,n),bin(dp2,n),ain(dp1,n),lam,kernl(102),sigma2(n),
     2 dmat(dp1,dp1),dmati(dp1,dp1),thij(dp1),thji(dp1),cb(dp1,dp1),
     3 bi0(dp2,n),dmat0(dp1,dp1),pij,epij,tau,lambda,tauakt,h,si0(n),
     4 sij,gammaij,wij,z,az,xij,gamma0,ha,ha2,eps,xi,
     5 zj
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      eps=.2*kernl(100)
      do 1 i=1,n
C        loop over design points
         lambda=lam*1.d-2
         if(sym) lambda=2.*lambda
         tauakt=tau*1.d-2
         if(bi(1,i).lt.5.d0*dp1) tauakt=1d10
C  disable extension penalty if estimate is very unstable
         ha=h
         xi=x(i)
C     first search for ja and je
1099     do 1000 j=i,1,-1
            if(x(j).le.xi-ha) goto 1001
            ja=j
1000     continue
1001     do 1002 j=i,n
            if(x(j).gt.xi+ha) goto 1003
            je=j
1002     continue
1003     if(je-ja-1.gt.dp1) goto 1090
             ha=ha*1.25
C            not enough points in neighborhood to estimate parameters
C            increase ha
             goto 1099
1090     continue
         ha2=ha*ha/1d2
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmati and dmat0
C    expand Bi in dmat
         l=1
         do 101 j=1,dp1
            do 101 k=1,dp1
               dmat(j,k)=bi(j+k-1,i)
               if(j.gt.k) then 
                  dmati(j,k)=0.0d0
                  dmat0(j,k)=0.0d0
               else
                  dmati(j,k)=bi(j+k-1,i)
                  dmat0(j,k)=bi0(j+k-1,i)
               end if
101      continue
C     now dmati contains inverse of B_i, dmat0 the inverse of B^(0), 
         if(dp1.eq.1) goto 8999
         call invers(dmat0,dp1,info)
         call invers(dmati,dp1,info)
C     now dmati contains inverse of B_i, dmat0 the inverse of B^(0), 
C     not needed for dp1=1
8999     nwij=0
         do 108 l=1,dp2
            bin(l,i)=0.d0
            bi0(l,i)=0.d0
            if(l.gt.dp1) goto 108
            ain(l,i)=0.d0
108      continue
C      if not enough points with positive weights (counted in nwij)
C      lambda and tau will be increased 
C      (to reduce panalization by stochastic and influence term)
         do 11 j=ja,je
            xij=(x(j)-x(i))
            pij=xij*xij/ha2
            iz=pij
            az=pij-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
              do 111 l=1,dp1
                 thij(l)=theta(l,j)
                 if(sym) thji(l)=theta(l,i)
111              continue
                 if(dp1.eq.1) goto 1111
                 z=1.d0
                 zj=1.d0
                 do 112 l=2,dp1
                 z=-z*xij
                 if(sym) zj=zj*xij
                 do 112 k=1,dp1-l+1
                    thij(k)=thij(k)+cb(k+l-1,k)*z*theta(k+l-1,j)
                 if(sym) thji(k)=thji(k)+cb(k+l-1,k)*zj*theta(k+l-1,i)
112              continue
C     thij contains theta_{j} in the model centered at xi
C     thji contains theta_{i} in the model centered at xj (if needed)
C
C    now get sij
C
1111        do 1210 l=1,dp1
               thij(l)=theta(l,i)-thij(l)
               if(sym) thji(l)=theta(l,j)-thji(l)
1210           continue
C
C   thats the difference between thetai and thetaij
C
            sij=0.0d0
            do 121 l=1,dp1
               do 122 k=1,dp1
                  sij=sij+dmat(k,l)*thij(l)*thij(k)
                  if(sym) sij=sij+bi(l+k-1,j)*thji(l)*thji(k)
122            continue
121         continue
C    now get gamma_ij
            z=1.d0
            do 1171 k=1,dp2
               psix(k)=z
               if(k.gt.dp1) goto 1172
               psiy(k)=z*y(j)
1172           z=z*xij
1171        continue
            gammaij=0.d0
            if(dp1.eq.1.or.tau.le.1.d-8) goto 1202
            gamma0=0.d0
            do 119 l=1,dp1
               do 119 k=1,dp1
                  z=psix(l)*psix(k)
                  gammaij=gammaij+dmati(l,k)*z
                  gamma0=gamma0+dmat0(l,k)*z
119         continue
            gammaij=gammaij*bi(1,i)/(gamma0*si0(i))
            gammaij=dmax1(gammaij-1,0.d0)
C
C     now we have everything to compute  w_{ij}
C
1202        do 1201 k=1,dp2
               bi0(k,i)=bi0(k,i)+psix(k)*epij
1201        continue   
C
C   thats needed to compute gamma0 in the next iteration
C 
            z=sij/lambda+gammaij/tauakt
            if(z.ge.1.d2) goto 11
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
            if(wij.gt.eps) nwij=nwij+1
            wij=wij/sigma2(j)
C        now compute contributions to bi(i),ai(i)  
            do 123 k=1,dp2
               bin(k,i)=bin(k,i)+wij*psix(k)
               if(k.gt.dp1) goto 123
               ain(k,i)=ain(k,i)+wij*psiy(k)
123         continue    
11       continue    
C    this was the j - loop
         if(nwij.ge.dp1)  goto 1
         lambda=lambda*1.25
         tauakt=tauakt*1.25
C    increase lambda and tau to weaken stochastic and influence penalization
         goto 8999
1     continue      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate polynomial aws)
C           and inverse of bi
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mphawsun(n,dp1,dp2,ai,bi,theta,dmat,biinv)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     dpm        number of components in di  (dp1+1)*dp1/2
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi
C     di         inverse of bi     
C     di0         inverse of bi0
C     theta      new parameter estimate
C     dmat       working array
C
C      implicit logical(a-z)
      integer n,dp1,dp2
      real*8 ai(dp1,n),bi(dp2,n),theta(dp1,n),dmat(dp1,dp1),
     1       biinv(dp1,n)
      integer i,j,k,info
      real*8 d
      do 1 i=1,n
         do 11 j=1,dp1
            do 11 k=1,dp1
               if(j.gt.k) then 
                  dmat(j,k)=0.0d0
               else
                  dmat(j,k)=bi(j+k-1,i)
               end if
11       continue
         if(dp1.eq.1) goto 12
         call invers(dmat,dp1,info)
         if(info.ne.0) goto 99 
         goto 13 
12       dmat(1,1)=1.d0/dmat(1,1)  
C      now dmati contains inverse of B_i 
C      now calculate theta as B_i^{-1} A_i
13       do 25 j=1,dp1
            d=0.0d0
            do 26 k=1,dp1
               d=d+dmat(j,k)*ai(k,i)
26          continue
            biinv(j,i)=dmat(j,j)
            theta(j,i)=d
25       continue
99    continue
C     just keep the old estimate
1     continue
      return
      end
