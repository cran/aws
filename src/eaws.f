C
C    Copyright (C) 2002 Weierstrass-Institut für 
C                       Angewandte Analysis und Stochastik (WIAS)
C
C    Author:  Jörg Polzehl
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
C  USA.
C
C  The following routines are part of the aws package and contain  
C  FORTRAN 77 code needed in R function laws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant initialisation for irregular design
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize estimates in multivariate local constant aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine iawsmul(n,y,nn,distm,ih,hinit,bi,ai,kernl)
C    
C     n          number of design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     dist       distances ordered by nearest neighbors
C     ih         initial number of nearest neighbors
C     hinit      initial bandwidth (scaled as number of nearest neighbors)
C     bi         \sum \Psi^T Wi \Psi  (output)
C     ai         \sum \Psi^T Wi Y     (output)
C     kernl      discretized localization kernel 
C
      integer n,i,j,ja,iz,ih,nn(ih,n)
      real*8 bi(n),ai(n),ha2,az,z,epij,y(n),hinit,ha,kernl(102),
     1       distm(ih,n),aii,bii
C     loop over i=1,n
      do 1 i=1,n
         ha=hinit
         ha2=ha*ha/1.d2
         aii=0.d0
         bii=0.d0
         do 11 j=1,ih
            z=distm(j,i)
            if(z.gt.ha) goto 11
            z=z*z/ha2
            iz=z
            az=z-iz
            epij=(1.-az)*kernl(iz+1)+az*kernl(iz+2)
C
C         now calculate contributions to the estimates
C
            ja=nn(j,i)
            aii=aii+y(ja)*epij
            bii=bii+epij
11       continue
         ai(i)=aii
         bi(i)=bii 
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize estimates in multivariate local constant aws (NN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine iawsmnn(n,y,nn,ih,hinit,bi,ai,kernl)
C    
C     n          number of design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     ih         initial number of nearest neighbors
C     hinit      initial bandwidth (scaled as number of nearest neighbors)
C     bi         \sum \Psi^T Wi \Psi  (output)
C     ai         \sum \Psi^T Wi Y     (output)
C     kernl      discretized localization kernel 
C
      integer n,i,j,ja,iz,ih,nn(ih,n),ij
      real*8 bi(n),ai(n),ha2,az,z,epij,y(n),hinit,ha,kernl(102),aii,bii
C     loop over i=1,n
C     use points within (xi-hinit,xi+hinit) but at least dp1 points
      do 1 i=1,n
         ha=hinit-0.9d0
C        not hinit -1 to avoid devision by 0 in case of hinit==ihinit
         ha2=ha*ha
         aii=0.d0
         bii=0.d0
         do 11 j=1,ih
            ja=nn(j,i)
            ij=(j-1)
            z=1.d2*ij*ij/ha2
            iz=z
            az=z-iz
            epij=(1.-az)*kernl(iz+1)+az*kernl(iz+2)
            aii=aii+y(ja)*epij
            bii=bii+epij
11       continue
         ai(i)=aii
         bi(i)=bii
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws for gaussian variables for irregular design
C
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws (bernoulli)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lawsmul(n,y,nn,distm,ihakt,theta,
     1                    bi,bin,ai,ain,lam,h,kernl,kerns,sym)
C    
C     n          number of design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     distm      distances ordered by nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ain        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic  (exponential)
C     sym      asymmetric or symmetric test (logical)
C     
      implicit logical (a-z)
      integer n,i,j,iz,ja,o,nwij,ij,ihakt,nn(ihakt,n)      
      logical sym
      real*8 y(n),theta(n),kerns(102),bi(n),bin(n),ai(n),ain(n),lam,
     1     kernl(102),distm(ihakt,n),epij,lambda,h,wij,z,az,ha2,bii,
     2     thetai

C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      ha2=h*h/1.d2
      lambda=lam*1.d-2
      if(sym) lambda=2*lambda
      do 1 i=1,n
C        loop over design points
         thetai=theta(i)
         bin(i)=0.d0
         ain(i)=0.d0
         do 11 j=1,ihakt
            z=distm(j,i)
            if(z.ge.h) goto 11
            z=z*z/ha2
            iz=z
            az=z-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            ja=nn(j,i)
            bii=bi(i)
            if(sym) bii=bii+bi(ja)
            z=thetai-theta(ja)
            z=bii/lambda*z*z
            if(z.ge.1.d2) goto 11
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
C        now compute contributions to bi(i),ai(i)  
            ain(i)=ain(i)+wij*y(ja)
            bin(i)=bin(i)+wij
11       continue 
C    this was the j - loop
1     continue      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws (NN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lawsmnn(n,y,nn,ihakt,theta,bi,bin,
     1                    ai,ain,lam,h,kernl,kerns,sym)
C    
C     n          number of design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ain        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic  (exponential)
C     sym        asymmetric or symmetric test (logical)
C     
      implicit logical (a-z)
      integer n,i,j,iz,ja,o,nwij,ij,ihakt,nn(ihakt,n)      
      logical sym
      real*8 y(n),theta(n),kerns(102),bi(n),bin(n),ai(n),ain(n),lam,
     1    kernl(102),epij,lambda,h,ha,wij,z,az,ha2,bii,thetai
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      ha=h-.9d0
      ha2=ha*ha/1.d2
      lambda=lam*1.d-2
      if(sym) lambda=2*lambda
      do 1 i=1,n
C        loop over design points
         thetai=theta(i)
         bin(i)=0.d0
         ain(i)=0.d0
         do 11 j=1,ihakt
            ja=nn(j,i)
            bii=bi(i)
            if(sym) bii=bii+bi(ja)
            ij=j-1
            z=ij*ij/ha2
            iz=z
            az=z-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            z=thetai-theta(ja)
            z=bii/lambda*z*z
            if(z.ge.1.d2) goto 11
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
C        now compute contributions to bi(i),ai(i)  
            ain(i)=ain(i)+wij*y(ja)
            bin(i)=bin(i)+wij
11       continue 
C    this was the j - loop
1     continue      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws for bernoulli variables 
C
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant univariate aws for bernoulli variables 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lberuni(y,n,hakt,lamakt,theta,ltheta,lctheta,bi,ai,
     1                   kernl,kerns,sym)
C   
C   y        observed values of regression function
C   n        number of observations
C   hakt     actual bandwidth
C   lamakt   lambda*sigma2 
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   lctheta   log of (1-estimates)    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n
      logical sym
      real*8 y(n),hakt,theta(n),bi(n),ai(n),kernl(102),kerns(102),
     1      lamakt,lam,ltheta(n),lctheta(n)
      integer i,j,ja,je,ih,iz
      real*8 z,z1,wj,az,swj,swjy,bii,thetai,lthetai,lcthetai,cthetai
      ih=hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i=1,n
         bii=bi(i)/lam
         thetai=theta(i)
         cthetai=1.d0-thetai
         lthetai=ltheta(i)
         lcthetai=lctheta(i)
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         swj=0.d0
         swjy=0.d0
         do 2 j=ja,je
C  first stochastic term
            z=bii*(thetai*(lthetai-ltheta(j))+
     1         cthetai*(lcthetai-lctheta(j)))
            if(sym) z=z+bi(j)*(theta(j)*(ltheta(j)-lthetai)+
     1         (1.-theta(j))*(lctheta(j)-lcthetai))/lam
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
            swj=swj+wj
            swjy = swjy+wj*y(j)
2        continue
         ai(i)=swjy
         bi(i)=swj
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform all iterations in local constant univariate aws for bernoulli variables 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gberuni(y,n,hinit,hincr,hmax,lamakt,eta,theta,ltheta, 
     1                   lctheta,bi,ai,kernl,kerns,biold,aiold,sym)
C   
C   y        observed values of regression function
C   n        number of observations
C   hinit    initial bandwidth
C   hincr    factor used to increase bandwiths
C   hmax     maximal bandwidth
C   lamakt   lambda*sigma2 
C   eta      memory parameter
C   theta    estimates    (output)
C   lctheta  log of estimates    
C   lctheta  log of 1-estimates   
C   bi       \sum \Psi^T Wi \Psi  (output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   biold    working array to store old values of bi
C   aiold    working array to store old values of ai
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n
      logical sym
      real*8 y(n),hinit,hincr,hmax,theta(n),bi(n),ai(n),ltheta(n),
     1       kernl(102),kerns(102),lamakt,eta,biold(n),aiold(n),
     2       eps,lctheta(n)
      integer i
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1.d0-eta
      eps=1.d-50
1     do 21 i=1,n
         ltheta(i)=dlog(theta(i)+eps)     
         lctheta(i)=dlog(1.d0-theta(i)+eps)     
21    continue
      call lberuni(y,n,hakt,lamakt,theta,ltheta,lctheta,bi,ai,
     1             kernl,kerns,sym)
      do 31 i=1,n
         ai(i)=onemeta*ai(i)+eta*aiold(i)
         bi(i)=onemeta*bi(i)+eta*biold(i)
         theta(i)=ai(i)/bi(i)
         biold(i)=bi(i)
         aiold(i)=ai(i)
31    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant bivariate aws for bernoulli variables (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lberbi(y,n1,n2,hakt,lamakt,theta,ltheta,lctheta,
     1                  bi,ai,kernl,kerns,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   hakt     actual bandwidth
C   lamakt   lambda*sigma2 
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   lctheta  log of 1-estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2
      logical sym
      real*8 y(n1,n2),hakt,theta(n1,n2),bi(n1,n2),ai(n1,n2),lam,
     1      kernl(102),kerns(102),lamakt,ltheta(n1,n2),lctheta(n1,n2)
      integer i1,i2,j1,j2,ja1,ja2,je1,je2,ih1,ih2,iz
      real*8 z,z1,z2,wj,az,swj,swjy,bii,thetai,hakt2,lthetai,
     1       lcthetai,cthetai
      ih1=hakt
      hakt2=hakt*hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i1=1,n1
         do 1 i2=1,n2
            bii=bi(i1,i2)/lam
            thetai=theta(i1,i2)
            cthetai=1.d0-thetai
            lthetai=ltheta(i1,i2)
            lcthetai=lctheta(i1,i2)
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            swj=0.d0
            swjy=0.d0
            do 2 j1=ja1,je1
               z1=(i1-j1)
               z1=z1*z1
               ih2=dsqrt(hakt2-z1)
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               do 2 j2=ja2,je2
C  first stochastic term
                  z=bii*(thetai*(lthetai-
     1              ltheta(j1,j2))+cthetai*(lcthetai-lctheta(j1,j2)))
                  if(sym) z=z+
     1              bi(j1,j2)*(theta(j1,j2)*(ltheta(j1,j2)-lthetai)+
     2              (1.-theta(j1,j2))*(lctheta(j1,j2)-lcthetai))/lam
                  if(z.ge.1.d2) goto 2
                  iz=z
                  az=z-iz
                  wj=kerns(iz+1)*(1-az)+kerns(iz+2)*az
                  z2=(i2-j2)
                  z=1.d2*(z1+z2*z2)/hakt2
                  if(z.ge.1.d2) goto 2
                  iz=z
                  az=z-iz
                  wj=wj*(kernl(iz+1)*(1-az)+kernl(iz+2)*az)
                  swj=swj+wj
                  swjy = swjy+wj*y(j1,j2)
2           continue
            ai(i1,i2)=swjy
            bi(i1,i2)=swj
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform all iterations in local constant bivariate aws for bernoulli variables (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gberbi(y,n1,n2,hinit,hincr,hmax,lamakt,eta,theta,
     1              ltheta,lctheta,bi,ai,kernl,kerns,biold,aiold,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   hinit    initial bandwidth
C   hincr    factor used to increase bandwiths
C   hmax     maximal bandwidth
C   lamakt   lambda*sigma2 
C   eta      memory parameter
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   lctheta  log of 1-estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   biold    working array to store old values of bi
C   aiold    working array to store old values of ai
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2
      logical sym
      real*8 y(n1,n2),hinit,hincr,hmax,theta(n1,n2),bi(n1,n2),
     1       ai(n1,n2),kernl(102),kerns(102),lamakt,eta,biold(n1,n2),
     2       aiold(n1,n2),ltheta(n1,n2),lctheta(n1,n2),eps
      integer i1,i2
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1-eta
      eps=1.d-50
1     do 21 i1=1,n1
         do 21 i2=1,n2
            ltheta(i1,i2)=dlog(theta(i1,i2)+eps)     
            lctheta(i1,i2)=dlog(1.d0-theta(i1,i2)+eps)     
21    continue
      call lberbi(y,n1,n2,hakt,lamakt,theta,ltheta,lctheta,bi,ai,
     1            kernl,kerns,sym)
      do 31 i1=1,n1
         do 31 i2=1,n2
            ai(i1,i2)=onemeta*ai(i1,i2)+eta*aiold(i1,i2)
            bi(i1,i2)=onemeta*bi(i1,i2)+eta*biold(i1,i2)
            theta(i1,i2)=ai(i1,i2)/bi(i1,i2)
            biold(i1,i2)=bi(i1,i2)
            aiold(i1,i2)=ai(i1,i2)
31    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant trivariate aws for bernoulli variables (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lbertri(y,n1,n2,n3,hakt,lamakt,theta,ltheta,lctheta,
     1                  bi,ai,kernl,kerns,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   n3       number of grid-points in third dimension
C   hakt     actual bandwidth
C   lamakt   lambda*sigma2 
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   lctheta  log of 1-estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2,n3
      logical sym
      real*8 y(n1,n2,n3),hakt,theta(n1,n2,n3),bi(n1,n2,n3),kernl(102),
     1      ai(n1,n2,n3),kerns(102),lamakt,lam,ltheta(n1,n2,n3),
     2      lctheta(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,ja1,ja2,ja3,je1,je2,je3,ih1,ih2,ih3,iz
      real*8 z,z1,z2,z3,wj,az,swj,swjy,bii,thetai,hakt2,lthetai,
     1       lcthetai,cthetai
      ih1=hakt
      hakt2=hakt*hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i1=1,n1
         do 1 i2=1,n2
            do 1 i3=1,n3
               bii=bi(i1,i2,i3)/lam
               thetai=theta(i1,i2,i3)
               cthetai=1.d0-thetai
               lthetai=ltheta(i1,i2,i3)
               lcthetai=lctheta(i1,i2,i3)
               ja1=max0(1,i1-ih1)
               je1=min0(n1,i1+ih1)
               swj=0.d0
               swjy=0.d0
               do 2 j1=ja1,je1
                  z1=(i1-j1)
                  z1=z1*z1
                  ih2=dsqrt(hakt2-z1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
                  do 2 j2=ja2,je2
                     z2=(i2-j2)
                     z2=z2*z2
                     ih3=dsqrt(hakt2-z1-z2)
                     ja3=max0(1,i3-ih3)
                     je3=min0(n3,i3+ih3)
                     do 2 j3=ja3,je3
C  first stochastic term
                        z=bii*(thetai*(lthetai-ltheta(j1,j2,j3))+
     1                     cthetai*(lcthetai-lctheta(j1,j2,j3)))
                        if(sym) z=z+
     1                     bi(j1,j2,j3)*(theta(j1,j2,j3)*
     2                              (ltheta(j1,j2,j3)-lthetai)+
     3                              (1.-theta(j1,j2,j3))*
     4                              (lctheta(j1,j2,j3)-lcthetai))/lam
                        if(z.ge.1.d2) goto 2
                        iz=z
                        az=z-iz
                        wj=kerns(iz+1)*(1-az)+kerns(iz+2)*az
                        z3=(i3-j3)
                        z3=z3*z3
                        z=1.d2*(z1+z2+z3)/hakt2
                        if(z.ge.1.d2) goto 2
                        iz=z
                        az=z-iz
                        wj=wj*(kernl(iz+1)*(1-az)+kernl(iz+2)*az)
                        swj=swj+wj
                        swjy = swjy+wj*y(j1,j2,j3)
2              continue
               ai(i1,i2,i3)=swjy
               bi(i1,i2,i3)=swj
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform all iterations in local constant bivariate aws for bernoulli variables (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gbertri(y,n1,n2,n3,hinit,hincr,hmax,lamakt,eta,theta,
     1               ltheta,lctheta,bi,ai,kernl,kerns,biold,aiold,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   n3       number of grid-points in second dimension
C   hinit    initial bandwidth
C   hincr    factor used to increase bandwiths
C   hmax     maximal bandwidth
C   lamakt   lambda*sigma2 
C   eta      memory parameter
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   lctheta  log of 1-estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   biold    working array to store old values of bi
C   aiold    working array to store old values of ai
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2,n3
      logical sym
      real*8 y(n1,n2,n3),hinit,hincr,hmax,theta(n1,n2,n3),bi(n1,n2,n3),
     1       ai(n1,n2,n3),kernl(102),kerns(102),lamakt,eta,
     2       biold(n1,n2,n3),aiold(n1,n2,n3),ltheta(n1,n2,n3),
     3       lctheta(n1,n2,n3),eps
      integer i1,i2,i3
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1-eta
      eps=1.d-50
1     do 21 i1=1,n1
         do 21 i2=1,n2
            do 21 i3=1,n3
               ltheta(i1,i2,i3)=dlog(theta(i1,i2,i3)+eps)     
               lctheta(i1,i2,i3)=dlog(1.d0-theta(i1,i2,i3)+eps)     
21    continue
      call lbertri(y,n1,n2,n3,hakt,lamakt,theta,ltheta,lctheta,
     1             bi,ai,kernl,kerns,sym)
      do 31 i1=1,n1
         do 31 i2=1,n2
            do 31 i3=1,n3
               ai(i1,i2,i3)=onemeta*ai(i1,i2,i3)+eta*aiold(i1,i2,i3)
               bi(i1,i2,i3)=onemeta*bi(i1,i2,i3)+eta*biold(i1,i2,i3)
               theta(i1,i2,i3)=ai(i1,i2,i3)/bi(i1,i2,i3)
               biold(i1,i2,i3)=bi(i1,i2,i3)
               aiold(i1,i2,i3)=ai(i1,i2,i3)
31    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws (bernoulli)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lbermul(n,y,nn,distm,ihakt,theta,ltheta,lctheta,
     1                    bi,bin,ai,ain,lam,h,kernl,kerns,sym)
C    
C     n          number of design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     distm      distances ordered by nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     ltheta     log of estimates    (input)
C     lctheta    log of (1-estimates)    (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ain        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic  (exponential)
C     sym        asymmetric or symmetric test (logical)
C     
      implicit logical (a-z)
      integer n,i,j,iz,ja,o,nwij,ij,ihakt,nn(ihakt,n)      
      logical sym
      real*8 y(n),theta(n),kerns(102),bi(n),bin(n),ai(n),ain(n),lam,
     1     kernl(102),distm(ihakt,n),epij,lambda,h,wij,z,az,ha2,bii,
     2     ltheta(n),lctheta(n),thetai,cthetai,lthetai,lcthetai

C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      ha2=h*h/1.d2
      lambda=lam*1.d-2
      if(sym) lambda=2*lambda
      do 1 i=1,n
C        loop over design points
         bii=bi(i)/lambda
         thetai=theta(i)
         cthetai=1.d0-thetai
         lthetai=ltheta(i)
         lcthetai=lctheta(i)
         bin(i)=0.d0
         ain(i)=0.d0
         do 11 j=1,ihakt
            z=distm(j,i)
            if(z.ge.h) goto 11
            z=z*z/ha2
            iz=z
            az=z-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            ja=nn(j,i)
            z=bii*(thetai*(lthetai-ltheta(ja))+
     1         cthetai*(lcthetai-lctheta(ja)))
            if(sym) z=z+bi(ja)*(theta(ja)*(ltheta(ja)-lthetai)+
     2                (1.-theta(ja))*(lctheta(ja)-lcthetai))/lambda
            if(z.ge.1.d2) goto 11
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
C        now compute contributions to bi(i),ai(i)  
            ain(i)=ain(i)+wij*y(ja)
            bin(i)=bin(i)+wij
11       continue 
C    this was the j - loop
1     continue      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws (NN)(bernoulli)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lbermnn(n,y,nn,ihakt,theta,ltheta,lctheta,bi,bin,
     1                    ai,ain,lam,h,kernl,kerns,sym)
C    
C     n          number of design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     ltheta     log of estimates    (input)
C     lctheta    log of (1-estimates)    (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ain        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic (exponential)
C     sym        asymmetric or symmetric test (logical)
C     
      implicit logical (a-z)
      integer n,i,j,iz,ja,o,nwij,ij,ihakt,nn(ihakt,n)      
      logical sym
      real*8 y(n),theta(n),kerns(102),bi(n),bin(n),ai(n),ain(n),lam,
     1    kernl(102),epij,lambda,h,ha,wij,z,az,ha2,bii,
     2    ltheta(n),lctheta(n),thetai,cthetai,lthetai,lcthetai
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      ha=h-.9d0
      ha2=ha*ha/1.d2
      lambda=lam*1.d-2
      if(sym) lambda=2*lambda
      do 1 i=1,n
C        loop over design points
         bii=bi(i)/lambda
         thetai=theta(i)
         cthetai=1.d0-thetai
         lthetai=ltheta(i)
         lcthetai=lctheta(i)
         bin(i)=0.d0
         ain(i)=0.d0
         do 11 j=1,ihakt
            ja=nn(j,i)
            ij=j-1
            z=ij*ij/ha2
            iz=z
            az=z-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            z=bii*(thetai*(lthetai-ltheta(ja))+
     1         cthetai*(lcthetai-lctheta(ja)))
            if(sym) z=z+bi(ja)*(theta(ja)*(ltheta(ja)-lthetai)+
     2                (1.-theta(ja))*(lctheta(ja)-lcthetai))/lambda
            if(z.ge.1.d2) goto 11
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
C        now compute contributions to bi(i),ai(i)  
            ain(i)=ain(i)+wij*y(ja)
            bin(i)=bin(i)+wij
11       continue 
C    this was the j - loop
1     continue      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws for poisson counts 
C
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant univariate aws for poisson counts 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpoiuni(y,n,hakt,lamakt,theta,ltheta,bi,ai,
     1                   kernl,kerns,sym)
C   
C   y        observed values of regression function
C   n        number of observations
C   hakt     actual bandwidth
C   lamakt   lambda*sigma2 
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n
      logical sym
      real*8 y(n),hakt,theta(n),bi(n),ai(n),kernl(102),kerns(102),
     1      lamakt,lam,ltheta(n)
      integer i,j,ja,je,ih,iz
      real*8 z,z1,wj,az,swj,swjy,bii,thetai,lthetai
      ih=hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i=1,n
         bii=bi(i)/lam
         thetai=theta(i)
         lthetai=ltheta(i)
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         swj=0.d0
         swjy=0.d0
         do 2 j=ja,je
C  first stochastic term
            z=(theta(j)-thetai*(1.d0-lthetai+ltheta(j)))*bii
            if(sym) z=z+(thetai-theta(j)*(1.d0-ltheta(j)+lthetai))*
     1                  bi(j)/lam
            if(z.ge.1.d2) goto 2
            iz=z
            az=z-iz
            wj=kerns(iz+1)*(1-az)+kerns(iz+2)*az
            z=(i-j)/hakt
            z=z*z*100
            if(z.ge.1.d2) goto 2
            iz=z
            az=z-iz
            wj=wj*(kernl(iz+1)*(1-az)+kernl(iz+2)*az)
            swj=swj+wj
            swjy = swjy+wj*y(j)
2        continue
         ai(i)=swjy
         bi(i)=swj
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform all iterations in local constant univariate aws for poisson counts 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gpoiuni(y,n,hinit,hincr,hmax,lamakt,eta,theta,ltheta, 
     1                   bi,ai,kernl,kerns,biold,aiold,sym)
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
C   aiold    working array to store old values of ai
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n
      logical sym
      real*8 y(n),hinit,hincr,hmax,theta(n),bi(n),ai(n),ltheta(n),
     1       kernl(102),kerns(102),lamakt,eta,biold(n),aiold(n),
     2       eps
      integer i
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1.d0-eta
1     eps=1.d-50
      do 21 i=1,n
         ltheta(i)=dlog(theta(i)+eps)     
21    continue
      call lpoiuni(y,n,hakt,lamakt,theta,ltheta,bi,ai,kernl,kerns,sym)
      do 31 i=1,n
         ai(i)=onemeta*ai(i)+eta*aiold(i)
         bi(i)=onemeta*bi(i)+eta*biold(i)
         theta(i)=ai(i)/bi(i)
         biold(i)=bi(i)
         aiold(i)=ai(i)
31    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant bivariate aws for poisson counts (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpoibi(y,n1,n2,hakt,lamakt,theta,ltheta,bi,ai,
     1                  kernl,kerns,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   hakt     actual bandwidth
C   lamakt   lambda*sigma2 
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2
      logical sym
      real*8 y(n1,n2),hakt,theta(n1,n2),bi(n1,n2),ai(n1,n2),
     1      kernl(102),kerns(102),lamakt,lam,ltheta(n1,n2)
      integer i1,i2,j1,j2,ja1,ja2,je1,je2,ih1,ih2,iz
      real*8 z,z1,z2,wj,az,swj,swjy,bii,thetai,hakt2,lthetai
      ih1=hakt
      hakt2=hakt*hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i1=1,n1
         do 1 i2=1,n2
            bii=bi(i1,i2)/lam
            thetai=theta(i1,i2)
            lthetai=ltheta(i1,i2)
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            swj=0.d0
            swjy=0.d0
            do 2 j1=ja1,je1
               z1=(i1-j1)
               z1=z1*z1
               ih2=dsqrt(hakt2-z1)
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               do 2 j2=ja2,je2
C  first stochastic term
                  z=(theta(j1,j2)-
     1               thetai*(1.d0-lthetai+ltheta(j1,j2)))*bii
                  if(sym) z=z+(thetai-theta(j1,j2)*
     1                     (1.d0-ltheta(j1,j2)+lthetai))*bi(j1,j2)/lam
                  if(z.ge.1.d2) goto 2
                  iz=z
                  az=z-iz
                  wj=kerns(iz+1)*(1-az)+kerns(iz+2)*az
                  z2=(i2-j2)
                  z=1.d2*(z1+z2*z2)/hakt2
                  if(z.ge.1.d2) goto 2
                  iz=z
                  az=z-iz
                  wj=wj*(kernl(iz+1)*(1-az)+kernl(iz+2)*az)
                  swj=swj+wj
                  swjy = swjy+wj*y(j1,j2)
2           continue
            ai(i1,i2)=swjy
            bi(i1,i2)=swj
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform all iterations in local constant bivariate aws for poisson counts (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gpoibi(y,n1,n2,hinit,hincr,hmax,lamakt,eta,theta,
     1                  ltheta,bi,ai,kernl,kerns,biold,aiold,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   hinit    initial bandwidth
C   hincr    factor used to increase bandwiths
C   hmax     maximal bandwidth
C   lamakt   lambda*sigma2 
C   eta      memory parameter
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   biold    working array to store old values of bi
C   aiold    working array to store old values of ai
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2
      logical sym
      real*8 y(n1,n2),hinit,hincr,hmax,theta(n1,n2),bi(n1,n2),
     1       ai(n1,n2),kernl(102),kerns(102),lamakt,eta,biold(n1,n2),
     2       aiold(n1,n2),ltheta(n1,n2),eps
      integer i1,i2
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1-eta
1     eps=1.d-50
      do 21 i1=1,n1
         do 21 i2=1,n2
            ltheta(i1,i2)=dlog(theta(i1,i2)+eps)     
21    continue
      call lpoibi(y,n1,n2,hakt,lamakt,theta,ltheta,bi,ai,
     1            kernl,kerns,sym)
      do 31 i1=1,n1
         do 31 i2=1,n2
            ai(i1,i2)=onemeta*ai(i1,i2)+eta*aiold(i1,i2)
            bi(i1,i2)=onemeta*bi(i1,i2)+eta*biold(i1,i2)
            theta(i1,i2)=ai(i1,i2)/bi(i1,i2)
            biold(i1,i2)=bi(i1,i2)
            aiold(i1,i2)=ai(i1,i2)
31    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant trivariate aws for poisson counts (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpoitri(y,n1,n2,n3,hakt,lamakt,theta,ltheta,bi,ai,
     1                  kernl,kerns,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   n3       number of grid-points in third dimension
C   hakt     actual bandwidth
C   lamakt   lambda*sigma2 
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2,n3
      logical sym
      real*8 y(n1,n2,n3),hakt,theta(n1,n2,n3),bi(n1,n2,n3),kernl(102),
     1      ai(n1,n2,n3),kerns(102),lamakt,lam,ltheta(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,ja1,ja2,ja3,je1,je2,je3,ih1,ih2,ih3,iz
      real*8 z,z1,z2,z3,wj,az,swj,swjy,bii,thetai,hakt2,lthetai
      ih1=hakt
      hakt2=hakt*hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i1=1,n1
         do 1 i2=1,n2
            do 1 i3=1,n3
               bii=bi(i1,i2,i3)/lam
               thetai=theta(i1,i2,i3)
               lthetai=ltheta(i1,i2,i3)
               ja1=max0(1,i1-ih1)
               je1=min0(n1,i1+ih1)
               swj=0.d0
               swjy=0.d0
               do 2 j1=ja1,je1
                  z1=(i1-j1)
                  z1=z1*z1
                  ih2=dsqrt(hakt2-z1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
                  do 2 j2=ja2,je2
                     z2=(i2-j2)
                     z2=z2*z2
                     ih3=dsqrt(hakt2-z1-z2)
                     ja3=max0(1,i3-ih3)
                     je3=min0(n3,i3+ih3)
                     do 2 j3=ja3,je3
C  first stochastic term
                        z=(theta(j1,j2,j3)-thetai*
     1                      (1.d0-lthetai+ltheta(j1,j2,j3)))*bii
                        if(sym) z=z+(thetai-theta(j1,j2,j3)*
     1                            (1.d0-ltheta(j1,j2,j3)+lthetai))*
     2                                bi(j1,j2,j3)/lam
                        if(z.ge.1.d2) goto 2
                        iz=z
                        az=z-iz
                        wj=kerns(iz+1)*(1-az)+kerns(iz+2)*az
                        z3=(i3-j3)
                        z3=z3*z3
                        z=1.d2*(z1+z2+z3)/hakt2
                        if(z.ge.1.d2) goto 2
                        iz=z
                        az=z-iz
                        wj=wj*(kernl(iz+1)*(1-az)+kernl(iz+2)*az)
                        swj=swj+wj
                        swjy = swjy+wj*y(j1,j2,j3)
2              continue
               ai(i1,i2,i3)=swjy
               bi(i1,i2,i3)=swj
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform all iterations in local constant bivariate aws for poisson counts (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gpoitri(y,n1,n2,n3,hinit,hincr,hmax,lamakt,eta,theta,
     1                   ltheta,bi,ai,kernl,kerns,biold,aiold,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   n3       number of grid-points in second dimension
C   hinit    initial bandwidth
C   hincr    factor used to increase bandwiths
C   hmax     maximal bandwidth
C   lamakt   lambda*sigma2 
C   eta      memory parameter
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   biold    working array to store old values of bi
C   aiold    working array to store old values of ai
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2,n3
      logical sym
      real*8 y(n1,n2,n3),hinit,hincr,hmax,theta(n1,n2,n3),bi(n1,n2,n3),
     1       ai(n1,n2,n3),kernl(102),kerns(102),lamakt,eta,
     2       biold(n1,n2,n3),aiold(n1,n2,n3),ltheta(n1,n2,n3),eps
      integer i1,i2,i3
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1-eta
1     eps=1.d-50
      do 21 i1=1,n1
         do 21 i2=1,n2
            do 21 i3=1,n3
               ltheta(i1,i2,i3)=dlog(theta(i1,i2,i3)+eps)     
21    continue
      call lpoitri(y,n1,n2,n3,hakt,lamakt,theta,ltheta,bi,ai,kernl,
     1             kerns,sym)
      do 31 i1=1,n1
         do 31 i2=1,n2
            do 31 i3=1,n3
               ai(i1,i2,i3)=onemeta*ai(i1,i2,i3)+eta*aiold(i1,i2,i3)
               bi(i1,i2,i3)=onemeta*bi(i1,i2,i3)+eta*biold(i1,i2,i3)
               theta(i1,i2,i3)=ai(i1,i2,i3)/bi(i1,i2,i3)
               biold(i1,i2,i3)=bi(i1,i2,i3)
               aiold(i1,i2,i3)=ai(i1,i2,i3)
31    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws (bernoulli)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpoimul(n,y,nn,distm,ihakt,theta,ltheta,
     1                    bi,bin,ai,ain,lam,h,kernl,kerns,sym)
C    
C     n          number of design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     distm      distances ordered by nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     ltheta   log of estimates    (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ain        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic  (exponential)
C     sym        asymmetric or symmetric test (logical)
C     
      implicit logical (a-z)
      integer n,i,j,iz,ja,o,nwij,ij,ihakt,nn(ihakt,n)      
      logical sym
      real*8 y(n),theta(n),kerns(102),bi(n),bin(n),ai(n),ain(n),lam,
     1     kernl(102),distm(ihakt,n),epij,lambda,h,wij,z,az,ha2,bii,
     2     ltheta(n),thetai,lthetai

C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      ha2=h*h/1.d2
      lambda=lam*1.d-2
      if(sym) lambda=2*lambda
      do 1 i=1,n
C        loop over design points
         bii=bi(i)/lambda
         thetai=theta(i)
         lthetai=ltheta(i)
         bin(i)=0.d0
         ain(i)=0.d0
         do 11 j=1,ihakt
            z=distm(j,i)
            if(z.ge.h) goto 11
            z=z*z/ha2
            iz=z
            az=z-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            ja=nn(j,i)
            z=(theta(ja)-thetai*(1.d0-lthetai+ltheta(ja)))*bii
            if(sym) z=z+(thetai-theta(ja)*(1.d0-ltheta(ja)+lthetai))*
     1                                bi(ja)/lambda
            if(z.ge.1.d2) goto 11
C    this is just to avoid numerical problems
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
C        now compute contributions to bi(i),ai(i)  
            ain(i)=ain(i)+wij*y(ja)
            bin(i)=bin(i)+wij
11       continue 
C    this was the j - loop
1     continue      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws (NN)(bernoulli)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpoimnn(n,y,nn,ihakt,theta,ltheta,bi,bin,
     1                    ai,ain,lam,h,kernl,kerns,sym)
C    
C     n          number of design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     ltheta   log of estimates    (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ain        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic  (exponential)
C     sym        asymmetric or symmetric test (logical)
C     
      implicit logical (a-z)
      integer n,i,j,iz,ja,o,nwij,ij,ihakt,nn(ihakt,n)   
      logical sym   
      real*8 y(n),theta(n),kerns(102),bi(n),bin(n),ai(n),ain(n),lam,
     1    kernl(102),epij,lambda,h,ha,wij,z,az,ha2,bii,
     2    ltheta(n),thetai,lthetai
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      ha=h-.9d0
      ha2=ha*ha/1.d2
      lambda=lam*1.d-2
      if(sym) lambda=2*lambda
      do 1 i=1,n
C        loop over design points
         bii=bi(i)/lambda
         thetai=theta(i)
         lthetai=ltheta(i)
         bin(i)=0.d0
         ain(i)=0.d0
         do 11 j=1,ihakt
            ja=nn(j,i)
            ij=j-1
            z=ij*ij/ha2
            iz=z
            az=z-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            z=(theta(ja)-thetai*(1.d0-lthetai+ltheta(ja)))*bii
            if(sym) z=z+(thetai-theta(ja)*(1.d0-ltheta(ja)+lthetai))*
     1                                bi(ja)/lambda
            if(z.ge.1.d2) goto 11
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
C        now compute contributions to bi(i),ai(i)  
            ain(i)=ain(i)+wij*y(ja)
            bin(i)=bin(i)+wij
11       continue 
C    this was the j - loop
1     continue      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws for exponential variates
C
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant univariate aws for exponential variates 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lexpuni(y,n,hakt,lamakt,theta,ltheta,bi,ai,
     1                   kernl,kerns,sym)
C   
C   y        observed values of regression function
C   n        number of observations
C   hakt     actual bandwidth
C   lamakt   lambda*sigma2 
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n
      logical sym
      real*8 y(n),hakt,theta(n),bi(n),ai(n),kernl(102),kerns(102),
     1      lamakt,lam,ltheta(n)
      integer i,j,ja,je,ih,iz
      real*8 z,z1,wj,az,swj,swjy,bii,thetai,lthetai,hakt0
      ih=hakt
      hakt0=hakt*1.d-1
      lam=lamakt*1.d-2
      if(sym) lam=2*lam
      do 1 i=1,n
         bii=bi(i)/lam
         thetai=theta(i)
         lthetai=ltheta(i)
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         swj=0.d0
         swjy=0.d0
         do 2 j=ja,je
C  first stochastic term
            z=(ltheta(j)-lthetai+thetai/theta(j)-1.d0)*bii
            if(sym) z=z+(lthetai-ltheta(j)+theta(j)/thetai-1.d0)*
     1                  bi(j)/lam
            if(z.ge.1.d2) goto 2
            iz=z
            az=z-iz
            wj=kerns(iz+1)*(1.d0-az)+kerns(iz+2)*az
            z=(i-j)/hakt0
            z=z*z
            if(z.ge.1.d2) goto 2
            iz=z
            az=z-iz
            wj=wj*(kernl(iz+1)*(1-az)+kernl(iz+2)*az)
            swj=swj+wj
            swjy = swjy+wj*y(j)
2        continue
         ai(i)=swjy
         bi(i)=swj
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform all iterations in local constant univariate aws for exponential variates 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gexpuni(y,n,hinit,hincr,hmax,lamakt,eta,theta,ltheta, 
     1                   bi,ai,kernl,kerns,biold,aiold,sym)
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
C   aiold    working array to store old values of ai
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n
      logical sym
      real*8 y(n),hinit,hincr,hmax,theta(n),bi(n),ai(n),ltheta(n),
     1       kernl(102),kerns(102),lamakt,eta,biold(n),aiold(n),
     2       eps
      integer i
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1.d0-eta
      eps=1.d-20
1     do 21 i=1,n
         ltheta(i)=dlog(theta(i)+eps)     
21    continue
      call lexpuni(y,n,hakt,lamakt,theta,ltheta,bi,ai,kernl,kerns,sym)
      do 31 i=1,n
         ai(i)=onemeta*ai(i)+eta*aiold(i)
         bi(i)=onemeta*bi(i)+eta*biold(i)
         theta(i)=ai(i)/bi(i)
         biold(i)=bi(i)
         aiold(i)=ai(i)
31    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant bivariate aws for exponential variates (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lexpbi(y,n1,n2,hakt,lamakt,theta,ltheta,bi,ai,
     1                  kernl,kerns,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   hakt     actual bandwidth
C   lamakt   lambda*sigma2 
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2
      logical sym
      real*8 y(n1,n2),hakt,theta(n1,n2),bi(n1,n2),ai(n1,n2),
     1      kernl(102),kerns(102),lamakt,lam,ltheta(n1,n2)
      integer i1,i2,j1,j2,ja1,ja2,je1,je2,ih1,ih2,iz
      real*8 z,z1,z2,wj,az,swj,swjy,bii,thetai,hakt2,lthetai
      ih1=hakt
      hakt2=hakt*hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i1=1,n1
         do 1 i2=1,n2
            bii=bi(i1,i2)/lam
            thetai=theta(i1,i2)
            lthetai=ltheta(i1,i2)
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            swj=0.d0
            swjy=0.d0
            do 2 j1=ja1,je1
               z1=(i1-j1)
               z1=z1*z1
               ih2=dsqrt(hakt2-z1)
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               do 2 j2=ja2,je2
C  first stochastic term
                  z=(ltheta(j1,j2)-lthetai+thetai/theta(j1,j2)-1.d0)*
     1               bii
                  if(sym) z=z+(lthetai-ltheta(j1,j2)+
     1                      theta(j1,j2)/thetai-1)*bi(j1,j2)/lam
                  if(z.ge.1.d2) goto 2
                  iz=z
                  az=z-iz
                  wj=kerns(iz+1)*(1-az)+kerns(iz+2)*az
                  z2=(i2-j2)
                  z=1.d2*(z1+z2*z2)/hakt2
                  if(z.ge.1.d2) goto 2
                  iz=z
                  az=z-iz
                  wj=wj*(kernl(iz+1)*(1-az)+kernl(iz+2)*az)
                  swj=swj+wj
                  swjy = swjy+wj*y(j1,j2)
2           continue
            ai(i1,i2)=swjy
            bi(i1,i2)=swj
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform all iterations in local constant bivariate aws for exponential variates (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gexpbi(y,n1,n2,hinit,hincr,hmax,lamakt,eta,theta,
     1                  ltheta,bi,ai,kernl,kerns,biold,aiold,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   hinit    initial bandwidth
C   hincr    factor used to increase bandwiths
C   hmax     maximal bandwidth
C   lamakt   lambda*sigma2 
C   eta      memory parameter
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   biold    working array to store old values of bi
C   aiold    working array to store old values of ai
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2
      logical sym
      real*8 y(n1,n2),hinit,hincr,hmax,theta(n1,n2),bi(n1,n2),
     1       ai(n1,n2),kernl(102),kerns(102),lamakt,eta,biold(n1,n2),
     2       aiold(n1,n2),ltheta(n1,n2),eps
      integer i1,i2
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1-eta
1     eps=1.d-50
      do 21 i1=1,n1
         do 21 i2=1,n2
            ltheta(i1,i2)=dlog(theta(i1,i2)+eps)     
21    continue
      call lexpbi(y,n1,n2,hakt,lamakt,theta,ltheta,bi,ai,kernl,
     1            kerns,sym)
      do 31 i1=1,n1
         do 31 i2=1,n2
            ai(i1,i2)=onemeta*ai(i1,i2)+eta*aiold(i1,i2)
            bi(i1,i2)=onemeta*bi(i1,i2)+eta*biold(i1,i2)
            theta(i1,i2)=ai(i1,i2)/bi(i1,i2)
            biold(i1,i2)=bi(i1,i2)
            aiold(i1,i2)=ai(i1,i2)
31    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant trivariate aws for exponential variates (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lexptri(y,n1,n2,n3,hakt,lamakt,theta,ltheta,bi,ai,
     1                  kernl,kerns,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   n3       number of grid-points in third dimension
C   hakt     actual bandwidth
C   lamakt   lambda*sigma2 
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2,n3
      logical sym
      real*8 y(n1,n2,n3),hakt,theta(n1,n2,n3),bi(n1,n2,n3),kernl(102),
     1      ai(n1,n2,n3),kerns(102),lamakt,lam,ltheta(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,ja1,ja2,ja3,je1,je2,je3,ih1,ih2,ih3,iz
      real*8 z,z1,z2,z3,wj,az,swj,swjy,bii,thetai,hakt2,lthetai
      ih1=hakt
      hakt2=hakt*hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i1=1,n1
         do 1 i2=1,n2
            do 1 i3=1,n3
               bii=bi(i1,i2,i3)/lam
               thetai=theta(i1,i2,i3)
               lthetai=ltheta(i1,i2,i3)
               ja1=max0(1,i1-ih1)
               je1=min0(n1,i1+ih1)
               swj=0.d0
               swjy=0.d0
               do 2 j1=ja1,je1
                  z1=(i1-j1)
                  z1=z1*z1
                  ih2=dsqrt(hakt2-z1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
                  do 2 j2=ja2,je2
                     z2=(i2-j2)
                     z2=z2*z2
                     ih3=dsqrt(hakt2-z1-z2)
                     ja3=max0(1,i3-ih3)
                     je3=min0(n3,i3+ih3)
                     do 2 j3=ja3,je3
C  first stochastic term
                        z=(ltheta(j1,j2,j3)-lthetai+
     1                     thetai/theta(j1,j2,j3)-1.d0)*bii
                        if(sym) z=z+(lthetai-ltheta(j1,j2,j3)+
     1                     theta(j1,j2,j3)/thetai-1)*bi(j1,j2,j3)/lam
                        if(z.ge.1.d2) goto 2
                        iz=z
                        az=z-iz
                        wj=kerns(iz+1)*(1-az)+kerns(iz+2)*az
                        z3=(i3-j3)
                        z3=z3*z3
                        z=1.d2*(z1+z2+z3)/hakt2
                        if(z.ge.1.d2) goto 2
                        iz=z
                        az=z-iz
                        wj=wj*(kernl(iz+1)*(1-az)+kernl(iz+2)*az)
                        swj=swj+wj
                        swjy = swjy+wj*y(j1,j2,j3)
2              continue
               ai(i1,i2,i3)=swjy
               bi(i1,i2,i3)=swj
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform all iterations in local constant bivariate aws for exponential variates (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gexptri(y,n1,n2,n3,hinit,hincr,hmax,lamakt,eta,theta,
     1                   ltheta,bi,ai,kernl,kerns,biold,aiold,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   n3       number of grid-points in second dimension
C   hinit    initial bandwidth
C   hincr    factor used to increase bandwiths
C   hmax     maximal bandwidth
C   lamakt   lambda*sigma2 
C   eta      memory parameter
C   theta    estimates            (input)
C   ltheta   log of estimates    (input)
C   bi       \sum \Psi^T Wi \Psi  (input/output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kernl    discretized localization kernel 
C   kerns    discretized stochastic kernel 
C   biold    working array to store old values of bi
C   aiold    working array to store old values of ai
C   sym      asymmetric or symmetric test (logical)
C
      implicit logical (a-z)
      integer n1,n2,n3
      logical sym
      real*8 y(n1,n2,n3),hinit,hincr,hmax,theta(n1,n2,n3),bi(n1,n2,n3),
     1       ai(n1,n2,n3),kernl(102),kerns(102),lamakt,eta,
     2       biold(n1,n2,n3),aiold(n1,n2,n3),ltheta(n1,n2,n3),eps
      integer i1,i2,i3
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1-eta
1     eps=1.d-50
      do 21 i1=1,n1
         do 21 i2=1,n2
            do 21 i3=1,n3
               ltheta(i1,i2,i3)=dlog(theta(i1,i2,i3)+eps)     
21    continue
      call lexptri(y,n1,n2,n3,hakt,lamakt,theta,ltheta,bi,ai,kernl,
     1             kerns,sym)
      do 31 i1=1,n1
         do 31 i2=1,n2
            do 31 i3=1,n3
               ai(i1,i2,i3)=onemeta*ai(i1,i2,i3)+eta*aiold(i1,i2,i3)
               bi(i1,i2,i3)=onemeta*bi(i1,i2,i3)+eta*biold(i1,i2,i3)
               theta(i1,i2,i3)=ai(i1,i2,i3)/bi(i1,i2,i3)
               biold(i1,i2,i3)=bi(i1,i2,i3)
               aiold(i1,i2,i3)=ai(i1,i2,i3)
31    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws (bernoulli)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lexpmul(n,y,nn,distm,ihakt,theta,ltheta,
     1                    bi,bin,ai,ain,lam,h,kernl,kerns,sym)
C    
C     n          number of design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     distm      distances ordered by nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     ltheta     log of estimates    (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ain        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic  (exponential)
C     sym        asymmetric or symmetric test (logical)
C     
      implicit logical (a-z)
      integer n,i,j,iz,ja,o,nwij,ij,ihakt,nn(ihakt,n) 
      logical sym     
      real*8 y(n),theta(n),kerns(102),bi(n),bin(n),ai(n),ain(n),lam,
     1     kernl(102),distm(ihakt,n),epij,lambda,h,wij,z,az,ha2,bii,
     2     ltheta(n),thetai,lthetai
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      ha2=h*h/1.d2
      lambda=lam*1.d-2
      if(sym) lambda=2*lambda
      do 1 i=1,n
C        loop over design points
         bii=bi(i)/lambda
         thetai=theta(i)
         lthetai=ltheta(i)
         bin(i)=0.d0
         ain(i)=0.d0
         do 11 j=1,ihakt
            z=distm(j,i)
            if(z.ge.h) goto 11
            z=z*z/ha2
            iz=z
            az=z-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            ja=nn(j,i)
            z=(ltheta(ja)-lthetai+thetai/theta(ja)-1)*bii
            if(sym) z=z+(lthetai-ltheta(ja)+
     1                     theta(ja)/thetai-1)*bi(ja)/lambda
            if(z.ge.1.d2) goto 11
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
C        now compute contributions to bi(i),ai(i)  
            ain(i)=ain(i)+wij*y(ja)
            bin(i)=bin(i)+wij
11       continue 
C    this was the j - loop
1     continue      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws (NN)(bernoulli)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lexpmnn(n,y,nn,ihakt,theta,ltheta,bi,bin,
     1                    ai,ain,lam,h,kernl,kerns,sym)
C    
C     n          number of design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     ltheta   log of estimates    (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ain        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic  (exponential)
C     sym        asymmetric or symmetric test (logical)
C     
      implicit logical (a-z)
      integer n,i,j,iz,ja,o,nwij,ij,ihakt,nn(ihakt,n)      
      logical sym
      real*8 y(n),theta(n),kerns(102),bi(n),bin(n),ai(n),ain(n),lam,
     1    kernl(102),epij,lambda,h,ha,wij,z,az,ha2,bii,
     2    ltheta(n),thetai,lthetai
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      ha=h-.9d0
      ha2=ha*ha/1.d2
      lambda=lam*1.d-2
      if(sym) lambda=2*lambda
      do 1 i=1,n
C        loop over design points
         bii=bi(i)/lambda
         thetai=theta(i)
         lthetai=ltheta(i)
         bin(i)=0.d0
         ain(i)=0.d0
         do 11 j=1,ihakt
            ja=nn(j,i)
            ij=j-1
            z=ij*ij/ha2
            iz=z
            az=z-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            z=(ltheta(ja)-lthetai+thetai/theta(ja)-1)*bii
            if(sym) z=z+(lthetai-ltheta(ja)+
     1                     theta(ja)/thetai-1)*bi(ja)/lambda
            if(z.ge.1.d2) goto 11
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
C        now compute contributions to bi(i),ai(i)  
            ain(i)=ain(i)+wij*y(ja)
            bin(i)=bin(i)+wij
11       continue 
C    this was the j - loop
1     continue      
      return
      end
