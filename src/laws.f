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
C  FORTRAN 77 code needed in R function aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws on a grid
C
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize estimates in local constant univariate aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine iawsuni(y,n,hinit,bi,ai,kern)
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
      real*8 y(n),hinit,bi(n),ai(n),kern(102)
      integer i,j,ja,je,ih,iz
      real*8 z,wj,az,swj,swjy
      ih=hinit
      do 1 i=1,n
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         swj=0.d0
         swjy=0.d0
         do 2 j=ja,je
            z=(i-j)/hinit
            z=z*z*1.d2
            if(z.ge.1.d2) goto 2
            iz=z
            az=z-iz
            wj=kern(iz+1)*(1-az)+kern(iz+2)*az
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
C   Perform one iteration in local constant univariate aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lawsuni(y,n,hakt,lamakt,theta,bi,ai,kernl,kerns,sym)
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
      real*8 y(n),hakt,theta(n),bi(n),ai(n),kernl(102),kerns(102),
     1      lamakt,lam
      integer i,j,ja,je,ih,iz
      real*8 z,z1,wj,az,swj,swjy,bii,thetai
      ih=hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i=1,n
         thetai=theta(i)
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         swj=0.d0
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
C   Perform all iterations in local constant univariate aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gawsuni(y,n,hinit,hincr,hmax,lamakt,eta,theta, 
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
      real*8 y(n),hinit,hincr,hmax,theta(n),bi(n),ai(n),
     1       kernl(102),kerns(102),lamakt,eta,biold(n),aiold(n)
      integer i
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1.d0-eta
1     call lawsuni(y,n,hakt,lamakt,theta,bi,ai,kernl,kerns,sym)
      do 11 i=1,n
         ai(i)=onemeta*ai(i)+eta*aiold(i)
         bi(i)=onemeta*bi(i)+eta*biold(i)
         theta(i)=ai(i)/bi(i)
         biold(i)=bi(i)
         aiold(i)=ai(i)
11    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize estimates in local constant bivariate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine iawsbi(y,n1,n2,hinit,bi,ai,kern)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   hinit    initial bandwidth
C   bi       \sum \Psi^T Wi \Psi  (output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kern     discretized localization kernel 
C
      implicit logical (a-z)
      integer n1,n2
      real*8 y(n1,n2),hinit,bi(n1,n2),ai(n1,n2),kern(102)
      integer i1,j1,ja1,je1,i2,j2,ja2,je2,ih1,ih2,iz
      real*8 z,z1,z2,wj,az,swj,swjy,hinit2
      ih1=hinit
      hinit2=hinit*hinit
      do 1 i1=1,n1
         do 1 i2=1,n2
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            swj=0.d0
            swjy=0.d0
            do 2 j1=ja1,je1
                z1=(i1-j1)
                z1=z1*z1
                ih2=dsqrt(hinit2-z1)
                ja2=max0(1,i2-ih2)
                je2=min0(n2,i2+ih2)
                do 2 j2=ja2,je2
                   z2=(i2-j2)
                   z=1.d2*(z1+z2*z2)/hinit2
                   if(z.ge.1.d2) goto 2
                   iz=z
                   az=z-iz
                   wj=kern(iz+1)*(1-az)+kern(iz+2)*az
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
C   Perform one iteration in local constant bivariate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lawsbi(y,n1,n2,hakt,lamakt,theta,bi,ai,kernl,kerns,
     1                  sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
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
      integer n1,n2
      logical sym
      real*8 y(n1,n2),hakt,theta(n1,n2),bi(n1,n2),ai(n1,n2),
     1      kernl(102),kerns(102),lamakt,lam
      integer i1,i2,j1,j2,ja1,ja2,je1,je2,ih1,ih2,iz
      real*8 z,z1,z2,wj,az,swj,swjy,bii,thetai,hakt2
      ih1=hakt
      hakt2=hakt*hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i1=1,n1
         do 1 i2=1,n2
            thetai=theta(i1,i2)
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
                  bii=bi(i1,i2)
                  if(sym) bii=bii+bi(j1,j2)
                  z=(thetai-theta(j1,j2))
                  z=z*z*bii/lam
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
C   Perform all iterations in local constant bivariate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gawsbi(y,n1,n2,hinit,hincr,hmax,lamakt,eta,theta,
     1                   bi,ai,kernl,kerns,biold,aiold,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
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
      integer n1,n2
      logical sym
      real*8 y(n1,n2),hinit,hincr,hmax,theta(n1,n2),bi(n1,n2),
     1       ai(n1,n2),kernl(102),kerns(102),lamakt,eta,biold(n1,n2),
     2       aiold(n1,n2)
      integer i1,i2
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1-eta
1     call lawsbi(y,n1,n2,hakt,lamakt,theta,bi,ai,kernl,kerns,sym)
      do 11 i1=1,n1
         do 11 i2=1,n2
            ai(i1,i2)=onemeta*ai(i1,i2)+eta*aiold(i1,i2)
            bi(i1,i2)=onemeta*bi(i1,i2)+eta*biold(i1,i2)
            theta(i1,i2)=ai(i1,i2)/bi(i1,i2)
            biold(i1,i2)=bi(i1,i2)
            aiold(i1,i2)=ai(i1,i2)
11    continue
      hakt=hakt*hincr
      if(hakt.le.hmax) goto 1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize estimates in local constant trivariate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine iawstri(y,n1,n2,n3,hinit,bi,ai,kern)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   n3       number of grid-points in third dimension
C   hinit    initial bandwidth
C   bi       \sum \Psi^T Wi \Psi  (output)
C   ai       \sum \Psi^T Wi Y     (output)
C   kern     discretized localization kernel 
C
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 y(n1,n2,n3),hinit,bi(n1,n2,n3),ai(n1,n2,n3),kern(102)
      integer i1,j1,ja1,je1,i2,j2,ja2,je2,i3,j3,ja3,je3,
     1       ih1,ih2,ih3,iz
      real*8 z,z1,z2,z3,wj,az,swj,swjy,hinit2
      ih1=hinit
      hinit2=hinit*hinit
      do 1 i1=1,n1
         do 1 i2=1,n2
            do 1 i3=1,n3
               ja1=max0(1,i1-ih1)
               je1=min0(n1,i1+ih1)
               swj=0.d0
               swjy=0.d0
               do 2 j1=ja1,je1
                  z1=(i1-j1)
                  z1=z1*z1
                  ih2=dsqrt(hinit2-z1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
                  do 2 j2=ja2,je2
                     z2=(i2-j2)
                     z2=z2*z2
                     ih3=dsqrt(hinit2-z1-z2)
                     ja3=max0(1,i3-ih3)
                     je3=min0(n3,i3+ih3)
                     do 2 j3=ja3,je3
                        z3=(i3-j3)
                        z3=z3*z3
                        z=1.d2*(z1+z2+z3)/hinit2
                        if(z.ge.1.d2) goto 2
                        iz=z
                        az=z-iz
                        wj=kern(iz+1)*(1-az)+kern(iz+2)*az
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
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lawstri(y,n1,n2,n3,hakt,lamakt,theta,bi,ai,kernl,
     1                  kerns,sym)
C   
C   y        observed values of regression function
C   n1       number of grid-points in first dimension
C   n2       number of grid-points in second dimension
C   n3       number of grid-points in third dimension
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
      integer n1,n2,n3
      logical sym
      real*8 y(n1,n2,n3),hakt,theta(n1,n2,n3),bi(n1,n2,n3),
     1      ai(n1,n2,n3),kernl(102),kerns(102),lamakt,lam
      integer i1,i2,i3,j1,j2,j3,ja1,ja2,ja3,je1,je2,je3,ih1,ih2,ih3,iz
      real*8 z,z1,z2,z3,wj,az,swj,swjy,bii,thetai,hakt2
      ih1=hakt
      hakt2=hakt*hakt
      lam=lamakt*1d-2
      if(sym) lam=2*lam
      do 1 i1=1,n1
         do 1 i2=1,n2
            do 1 i3=1,n3
               thetai=theta(i1,i2,i3)
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
                        bii=bi(i1,i2,i3)
                        if(sym) bii=bii+bi(j1,j2,j3)
                        z=(thetai-theta(j1,j2,j3))
                        z=z*z*bii/lam
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
C   Perform all iterations in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gawstri(y,n1,n2,n3,hinit,hincr,hmax,lamakt,eta,theta,
     1                   bi,ai,kernl,kerns,biold,aiold,sym)
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
      integer n1,n2,n3
      logical sym
      real*8 y(n1,n2,n3),hinit,hincr,hmax,theta(n1,n2,n3),bi(n1,n2,n3),
     1       ai(n1,n2,n3),kernl(102),kerns(102),lamakt,eta,
     2       biold(n1,n2,n3),aiold(n1,n2,n3)
      integer i1,i2,i3
      real*8 hakt,onemeta
      hakt=hinit*hincr
      onemeta=1-eta
1     call lawstri(y,n1,n2,n3,hakt,lamakt,theta,bi,ai,kernl,kerns,sym)
      do 11 i1=1,n1
         do 11 i2=1,n2
            do 11 i3=1,n3
            ai(i1,i2,i3)=onemeta*ai(i1,i2,i3)+eta*aiold(i1,i2,i3)
            bi(i1,i2,i3)=onemeta*bi(i1,i2,i3)+eta*biold(i1,i2,i3)
            theta(i1,i2,i3)=ai(i1,i2,i3)/bi(i1,i2,i3)
            biold(i1,i2,i3)=bi(i1,i2,i3)
            aiold(i1,i2,i3)=ai(i1,i2,i3)
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
      subroutine ipawsuni(n,dp1,dp2,x,y,hinit,bi,ai,theta,
     1                    kernl,dmat)
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
      integer n,dp1,dp2,i,j,k,info,je,ja,m,o,iz
      real*8 bi(dp2,n),ai(dp1,n),theta(dp1,n),dmat(dp1,dp1),ha2,d
      real*8 x(n),xij,xi,z,epij,y(n),hinit,ha,kernl(102),az
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
      subroutine lpawsuni(n,dp1,dp2,x,y,theta,bi,bin,si0,bi0,
     1    ain,lam,tau,h,kernl,kerns,cb,dmat,dmati,dmat0,
     2    thij,thji,psix,psiy,sym)
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
      integer n,dp1,dp2,i,j,k,l,info,iz,je,ja,o,nwij,pm  
      logical sym     
      real*8 x(n),y(n),psix(dp2),psiy(dp1),theta(dp1,n),kerns(102),
     1 bi(dp2,n),bin(dp2,n),ain(dp1,n),lam,kernl(102),
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
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpawsuni(n,dp1,dp2,ai,bi,theta,dmat)
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
      real*8 ai(dp1,n),bi(dp2,n),theta(dp1,n),dmat(dp1,dp1)
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
            theta(j,i)=d
25       continue
99    continue
C     just keep the old estimate
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     bivariate local polynomial aws 
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize estimates in bivariate local polynomial aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ipawsbi(n1,n2,dp1,dp2,y,hinit,bi,ai,theta,
     1                    kernl,dmat,si,siy,ind)
C    
C     n1         number of points in first dimension
C     n2         number of points in second dimension
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     y          observed values at design points
C     hinit      initial bandwidth
C     bi         \sum \Psi^T Wi \Psi  (output)
C     ai         \sum \Psi^T Wi Y     (output)
C     theta      initial estimates    (output)
C     kernl      discretized localization kernel 
C     dmat       working array
C     si         working array  for \sum \Psi^T Wi \Psi
C     siy        working array  for \sum \Psi^T Wi Y 
C
      implicit logical (a-z)
      integer n1,n2,dp1,dp2,i1,i2,j1,j2,k,info,je1,ja1,je2,ja2,iz,
     1        ih1,ih2,j,l,ind(dp1,dp1)
      real*8 y(n1,n2),bi(dp2,n1,n2),ai(dp1,n1,n2),
     1      si(dp2),siy(dp1),kernl(102),dmat(dp1,dp1),theta(dp1,n1,n2),
     2      z11,z12,z22,d,z,epij,hinit,ha,ha2,az,z1,z2,epijy
C     loop over i1=1,n1 and i2=1,n2
      do 1 i1=1,n1
         do 1 i2=1,n2
            ha=hinit
C  first fill si and siy with 0's
1099        ha2=ha*ha
            ih1=ha
            do 1001 j=1,dp2
            si(j)=0.d0
            if(j.le.dp1) siy(j)=0.d0
1001        continue
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            do 11 j1=ja1,je1
                z1=(i1-j1)
                z11=z1*z1
                ih2=dsqrt(ha2-z11)
                ja2=max0(1,i2-ih2)
                je2=min0(n2,i2+ih2)
                do 11 j2=ja2,je2
                   z2=(i2-j2)
                   z22=z2*z2
                   z12=z1*z2
                   z=1.d2*(z11+z22)/ha2
                   iz=z
                   az=z-iz
                   epij=(1.-az)*kernl(iz+1)+az*kernl(iz+2)
                   epijy=y(j1,j2)*epij
                   si(1)=si(1)+epij
                   siy(1)=siy(1)+epijy
                   if(dp1.le.1) goto 11
                   si(2)=si(2)-z1*epij
                   si(3)=si(3)-z2*epij
                   si(4)=si(4)+z11*epij
                   si(5)=si(5)+z12*epij
                   si(6)=si(6)+z22*epij
                   siy(2)=siy(2)-z1*epijy
                   siy(3)=siy(3)-z2*epijy
                   if(dp1.le.3) goto 11
                   si(7)=si(7)-z11*z1*epij
                   si(8)=si(8)-z11*z2*epij
                   si(9)=si(9)-z1*z22*epij
                   si(10)=si(10)-z2*z22*epij
                   si(11)=si(11)+z11*z11*epij
                   si(12)=si(12)+z11*z12*epij
                   si(13)=si(13)+z11*z22*epij
                   si(14)=si(14)+z12*z22*epij
                   si(15)=si(15)+z22*z22*epij
                   siy(4)=siy(4)+z11*epijy
                   siy(5)=siy(5)+z12*epijy
                   siy(6)=siy(6)+z22*epijy
11          continue
C        prepare matrix 
         do 13 k=1,dp1
            do 13 j=1,dp1
            if(j.gt.k) then
               dmat(j,k)=0.d0
            else
               dmat(j,k)=si(ind(j,k))
            endif
13       continue
C     compute choleski decomposition
         call invers(dmat,dp1,info)
         if(info.eq.0)  goto 14
            ha=ha*1.25
            goto 1099
C     compute inverse 
C     now dmat contains inverse of B_i 
C     now calculate theta as B_i^{-1} Z_i
14       l=1   
         do 15 j=1,dp1
            d=0.0d0
            do 16 k=1,dp1
               d=d+dmat(j,k)*siy(k)  
16          continue
            theta(j,i1,i2)=d
            ai(j,i1,i2)=siy(j)
15       continue
         do 20 j=1,dp2
            bi(j,i1,i2)=si(j)
20       continue
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in bivariate local polynomial aws (gridded) 
C
C   p > 0   only  !!!! 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpawsbi(n1,n2,dp1,dp2,y,theta,bi,bin,bi0,ain,lam,tau,
     1    h,kernl,kerns,dmat,dmati,dmat0,thij,thji,psix,
     2    si,si0,siy,sym,ind)
C
C     n1         number of points in first dimension
C     n2         number of points in second dimension
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     y          observed values at design points
C     theta      estimates from step (k-1)   (input)
C     bi         \sum \Psi^T Wi \Psi  from step (k-1)
C     bin        \sum \Psi^T Wi \Psi  (output)
C     bi0        \sum \Psi^T Wi0 \Psi  from step (k-1)
C     ain        \sum \Psi^T Wi Y     (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     tau        tau                      (extension penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic and extension kernel (exponential)
C     dmat, dmati, dmat0       working arrays  dp1 times dp1
C     thij, thji               vector for parameter differences
C     psix                     working memory for Psij
C     si         working array  for \sum \Psi^T Wi \Psi
C     si0        working array  for \sum \Psi^T Wi0 \Psi
C     siy        working array  for \sum \Psi^T Wi Y 
C     sym        asymmetric or symmetric test (logical)
C
C      implicit logical (a-z)
      integer n1,n2,dp1,dp2,i1,i2,j1,j2,k,l,info,je1,ja1,je2,ja2,
     1        iz,j,ih1,ih2,ind(dp1,dp1),m
      logical sym
      real*8 bi(dp2,n1,n2),bi0(dp2,n1,n2),bin(dp2,n1,n2),
     1       ain(dp1,n1,n2),kerns(102),
     2       theta(dp1,n1,n2),dmat(dp1,dp1),siy(dp1),dmati(dp1,dp1),
     3       dmat0(dp1,dp1),thij(dp1),
     4       thji(dp1),y(n1,n2),h,lam,tau,kernl(102),psix(dp2),
     5       si(dp2),si0(dp2),sii,wijy,d,z,epij,z11,z12,z22,ha2,pm,
     6       lambda,tauakt,gammaij,gammaji,wij,gamma0,gamma0j,s0i,xi,
     7       az,z1,z2
C     
C     in case of dp1==1  lawsbi should be called  (p=0)
C
      ha2=h*h
      ih1=h
      do 1 i1=1,n1
         do 1 i2=1,n2
            lambda=lam*1.d-2
            if(sym) lambda=2*lambda
            tauakt=tau*1.d-2
C  disable extension penalty if estimate is very unstable
            if(bi(1,i1,i2).lt.5.d0*dp1) tauakt=1d10
C  first fill si and siy with 0's
C  fields are used to sum components of ain, bin and bi0
8999        do 1001 j=1,dp2
               si(j)=0.d0
               si0(j)=0.d0
               if(j.le.dp1) siy(j)=0.d0
1001        continue
         do 1002 j=1,dp1
            do 1002 k=1,dp1
               m=ind(j,k)
               dmati(j,k)=bi(m,i1,i2)
               if(j.gt.k) then 
                  dmat(j,k)=0.0d0
                  dmat0(j,k)=0.0d0
               else
                  dmat(j,k)=bi(m,i1,i2)
                  dmat0(j,k)=bi0(m,i1,i2)
               end if
1002     continue
C        store     Bi in dmati
C        generate  Bi0^{-1} in dmat0  
C        generate  Bi^{-1} in dmat
         if(dp1.eq.1) goto 1003
         call invers(dmat0,dp1,info)
         call invers(dmat,dp1,info)
C        thats  Bi0^{-1} in dmat0  
C        thats  Bi^{-1} in dmat
1003     s0i=bi0(1,i1,i2)
         sii=bi(1,i1,i2)
C
C     Prepare for loop over j's
C
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            do 101 j1=ja1,je1
                z1=(i1-j1)
                z11=z1*z1
                ih2=dsqrt(ha2-z11)
                ja2=max0(1,i2-ih2)
                je2=min0(n2,i2+ih2)
                do 102 j2=ja2,je2
                   z2=(i2-j2)
                   z22=z2*z2
                   z12=z1*z2
C          first compute location part
                   z=1.d2*(z11+z22)/ha2
                   iz=z
                   az=z-iz
                   epij=(1.-az)*kernl(iz+1)+az*kernl(iz+2)
C          this is the location penalty, now fill si0
                   si0(1)=si0(1)+epij
                   si0(2)=si0(2)-z1*epij
                   si0(3)=si0(3)-z2*epij
                   si0(4)=si0(4)+z11*epij
                   si0(5)=si0(5)+z12*epij
                   si0(6)=si0(6)+z22*epij
                   if(dp1.le.3) goto 10201
                   si0(7)=si0(7)-z11*z1*epij
                   si0(8)=si0(8)-z11*z2*epij
                   si0(9)=si0(9)-z1*z22*epij
                   si0(10)=si0(10)-z2*z22*epij
                   si0(11)=si0(11)+z11*z11*epij
                   si0(12)=si0(12)+z11*z12*epij
                   si0(13)=si0(13)+z11*z22*epij
                   si0(14)=si0(14)+z12*z22*epij
                   si0(15)=si0(15)+z22*z22*epij
10201              continue
C          now fill psix 
                   psix(1)=1
                   psix(2)=-z1
                   psix(3)=-z2
                   psix(4)=z11
                   psix(5)=z12
                   psix(6)=z22
                   if(dp1.le.3) goto 10202
                   psix(7)=-z11*z1
                   psix(8)=-z11*z2
                   psix(9)=-z1*z22
                   psix(10)=-z2*z22
                   psix(11)=z11*z11
                   psix(12)=z11*z12
                   psix(13)=z11*z22
                   psix(14)=z12*z22
                   psix(15)=z22*z22
10202              continue
C         now get gamma
                   gammaij=0.d0
                   if(tau.le.1.d-8) goto 1202
                   gamma0=0.d0
                   do 1021 l=1,dp1
                      do 1021 k=1,dp1
                         z=psix(l)*psix(k)
                         gammaij=gammaij+dmat(l,k)*z
                         gamma0=gamma0+dmat0(l,k)*z
1021               continue
                   gammaij=gammaij*sii/(gamma0*s0i)
                   gammaij=dmax1(gammaij-1.d0,0.d0)
C           thats the main part of extension penalty
C
C           now translate thetaj into model centered in xi
C
1202               thij(1)=theta(1,j1,j2)
                   thij(1)=thij(1)+theta(2,j1,j2)*z1+theta(3,j1,j2)*z2
                   thij(2)=theta(2,j1,j2)
                   thij(3)=theta(3,j1,j2)
                   if(dp1.eq.3) goto 10211
                   thij(1)=thij(1)+theta(4,j1,j2)*z11+
     +                        theta(5,j1,j2)*z12+theta(6,j1,j2)*z22
                   thij(2)=thij(2)+theta(5,j1,j2)*z2+
     +                             2.d0*theta(4,j1,j2)*z1
                   thij(3)=thij(3)+theta(5,j1,j2)*z1+
     +                             2.d0*theta(6,j1,j2)*z2
                   thij(4)=theta(4,j1,j2)
                   thij(5)=theta(5,j1,j2)
                   thij(6)=theta(6,j1,j2)
10211              continue
                   if(sym) then
                      thji(1)=theta(1,i1,i2)
                      thji(1)=thji(1)-theta(2,i1,i2)*z1-
     +                        theta(3,i1,i2)*z2
                      thji(2)=theta(2,i1,i2)
                      thji(3)=theta(3,i1,i2)
                      if(dp1.eq.3) goto 10212
                      thji(1)=thji(1)+theta(4,i1,i2)*z11+
     +                        theta(5,i1,i2)*z12+theta(6,i1,i2)*z22
                      thji(2)=thji(2)-theta(5,i1,i2)*z2-
     +                             2.d0*theta(4,i1,i2)*z1
                      thji(3)=thji(3)-theta(5,i1,i2)*z1-
     +                             2.d0*theta(6,i1,i2)*z2
                      thji(4)=theta(4,i1,i2)
                      thji(5)=theta(5,i1,i2)
                      thji(6)=theta(6,i1,i2)
10212                 continue
                   endif
C  
C           get difference of thetas
C
                   do 1022 l=1,dp1
                      thij(l)=theta(l,i1,i2)-thij(l)
                      if(sym) thji(l)=theta(l,j1,j2)-thji(l)
1022               continue
C
C           get stochastic penalty
C
                   d=0.d0
                   do 1023 l=1,dp1
                      do 1023 k=1,dp1
                        d=d+dmati(k,l)*thij(l)*thij(k)
                        if(sym) d=d+bi(ind(k,l),j1,j2)*thji(l)*thji(k)
1023               continue
C
C           now compute weights 
C           stochastic and extension penalty can be added because kernel is exp
C
                   z=d/lambda+gammaij/tauakt
                   if(z.ge.1.d2) goto 102
                   iz=z
                   az=z-iz
                   wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
C           now compute contributions to bi(i),ai(i)  
                   wijy=y(j1,j2)*wij
                   si(1)=si(1)+wij
                   siy(1)=siy(1)+wijy
                   si(2)=si(2)-z1*wij
                   si(3)=si(3)-z2*wij
                   si(4)=si(4)+z11*wij
                   si(5)=si(5)+z12*wij
                   si(6)=si(6)+z22*wij
                   siy(2)=siy(2)-z1*wijy
                   siy(3)=siy(3)-z2*wijy
                   if(dp1.le.3) goto 102
                   si(7)=si(7)-z11*z1*wij
                   si(8)=si(8)-z11*z2*wij
                   si(9)=si(9)-z1*z22*wij
                   si(10)=si(10)-z2*z22*wij
                   si(11)=si(11)+z11*z11*wij
                   si(12)=si(12)+z11*z12*wij
                   si(13)=si(13)+z11*z22*wij
                   si(14)=si(14)+z12*z22*wij
                   si(15)=si(15)+z22*z22*wij
                   siy(4)=siy(4)+z11*wijy
                   siy(5)=siy(5)+z12*wijy
                   siy(6)=siy(6)+z22*wijy
102             continue
101          continue
C        prepare matrix to test for singularity of Bi 
C        this should be changed to SVD at some point
12           do 121 k=1,dp1
                do 121 j=1,dp1
                   if(j.gt.k) then
                      dmat(j,k)=0.d0
                   else
                      dmat(j,k)=si(ind(j,k))
                   endif
121          continue
C
C     compute choleski decomposition
C
             call invers(dmat,dp1,info)
             if(info.eq.0) goto 20
C
C          if singular relax stochastic and extension penalty 
C
                lambda=1.5*lambda
                tauakt=1.5*tauakt
                goto 8999
C     
C     now fill ain, bin and bi0
C
20       do 201 j=1,dp1
            ain(j,i1,i2)=siy(j)
201      continue
         do 202 j=1,dp2
            bin(j,i1,i2)=si(j)
            bi0(j,i1,i2)=si0(j)
202      continue
1     continue      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (bivariate case)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpawsbi(n,dp1,dp2,ai,bi,theta,dmat,ind)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi    
C     theta      new parameter estimate
C     dmat       working arrays
C
C      implicit logical (a-z)
      integer n,dp1,dp2
      logical sym
      real*8 ai(dp1,n),bi(dp2,n),theta(dp1,n),dmat(dp1,dp1)
      integer i,j,k,info,l,ind(dp1,dp1)
      real*8 d
      do 1 i=1,n
         do 11 k=1,dp1
            do 11 j=1,dp1
               if(j.gt.k) then
                  dmat(j,k)=0.d0
               else
                  dmat(j,k)=bi(ind(j,k),i)
               endif
11          continue
         call invers(dmat,dp1,info)
         if(info.gt.0) goto 99  
C     now dmat contains inverse of B_i 
C     now calculate theta as B_i^{-1} A_i
         do 20 j=1,dp1
            d=0.0d0
            do 21 k=1,dp1
               d=d+dmat(j,k)*ai(k,i)  
21          continue
            theta(j,i)=d
20       continue
99    continue
C     just keep the old estimate
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Multivariate local polynomial ( p=0 or 1 ) aws ( p=0 or 1 )
C
C   Nongridded design and Nearest Neighbor approach
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize estimates in multivariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ipawsmul(n,px,dp1,dp2,x,y,nn,distm,ihinit,hinit,
     1                    bi,ai,theta,kernl,dmat,xij,info)
C    
C     n          number of design points
C     px         dimension of xi
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     x          design points (px,n)
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     dist       distances ordered by nearest neighbors
C     ihinit     initial number of nearest neighbors
C     hinit      initial bandwidth (scaled as number of nearest neighbors)
C     bi         \sum \Psi^T Wi \Psi  (output)
C     ai         \sum \Psi^T Wi Y     (output)
C     theta      initial estimates    (output)
C     kernl      discretized localization kernel 
C     dmat       working array
C     xij        working vector 
C     info       error indicator
C
      integer n,dp1,dp2,i,j,k,l,m,info,ja,o,iz,px,ihinit,
     1        nn(ihinit,n),ij
      real*8 bi(dp2,n),ai(dp1,n),theta(dp1,n),dmat(dp1,dp1),ha2,az,d,
     1    x(px,n),xij(dp1),xi,z,epij,y(n),hinit,ha,kernl(102),
     2    distm(ihinit,n)
C     loop over i=1,n
C     use points within (xi-hinit,xi+hinit) but at least dp1 points
      do 1 i=1,n
         ha=hinit
C        first fill bi(,i) and ai(,i)  
         ha2=ha*ha/1.d2
         do 11 j=1,ihinit
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
            z=1.d0
            xij(1)=1.d0
            ai(1,i)=ai(1,i)+y(ja)*epij
            if(dp1.eq.1) goto 112
            do 111 k=1,px
               z=x(k,i)-x(k,ja)
               xij(k+1)=-z
               ai(k+1,i)=ai(k+1,i)-z*y(ja)*epij
111            continue
112         m=1
            do 12 k=1,dp1
               do 12 l=1,k
                  bi(m,i)=bi(m,i)+xij(k)*xij(l)*epij
                  m=m+1
12         continue
11      continue
C       expand bi as p times p matrix
        m=1
        do 13 k=1,dp1
           do 13 l=1,k
              dmat(l,k)=bi(m,i)
              m=m+1
13      continue
C     compute inverse of bi by choleski decomposition
         if(dp1.eq.1) goto 15
         call invers(dmat,dp1,info)
         if(info.gt.0) goto 9999
C        ihinit to samll
         goto 16
15       dmat(1,1)=1.d0/dmat(1,1)
         info=0
C     now dmat contains inverse of B_i 
C     now calculate theta as B_i^{-1} Z_i
16      do 20 j=1,dp1
           d=0.0d0
           do 21 k=1,dp1
              d=d+dmat(j,k)*ai(k,i)  
21         continue
           theta(j,i)=d
20      continue
1     continue
9999  return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local polynomial aws 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpawsmul(n,px,dp1,dp2,x,y,nn,distm,ihakt,theta,
     1                    bi,bin,bi0,ai,ain,lam,tau,h,kernl,kerns,
     2                    dmat,dmati,dmat0,thij,psix,psiy,sym)
C    
C     n          number of design points
C     px         dimension of x
C     dp1        number of parameters  (1+p*px)
C     dp2        number of components in bi  (dp1*(dp1+1)/2)
C     x          design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     distm      distances ordered by nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ai         \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     tau        tau                      (extension penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic and extension kernel (exponential)
C     dmat, dmati, dmat0       working arrays  dp1 times dp1
C     thij                     vector for parameter differences
C     psix,psiy,xij                working memory 
C     sym        asymmetric or symmetric test (logical)
C     
      implicit logical (a-z)
      integer n,dp1,dp2,px,i,j,k,l,info,iz,je,ja,o,nwij,ij,
     1        ihakt,m,nn(ihakt,n)      
      logical sym
      real*8 x(px,n),y(n),psix(dp1),psiy(dp1),theta(dp1,n),kerns(102),
     1 bi(dp2,n),bin(dp2,n),ai(dp1,n),ain(dp1,n),lam,kernl(102),
     2 dmat(dp1,dp1),dmati(dp1,dp1),thij(dp1),distm(ihakt,n),
     3 bi0(dp2,n),dmat0(dp1,dp1),pij,epij,tau,lambda,tauakt,h,
     4 d,gammaij,wij,z,az,xij,gamma0,s0i,ha,ha2,eps
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      eps=.2*kernl(100)
      ha=h
      ha2=ha*ha/1.d2
      do 1 i=1,n
C        loop over design points
         lambda=lam*1.d-2
         if(sym) lambda=2*lambda
         tauakt=tau*1.d-2
         if(bi(1,i).lt.5.d0*dp1) tauakt=1d10
C  disable extension penalty if estimate is very unstable
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmati and dmat0
C    expand Bi in dmat
         s0i=bi0(1,i)
         m=1
         do 101 k=1,dp1
            do 101 l=1,k
               dmat(l,k)=bi(m,i)
               dmat(k,l)=bi(m,i)
               dmati(l,k)=bi(m,i)
               dmat0(l,k)=bi0(m,i)
               m=m+1
101      continue
         do 105 l=1,dp2
         bi0(l,i)=0.d0
105      continue
C     compute inverse of dmat0 and dmati
         if(dp1.eq.1) goto 8999
         call invers(dmat0,dp1,info)
         call invers(dmati,dp1,info)
C     now dmati contains inverse of B_i, dmat0 the inverse of B^(0), 
C     not needed for dp1=1
8999     nwij=0
         do 108 l=1,dp2
            bin(l,i)=0.d0
            if(l.gt.dp1) goto 108
            ain(l,i)=0.d0
108      continue
C      if not enough points with positive weights (counted in nwij)
C      lambda and tau will be increased 
C      (to reduce panalization by stochastic and influence term)
         do 11 j=1,ihakt
            z=distm(j,i)
            if(z.ge.ha) goto 11
            pij=z*z/ha2
            iz=pij
        if(iz.gt.1.d2) goto 11
            az=pij-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
              ja=nn(j,i)
              do 111 l=1,dp1
                 thij(l)=theta(l,ja)
111              continue
                 psix(1)=1
                 psiy(1)=y(ja)
                 gammaij=0.d0
                 if(dp1.eq.1) goto 1202
                 z=1.d0
                 do 112 l=2,dp1
                    xij=(x(l-1,i)-x(l-1,ja))
                    thij(1)=thij(1)+xij*theta(l,ja)
                    psix(l)=-xij
                    psiy(l)=-xij*y(ja)
112              continue
C     thij contains theta_{j} in the model centered at xi
C    now get gamma_ij
            gamma0=0.d0
            do 119 l=1,dp1
               do 119 k=1,dp1
                  z=psix(l)*psix(k)
                  gammaij=gammaij+dmati(l,k)*z
                  gamma0=gamma0+dmat0(l,k)*z
119         continue
            gammaij=gammaij*bi(1,i)/(gamma0*s0i)
            gammaij=dmax1(gammaij-1,0.d0)/tauakt
C
C     now we have everything to compute  w_{ij}
C
            m=1       
            do 1201 k=1,dp1
               do 1201 l=1,k
                  bi0(m,i)=bi0(m,i)+psix(k)*psix(l)*epij
                  m=m+1
1201        continue   
C
C   thats needed to compute gamma0 in the next iteration
C
1202        d=0.0d0
            do 1210 l=1,dp1
               thij(l)=theta(l,i)-thij(l)
1210           continue
C
C   thats the difference between thetai and thetaij
C
            do 121 l=1,dp1
               do 122 k=1,dp1
                  d=d+dmat(k,l)*thij(l)*thij(k)
122            continue
121         continue
            z=d/lambda+gammaij
            if(z.ge.1.d2) goto 11
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
            if(wij.gt.eps) nwij=nwij+1
C        now compute contributions to bi(i),ai(i)  
            m=1
            do 123 k=1,dp1
               ain(k,i)=ain(k,i)+wij*psiy(k)
               do 123 l=1,k
                  bin(m,i)=bin(m,i)+wij*psix(k)*psix(l)
                  m=m+1
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
C      Generate estimates from ai and bi (multivariate polynomial aws)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpawsmul(n,dp1,dp2,ai,bi,theta,dmat)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi    
C     theta      new parameter estimate
C     dmat       working array
C
      integer n,dp1,dp2
      real*8 ai(dp1,n),bi(dp2,n),theta(dp1,n),dmat(dp1,dp1)
      integer i,j,l,info
      real*8 d
      do 1 i=1,n
         l=1
         m=1
         do 11 j=1,dp1
            do 11 k=1,j
               dmat(k,j)=bi(m,i)
               m=m+1
11       continue
         if(dp1.eq.1) goto 12
         call invers(dmat,dp1,info)
         if(info.ne.0) goto 99 
         goto 13 
12       dmat(1,1)=1.d0/dmat(1,1)  
C     now dmati contains inverse of B_i 
C     now calculate theta as B_i^{-1} A_i
13       do 25 j=1,dp1
            d=0.0d0
            do 26 k=1,dp1
               d=d+dmat(j,k)*ai(k,i)  
26          continue
            theta(j,i)=d
25       continue
99    continue
C     just keep the old estimate
1     continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Multivariate local polynomial ( p=0 or 1 ) aws ( p=0 or 1 )
C
C   Nongridded design and Nearest Neighbor approach
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize estimates in multivariate local polynomial aws (NN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ipawsmnn(n,px,dp1,dp2,x,y,nn,ihinit,hinit,bi,ai,
     1                    theta,kernl,dmat,xij,info)
C    
C     n          number of design points
C     px         dimension of xi
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     x          design points (px,n)
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     ihinit     initial number of nearest neighbors
C     hinit      initial bandwidth (scaled as number of nearest neighbors)
C     bi         \sum \Psi^T Wi \Psi  (output)
C     ai         \sum \Psi^T Wi Y     (output)
C     theta      initial estimates    (output)
C     kernl      discretized localization kernel 
C     dmat       working array
C     xij        working vector 
C     info       error indicator
C
      integer n,dp1,dp2,i,j,k,l,m,info,ja,o,iz,px,ihinit,
     1        nn(ihinit,n),ij
      real*8 bi(dp2,n),ai(dp1,n),theta(dp1,n),dmat(dp1,dp1),ha2,az,d,
     1    x(px,n),xij(dp1),xi,z,epij,y(n),hinit,ha,kernl(102)
C     loop over i=1,n
C     use points within (xi-hinit,xi+hinit) but at least dp1 points
      do 1 i=1,n
         ha=hinit-0.9d0
C        not hinit -1 to avoid devision by 0 in case of hinit==ihinit
C        first fill bi(,i) and ai(,i)  
         ha2=ha*ha
         do 11 j=1,ihinit
            ja=nn(j,i)
            ij=(j-1)
            z=1.d2*ij*ij/ha2
            iz=z
            az=z-iz
            epij=(1.-az)*kernl(iz+1)+az*kernl(iz+2)
            z=1.d0
            xij(1)=1.d0
            ai(1,i)=ai(1,i)+y(ja)*epij
            if(dp1.eq.1) goto 112
            do 111 k=1,px
               z=x(k,i)-x(k,ja)
               xij(k+1)=-z
               ai(k+1,i)=ai(k+1,i)-z*y(ja)*epij
111            continue
112         m=1
            do 12 k=1,dp1
               do 12 l=1,k
                  bi(m,i)=bi(m,i)+xij(k)*xij(l)*epij
                  m=m+1
12         continue
11      continue
C       expand bi as p times p matrix
        m=1
        do 13 k=1,dp1
           do 13 l=1,k
              dmat(l,k)=bi(m,i)
              m=m+1
13      continue
C     compute inverse of bi by choleski decomposition
         if(dp1.eq.1) goto 15
         call invers(dmat,dp1,info)
         if(info.gt.0) goto 9999
C        ihinit to samll
         goto 16
15       dmat(1,1)=1.d0/dmat(1,1)
         info=0
C     now dmat contains inverse of B_i 
C     now calculate theta as B_i^{-1} Z_i
16      do 20 j=1,dp1
           d=0.0d0
           do 21 k=1,dp1
              d=d+dmat(j,k)*ai(k,i)  
21         continue
           theta(j,i)=d
20      continue
1     continue
9999  return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local polynomial aws (NN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpawsmnn(n,px,dp1,dp2,x,y,nn,ihakt,theta,bi,bin,bi0,
     1  ai,ain,lam,tau,h,kernl,kerns,dmat,dmati,dmat0,thij,psix,psiy,
     2  sym)
C    
C     n          number of design points
C     px         dimension of x
C     dp1        number of parameters  (1+p*px)
C     dp2        number of components in bi  (dp1*(dp1+1)/2)
C     x          design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ai         \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     tau        tau                      (extension penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     kernl      discretized localization kernel 
C     kerns      discretized stochastic and extension kernel (exponential)
C     dmat, dmati, dmat0       working arrays  dp1 times dp1
C     thij                     vector for parameter differences
C     psix,psiy,xij                working memory 
C     sym        asymmetric or symmetric test (logical)
C     
      implicit logical (a-z)
      integer n,dp1,dp2,px,i,j,k,l,info,iz,je,ja,o,nwij,ij,
     1        ihakt,m,nn(ihakt,n)      
      logical sym
      real*8 x(px,n),y(n),psix(dp1),psiy(dp1),theta(dp1,n),kerns(102),
     1 bi(dp2,n),bin(dp2,n),ai(dp1,n),ain(dp1,n),lam,kernl(102),
     2 dmat(dp1,dp1),dmati(dp1,dp1),thij(dp1),
     3 bi0(dp2,n),dmat0(dp1,dp1),pij,epij,tau,lambda,tauakt,h,
     4 d,gammaij,wij,z,az,xij,gamma0,s0i,ha,ha2,eps
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      eps=.2*kernl(100)
      ha=h-.9d0
      ha2=ha*ha/1.d2
      do 1 i=1,n
C        loop over design points
         lambda=lam*1.d-2
         if(sym) lambda=2*lambda
         tauakt=tau*1.d-2
         if(bi(1,i).lt.5.d0*dp1) tauakt=1d10
C  disable extension penalty if estimate is very unstable
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmati and dmat0
C    expand Bi in dmat
         s0i=bi0(1,i)
         m=1
         do 101 k=1,dp1
            do 101 l=1,k
               dmat(l,k)=bi(m,i)
               dmat(k,l)=bi(m,i)
               dmati(l,k)=bi(m,i)
               dmat0(l,k)=bi0(m,i)
               m=m+1
101      continue
         do 105 l=1,dp2
         bi0(l,i)=0.d0
105      continue
C     compute inverse of dmat0 and dmati
         if(dp1.eq.1) goto 8999
         call invers(dmat0,dp1,info)
         call invers(dmati,dp1,info)
C     now dmati contains inverse of B_i, dmat0 the inverse of B^(0), 
C     not needed for dp1=1
8999     nwij=0
         do 108 l=1,dp2
            bin(l,i)=0.d0
            if(l.gt.dp1) goto 108
            ain(l,i)=0.d0
108      continue
C      if not enough points with positive weights (counted in nwij)
C      lambda and tau will be increased 
C      (to reduce panalization by stochastic and influence term)
         do 11 j=1,ihakt
            ja=nn(j,i)
            ij=j-1
            pij=ij*ij/ha2
            iz=pij
            az=pij-iz
            epij=kernl(iz+1)*(1.d0-az)+kernl(iz+2)*az
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
              do 111 l=1,dp1
                 thij(l)=theta(l,ja)
111              continue
                 psix(1)=1
                 psiy(1)=y(ja)
                 gammaij=0.d0
                 if(dp1.eq.1) goto 1202
                 z=1.d0
                 do 112 l=2,dp1
                    xij=(x(l-1,i)-x(l-1,ja))
                    thij(1)=thij(1)+xij*theta(l,ja)
                    psix(l)=-xij
                    psiy(l)=-xij*y(ja)
112              continue
C     thij contains theta_{j} in the model centered at xi
C    now get gamma_ij
            gamma0=0.d0
            do 119 l=1,dp1
               do 119 k=1,dp1
                  z=psix(l)*psix(k)
                  gammaij=gammaij+dmati(l,k)*z
                  gamma0=gamma0+dmat0(l,k)*z
119         continue
            gammaij=gammaij*bi(1,i)/(gamma0*s0i)
            gammaij=dmax1(gammaij-1,0.d0)/tauakt
C
C     now we have everything to compute  w_{ij}
C
            m=1       
            do 1201 k=1,dp1
               do 1201 l=1,k
                  bi0(m,i)=bi0(m,i)+psix(k)*psix(l)*epij
                  m=m+1
1201        continue   
C
C   thats needed to compute gamma0 in the next iteration
C
1202        d=0.0d0
            do 1210 l=1,dp1
               thij(l)=theta(l,i)-thij(l)
1210           continue
C
C   thats the difference between thetai and thetaij
C
            do 121 l=1,dp1
               do 122 k=1,dp1
                  d=d+dmat(k,l)*thij(l)*thij(k)
122            continue
121         continue
            z=d/lambda+gammaij
            if(z.ge.1.d2) goto 11
            iz=z
            az=z-iz
            wij=epij*(kerns(iz+1)*(1-az)+kerns(iz+2)*az)
            if(wij.gt.eps) nwij=nwij+1
C        now compute contributions to bi(i),ai(i)  
            m=1
            do 123 k=1,dp1
               ain(k,i)=ain(k,i)+wij*psiy(k)
               do 123 l=1,k
                  bin(m,i)=bin(m,i)+wij*psix(k)*psix(l)
                  m=m+1
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
c
c     invers computes the inverse of a certain
c     double precision symmetric positive definite matrix (see below)
c     
c     this code is based on choleski decomposition and
c     integrates code from linpack (dpofa.f and dpodi.f, version of 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.)
c     and BLAS
C
c     on entry
c
c        a       double precision(n, n)
c
c        n       integer
c                the order of the matrix  a .
c
c
c     on return
c
c        a       invers produces the upper half of inverse(a) .
c
c     subroutines and functions
c
c     fortran dsqrt
c
      subroutine invers(a, n, info)
      integer n,info
      double precision a(n,n)
c
c     internal variables
c
      double precision ddot,t
      double precision s
      integer j,jm1,k,l,i,kp1,im1
c     begin block with ...exits to 940 if singular (info.ne.0)
c
c
         do 30 j = 1, n
            info = j
            s = 0.0d0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
           ddot=0.0d0
           do 11 l=1,k-1
              ddot=ddot+a(l,k)*a(l,j)
   11          continue           
               t = a(k,j) - ddot
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
   10       continue
   20       continue
            s = a(j,j) - s
c     ......exit
            if (s .le. 1.d-100) go to 940
            a(j,j) = dsqrt(s)
   30    continue
         info = 0
c
c     now we have the choeski decomposition in a     
c
c     next code from dpodi to compute inverse
c
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
        do 101 l=1,k-1
           a(l,k)=t*a(l,k)
  101       continue          
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
           do 81 l=1,k
              a(l,j)=a(l,j)+t*a(l,k)
   81          continue       
   80       continue
   90       continue
  100    continue
c
c        form  inverse(r) * trans(inverse(r))
c
         do 130 j = 1, n
            jm1 = j - 1
            if (jm1 .lt. 1) go to 120
            do 110 k = 1, jm1
               t = a(k,j)
           do 111 l=1,k
              a(l,k)=a(l,k)+t*a(l,j)
  111          continue       
  110       continue
  120       continue
            t = a(j,j)
        do 121 l=1,j
           a(l,j)=t*a(l,j)
  121       continue          
  130    continue
c
c     now fill lower triangle       
c  
      do 240 i=1,n
        im1 = i-1
        do 230 j=1,im1
          a(i,j) = a(j,i)
  230   continue
  240 continue
  940 continue     
      return
      end
