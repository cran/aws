C
C    Copyright (C) 2019 Weierstrass-Institut fuer
C                       Angewandte Analysis und Stochastik (WIAS)
C
C    Author:  Joerg Polzehl
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
C  FORTRAN 77 code needed in R functions aws, vaws,
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Patch based aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pcaws(y,n1,n2,n3,hakt,lambda,theta,bi,bi2,
     1                bi0,bin,ai,model,kern,spmin,lwght,wght,npsize)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   npsize       patch size
C
      implicit none

      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,model,kern,npsize
      logical aws
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,wght(2),
     1       bi2(*),bin(*),hakt,lwght(*),spmin,spf
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,ip1,ip2,ip3,nph1,nph2,nph3,
     3        ipind,jp1,jp2,jp3,jpind,np1,np2,np3,ipi0,jpi0
      double precision sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,
     1       hakt2,hmax2,w1,w2
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      np1=2*npsize+1
      np2=min(n2,2*npsize+1)
      np3=min(n3,2*npsize+1)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      nph1=(np1-1)/2
      nph2=(np2-1)/2
      nph3=(np3-1)/2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
C   rescale bi with 1/lambda
      DO iind=1,n1*n2*n3
         bi(iind) = bi(iind)/lambda
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,bin,n1,n2,n3,hakt2,hmax2,theta,
C$OMP& ih3,lwght,wght,y,nph1,nph2,nph3,npsize)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& model,spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,swj,swj2,
C$OMP& swj0,swjy,sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,ip1,ip2,ip3,ipind,ipi0,
C$OMP& jp1,jp2,jp3,jpind,jpi0)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
C   scaling of sij outside the loop
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  wj=lwght(jw1+clw1+1+jwind2)
                  if(wj.le.1e-10) CYCLE
                  swj0=swj0+wj
                  IF (aws) THEN
                     sij=0.d0
                      DO ip1=i1-nph1,i1+nph1
                         if(ip1.le.0.or.ip1.gt.n1) CYCLE
                         jp1=ip1+jw1
                         if(jp1.le.0.or.jp1.gt.n1) CYCLE
                         DO ip2=i2-nph2,i2+nph2
                            if(ip2.le.0.or.ip2.gt.n2) CYCLE
                            ipi0=ip1+(ip2-1)*n1
                            jp2=ip2+jw2
                            if(jp2.le.0.or.jp2.gt.n2) CYCLE
                            jpi0=jp1+(jp2-1)*n1
                            DO ip3=i3-nph3,i3+nph3
                               if(sij.gt.1.d0) CYCLE
                               if(ip3.le.0.or.ip3.gt.n3) CYCLE
                               jp3=ip3+jw3
                               if(jp3.le.0.or.jp3.gt.n3) CYCLE
                               jpind=jpi0+(jp3-1)*n12
                               ipind=ipi0+(ip3-1)*n12
                               sij=max(sij,bi(ipind)*
     1                      kldist(model,theta(ipind),theta(jpind)))
                            END DO
                         END DO
                      END DO
                     IF(sij.gt.1.d0) CYCLE
                     IF (sij.gt.spmin) THEN
                        wj=wj*(1.d0-spf*(sij-spmin))
                     END IF
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  swjy=swjy+wj*y(jind)
               END DO
            END DO
         END DO
         ai(iind)=swjy
         bin(iind)=swj
         bi2(iind)=swj2
         bi0(iind)=swj0
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bin,bi0,bi2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Patch based aws using mask
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pcawsm(y,pos,n1,n2,n3,hakt,lambda,theta,bi,bi2,
     1                bin,thnew,model,kern,spmin,lwght,wght,npsize)
C
C   y        observed values of regression function
C   mask     image mask
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew       \sum  Wi Y / \sum Wi    (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   npsize       patch size
C
      implicit none

      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,model,kern,np,npsize,pos(*)
      logical aws
      double precision y(*),theta(*),bi(*),thnew(*),lambda,wght(2),
     1       bi2(*),bin(*),hakt,lwght(*),spmin,spf
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,ip1,ip2,ip3,nph1,nph2,nph3,
     3        jp1,jp2,jp3,np1,np2,np3,iindp,ipindp,jindp,jpindp
      double precision sij,swj,swj2,swjy,z1,z2,z3,wj,
     1       hakt2,hmax2,w1,w2
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      np1=2*npsize+1
      np2=min(n2,2*npsize+1)
      np3=min(n3,2*npsize+1)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      nph1=(np1-1)/2
      nph2=(np2-1)/2
      nph3=(np3-1)/2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
          ELSE
            jind3=0
          END IF
          DO j2=-ih2,ih2
            if(n2.gt.1) THEN
                z2=j2*w1
                z2=z3+z2*z2
                ih1=FLOOR(sqrt(hakt2-z2))
                jind2=jind3+(j2+clw2)*dlw1
            ELSE
                jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
                jind=j1+clw1+1+jind2
                z1=j1
                lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
                if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
          END DO
      END DO
C   rescale bi with 1/lambda
      DO iind=1,n1*n2*n3
          bi(iind) = bi(iind)/lambda
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,bi2,bin,n1,n2,n3,hakt2,hmax2,theta,pos,
C$OMP& ih3,lwght,wght,y,nph1,nph2,nph3,np,npsize)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& model,spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,swj,swj2,
C$OMP& swjy,sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,ip1,ip2,ip3,
C$OMP& jp1,jp2,jp3,iindp,ipindp,jindp,jpindp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
        iindp = pos(iind)
        if(iindp.eq.0) CYCLE
C voxel not in mask
          i1=mod(iind,n1)
          if(i1.eq.0) i1=n1
          i2=mod((iind-i1)/n1+1,n2)
          if(i2.eq.0) i2=n2
          i3=(iind-i1-(i2-1)*n1)/n12+1
C   scaling of sij outside the loop
          swj=0.d0
          swj2=0.d0
          swjy=0.d0
          DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
                j2=jw2+i2
                if(j2.lt.1.or.j2.gt.n2) CYCLE
                jwind2=jwind3+(jw2+clw2)*dlw1
                jind2=(j2-1)*n1+jind3
                z2=jw2*w1
                z2=z3+z2*z2
                ih1=FLOOR(sqrt(hakt2-z2))
                DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jindp=pos(jind)
                  if(jindp.eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                      sij=0.d0
                      DO ip1=i1-nph1,i1+nph1
                          if(ip1.le.0.or.ip1.gt.n1) CYCLE
                          jp1=ip1+jw1
                          DO ip2=i2-nph2,i2+nph2
                            if(ip2.le.0.or.ip2.gt.n2) CYCLE
                            jp2=ip2+jw2
                            DO ip3=i3-nph3,i3+nph3
                                if(sij.gt.1.d0) CYCLE
                                if(ip3.le.0.or.ip3.gt.n3) CYCLE
                                ipindp=pos(ip1+(ip2-1)*n1+(ip3-1)*n12)
                                if(ipindp.eq.0) CYCLE
                                if(jp1.le.0.or.jp1.gt.n1) CYCLE
                                if(jp2.le.0.or.jp2.gt.n2) CYCLE
                                jp3=ip3+jw3
                                if(jp3.le.0.or.jp3.gt.n3) CYCLE
                                jpindp=pos(jp1+(jp2-1)*n1+(jp3-1)*n12)
                                if(jpindp.eq.0) CYCLE
                                sij=max(sij,bi(ipindp)*
     1                kldist(model,theta(ipindp),theta(jpindp)))
                            END DO
                          END DO
                      END DO
                      IF(sij.gt.1.d0) CYCLE
                      IF (sij.gt.spmin) THEN
                        wj=wj*(1.d0-spf*(sij-spmin))
                      END IF
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  swjy=swjy+wj*y(jindp)
                END DO
            END DO
          END DO
          thnew(iindp)=swjy/swj
          bin(iindp)=swj
          bi2(iindp)=swj2
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bin,bi2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pvaws(y,mask,nv,n1,n2,n3,hakt,lambda,theta,bi,
     1                bin,thnew,ncores,spmin,lwght,wght,swjy,
     2                np1,np2,np3)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none

      integer nv,n1,n2,n3,ncores,mask(*)
      logical aws
      double precision y(nv,*),theta(nv,*),bi(*),thnew(nv,*),lambda,
     1  wght(2),hakt,lwght(*),spmin,spf,swjy(nv,ncores),bin(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr
      double precision biinv,sij,swj,z,z1,z2,z3,wj,hakt2,hmax2,
     1        w1,w2,spmb,sijp
      integer np1,np2,np3
      integer ip1,ip2,ip3,nph1,nph2,nph3,ipind,jp1,jp2,jp3,jpind
      external lkern
      double precision lkern
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=FLOOR(hakt)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      nph1=(np1-1)/2
      nph2=(np2-1)/2
      nph3=(np3-1)/2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(2,(z1*z1+z2)/hakt2)
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
C   rescale bi with 1/lambda
      DO iind=1,n1*n2*n3
          bi(iind) = bi(iind)/lambda
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,nv,n1,n2,n3,hakt2,hmax2,theta,bin,
C$OMP& ih3,lwght,wght,y,swjy,mask,nph1,nph2,nph3)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,biinv,swj,spmb,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z,thrednr,ip1,ip2,ip3,ipind,
C$OMP& jp1,jp2,jp3,jpind,sijp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         if(mask(iind).eq.0) CYCLE
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
C   scaling of sij outside the loop
         swj=0.d0
         DO k=1,nv
            swjy(k,thrednr)=0.d0
         END DO
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  if(mask(jind).eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     sij=0.d0
                     DO ip1=i1-nph1,i1+nph1
                        if(ip1.le.0.or.ip1.gt.n1) CYCLE
                        jp1=ip1+jw1
                        if(jp1.le.0.or.jp1.gt.n1) CYCLE
                        DO ip2=i2-nph2,i2+nph2
                           if(ip2.le.0.or.ip2.gt.n2) CYCLE
                           jp2=ip2+jw2
                           if(jp2.le.0.or.jp2.gt.n2) CYCLE
                           DO ip3=i3-nph3,i3+nph3
                              if(sij.gt.1.d0) CYCLE
                              if(ip3.le.0.or.ip3.gt.n3) CYCLE
                              jp3=ip3+jw3
                              if(jp3.le.0.or.jp3.gt.n3) CYCLE
                              jpind=jp1+(jp2-1)*n1+(jp3-1)*n12
                              ipind=jp1+(jp2-1)*n1+(jp3-1)*n12
                              if(mask(jpind).eq.0) CYCLE
                              sijp=0.d0
                              DO k=1,nv
                                 z=theta(k,ipind)-theta(k,jpind)
                                 sijp=sijp+z*z
                              END DO
                              sij=max(sij,bi(ipind)*sijp)
                           END DO
                        END DO
                     END DO
                     IF (sij.ge.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  swj=swj+wj
                  DO k=1,nv
                     swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,jind)
                  END DO
               END DO
            END DO
         END DO
         DO k=1,nv
            thnew(k,iind)=swjy(k,thrednr)/swj
         END DO
         bin(iind)=swj
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bin)
      RETURN
      END
      subroutine pvaws2(y,mask,nv,nvd,n1,n2,n3,hakt,lambda,theta,bi,
     1                bin,thnew,invcov,ncores,spmin,lwght,wght,swjy,
     2                np1,np2,np3)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none

      integer nv,n1,n2,n3,ncores,nvd,mask(*)
      logical aws
      double precision y(nv,*),theta(nv,*),bi(*),thnew(nv,*),lambda,
     1  wght(2),hakt,lwght(*),spmin,spf,swjy(nv,ncores),invcov(nvd,*),
     2  bin(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr
      double precision biinv,sij,swj,z,z1,z2,z3,wj,hakt2,hmax2,
     1        w1,w2,spmb,sijp
      integer np1,np2,np3,l,m
      integer ip1,ip2,ip3,nph1,nph2,nph3,ipind,jp1,jp2,jp3,jpind
      external lkern, KLdistsi
      double precision lkern, KLdistsi
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=FLOOR(hakt)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      nph1=(np1-1)/2
      nph2=(np2-1)/2
      nph3=(np3-1)/2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(2,(z1*z1+z2)/hakt2)
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
C   rescale bi with 1/lambda
      DO iind=1,n1*n2*n3
        bi(iind) = bi(iind)/lambda
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,nv,nvd,n1,n2,n3,hakt2,hmax2,theta,invcov,
C$OMP& ih3,lwght,wght,y,swjy,mask,nph1,nph2,nph3,bin)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,biinv,swj,spmb,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z,thrednr,ip1,ip2,ip3,ipind,
C$OMP& jp1,jp2,jp3,jpind,sijp,l,m)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         if(mask(iind).eq.0) CYCLE
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
C   scaling of sij outside the loop
         swj=0.d0
         DO k=1,nv
            swjy(k,thrednr)=0.d0
         END DO
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  if(mask(jind).eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     sij=0.d0
                     DO ip1=i1-nph1,i1+nph1
                        if(ip1.le.0.or.ip1.gt.n1) CYCLE
                        jp1=ip1+jw1
                        if(jp1.le.0.or.jp1.gt.n1) CYCLE
                        DO ip2=i2-nph2,i2+nph2
                           if(ip2.le.0.or.ip2.gt.n2) CYCLE
                           jp2=ip2+jw2
                           if(jp2.le.0.or.jp2.gt.n2) CYCLE
                           DO ip3=i3-nph3,i3+nph3
                              if(sij.gt.1.d0) CYCLE
                              if(ip3.le.0.or.ip3.gt.n3) CYCLE
                              ipind=ip1+(ip2-1)*n1+(ip3-1)*n12
                              if(mask(ipind).eq.0) CYCLE
                              jp3=ip3+jw3
                              if(jp3.le.0.or.jp3.gt.n3) CYCLE
                              jpind=jp1+(jp2-1)*n1+(jp3-1)*n12
                              if(mask(jpind).eq.0) CYCLE
C   need both ipind and jpind in mask,
                              sijp=KLdistsi(theta(1,jpind),
     1                                theta(1,ipind),
     2                                invcov(1,ipind),nv)
                              sij=max(sij,bi(ipind)*sijp)
                           END DO
                        END DO
                     END DO
                     IF (sij.ge.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  swj=swj+wj
                  DO k=1,nv
                     swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,jind)
                  END DO
               END DO
            END DO
         END DO
         DO k=1,nv
            thnew(k,iind)=swjy(k,thrednr)/swj
         END DO
         bin(iind)=swj
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bin)
      RETURN
      END

      double precision function KLdistsi(thi,thj,si2,nv)
      implicit logical (a-z)
      integer nv
      double precision thi(nv), thj(nv), si2(*)
      integer k,l,m
      double precision z,zdk
      z=0.d0
      m=1
      DO k=1,nv
         zdk=thi(k)-thj(k)
         if(k.gt.1) THEN
           DO l=1,k-1
              z=z+2.d0*(thi(l)-thj(l))*zdk*si2(m)
              m=m+1
           END DO
         ENDIF
         z=z+zdk*zdk*si2(m)
         m=m+1
      END DO
      KLdistsi=z
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Patch based aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pawswght(n1,n2,n3,i1,i2,i3,hakt,lambda,theta,bi,
     1              model,kern,spmin,lwght,wght,npsize,wi)
C
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (input)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   npsize       patch size
C
      implicit none

      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,i1,i2,i3,model,kern,npsize
      double precision theta(*),bi(*),lambda,wght(2),
     1       wi(*),hakt,lwght(*),spmin,spf
      integer ih1,ih2,ih3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,ip1,ip2,ip3,nph1,nph2,nph3,
     3        ipind,jp1,jp2,jp3,jpind,np1,np2,np3
      double precision sij,z1,z2,z3,wj,
     1       hakt2,hmax2,w1,w2
       hakt2=hakt*hakt
       spf=1.d0/(1.d0-spmin)
C
C   first calculate location weights
C
       w1=wght(1)
       w2=wght(2)
       ih3=FLOOR(hakt/w2)
       ih2=FLOOR(hakt/w1)
       ih1=FLOOR(hakt)
       np1=2*npsize+1
       np2=min(n2,2*npsize+1)
       np3=min(n3,2*npsize+1)
       if(n3.eq.1) ih3=0
       if(n2.eq.1) ih2=0
       clw1=ih1
       clw2=ih2
       clw3=ih3
       dlw1=ih1+clw1+1
       dlw2=ih2+clw2+1
       dlw3=ih3+clw3+1
       dlw12=dlw1*dlw2
       nph1=(np1-1)/2
       nph2=(np2-1)/2
       nph3=(np3-1)/2
       n12=n1*n2
       z2=0.d0
       z3=0.d0
       hmax2=0.d0
       DO j3=-clw3,clw3
          if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
          ELSE
            jind3=0
          END IF
          DO j2=-ih2,ih2
            if(n2.gt.1) THEN
                z2=j2*w1
                z2=z3+z2*z2
                ih1=FLOOR(sqrt(hakt2-z2))
                jind2=jind3+(j2+clw2)*dlw1
            ELSE
                jind2=0
            END IF
            DO j1=-ih1,ih1
                jind=j1+clw1+1+jind2
                z1=j1
                lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
                if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
          END DO
      END DO
      call rchkusr()
      iind=i1+(i2-1)*n1+(i3-1)*n12
      DO jw3=-clw3,clw3
        j3=jw3+i3
        if(j3.lt.1.or.j3.gt.n3) CYCLE
        jwind3=(jw3+clw3)*dlw12
        jind3=(j3-1)*n12
        z3=jw3*w2
        z3=z3*z3
        if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
        DO jw2=-ih2,ih2
            j2=jw2+i2
            if(j2.lt.1.or.j2.gt.n2) CYCLE
            jwind2=jwind3+(jw2+clw2)*dlw1
            jind2=(j2-1)*n1+jind3
            z2=jw2*w1
            z2=z3+z2*z2
            ih1=FLOOR(sqrt(hakt2-z2))
            DO jw1=-ih1,ih1
C  first stochastic term
                j1=jw1+i1
                if(j1.lt.1.or.j1.gt.n1) CYCLE
                jind=j1+jind2
                wj=lwght(jw1+clw1+1+jwind2)
                sij=0.d0
                DO ip1=i1-nph1,i1+nph1
                    if(sij.gt.1.d0) CYCLE
                    if(ip1.le.0.or.ip1.gt.n1) CYCLE
                    jp1=ip1+jw1
                    if(jp1.le.0.or.jp1.gt.n1) CYCLE
                    DO ip2=i2-nph2,i2+nph2
                      if(sij.gt.1.d0) CYCLE
                      if(ip2.le.0.or.ip2.gt.n2) CYCLE
                      jp2=ip2+jw2
                      if(jp2.le.0.or.jp2.gt.n2) CYCLE
                      DO ip3=i3-nph3,i3+nph3
                          if(sij.gt.1.d0) CYCLE
                          if(ip3.le.0.or.ip3.gt.n3) CYCLE
                          jp3=ip3+jw3
                          if(jp3.le.0.or.jp3.gt.n3) CYCLE
                          jpind=jp1+(jp2-1)*n1+(jp3-1)*n12
                          ipind=ip1+(ip2-1)*n1+(ip3-1)*n12
                          sij=max(sij,bi(ipind)/lambda*
     1                kldist(model,theta(ipind),theta(jpind)))
                      END DO
                    END DO
                END DO
                IF(sij.gt.1.d0) CYCLE
                IF (sij.gt.spmin) THEN
                    wj=wj*(1.d0-spf*(sij-spmin))
                END IF
                wi(jind) = wj
            END DO
        END DO
      END DO
      RETURN
      END
