CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   segmentation for Gaussian Models with fixed variance
C   access kritical values
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine segmgkrv(y,si2,n1,n2,n3,hakt,lambda,th,thn,bi,
     1                   lwght,wght,fov,maxvalue,minvalue)
C
C   y        observed values of regression function
C   si2      inverse of variances
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   th       estimates from last step   (input)
C   thn      new estimates (output)
C   bi       \sum  Wi   (input/output)
C   lwght    array for location weights
C   wght     scaling factor for second and third dimension (larger values shrink)
C   fov      field of view (may differ from n1*n2*n3 in case of spatial correlation)
C   maxvalue maximum observed value of test statistics
C   minvalue minimum observed value of test statistics

      implicit logical (a-z)
      integer n1,n2,n3
      logical aws
      real*8 y(1),si2(1),hakt,lambda,th(1),thn(1),bi(1),wght(2),
     1       lwght(1),fov,maxvalue,minvalue
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3
      real*8 a,b,spf,thi,bii,s2i,hakt2,z,z1,z2,z3,swj,swj2,swjy,wj,
     1       sij,fovh,lkern,biin
      external lkern
      hakt2=hakt*hakt
      spf=4.d0/3.d0
      ih1=hakt
      aws=lambda.lt.1d40
      maxvalue=-1.d4
      minvalue=1.d4
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1+1
      clw2=ih2+1
      clw3=ih3+1
      dlw1=ih1+clw1
      dlw2=ih2+clw2
      dlw3=ih3+clw3
      z2=0.d0
      z3=0.d0
      DO j3=1,dlw3
         if(n3.gt.1) THEN
            z3=(clw3-j3)*wght(2)
            z3=z3*z3
            ih2=sqrt(hakt2-z3)/wght(1)
            jind3=(j3-1)*dlw1*dlw2
         ELSE
            jind3=0
         END IF
         DO j2=clw2-ih2,clw2+ih2
            if(n2.gt.1) THEN
               z2=(clw2-j2)*wght(1)
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jind2=jind3+(j2-1)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=clw1-ih1,clw1+ih1
C  first location weights
               jind=j1+jind2
               z1=clw1-j1
               wj=lkern(1,(z1*z1+z2)/hakt2)
C  this is the plateau kernel
               lwght(jind)=wj
            END DO
         END DO
      END DO
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               thi = th(iind)
               s2i = si2(iind)
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               if(n3.gt.1) ih3=hakt/wght(2)
               swj=0.d0
               swj2=0.d0
               swjy=0.d0
               DO jw3=1,dlw3
                  j3=jw3-clw3+i3
                  if(j3.lt.1.or.j3.gt.n3) CYCLE
                  jind3=(j3-1)*n1*n2
                  z3=(clw3-jw3)*wght(2)
                  z3=z3*z3
                  if(n2.gt.1) ih2=sqrt(hakt2-z3)/wght(1)
                  jwind3=(jw3-1)*dlw1*dlw2
                  DO jw2=clw2-ih2,clw2+ih2
                     j2=jw2-clw2+i2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     jind2=(j2-1)*n1+jind3
                     z2=(clw2-jw2)*wght(1)
                     z2=z3+z2*z2
                     ih1=sqrt(hakt2-z2)
                     jwind2=jwind3+(jw2-1)*dlw1
                     DO jw1=clw1-ih1,clw1+ih1
C  first stochastic term
                        j1=jw1-clw1+i1
                        if(j1.lt.1.or.j1.gt.n1) CYCLE
                        jind=j1+jind2
                        wj=lwght(jw1+jwind2)
                        z1=(clw1-jw1)
                        z1=z2+z1*z1
                        IF (aws) THEN
C
C      gaussian case only
C
                           z=(thi-th(jind))
                           sij=bii*z*z
                           IF (sij.gt.1.d0) CYCLE
                           IF (sij.gt.0.25d0) THEN
                               wj=wj*(1.d0-spf*(sij-0.25d0))
                           END IF
                        END IF
                        swj=swj+wj*si2(jind)
                        swj2=swj2+wj*wj*si2(jind)
                        swjy=swjy+wj*si2(jind)*y(jind)
                     END DO
                  END DO
               END DO
               thi=swjy/swj
               thn(iind)=thi
               biin=swj*swj/swj2
               bi(iind)=biin
               fovh=fov/biin*s2i
               call getnab(fovh,a,b)
               biin=sqrt(biin)
C    both are equivalent for  homogeneous si2
               maxvalue=max(maxvalue,a*thi*biin-b)
               minvalue=min(minvalue,a*thi*biin+b)
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine getnab(n,a,b)
C
C   this function computes constants a=1/a_n and b=b_n/a_n
C   such that for the maximum T_n of Gaussian R.V. 
C   a_n T_n +b_n  ~ \Lambda   (Gumbel distribution)
C
      implicit logical (a-z)
      real*8 a,b,n
      a=sqrt(2.d0*log(n))
      b=a*a-(log(log(n))+log(1.256637d1))/2.d0
      RETURN
      END
      subroutine getnabp(n,a,b,lambda)
C
C   this function computes constants a=1/a_n and b=b_n/a_n
C   such that for the maximum T_n of standardized Poisson R.V.
C   a_n T_n +b_n  ~ \Lambda   (Gumbel distribution)
C   see C. W. Anderson, S. G. Coles and J. Husler (1997), 
C       Annals of Applied Probability, 7, 953-971
C
      implicit logical (a-z)
      real*8 a,b,n,lambda
      real*8 a2
      a=sqrt(2.d0*log(n))
      a2=a*a
      b=a2-(log(log(n))+log(1.256637d1))/2.d0+
     1  a*a2/6.d0/sqrt(lambda)-a2*a2/24.d0/lambda
      RETURN
      END
