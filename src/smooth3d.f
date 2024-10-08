CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform 3D nonadaptive smoothing on a grid (gridded)
C   mask is coded via position
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smooth3d(y,si2,pos,wlse,nvox,n1,n2,n3,dv,hakt,
     1                    thn,bi,kern,lwght,wght,swjy)
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
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none
      integer nvox,n1,n2,n3,kern,dv,wlse,pos(n1,n2,n3)
      double precision y(nvox,dv),thn(nvox,dv),bi(nvox),wght(2),
     1       si2(nvox),hakt,lwght(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n,dlw12,iindp,jindp
      double precision swj,swjy(dv),wj,hakt2
      hakt2=hakt*hakt
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/wght(2))
      ih2=FLOOR(hakt/wght(1))
      ih1=FLOOR(hakt)
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      dlw12=dlw1*dlw2
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               iindp = pos(i1,i2,i3)
               if(iindp.eq.0) CYCLE
               swj=0.d0
               DO k=1,dv
                  swjy(k)=0.d0
               END DO
               DO jw3=1,dlw3
                  j3=jw3-clw3+i3
                  if(j3.lt.1.or.j3.gt.n3) CYCLE
                  DO jw2=1,dlw2
                     j2=jw2-clw2+i2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     DO jw1=1,dlw1
C  first stochastic term
                        j1=jw1-clw1+i1
                        if(j1.lt.1.or.j1.gt.n1) CYCLE
                        jindp = pos(j1,j2,j3)
                        IF(jindp.eq.0) CYCLE
                        wj=lwght(jw1+(jw2-1)*dlw1+(jw3-1)*dlw12)
                        if(wj.le.0.d0) CYCLE
                        if(wlse.ne.0) THEN
                           wj=wj*si2(jindp)
                        END IF
                        swj=swj+wj
                        DO k=1,dv
                           swjy(k)=swjy(k)+wj*y(jindp,k)
                        END DO
                     END DO
                  END DO
               END DO
               DO k=1,dv
                  thn(iindp,k)=swjy(k)/swj
               END DO
               bi(iindp) = swj
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END

C
C   calculate location weights in lwght
C
      subroutine locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      implicit none
      integer dlw1,dlw2,dlw3,kern
      double precision wght(2),hakt2,lwght(dlw1,dlw2,dlw3),lkern
      external lkern
      double precision z1,z2,z3
      integer j1,j2,j3,clw1,clw2,clw3,ih1,ih2
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      DO j3=1,dlw3
         Do j2=1,dlw2
            DO j1=1,dlw1
               lwght(j1,j2,j3)=0.d0
            END DO
         END DO
         z3=(clw3-j3)*wght(2)
         z3=z3*z3
         ih2=FLOOR(sqrt(hakt2-z3)/wght(1))
         DO j2=clw2-ih2,clw2+ih2
            IF(j2.lt.1.or.j2.gt.dlw2) CYCLE
            z2=(clw2-j2)*wght(1)
            z2=z3+z2*z2
            ih1=FLOOR(sqrt(hakt2-z2))
            DO j1=clw1-ih1,clw1+ih1
               IF(j1.lt.1.or.j1.gt.dlw1) CYCLE
               z1=clw1-j1
               lwght(j1,j2,j3)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          Compute aws-weights  w_{ij} Gaussian case
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awswght3(thi,theta,skern,spf,spmin,spmax,bii,wj)
      implicit none
      integer skern
      double precision thi,theta,spf,spmin,spmax,bii,wj,wjin
      double precision sij,z
      wjin=wj
C  compute distance in sij
      z=thi-theta
      sij=bii*z*z
      IF (sij.gt.spmax) THEN
         wj=0.d0
      ELSE IF (skern.eq.1) THEN
C  skern == "Plateau"
         wj=wj*min(1.d0,1.d0-spf*(sij-spmin))
      ELSE IF (skern.eq.2) THEN
C  skern == "Triangle"
         wj=wj*(1.d0-sij)
      ELSE
C  skern == "Exp"
         IF (sij.gt.spmin) wj=wj*exp(-spf*(sij-spmin))
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   local constant 3D heteroskedastic aws for y
C   si2   inverse of variance
C   mask is coded via position
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chaws2(y,si2,pos,wlse,n1,n2,n3,hakt,lambda,
     1                  theta,bi,thn,kern,skern,spmin,spmax,
     2                  lwght,wght)
C
C   y        3D parameter map
C   si2      variance estimate^{-1} (1/arg to function segm3D)
C   mask     brain mask
C   wlse     (arg weighted to function segm3D)
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda
C   theta    estimates from last step   (input)
C   bi       \sum  Wi (wi depends on si2 in case of wlse)  (output)
C   thn      new estimates (output)
C   kern     indicator for K_{loc}
C   skern    indicator for K_{st}
C   spmin    shape parameter for plateo kernel  (skern=1)
C   spmax    shape parameter for kernel K_{st}  (skern=1,2)
C   lwght    array for non-adaptive weights (auxiliary)
C   wght     rations of voxel dimensions
C
      implicit none
      integer n1,n2,n3,kern,skern,wlse,pos(*)
      logical aws
      double precision y(*),theta(*),bi(*),thn(*),lambda,spmax,
     1       wght(2),si2(*),hakt,lwght(*),spmin
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,clw1,
     1        clw2,clw3,dlw1,dlw2,dlw3,n,iind,jind,n12,jind2,jind3,
     2        a1,e1,a2,e2,a3,e3,dlw12,iindp,jindp
      double precision bii,swj,swjy,thi,wj,hakt2,spf,si2i,swjv
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/wght(2))
      ih2=FLOOR(hakt/wght(1))
      ih1=FLOOR(hakt)
      n=n1*n2*n3
      n12=n1*n2
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      dlw12=dlw1*dlw2
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      a1=1
      e1=n1
      a2=1
      e2=n2
      a3=1
      e3=n3
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,kern,skern,aws,wlse,pos,y,theta,bi,thn,
C$OMP& lambda,spmax,wght,si2,hakt,lwght,spmin,spf,dlw12,
C$OMP& ih1,ih2,ih3,clw1,clw2,clw3,dlw1,dlw2,dlw3,n,hakt2,n12)
C$OMP& FIRSTPRIVATE(a1,e1,a2,e2,a3,e3)
C$OMP& PRIVATE(iind,i1,i2,i3,si2i,bii,swjv,swj,thi,swjy,
C$OMP& j1,j2,j3,jw1,jw2,jw3,jind,jind2,jind3,wj,iindp,jindp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n
         iindp = pos(iind)
         if(iindp.eq.0) CYCLE
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
C  precompute range of jw1,jw2,jw3
         if(i3.ge.clw3) THEN
            a3=1
         ELSE
            a3=clw3-i3+1
         END IF
         if(i3.le.n3-clw3+1) THEN
            e3=dlw3
         ELSE
            e3=n3+clw3-i3
         END IF
         if(i2.ge.clw2) THEN
            a2=1
         ELSE
            a2=clw2-i2+1
         END IF
         if(i2.le.n2-clw2+1) THEN
            e2=dlw2
         ELSE
            e2=n2+clw2-i2
         END IF
         if(i1.ge.clw1) THEN
            a1=1
         ELSE
            a1=clw1-i1+1
         END IF
         if(i1.le.n1-clw1+1) THEN
            e1=dlw1
         ELSE
            e1=n1+clw1-i1
         END IF
         si2i=si2(iindp)
         bii=bi(iindp)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         swjv=0.d0
         swjy=0.d0
         thi=theta(iindp)
         DO jw3=a3,e3
            j3=jw3-clw3+i3
            jind3=n12*(j3-1)
            DO jw2=a2,e2
               j2=jw2-clw2+i2
               jind2=jind3+n1*(j2-1)
               DO jw1=a1,e1
C  first stochastic term
                  j1=jw1-clw1+i1
                  jind=j1+jind2
                  jindp = pos(jind)
                  IF(jindp.eq.0) CYCLE
                  wj=lwght(jw1+(jw2-1)*dlw1+(jw3-1)*dlw12)
                  if(wj.le.0.d0) CYCLE
                  IF (aws) THEN
                     call awswght3(thi,theta(jindp),
     1                             skern,spf,spmin,spmax,bii,wj)
                     if(wj.le.0.d0) CYCLE
                  END IF
                  if(wlse.ne.0) THEN
                     wj=wj*si2(jindp)
                  ELSE
                     swjv=swjv+wj/si2(jindp)
                  END IF
                  swj=swj+wj
                  swjy=swjy+wj*y(jindp)
               END DO
            END DO
         END DO
         thn(iindp)=swjy/swj
         IF(wlse.ne.0) THEN
            bi(iindp)=swj
         ELSE
            bi(iindp)=swj*swj/swjv
         END IF
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,bi)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   local constant 3D heteroskedastic aws for y
C   si2   inverse of variance
C   mask is coded via position
C   vectorial data in res is smoothed using same weights
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chawsv(y,res,si2,pos,wlse,n1,n2,n3,n4,hakt,
     1                  lambda,theta,bi,resnew,thn,kern,skern,
     2                  spmin,spmax,lwght,wght,resi)
C
C   y        3D array (estimated spm)
C   res      residual array (4D)
C   si2      variance estimate^{-1} (1/arg to function segm3D)
C   mask     brain mask
C   wlse     (arg weighted to function segm3D)
C   n1,n2,n3,n4    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda
C   theta    estimates from last step   (input)
C   bi       \sum  Wi (wi depends on si2 in case of wlse)  (output)
C   resnew   smoothed residuals (output)
C   thn      new estimates (output)
C   kern     indicator for K_{loc}
C   skern    indicator for K_{st}
C   spmin    shape parameter for plateo kernel  (skern=1)
C   spmax    shape parameter for kernel K_{st}  (skern=1,2)
C   lwght    array for non-adaptive weights (auxiliary)
C   wght     rations of voxel dimensions
C   resi     auxilary array
C
!$      use omp_lib
      implicit none
      integer n1,n2,n3,n4,kern,skern,wlse,pos(*)
      logical aws
      double precision res(n4,*),y(*),theta(*),bi(*),thn(*),lambda,
     1   spmax,wght(2),si2(*),hakt,lwght(*),spmin,resi(*),resnew(n4,*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,clw1,clw2,clw3,
     1  dlw1,dlw2,dlw3,k,n,thrednr,iind,jind,rthrednr,dlw12,jindp,iindp
      double precision bii,swj,swjy,thi,wj,hakt2,spf,si2i,sresisq,resik
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/wght(2))
      ih2=FLOOR(hakt/wght(1))
      ih1=FLOOR(hakt)
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      dlw12=dlw1*dlw2
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      thrednr = 1
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(kern,skern,pos,res,
C$OMP& y,theta,bi,thn,spmax,wght,si2,hakt,lwght,spmin,
C$OMP& resi,resnew,ih1,ih2,ih3,hakt2,spf,dlw12)
C$OMP& FIRSTPRIVATE(n,n1,n2,n3,n4,lambda,clw1,clw2,clw3,
C$OMP& dlw1,dlw2,dlw3,aws,wlse)
C$OMP& PRIVATE(iind,jind,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,k,si2i,
C$OMP& bii,swj,wj,sresisq,thrednr,rthrednr,resik,thi,swjy,iindp,jindp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n
!$         thrednr = omp_get_thread_num()+1
         iindp=pos(iind)
         if(iindp.eq.0) CYCLE
         rthrednr = (thrednr-1)*n4
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1
         si2i=si2(iindp)
         bii=bi(iindp)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         swjy=0.d0
         DO k=1,n4
            resi(k+rthrednr)=0.d0
         END DO
         thi=theta(iindp)
         DO jw3=1,dlw3
            j3=jw3-clw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            DO jw2=1,dlw2
               j2=jw2-clw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               DO jw1=1,dlw1
C  first stochastic term
                  j1=jw1-clw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+n1*(j2-1)+n1*n2*(j3-1)
                  jindp=pos(jind)
                  IF(jindp.eq.0) CYCLE
                  wj=lwght(jw1+(jw2-1)*dlw1+(jw3-1)*dlw12)
                  if(wj.le.0.d0) CYCLE
                  IF (aws) THEN
                     call awswght3(thi,theta(jindp),
     1                             skern,spf,spmin,spmax,bii,wj)
                     if(wj.le.0.d0) CYCLE
                  END IF
                  if(wlse.ne.0) THEN
                     wj=wj*si2(jindp)
                  END IF
                  swj=swj+wj
                  swjy=swjy+wj*y(jindp)
                  call daxpy(n4,wj,res(1,jindp),1,resi(1+rthrednr),1)
               END DO
            END DO
         END DO
         thn(iindp)=swjy/swj
         sresisq=0.d0
         DO k=1,n4
            resik=resi(k+rthrednr)
            resnew(k,iindp)=resik/swj
            sresisq=sresisq+resik*resik
         END DO
         bi(iindp)=swj*swj/sresisq*(n4-1.d0)
C  thats the inverse variance
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,bi,resnew)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   local constant 3D heteroskedastic aws for vector valued y
C   si2   inverse of variance (constant along vector y)
C   mask is coded via position
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ihaws2(y,si2,pos,wlse,n1,n2,n3,dv,hakt,lambda,
     1                  theta,ncores,bi,thn,kern,skern,spmin,spmax,
     2                  lwght,wght,swjy)
C
C   y        values of residuals
C   si2      variance estimate^{-1} (1/arg to function segm3D)
C   mask     brain mask
C   wlse     (arg weighted to function segm3D)
C   n1,n2,n3    design dimensions
C   dv       length of fMRI time serien
C   hakt     actual bandwidth
C   lambda   lambda
C   theta    estimates from last step   (input)
C   ncores   number of cores used in openMP
C   bi       \sum  Wi (wi depends on si2 in case of wlse)  (output)
C   thn      new estimates (output)
C   kern     indicator for K_{loc}
C   skern    indicator for K_{st}
C   spmin    shape parameter for plateo kernel  (skern=1)
C   spmax    shape parameter for kernel K_{st}  (skern=1,2)
C   lwght    array for non-adaptive weights (auxiliary)
C   wght     rations of voxel dimensions
C   swjy     auxilary array
C
!$      use omp_lib
      implicit none
      integer dv,n1,n2,n3,kern,skern,ncores,wlse,pos(*)
      logical aws
      double precision theta(*),bi(*),y(dv,*),
     1       lambda,spmax,wght(2),si2(*),thn(dv,*),
     1       hakt,lwght(*),spmin,swjy(dv,ncores)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,dlw12,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n,thrednr,iind,jind,
     2        iindp,jindp
      double precision bii,swj,wj,hakt2,spf,si2i,thi
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/wght(2))
      ih2=FLOOR(hakt/wght(1))
      ih1=FLOOR(hakt)
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      dlw12=dlw1*dlw2
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      thrednr = 1
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(dv,n1,n2,n3,kern,skern,ncores,aws,wlse,pos,theta,
C$OMP& bi,y,lambda,spmax,wght,si2,thn,hakt,lwght,spmin,dlw12,
C$OMP& ih1,ih2,ih3,clw1,clw2,clw3,dlw1,dlw2,dlw3,n,swjy,hakt2,spf)
C$OMP& PRIVATE(iind,jind,i1,i2,i3,k,si2i,bii,swj,j1,j2,j3,thi,
C$OMP& jw1,jw2,jw3,wj,thrednr,iindp,jindp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n
!$         thrednr = omp_get_thread_num()+1
         iindp=pos(iind)
         if(iindp.eq.0) CYCLE
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1
         si2i=si2(iindp)
         bii=bi(iindp)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         DO k=1,dv
            swjy(k,thrednr)=0.d0
         END DO
         thi=theta(iindp)
         DO jw3=1,dlw3
            j3=jw3-clw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            DO jw2=1,dlw2
               j2=jw2-clw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               DO jw1=1,dlw1
C  first stochastic term
                  j1=jw1-clw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+n1*(j2-1)+n1*n2*(j3-1)
                  jindp=pos(jind)
                  IF(jindp.eq.0) CYCLE
                  wj=lwght(jw1+(jw2-1)*dlw1+(jw3-1)*dlw12)
                  if(wj.le.0.d0) CYCLE
                  IF (aws) THEN
                     call awswght3(thi,theta(jindp),
     1                             skern,spf,spmin,spmax,bii,wj)
                     if(wj.le.0.d0) CYCLE
                  END IF
                  if(wlse.ne.0) THEN
                     wj=wj*si2(jindp)
                  END IF
                  swj=swj+wj
                  call daxpy(dv,wj,y(1,jindp),1,swjy(1,thrednr),1)
               END DO
            END DO
         END DO
         DO k=1,dv
            thn(k,iindp)=swjy(k,thrednr)/swj
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn)
      RETURN
      END
