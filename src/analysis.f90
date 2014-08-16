MODULE Analysis
 USE Params, ONLY : db, wflag, pi,lnPmin, Rkm, NucPh,rhos,EDtoKM, Love, MomI,iseed
 USE Data
 USE Parallel
 USE Chain
 IMPLICIT NONE
CONTAINS
 SUBROUTINE Bayesian
  INTEGER :: k,kk,i,j,kkk,cnt
  REAL(db) :: lnP,PP,dl,RR,MM,RRinf,zz
  INTEGER, ALLOCATABLE :: hist(:,:)
  lnPmax=lnPmin
  IF(mMaxAll<MmaxNS) RETURN
  PP=0
!$OMP PARALLEL DO PRIVATE(k,kk,lnP,dl,RR,MM,RRinf,zz) REDUCTION(+:PP)
  DO k=2,strmax
   IF(NSchain(k,1)<0) CYCLE
   kk=k-1
   DO WHILE(NSchain(kk,1)<0)
    kk=kk-1
    IF(kk==0) EXIT
   ENDDO
   IF(kk==0) CYCLE
   ! collapse to Black Hole check
   IF(NSchain(k,1)<NSchain(kk,1)) CYCLE
   MM=NSchain(k,1)
   RR=NSchain(k,2)
   RRinf=RR/SQRT(1-2*MM*Rkm/RR)
   zz=RRinf/RR-1
   lnP=0
   IF(M1maxNS>0) lnP=lnP-(MM-M1maxNS)**2*0.5D0/sM1maxNS**2
   IF(M2maxNS>0) lnP=lnP-(MM-M2maxNS)**2*0.5D0/sM2maxNS**2
   IF(M3maxNS>0) lnP=lnP-(MM-M3maxNS)**2*0.5D0/sM3maxNS**2
   IF(RsNS>0.AND.MM>=MmnNS) lnP=lnP-(RR-RsNS)**2*0.5D0/sRsNS**2
   IF(M1>0) lnP=lnP-(MM-M1)**2*0.5D0/sM1**2
   IF(R1>0) lnP=lnP-(RR-R1)**2*0.5D0/sR1**2
   IF(R1inf>0) lnP=lnP-(RRinf-R1inf)**2*0.5D0/sR1inf**2
   IF(R2inf>0) lnP=lnP-(RRinf-R2inf)**2*0.5D0/sR2inf**2
   IF(z2>0) lnP=lnP-(zz-z2)**2*0.5D0/sz2**2
   dl=SQRT(((NSchain(k,1)-NSchain(kk,1))*Rkm)**2+(NSchain(k,2)-NSchain(kk,2))**2)
   PP=PP+EXP(lnP)*dl
  END DO
!$OMP END PARALLEL DO
  IF(PP==0) THEN
   lnPmax=lnPmin
  ELSE
   lnPmax=LOG(PP)
  ENDIF
  IF(NucPh) THEN
   kk=1
   DO WHILE(NSchain(kk,3)<xlEDF)
    kk=kk+1
   ENDDO
   kkk=kk
   DO WHILE(NSchain(kkk+1,3)<xrEDF)
    kkk=kkk+1
   ENDDO
   cnt=0
   PP=0
   ALLOCATE(hist(kkk-kk+1,2))
! OMP PARALLEL DO PRIVATE(k,i,j) REDUCTION(+:PP,cnt)
   DO k=kk,kkk
    i=INT((NSchain(k,3)-xlEDF-dxEDF/2)/dxEDF)+1
    j=INT((NSchain(k,4)-ydEDF-dyEDF/2)/dyEDF)+1
    IF(j<1.OR.j>NyEDF) CYCLE
! OMP CRITICAL
    IF(.NOT.gridEDF(i,j)) THEN
     gridEDF(i,j)=.TRUE. 
     cnt=cnt+1
     hist(cnt,1)=i
     hist(cnt,2)=j
     PP=PP+pEDF(i,j)
    ENDIF
! OMP END CRITICAL
   END DO
! OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,cnt
    gridEDF(hist(i,1),hist(i,2))=.FALSE.
   ENDDO
!$OMP END PARALLEL DO
   DEALLOCATE(hist)
   IF(PP==0) THEN
    lnPmax=lnPmin
   ELSE
    lnPmax=lnPmax+LOG(PP)
   ENDIF
  ENDIF
  IF(Love) THEN
!$OMP PARALLEL DO PRIVATE(k,MM,RR) REDUCTION(MIN:lnPmax)
   DO k=1,strmax
    RR=NSchain(k,5)
    IF(NSchain(k,1)<Mlimit.OR.NSchain(k,2)>Rlimit.OR.(RR.NE.RR)) CYCLE
    IF(RR>LoveLimit) lnPmax=MIN(lnPmin,lnPmax)
   END DO
!$OMP END PARALLEL DO
  ENDIF
  IF(MomI.AND.Ntrn>0.AND.glitchAmax/=0) THEN
   IF(Icfmx<glitchAmax) lnPmax=lnPmin
  ENDIF
  IF(ISNAN(lnPmax)) lnPmax=lnPmin
 END SUBROUTINE Bayesian
END MODULE Analysis