MODULE Data
 USE Params, ONLY : db, scratch, wflag, gcm3,dyncm2, &
	    rhos, xMax, xEOSmx, yMax, NucPh,EDtoKM, Rkm
 IMPLICIT NONE
 SAVE
 CHARACTER(LEN=50), PUBLIC :: fileEOS='sly4.dat',fileEDF='UNEDF1.map',EDFname=''
 INTEGER, PUBLIC :: dotsEOS=152, jumps,units=1
 INTEGER :: iost, NxEDF,NyEDF
 REAL(db), ALLOCATABLE, PUBLIC :: xEOS(:),yEOS(:),xEDF(:),yEDF(:),pEDF(:,:),nEOS(:)
 LOGICAL, ALLOCATABLE, PUBLIC :: gridEDF(:,:)
 REAL(db) :: xMin,dxEDF,dyEDF,xlEDF,xrEDF,ydEDF,yuEDF
 REAL(db) :: RsNS=9.1,sRsNS=1.5,MmnNS=0.4
 REAL(db) :: M1maxNS=1.97,sM1maxNS=0.04, &
	     M2maxNS=2.01,sM2maxNS=0.04, &
	     M3maxNS=2.44,sM3maxNS=0.27, &
	     MmaxNS=2.0
 REAL(db) :: M1=1.33, sM1=0.09, R1=11.9, sR1=0.8, R1inf=14.6,sR1inf=1
 REAL(db) :: R2inf=16.1,sR2inf=1.8,z2=0.37,sz2=0.03
 REAL(db) :: LoveLimit=1000000, Mlimit=0.5, Rlimit=25
 REAL(db) :: glitchAmax=0.07, glitchRmax=25.0, glitchPmx=1
CONTAINS
!SUBROUTINE loadSTARS
! 
!END SUBROUTINE loadSTARS

SUBROUTINE loadEOS
  INTEGER :: i,indx
  REAL(db) :: nn,rr,pp, kr, kp, lns
  IF(wflag) WRITE(*,*) 'Loading ',dotsEOS,' EOS data points from ',fileEOS
  OPEN(UNIT=scratch,FILE=TRIM('DAT/'//fileEOS),FORM='FORMATTED',IOSTAT=iost)
  IF(iost/=0) STOP 'Error when reading EOS data file'
  ALLOCATE(xEOS(dotsEOS),yEOS(dotsEOS),nEOS(dotsEOS))
  lns=LOG(rhos)
  IF(units==1) THEN
   IF(wflag) WRITE(*,*) ' units=',units,' pressure in dyn/cm^2 and density in g/cm^3'
   kp=LOG(dyncm2)-lns ! SLy4, FPS, etc.
   kr=LOG(gcm3)-lns   ! SLy4, FPS
  ELSEIF(units==0) THEN
   IF(wflag) WRITE(*,*) ' units=',units,' pressure and density in M_Sun km^-3'
   kp=-lns+LOG(Rkm)  ! AP4, PAL32 etc.
   kr=-lns+LOG(Rkm)  ! AP4, PAL32
  ELSEIF(units==2) THEN
   IF(wflag) WRITE(*,*) ' units=',units,' pressure and density in km^-2'
   kp=-lns ! SYN0
   kr=-lns ! SYN0
  ELSEIF(units==3) THEN
   IF(wflag) WRITE(*,*) ' units=',units,' pressure in MeV/fm^3 and density in MeV/fm^3'
   kp=LOG(EDtoKM)-lns
   kr=LOG(EDtoKM)-lns
  ELSE
   STOP 'Specify units of input seed EOS or edit data.f90 file to add new ones'
  ENDIF
  DO i=1,dotsEOS
   READ(scratch,*,END=11) indx,nn,rr,pp
   xEOS(i)=LOG(pp)+kp
   yEOS(i)=LOG(rr)+kr
   nEOS(i)=nn
  END DO
11  CLOSE(UNIT=scratch)
  IF(wflag) WRITE(*,*) 'Sorting input EOS according to pressure'
  CALL sort_it(xEOS,yEOS,dotsEOS)
  xMin=xEOS(1)
  IF(wflag) WRITE(*,*) 'xMin=',xMin
  xEOSmx=xEOS(dotsEOS)
  IF(wflag) WRITE(*,*) 'xEOSmx=',xEOSmx
  IF(xMax>xEOS(dotsEOS).AND.wflag) WRITE(*,*) 'Warning: xMax=',xMax,' > ',xEOSmx, ', causual extrapolation is used'
  ! Extrapolation with light speed for sound
  IF(xMax>xEOSmx) THEN
   yMax=yEOS(dotsEOS)+(xMax-xEOSmx)
  ENDIF
END SUBROUTINE loadEOS
SUBROUTINE loadEDF
  INTEGER :: i,j
  REAL(db) :: xxx,yyy,ppp, lns
  IF(wflag) WRITE(*,*) 'Loading data from ',fileEDF
  OPEN(UNIT=scratch,FILE=TRIM('DAT/'//fileEDF),FORM='FORMATTED',IOSTAT=iost)
  IF(iost/=0) STOP 'Error when reading EDF data file'
  READ(scratch,*) EDFname
  IF(wflag) WRITE(*,*) 'EDF used is ',TRIM(EDFname)
  READ(scratch,*) dyEDF,dxEDF,NyEDF,NxEDF
  IF(wflag) WRITE(*,*) 'EDF map resolution is ',dxEDF,dyEDF,NxEDF,NyEDF
  ALLOCATE(xEDF(NxEDF),yEDF(NyEDF),pEDF(NxEDF,NyEDF),gridEDF(NxEDF,NyEDF))
  READ(scratch,*) ydEDF,yuEDF,xlEDF,xrEDF
  lns=LOG(rhos/EDtoKM)
  xlEDF=xlEDF-lns
  xrEDF=xrEDF-lns
  ydEDF=ydEDF-lns
  yuEDF=yuEDF-lns
  IF(wflag) WRITE(*,*) 'EDF map boundaries are ',xlEDF,xrEDF,ydEDF,yuEDF
  IF(wflag) WRITE(*,*) 'EDF: assumed units are ln(MeV/fm^3)'
  DO i=1,NxEDF
   DO j=1,NyEDF
    READ(scratch,*,END=12) xxx,yyy,ppp
    pEDF(i,j)=ppp
    gridEDF(i,j)=.FALSE.
    IF(i==1) yEDF(j)=xxx-lns
   ENDDO
   xEDF(i)=yyy-lns
  END DO
12  CLOSE(UNIT=scratch)
END SUBROUTINE loadEDF

SUBROUTINE kill_data
  IF(dotsEOS>0) DEALLOCATE(xEOS,yEOS)
  IF(NucPh) DEALLOCATE(xEDF,yEDF,pEDF,gridEDF)
END SUBROUTINE kill_data

SUBROUTINE sort_it(a,b,n)
  REAL(db), INTENT(INOUT):: a(:),b(:)
  INTEGER,INTENT(IN) :: n
  REAL(db) :: t(n,2), tt(2)
  INTEGER :: i,j
!$OMP PARALLEL DO PRIVATE(i)
  DO i=1,n
   t(i,1)=a(i)
   t(i,2)=b(i)
  END DO
!$OMP END PARALLEL DO
  i=1
  DO WHILE(i<n)
   j=i+1
   DO WHILE(j<=n)
    IF(t(j,1)<t(i,1)) EXIT
    j=j+1
   ENDDO
   IF(j>n) THEN
    i=i+1
   ELSE
    tt(:)=t(i,:)
    t(i,:)=t(j,:)
    t(j,:)=tt(:)
   ENDIF
  ENDDO
!$OMP PARALLEL DO PRIVATE(i)
  DO i=1,n
   a(i)=t(i,1)
   b(i)=t(i,2)
  END DO
!$OMP END PARALLEL DO
  j=0
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:j)
  DO i=1,n-1
   IF(a(i)==a(i+1)) j=j+1
  END DO
!$OMP END PARALLEL DO
  jumps=j
  IF(wflag) WRITE(*,*) 'Number of density jumps: ',jumps
END SUBROUTINE sort_it

END MODULE Data