MODULE Chain
 USE Params, ONLY : db, do, xseed, yseed, iseed, CRUSTknots, COREknots, linkpoints, &
		    wflag,scratch,scratch2,rhos,EDtoKM, step, restart, lnPmin, &
		    Love, MomI, nDen, Rkm, NusrX,NusrY,Ugrid, &
		    BE, addUSR, step0
 USE Parallel, ONLY : mpi_myproc
 USE Data
 USE User
IMPLICIT NONE
 SAVE
 INTEGER, PUBLIC :: Nknots, writeevery=10,NNS,Nx,Ny,NM,NR, No, writeNS=10, MINknot=2, NLove
 INTEGER, PUBLIC :: Nfrc, Ntrn, strmax
 REAL(db) :: mMaxAll,mMax, lnPal, Icfmx
 REAL(db), PUBLIC :: lnPmax=-9.9D-99,xL,xR,yD,yU,ML,MR,RD,RU,dxx,dy,dM,dR,cut=3,djump=0.5,xCHmax,dyR=0
 REAL(db), PUBLIC :: LoveL,LoveR,dLove, frcU,frcD,dfrc, trnL,trnR,dtrn
 REAL(db),ALLOCATABLE,PUBLIC :: NSchain(:,:),ychain(:),xchain(:),fchain(:), &
				xall(:),yall(:),band(:,:), &
				Mall(:),Rall(:),pmap(:,:), &
				Loveall(:),lmap(:,:), &
				trnall(:), frcall(:), imap(:,:)
 LOGICAL, PUBLIC :: addBAND=.FALSE.,addMAP=.FALSE.,addLOVE=.FALSE.,addMOM=.FALSE.
 LOGICAL, ALLOCATABLE,PRIVATE :: grid(:,:),gridEOS(:,:),Lgrid(:,:),Igrid(:,:,:)
 REAL(db), ALLOCATABLE,PUBLIC :: profile(:,:,:)
! CHARACTER(LEN=100), PUBLIC :: fileBAND='band.all',fileMAP='pmap.all'
CONTAINS
SUBROUTINE alloc_chain
  Nknots=3+CRUSTknots+COREknots
  No=2+CRUSTknots
  ALLOCATE(xchain(Nknots),ychain(Nknots))
  NNS=(Nknots-1)*linkpoints+1
  ALLOCATE(fchain(NNS))
  IF(BE) THEN
   ALLOCATE(NSchain(NNS,8))
  ELSEIF(nDen) THEN
   ALLOCATE(NSchain(NNS,7))
  ELSEIF(MomI) THEN
   ALLOCATE(NSchain(NNS,6))
  ELSEIF(Love) THEN
    ALLOCATE(NSchain(NNS,5))
  ELSE
    ALLOCATE(NSchain(NNS,4))
  ENDIF
  IF(Nx>0) ALLOCATE(xall(Nx),yall(Ny),band(Nx,Ny))
  IF(NM>0) ALLOCATE(Mall(NM),Rall(NR),pmap(NM,NR))
  IF(Love) ALLOCATE(Loveall(NLove),lmap(NLove,NR))
  IF(MomI.AND.Ntrn>0) ALLOCATE(frcall(Nfrc),trnall(Ntrn),imap(Ntrn,Nfrc))
END SUBROUTINE alloc_chain

SUBROUTINE kill_chain
 IF(wflag) WRITE(*,*) 'Killing chain'
 DEALLOCATE(ychain,xchain,NSchain,fchain)
 IF(Nx>0) DEALLOCATE(xall,yall,band)
 IF(NM>0) DEALLOCATE(Mall,Rall,pmap)
 IF(Nfrc>0) DEALLOCATE(frcall,trnall,imap)
END SUBROUTINE kill_chain

 SUBROUTINE init_profile
  ALLOCATE(profile(iseed+NNS,5,NNS))
 END SUBROUTINE init_profile

 SUBROUTINE kill_profile
  DEALLOCATE(profile)
 END SUBROUTINE kill_profile

SUBROUTINE init_grid
  INTEGER :: i,j
  ALLOCATE(grid(NM,NR))
!$OMP PARALLEL DO PRIVATE(i,j)
  DO i=1,NM
   DO j=1,NR
    grid(i,j)=.FALSE.
   ENDDO
  END DO
!$OMP END PARALLEL DO
 IF(MomI.AND.Ntrn>0)THEN
   ALLOCATE(Igrid(Ntrn,Nfrc,NNS))
!$OMP PARALLEL DO PRIVATE(i,j)
   DO i=1,Ntrn
    DO j=1,Nfrc
     Igrid(i,j,:)=.FALSE.
    ENDDO
   END DO
!$OMP END PARALLEL DO
 ENDIF
 IF(.NOT.Love) RETURN
  ALLOCATE(Lgrid(NLove,NR))
!$OMP PARALLEL DO PRIVATE(i,j)
  DO i=1,NLove
   DO j=1,NR
    Lgrid(i,j)=.FALSE.
   ENDDO
  END DO
!$OMP END PARALLEL DO
 lnPal=lnPmin
 IF(addUSR) THEN
    CALL loadUSR(lnPal)
 ELSEIF(NusrX>0) THEN
    CALL init_user
 ENDIF
END SUBROUTINE init_grid
SUBROUTINE init_gridEOS
  INTEGER :: i,j
  ALLOCATE(gridEOS(Nx,Ny))
!$OMP PARALLEL DO PRIVATE(i,j)
  DO i=1,Nx
   DO j=1,Ny
    gridEOS(i,j)=.FALSE.
   ENDDO
  END DO
!$OMP END PARALLEL DO
END SUBROUTINE init_gridEOS

SUBROUTINE kill_grids
  IF(NM>0) DEALLOCATE(grid)
  IF(Nx>0) DEALLOCATE(gridEOS)
  IF(Love) DEALLOCATE(Lgrid)
  IF(MomI.AND.Ntrn>0) DEALLOCATE(Igrid)
END SUBROUTINE kill_grids

SUBROUTINE init_chain(lnPmap,lnPband)
 INTEGER :: i,j
 REAL(db), INTENT(OUT) :: lnPmap,lnPband
 REAL(db) :: yy,qw
 IF(wflag) WRITE(*,*) 'Init. chain'
 lnPmap=lnPmin
 lnPband=lnPmin
 IF(Nx>0) THEN
  IF(wflag) WRITE(*,*) 'Init. band'
  IF(addBAND) THEN
   CALL loadBAND(lnPband)
   lnPal=MAX(lnPal,lnPband)
  ELSE
   dxx=(xR-xL)/Nx
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,Nx
    xall(i)=xL+(i-0.5D0)*dxx
   ENDDO
!$OMP END PARALLEL DO
   dy=(yU-yD)/Ny
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,Ny
    yall(i)=yD+(i-0.5D0)*dy
   ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i,j)
   DO j=1,Ny
    DO i=1,Nx
     band(i,j)=lnPmin
    END DO
   END DO
!$OMP END PARALLEL DO
  ENDIF
 ENDIF
 IF(NM>0) THEN
  IF(wflag) WRITE(*,*) 'Init. map'
  IF(addMAP) THEN
   CALL loadMAP(lnPmap)
   lnPal=MAX(lnPal,lnPmap)
  ELSE
   dM=(MR-ML)/NM
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,NM
    Mall(i)=ML+(i-0.5D0)*dM
   ENDDO
!$OMP END PARALLEL DO
   dR=(RU-RD)/NR
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,NR
    Rall(i)=RD+(i-0.5D0)*dR
   ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i,j)
   DO j=1,NR
    DO i=1,NM
     pmap(i,j)=lnPmin
    END DO
   END DO
!$OMP END PARALLEL DO
  ENDIF
 ENDIF
 IF(NLove>0) THEN
  IF(wflag) WRITE(*,*) 'Init. Love map'
  IF(addLOVE) THEN
   CALL loadLOVE(lnPmap)
   lnPal=MAX(lnPal,lnPmap)
  ELSE
   dLove=(LoveR-LoveL)/NLove
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,NLove
    Loveall(i)=LoveL+(i-0.5D0)*dLove
   ENDDO
!$OMP END PARALLEL DO
  IF(NM==0) THEN
   dR=(RU-RD)/NR
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,NR
    Rall(i)=RD+(i-0.5D0)*dR
   ENDDO
!$OMP END PARALLEL DO
  ENDIF
!$OMP PARALLEL DO PRIVATE(i,j)
   DO j=1,NR
    DO i=1,NM
     lmap(i,j)=lnPmin
    END DO
   END DO
!$OMP END PARALLEL DO
  ENDIF
 ENDIF
 IF(Nfrc>0) THEN
  IF(wflag) WRITE(*,*) 'Init. I map'
  IF(addMOM) THEN
   CALL loadMOM(lnPmap)
   lnPal=MAX(lnPal,lnPmap)
  ELSE
   dfrc=(frcU-frcD)/Nfrc
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,Nfrc
    frcall(i)=frcD+(i-0.5D0)*dfrc
   ENDDO
!$OMP END PARALLEL DO
   dtrn=(trnR-trnL)/Ntrn
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,Ntrn
    trnall(i)=trnL+(i-0.5D0)*dtrn
   ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i,j)
   DO j=1,Nfrc
    DO i=1,Ntrn
     imap(i,j)=lnPmin
    END DO
   END DO
!$OMP END PARALLEL DO
  ENDIF
 ENDIF
 IF(.NOT.restart) THEN
  WRITE(*,*) mpi_myproc,': Init. chain from seed EoS'
  xchain(1)=xseed
  ychain(1)=yseed
  DO i=1,CRUSTknots
   xchain(i+1)=xseed-i*xseed/(1+CRUSTknots)
  ENDDO
  i=2
  ns: DO j=iseed+1,dotsEOS
   qw=(yEOS(j)-yEOS(j-1))/(xEOS(j)-xEOS(j-1))
   DO WHILE(xEOS(j)>=xchain(i))
    ychain(i)=yEOS(j-1)+qw*(xchain(i)-xEOS(j-1))
!    dchain(i)=yy-LOG(do+EXP(xchain(i)))
    i=i+1
    IF(i>CRUSTknots+1) EXIT ns
   ENDDO
  END DO ns
  
  DO WHILE(xEOS(j)<0)
   j=j+1
  ENDDO

  xchain(No)=0
  ychain(No)=yEOS(j-1)-(yEOS(j)-yEOS(j-1))*xEOS(j-1)/(xEOS(j)-xEOS(j-1))

  DO i=1,COREknots
   xchain(No+i)=i*xEOS(dotsEOS)/(1+COREknots)
  ENDDO
  i=No+1
  ns2: DO WHILE(j<=dotsEOS)
   qw=(yEOS(j)-yEOS(j-1))/(xEOS(j)-xEOS(j-1))
   DO WHILE(xEOS(j)>xchain(i))
    ychain(i)=yEOS(j-1)+qw*(xchain(i)-xEOS(j-1))
!    dchain(i)=yy-LOG(do+EXP(xchain(i)))
    i=i+1
    IF(i>Nknots-1) EXIT ns2
   ENDDO
   j=j+1
  END DO ns2
  IF(xMax>xEOS(dotsEOS)) THEN
   xchain(Nknots)=xMax
   ychain(Nknots)=yMax
   xCHmax=xMax
  ELSE
   xchain(Nknots)=xEOS(dotsEOS)
   ychain(Nknots)=yEOS(dotsEOS)
   xCHmax=xchain(Nknots)
  ENDIF
!$OMP PARALLEL DO PRIVATE(i,j,yy,qw)
   DO i=1,Nknots-1
    DO j=1,linkpoints
     qw=xchain(i)+j*(xchain(i+1)-xchain(i))/linkpoints
     NSchain((i-1)*linkpoints+j+1,3)=qw
     yy=ychain(i)+j*(ychain(i+1)-ychain(i))/linkpoints
     NSchain((i-1)*linkpoints+j+1,4)=yy
     fchain((i-1)*linkpoints+j+1)=1/(1+EXP(yy-qw))
    ENDDO
   ENDDO
!$OMP END PARALLEL DO
 ELSE
  lnPmax=MAX(lnPmap,lnPband)
  CALL loadNS
!$OMP PARALLEL DO PRIVATE(i,j)
   DO i=1,Nknots
    j=(i-1)*linkpoints+1
    xchain(i)=NSchain(j,3)
    ychain(i)=NSchain(j,4)
    fchain(i)=1.0D0/(1+EXP(ychain(i)-xchain(i)))
   ENDDO
!$OMP END PARALLEL DO
 ENDIF
END SUBROUTINE init_chain

SUBROUTINE saveNS(lnPa)
  REAL(db) :: lnPa
   INTENT(IN) :: lnPa
  INTEGER :: i, iost
  CHARACTER(13) :: names
  REAL(db) :: kk
  WRITE(names,'(I6.6,A1,I6.6)') step,'.',mpi_myproc
  WRITE(*,*) mpi_myproc,': writing ',names
  OPEN(UNIT=scratch,FILE=TRIM('NS/'//names),FORM='FORMATTED',IOSTAT=iost,STATUS='REPLACE')
  IF(iost/=0) STOP 'Error when saving NS data file'
  WRITE(scratch,'(A)') fileEOS
  WRITE(scratch,'(3E16.6)') lnPa, rhos, xseed
  WRITE(scratch,'(I6)') strmax
  kk=rhos/EDtoKM
  IF(BE) THEN
   DO i=1,strmax
    WRITE(scratch,'(8E16.6)') NSchain(i,1),NSchain(i,2),kk*EXP(NSchain(i,3)),&
    kk*EXP(NSchain(i,4)),NSchain(i,5),NSchain(i,6),NSchain(i,7),NSchain(i,8)
   ENDDO
  ELSEIF(nDen) THEN
   DO i=1,strmax
    WRITE(scratch,'(7E16.6)') NSchain(i,1),NSchain(i,2),kk*EXP(NSchain(i,3)),&
    kk*EXP(NSchain(i,4)),NSchain(i,5),NSchain(i,6),NSchain(i,7)
   ENDDO
  ELSEIF(MomI) THEN
   DO i=1,strmax
    WRITE(scratch,'(6E16.6)') NSchain(i,1),NSchain(i,2),kk*EXP(NSchain(i,3)),kk*EXP(NSchain(i,4)),NSchain(i,5),NSchain(i,6)
   ENDDO
  ELSEIF(Love) THEN
   DO i=1,strmax
    WRITE(scratch,'(5E16.6)') NSchain(i,1),NSchain(i,2),kk*EXP(NSchain(i,3)),kk*EXP(NSchain(i,4)),NSchain(i,5)
   ENDDO
  ELSE
   DO i=1,strmax
    WRITE(scratch,'(4E16.6)') NSchain(i,1),NSchain(i,2),kk*EXP(NSchain(i,3)),kk*EXP(NSchain(i,4))
   ENDDO
  ENDIF
  CLOSE(UNIT=scratch)
END SUBROUTINE saveNS
SUBROUTINE saveMAX(kmx)
  INTEGER, INTENT(IN) :: kmx
  INTEGER :: i, iost
  CHARACTER(13) :: names
  REAL(db) :: kk
  LOGICAL :: exists
  kk=rhos/EDtoKM
  WRITE(names,'(A4,I6.6)') 'max.',mpi_myproc
  WRITE(*,*) mpi_myproc,': writing ',names
  INQUIRE(FILE=TRIM('MAX/'//names),EXIST=exists)
  IF(exists.AND.step>0) THEN
   OPEN(unit=scratch,FILE=TRIM('MAX/'//names),POSITION='APPEND')
  ELSE
   OPEN(UNIT=scratch,FILE=TRIM('MAX/'//names),FORM='FORMATTED',STATUS='REPLACE')
  ENDIF
  IF(BE) THEN
    WRITE(scratch,'(I10,8E16.6)') step, lnPmax, NSchain(kmx,1), NSchain(kmx,2), kk*EXP(NSchain(kmx,3)),&
    kk*EXP(NSchain(kmx,4)),NSchain(kmx,5),NSchain(kmx,6),NSchain(kmx,7),NSchain(kmx,8)
  ELSEIF(nDen) THEN
    WRITE(scratch,'(I10,8E16.6)') step, lnPmax, NSchain(kmx,1), NSchain(kmx,2), kk*EXP(NSchain(kmx,3)),&
    kk*EXP(NSchain(kmx,4)),NSchain(kmx,5),NSchain(kmx,6),NSchain(kmx,7)
  ELSEIF(MomI) THEN
    WRITE(scratch,'(I10,7E16.6)') step, lnPmax, NSchain(kmx,1), NSchain(kmx,2), kk*EXP(NSchain(kmx,3)),&
    kk*EXP(NSchain(kmx,4)),NSchain(kmx,5),NSchain(kmx,6)
  ELSEIF(Love) THEN
    WRITE(scratch,'(I10,6E16.6)') step, lnPmax, NSchain(kmx,1), NSchain(kmx,2), kk*EXP(NSchain(kmx,3)),&
    kk*EXP(NSchain(kmx,4)),NSchain(kmx,5)
  ELSE
    WRITE(scratch,'(I10,5E16.6)') step, lnPmax, NSchain(kmx,1), NSchain(kmx,2), kk*EXP(NSchain(kmx,3)),&
    kk*EXP(NSchain(kmx,4))
  ENDIF
  CLOSE(UNIT=scratch)
END SUBROUTINE saveMAX
SUBROUTINE loadNS
  INTEGER :: i, iost
  CHARACTER(13) :: names
  REAL(db) :: kk, aa,bb,gg,hh,vv,iii,nnn
  WRITE(names,'(I6.6,A1,I6.6)') step,'.',mpi_myproc
  WRITE(*,*) mpi_myproc,': loading ',names
  OPEN(UNIT=scratch,FILE=TRIM('NSr/'//names),FORM='FORMATTED',IOSTAT=iost)
  IF(iost/=0) STOP 'Error when loading NS data file'
  READ(scratch,*)
  READ(scratch,*)
  READ(scratch,*) strmax
  kk=rhos/EDtoKM
  DO i=1,strmax
   IF(nDen) THEN
    READ(scratch,*,END=1) aa,bb,gg,hh,vv,iii,nnn
    NSchain(i,5)=vv
    NSchain(i,6)=iii
    NSchain(i,7)=nnn
   ELSEIF(MomI) THEN
    READ(scratch,*,END=1) aa,bb,gg,hh,vv,iii
    NSchain(i,5)=vv
    NSchain(i,6)=iii
   ELSEIF(Love) THEN
    READ(scratch,*,END=1) aa,bb,gg,hh,vv
    NSchain(i,5)=vv
   ELSE
    READ(scratch,*,END=1) aa,bb,gg,hh
   ENDIF
   NSchain(i,1)=aa
   NSchain(i,2)=bb
   NSchain(i,3)=LOG(gg/kk)
   NSchain(i,4)=LOG(hh/kk)
  END DO
1  CLOSE(UNIT=scratch)
! RETURN
END SUBROUTINE loadNS

SUBROUTINE updateBAND
  INTEGER :: i,j,k,cnt
  REAL(db) :: lnPh
  INTEGER :: hist(NNS,2)
  cnt=0
  DO k=1,NNS
   i=INT((NSchain(k,3)-xL-dxx/2)/dxx)+1
   IF(i<1.OR.i>Nx) CYCLE
   j=INT((NSchain(k,4)-yD-dy/2)/dy)+1
   IF(j<1.OR.j>Ny) CYCLE
   IF(gridEOS(i,j)) CYCLE
   cnt=cnt+1
   gridEOS(i,j)=.TRUE.
   hist(cnt,1)=i
   hist(cnt,2)=j
   lnPh=lnPmax-band(i,j)
   IF(lnPh<100) THEN
     band(i,j)=band(i,j)+LOG(1+EXP(lnPh))
   ELSE
     band(i,j)=band(i,j)+lnPh
   ENDIF
  END DO
!$OMP PARALLEL DO PRIVATE(i)
  DO i=1,cnt
   gridEOS(hist(i,1),hist(i,2))=.FALSE.
  ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE updateBAND
SUBROUTINE updateMOM
  INTEGER :: i,j,k,cnt,o
  REAL(db) :: lnPh,cfc,frc,ptr
  INTEGER :: hist(NNS+iseed,2)
  cnt=0
  cfc=rhos/EDtoKM
  Icfmx=0
!$OMP PARALLEL DO PRIVATE(o,cnt,k,i,j,hist,lnPh,frc,ptr) REDUCTION(MAX:Icfmx)
 DO o=1,strmax
  cnt=0
  IF(RU>0.AND.NSchain(o,2)>RU) CYCLE
  DO k=1,o+iseed
   ptr=cfc*EXP(profile(k,1,o))
   i=INT((ptr-trnL-dtrn/2)/dtrn)+1
   IF(i<1.OR.i>Ntrn) CYCLE
   frc=1-profile(k,5,o)/profile(1,5,o)
   j=INT((frc-frcD-dfrc/2)/dfrc)+1
   IF(j<1.OR.j>Nfrc) CYCLE
   IF(NSchain(o,1)>0.AND.NSchain(o,2)<=glitchRmax.AND.ptr<=glitchPmx) Icfmx=MAX(Icfmx,frc)
   IF(Igrid(i,j,o)) CYCLE
   cnt=cnt+1
   Igrid(i,j,o)=.TRUE.
   hist(cnt,1)=i
   hist(cnt,2)=j
!$OMP CRITICAL
   lnPh=lnPmax-imap(i,j)
   IF(lnPh<100) THEN
     imap(i,j)=imap(i,j)+LOG(1+EXP(lnPh))
   ELSE
     imap(i,j)=imap(i,j)+lnPh
   ENDIF
!$OMP END CRITICAL
  END DO
  DO i=1,cnt
   Igrid(hist(i,1),hist(i,2),o)=.FALSE.
  ENDDO
 ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE updateMOM

SUBROUTINE updateLOVE
  INTEGER :: i,j,k,cnt
  REAL(db) :: lnPh
  INTEGER :: hist(NNS,2)
  cnt=0
  DO k=1,strmax
   IF(k>1) THEN
    IF(NSchain(k,1)<0.OR.NSchain(k,2)>RU.OR.NSchain(k,1)<NSchain(k-1,1)) CYCLE
   ENDIF
   i=INT((NSchain(k,5)-LoveL-dLove/2)/dLove)+1
   IF(i<1.OR.i>NLove) CYCLE
   j=INT((NSchain(k,2)-RD-dR/2)/dR)+1
   IF(j<1.OR.j>NR) CYCLE
   IF(Lgrid(i,j)) CYCLE
   cnt=cnt+1
   Lgrid(i,j)=.TRUE.
   hist(cnt,1)=i
   hist(cnt,2)=j
   lnPh=lnPmax-lmap(i,j)
   IF(lnPh<100) THEN
     lmap(i,j)=lmap(i,j)+LOG(1+EXP(lnPh))
   ELSE
     lmap(i,j)=lmap(i,j)+lnPh
   ENDIF
  END DO
!$OMP PARALLEL DO PRIVATE(i)
  DO i=1,cnt
   Lgrid(hist(i,1),hist(i,2))=.FALSE.
  ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE updateLOVE

SUBROUTINE updateMAP
  INTEGER :: i,j,k,cnt,kmx
  REAL(db) :: lnPh,mx
  INTEGER :: hist(NNS,2)
  cnt=0
  mx=0
  DO k=1,NNS
   IF(k>1) THEN
    IF(NSchain(k,1)<0.OR.NSchain(k-1,1)<0.OR.NSchain(k,1)<NSchain(k-1,1)) CYCLE
   ENDIF
   i=INT((NSchain(k,1)-ML-dM/2)/dM)+1
   IF(i<1.OR.i>NM) CYCLE
   j=INT((NSchain(k,2)-RD-dR/2)/dR)+1
   IF(j<1.OR.j>NR) CYCLE
   IF(NSchain(k,1)>mx) THEN
    mx=NSchain(k,1)
    kmx=k
   ENDIF
   IF(grid(i,j)) CYCLE
   cnt=cnt+1
   grid(i,j)=.TRUE.
   hist(cnt,1)=i
   hist(cnt,2)=j
   lnPh=lnPmax-pmap(i,j)
   IF(lnPh<100) THEN
     pmap(i,j)=pmap(i,j)+LOG(1+EXP(lnPh))
   ELSE
     pmap(i,j)=pmap(i,j)+lnPh
   ENDIF
  END DO
  IF(mx>0) THEN
    mMax=mx
    CALL saveMAX(kmx)
  ENDIF
!$OMP PARALLEL DO PRIVATE(i)
  DO i=1,cnt
   grid(hist(i,1),hist(i,2))=.FALSE.
  ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE updateMAP

SUBROUTINE saveBAND(lnPa)
  INTEGER :: i,j, iost
  REAL(db) :: z,lnPa
   INTENT(IN) :: lnPa
  CHARACTER(13) :: names
  WRITE(names,'(A5,I6.6)') 'band.',mpi_myproc
  WRITE(*,*) mpi_myproc,': writing band'
  OPEN(UNIT=scratch,FILE=TRIM('BND/'//names),FORM='FORMATTED',IOSTAT=iost,STATUS='REPLACE')
  IF(iost/=0) STOP 'Error when saving band file'
  WRITE(scratch,'(A,E16.6)') fileEOS, rhos*EXP(xseed)/EDtoKM
  WRITE(scratch,'(3E16.6,2I16,E16.6)') lnPa,dxx,dy, Nx, Ny, rhos
  WRITE(scratch,'(4E16.6)') xL,xR,yD,yU
  DO i=1,Nx
   DO j=1,Ny
    WRITE(scratch,'(3E16.6)') yall(j),xall(i),band(i,j)
   ENDDO
  END DO
  CLOSE(UNIT=scratch)
END SUBROUTINE saveBAND
SUBROUTINE saveMOM(lnPa)
  INTEGER :: i,j, iost
  REAL(db) :: z,lnPa
   INTENT(IN) :: lnPa
  CHARACTER(13) :: names
  WRITE(names,'(A5,I6.6)') 'ivst.',mpi_myproc
  WRITE(*,*) mpi_myproc,': writing I map'
  OPEN(UNIT=scratch,FILE=TRIM('MOM/'//names),FORM='FORMATTED',IOSTAT=iost,STATUS='REPLACE')
  IF(iost/=0) STOP 'Error when saving I map file'
  WRITE(scratch,'(A,E16.6,I16)') fileEOS, rhos*EXP(xseed)/EDtoKM,step
  WRITE(scratch,'(3E16.6,2I16,E16.6)') lnPa,dtrn,dfrc, Ntrn, Nfrc, rhos
  WRITE(scratch,'(4E16.6)') trnL,trnR,frcD,frcU
  DO j=1,Nfrc
   DO i=1,Ntrn
    WRITE(scratch,'(3E16.6)') trnall(i),frcall(j),imap(i,j)
   ENDDO
  END DO
  CLOSE(UNIT=scratch)
END SUBROUTINE saveMOM
SUBROUTINE loadMOM(lnPa)
  INTEGER :: i,j, iost
  REAL(db) :: rhoss,yy,xx,bb,lnPa,lnPaa
   INTENT(OUT) :: lnPa
  CHARACTER(13) :: names
  WRITE(names,'(A5,I6.6)') 'ivst.',mpi_myproc
  WRITE(*,*) mpi_myproc,': loading I map'
  OPEN(UNIT=scratch,FILE=TRIM('MOM/'//names),FORM='FORMATTED',IOSTAT=iost)
  IF(iost/=0) STOP 'Error when loading I map file'
  READ(scratch,*)
  READ(scratch,*) lnPaa,dtrn,dfrc, Ntrn, Nfrc, rhoss
  lnPa=lnPaa
  IF(rhoss/=rhos) WRITE(*,*) ' WARNING rhoss/=rhos'
  rhos=rhoss
  READ(scratch,*) trnL,trnR,frcD,frcU
  DO j=1,Nfrc
   DO i=1,Ntrn
    READ(scratch,*,END=22) xx,yy,bb
    imap(i,j)=bb
    IF(j==1) trnall(i)=xx
   ENDDO
   frcall(j)=yy
  END DO
22 CLOSE(UNIT=scratch)
END SUBROUTINE loadMOM

SUBROUTINE loadBAND(lnPa)
  INTEGER :: i,j, iost
  REAL(db) :: rhoss,yy,xx,bb,lnPa,lnPaa
   INTENT(OUT) :: lnPa
  CHARACTER(13) :: names
  WRITE(names,'(A5,I6.6)') 'band.',mpi_myproc
  WRITE(*,*) mpi_myproc,': loading band'
  OPEN(UNIT=scratch,FILE=TRIM('BND/'//names),FORM='FORMATTED',IOSTAT=iost)
  IF(iost/=0) STOP 'Error when loading band file'
  READ(scratch,*)
  READ(scratch,*) lnPaa,dxx,dy, Nx, Ny, rhoss
  lnPa=lnPaa
  IF(rhoss/=rhos) WRITE(*,*) '...WARNING rhoss/=rhos'
  rhos=rhoss
  READ(scratch,*) xL,xR,yD,yU
  DO i=1,Nx
   DO j=1,Ny
    READ(scratch,*,END=2) yy,xx,bb
    band(i,j)=bb
    IF(i==1) yall(j)=yy
   ENDDO
   xall(i)=xx
  END DO
2  CLOSE(UNIT=scratch)
! RETURN
END SUBROUTINE loadBAND
SUBROUTINE saveMAP(lnPa)
  INTEGER :: i,j, iost
  REAL(db) :: lnPa
   INTENT(IN) :: lnPa
  CHARACTER(13) :: names
  WRITE(names,'(A5,I6.6)') 'pmap.',mpi_myproc
  WRITE(*,*) mpi_myproc,': writing map'
  OPEN(UNIT=scratch,FILE=TRIM('MAP/'//names),FORM='FORMATTED',IOSTAT=iost,STATUS='REPLACE')
  IF(iost/=0) STOP 'Error when saving map file'
  WRITE(scratch,'(A,E16.6,I16)') fileEOS, rhos*EXP(xseed)/EDtoKM,step
  WRITE(scratch,'(3E16.6,2I6)') lnPa, dM,dR, NM,NR
  WRITE(scratch,'(4E16.6)') ML,MR,RD,RU
  DO i=1,NM
   DO j=1,NR
    WRITE(scratch,'(3E16.6)') Rall(j),Mall(i),pmap(i,j)
   ENDDO
  END DO
  CLOSE(UNIT=scratch)
END SUBROUTINE saveMAP
SUBROUTINE loadMAP(lnPa)
  INTEGER :: i,j, iost,sss
  REAL(db) :: lnPa,rr,mm,pp,lnPaa
   INTENT(OUT) :: lnPa
  CHARACTER(13) :: names
  WRITE(names,'(A5,I6.6)') 'pmap.',mpi_myproc
  WRITE(*,*) mpi_myproc,': reading map'
  OPEN(UNIT=scratch2,FILE=TRIM('MAP/'//names),FORM='FORMATTED',IOSTAT=iost)
  IF(iost/=0) STOP 'Error when loading map file'
  READ(scratch2,*) names,pp,sss
  WRITE(*,*) mpi_myproc,': seed EOS ',names,' seed pressure ',pp,' MeV/fm^3'
  READ(scratch2,*) lnPaa, dM,dR, NM,NR
  lnPa=lnPaa
  READ(scratch2,*) ML,MR,RD,RU
  DO i=1,NM
   DO j=1,NR
    READ(scratch2,*,END=3) rr,mm,pp
    IF(i==1) Rall(j)=rr
    pmap(i,j)=pp
   ENDDO
   Mall(i)=mm
  END DO
3  CLOSE(UNIT=scratch2)
! RETURN
END SUBROUTINE loadMAP
SUBROUTINE saveLOVE(lnPa)
  INTEGER :: i,j, iost
  REAL(db) :: lnPa
   INTENT(IN) :: lnPa
  CHARACTER(13) :: names
  WRITE(names,'(A5,I6.6)') 'lmap.',mpi_myproc
  WRITE(*,*) mpi_myproc,': writing Love map'
  OPEN(UNIT=scratch,FILE=TRIM('LOV/'//names),FORM='FORMATTED',IOSTAT=iost,STATUS='REPLACE')
  IF(iost/=0) STOP 'Error when saving Love file'
  WRITE(scratch,'(A,E16.6,I16)') fileEOS, dyR,step
  WRITE(scratch,'(3E16.6,2I6)') lnPa, dLove,dR, NLove,NR
  WRITE(scratch,'(4E16.6)') LoveL,LoveR,RD,RU
  DO i=1,NLove
   DO j=1,NR
    WRITE(scratch,'(3E16.6)') Loveall(i),Rall(j),lmap(i,j)
   ENDDO
  END DO
  CLOSE(UNIT=scratch)
END SUBROUTINE saveLOVE
SUBROUTINE loadLOVE(lnPa)
  INTEGER :: i,j, iost,sss
  REAL(db) :: lnPa,rr,mm,pp,lnPaa
   INTENT(OUT) :: lnPa
  CHARACTER(13) :: names
  WRITE(names,'(A5,I6.6)') 'lmap.',mpi_myproc
  WRITE(*,*) mpi_myproc,': reading Love'
  OPEN(UNIT=scratch2,FILE=TRIM('LOV/'//names),FORM='FORMATTED',IOSTAT=iost)
  IF(iost/=0) STOP 'Error when loading Love file'
  READ(scratch2,*) names,pp,sss
  WRITE(*,*) mpi_myproc,': seed EOS ',names,' surface jump ',pp,' step=',sss
  READ(scratch2,*) lnPaa, dLove,dR, NLove,NR
  lnPa=lnPaa
  READ(scratch2,*) LoveL,LoveR,RD,RU
  DO i=1,NLove
   DO j=1,NR
    READ(scratch2,*,END=4) mm,rr,pp
    IF(i==1) Rall(j)=rr
    lmap(i,j)=pp
   ENDDO
   Loveall(i)=mm
  END DO
4  CLOSE(UNIT=scratch2)
END SUBROUTINE loadLOVE

END MODULE Chain