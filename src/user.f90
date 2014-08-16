MODULE User
  USE Params
  USE Data, ONLY : fileEOS
  USE Parallel, ONLY : mpi_myproc
  IMPLICIT NONE
  SAVE
CONTAINS
!***********************************************************************
!***********************************************************************
  SUBROUTINE finish_user
  END SUBROUTINE finish_user
!***********************************************************************
  SUBROUTINE updateUSR(lnPP,NSchn,strmx)
    REAL(db),INTENT(IN) :: NSchn(:,:),lnPP
    INTEGER,INTENT(IN) :: strmx
    INTEGER :: i, hist(strmx,2),cnt,j,k
    REAL(db) :: xusr,yusr,lnPh
  cnt=0
  DO k=1,strmx
   IF(k>1) THEN
    IF(NSchn(k,1)<0.OR.NSchn(k,1)<NSchn(k-1,1)) CYCLE
   ENDIF
!   xusr=NSchn(k,2)/(2*NSchn(k,1)*Rkm) ! user input for x
   xusr=NSchn(k,5)/(NSchn(k,1)*Rkm)**5 ! Love bar
   i=INT((xusr-usrXL-usrdX/2)/usrdX)+1
   IF(i<1.OR.i>NusrX) CYCLE
!   yusr=NSchn(k,6)/(NSchn(k,1)*NSchn(k,2)**2) ! user input for y
   yusr=NSchn(k,6)/(NSchn(k,1)**3*Rkm**2) ! I bar
   j=INT((yusr-usrYD-usrdY/2)/usrdY)+1
   IF(j<1.OR.j>NusrY) CYCLE
   IF(Ugrid(i,j)) CYCLE
   cnt=cnt+1
   Ugrid(i,j)=.TRUE.
   hist(cnt,1)=i
   hist(cnt,2)=j
   lnPh=lnPP-umap(i,j)
   IF(lnPh<100) THEN
     umap(i,j)=umap(i,j)+LOG(1+EXP(lnPh))
   ELSE
     umap(i,j)=umap(i,j)+lnPh
   ENDIF
  END DO
!$OMP PARALLEL DO PRIVATE(i)
  DO i=1,cnt
   Ugrid(hist(i,1),hist(i,2))=.FALSE.
  ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE updateUSR
SUBROUTINE saveUSR(lnPa)
  INTEGER :: i,j, iost
  REAL(db) :: lnPa
   INTENT(IN) :: lnPa
  CHARACTER(13) :: names
  WRITE(names,'(A5,I6.6)') 'umap.',mpi_myproc
  WRITE(*,*) mpi_myproc,': writing map'
  OPEN(UNIT=scratch,FILE=TRIM('USR/'//names),FORM='FORMATTED',IOSTAT=iost,STATUS='REPLACE')
  IF(iost/=0) STOP 'Error when saving map file'
  WRITE(scratch,'(A,E16.6,I16)') fileEOS, rhos*EXP(xseed)/EDtoKM,step
  WRITE(scratch,'(3E16.6,2I6)') lnPa, usrdX,usrdY, NusrX,NusrY
  WRITE(scratch,'(4E16.6)') usrXL,usrXR,usrYD,usrYU
  DO i=1,NusrY
   DO j=1,NusrX
    WRITE(scratch,'(3E16.6)') xusra(i),yusra(j),umap(i,j)
   ENDDO
  END DO
  CLOSE(UNIT=scratch)
END SUBROUTINE saveUSR
SUBROUTINE loadUSR(lnPa)
  INTEGER :: i,j, iost
  REAL(db),INTENT(OUT) :: lnPa
  CHARACTER(13) :: names
  WRITE(names,'(A5,I6.6)') 'umap.',mpi_myproc
  WRITE(*,*) mpi_myproc,': reading USR'
  OPEN(UNIT=scratch,FILE=TRIM('USR/'//names),FORM='FORMATTED',IOSTAT=iost)
  IF(iost/=0) STOP 'Error when reading file from USR'
  READ(scratch,*) 
  READ(scratch,*) lnPa, usrdX,usrdY, NusrX,NusrY
  READ(scratch,*) usrXL,usrXR,usrYD,usrYU
  CALL init_user
  DO i=1,NusrY
   DO j=1,NusrX
    READ(scratch,*) xusra(i),yusra(j),umap(i,j)
   ENDDO
  END DO
  CLOSE(UNIT=scratch)
  IF(wflag) WRITE(*,*) ' lnPa, usrdX,usrdY, NusrX,NusrY:',lnPa, usrdX,usrdY, NusrX,NusrY
  IF(wflag) WRITE(*,*) ' usrXL,usrXR,usrYD,usrYU:', usrXL,usrXR,usrYD,usrYU
END SUBROUTINE loadUSR

END MODULE User