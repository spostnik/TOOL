MODULE Params
  IMPLICIT NONE
  INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
  INTEGER,PARAMETER :: qp=SELECTED_REAL_KIND(12,100)
!  INTEGER, PARAMETER :: qp=SELECTED_REAL_KIND(32)
!  INTEGER,PARAMETER :: qp=selected_real_kind(33, 4931)
  ! input and output units for files
  INTEGER, PARAMETER :: scratch=11, scratch2=12
  !**********************************************************************
  !     constants                                                       *
  !**********************************************************************
  REAL(db),PARAMETER :: pi=3.14159265358979D0
  REAL(db),PARAMETER :: co=299792.458D0, hbarc=197.3269718D0
  REAL(db),PARAMETER :: mnc2=939.565D0, mAc2=939.565D0
  REAL(db),PARAMETER :: Rkm=1.477D0, EDtoKM=1.3237D-6, gcm3=5.6D-13*EDtoKM, dyncm2=6.23D-34*EDtoKM
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  SAVE
  LOGICAL, ALLOCATABLE :: Ugrid(:,:)
  REAL(db),ALLOCATABLE,PUBLIC :: xusra(:), yusra(:), umap(:,:)
  REAL(db) :: usrXL,usrXR,usrYD,usrYU,usrdX,usrdY
  INTEGER :: stepMax=10, linkdots=10,step0=0
  INTEGER :: step=0, CRUSTknots=10, COREknots=10, linkpoints=10, iseed=0,NusrX,NusrY
  REAL(db) :: rhos=0, refden=0.08 !fm^-3
  REAL(db) :: lnPmin=-10000, xMax=2, xseed=-2, xEOSmx, yseed=0, yMax=2,nseed=0.001
  REAL(db) :: do=0, rss, lnPall0=-10000
  REAL(8) :: wall_time_now, wall_time_start, wall_time_old
  LOGICAL :: wflag=.FALSE., restart=.FALSE., NucPh=.FALSE., withMPI=.TRUE., &
	     Love=.FALSE.,MomI=.FALSE., MomQ=.FALSE., nDen=.FALSE., BE=.FALSE., &
	     addUSR=.FALSE.
  CHARACTER(LEN=10) :: wall_time
CONTAINS
  SUBROUTINE init_user
    INTEGER :: i,j
    IF(wflag) WRITE(*,'(A)') 'User map init. '
    ALLOCATE(Ugrid(NusrX,NusrY),xusra(NusrX),yusra(NusrY),umap(NusrX,NusrY))
!$OMP PARALLEL DO PRIVATE(i,j)
    DO i=1,NusrX
     DO j=1,NusrY
      Ugrid(i,j)=.FALSE.
     ENDDO
    END DO
!$OMP END PARALLEL DO
    IF(addUSR) RETURN

    usrdX=(usrXR-usrXL)/NusrX
    usrdY=(usrYU-usrYD)/NusrY
!$OMP PARALLEL DO PRIVATE(i,j)
    DO i=1,NusrX
	DO j=1,NusrY
	    umap(i,j)=lnPmin
	ENDDO
    END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,NusrX
    xusra(i)=usrXL+(i-0.5D0)*usrdX
   ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i)
   DO i=1,NusrY
    yusra(i)=usrYD+(i-0.5D0)*usrdY
   ENDDO
!$OMP END PARALLEL DO
  END SUBROUTINE init_user
  SUBROUTINE kill_user
    INTEGER :: i,j
    IF(wflag) WRITE(*,'(A)') 'User map kill '
    DEALLOCATE(Ugrid,xusra,yusra,umap)
  END SUBROUTINE kill_user
  !timers
  SUBROUTINE timing_start(message)
    CHARACTER(*), INTENT(IN) :: message
    CALL DATE_AND_TIME(TIME=wall_time)
    READ(wall_time,*) wall_time_now
    wall_time_old=wall_time_now
    WRITE(*,*) message,':', wall_time
  END SUBROUTINE timing_start
  SUBROUTINE timing_end(message)
    INTEGER :: t,h,m,t0,h0,m0
    REAL(8) :: s,s0
    CHARACTER(*), INTENT(IN) :: message
    CALL DATE_AND_TIME(TIME=wall_time)
    READ(wall_time,*) wall_time_now
    t=INT(wall_time_now)
    h=t/10000
    m=(t-10000*h)/100
    s=wall_time_now-h*10000-m*100
    t0=INT(wall_time_old)
    h0=t0/10000
    m0=(t0-10000*h0)/100
    s0=wall_time_old-h0*10000-m0*100
    IF(s0>s) THEN
     m=m-1
     s=s+60
     IF(m<0) THEN
      h=h-1
      m=60+m
     END IF
    END IF
    IF(m0>m) THEN
     h=h-1
     m=m+60
    END IF
    WRITE(*,*) message, ':', wall_time, ', took ', h-h0,'h ',m-m0,'m ',s-s0,'s '
  END SUBROUTINE timing_end
END MODULE Params