MODULE Markov
 USE Params, ONLY : db, do, wflag, timing_start, timing_end, &
		     restart, step, Rkm, stepMax, &
		     scratch,lnPmin,scratch2, withMPI, Love,lnPmin,MomI, &
		     NusrX,step0,lnPall0,lnPmin
 USE Solver
 USE Parallel
 USE Chain
 USE Analysis
 USE User, ONLY : updateUSR,saveUSR,loadUSR
 IMPLICIT NONE
 SAVE
 REAL, PRIVATE :: rnd,drnd
 REAL(db), PRIVATE :: xold,dold, lnPold, lnPall,yold
 INTEGER, PRIVATE :: knot,err,ok
 LOGICAL, PRIVATE :: tst
CONTAINS

 SUBROUTINE Markov_start
  INTEGER :: i,stp
  LOGICAL :: fff
  IF(.NOT.restart) THEN
   step=step0

   CALL timing_start(' start NS ')
   CALL calc_NS(2,err)
   CALL timing_end(' end NS ')
   WRITE(*,*) mpi_myproc,': err=',err
   
   CALL timing_start(' Bayesian start ')
   CALL Bayesian
   CALL timing_end(' Bayesian end ')
   WRITE(*,*) mpi_myproc,': lnPmax=',lnPmax

   CALL saveNS(lnPmax)
   fff=.TRUE.
   IF(NM>0) CALL updateMAP
   IF(Nx>0) CALL updateBAND
   IF(Love) CALL updateLOVE
   IF(MomI.AND.Ntrn>0) CALL updateMOM
   IF(NusrX>0) CALL updateUSR(lnPmax,NSchain,strmax)
   IF(wflag) CALL startLOG
   IF(wflag) CALL saveLOG
  ENDIF

  IF(step0==0) THEN
    lnPall=lnPmax
  ELSEIF(lnPall0/=0) THEN
    lnPall=lnPall0
  ELSE
    lnPall=lnPal
  ENDIF

  drnd=1.0/Nknots
  stp=step
  DO WHILE(step<=stepMax)
   step=step+1
   IF(wflag) WRITE(*,*) '**********Step=',step,'******************************'
   knot=-1
   DO WHILE(knot<MINknot)
    CALL RANDOM_NUMBER(rnd)
    knot=INT(rnd/drnd)+1
   ENDDO
   
   CALL modify_chain(ok)
   IF(ok==0) THEN
    lnPall=lnPmin
    IF(withMPI) CALL mpi_allreduce(MPI_IN_PLACE,lnPall,1,mpi_double_PRECISION,MPI_MAX,mpi_comm_world,mpi_ierror)
    CYCLE
   ENDIF

   CALL calc_NS(knot,err)
   IF(err<=0) WRITE(*,*) mpi_myproc,': step=',step,' err=',err
   
   lnPold=lnPmax
   CALL Bayesian
   WRITE(*,*) mpi_myproc,': lnPmax=',lnPmax
   IF(lnPmax>=lnPall.AND.lnPmax/=lnPmin) THEN
    WRITE(*,*) mpi_myproc,': step=',step,' maximum case found and saved lnPmax=',lnPmax
    CALL saveNS(lnPmax)
   ENDIF
   lnPall=MAX(lnPall,lnPmax)
!   IF(lnPall.NE.lnPall) lnPall=lnPmin
!   WRITE(*,*) mpi_myproc,': waiting others...'
    IF(withMPI) CALL mpi_allreduce(MPI_IN_PLACE,lnPall,1,mpi_double_PRECISION,MPI_MAX,mpi_comm_world,mpi_ierror)
!   WRITE(*,*) mpi_myproc,': passed others...'
   tst=.TRUE.
   IF(lnPmax<lnPold) CALL test_chain
   IF(wflag) WRITE(*,*) ' lnPall=',lnPall
   IF(tst) THEN 
    IF(lnPall-cut**2*0.5<lnPmax) THEN
     IF(writeNS/=0.AND.MOD(step,writeNS)==0) CALL saveNS(lnPmax)
     IF(Nx>0) CALL updateBAND
     IF(NM>0) CALL updateMAP
     IF(Love) CALL updateLOVE
     IF(MomI.AND.Ntrn>0) CALL updateMOM
     IF(NusrX>0) CALL updateUSR(lnPmax,NSchain,strmax)
     IF(stp+writeevery<step) THEN 
      stp=step
      IF(Nx>0) CALL saveBAND(lnPall)
      IF(NM>0) CALL saveMAP(lnPall)
      IF(Love) CALL saveLOVE(lnPall)
      IF(MomI.AND.Ntrn>0) CALL saveMOM(lnPall)
      IF(NusrX>0) CALL saveUSR(lnPall)
      fff=.TRUE.
     ENDIF
    ENDIF
    IF(wflag) CALL saveLOG
   ENDIF
  END DO
  fff=.TRUE.
  IF(Nx>0) CALL saveBAND(lnPall)
  IF(NM>0) CALL saveMAP(lnPall)
  IF(Love) CALL saveLOVE(lnPall)
  IF(MomI.AND.Ntrn>0) CALL saveMOM(lnPall)
  IF(NusrX>0) CALL saveUSR(lnPall)
  CALL saveNS(lnPmax)
 END SUBROUTINE Markov_start

 SUBROUTINE modify_chain(ok)
  INTEGER, INTENT(OUT) :: ok
  INTEGER :: i,j,k
  REAL(db) :: ynew, xnew, dd,xx, aL, aR,yy,yp,yn,xp,xn
  LOGICAL :: chk
  xold=xchain(knot)
  yold=ychain(knot)
!  dold=yold-LOG(do+EXP(xold))
  xx=xold
  xnew=xold
!  WRITE(*,*) mpi_myproc,': start jump knot=',knot
  j=0
  IF(knot>1) THEN
   xp=xchain(knot-1)
   yp=ychain(knot-1)
  ENDIF
  IF(knot==Nknots) THEN
!    WRITE(*,*) mpi_myproc,': start jump 1'
!$OMP PARALLEL DO PRIVATE(i,k,dd,aL,rnd,yy,chk)
   DO i=1,320
    IF(j==0) k=0
    DO WHILE(k<100)
     CALL RANDOM_NUMBER(rnd)
     yy=yold+djump*(1-2*rnd)
     dd=yy-LOG(do+EXP(xx))
     aL=(yy-yp)/(xx-xp)
     IF(aL/=1.AND.aL>0) THEN
      aL=(LOG(aL)+yp)/(1-aL)
     ELSEIF(aL==1) THEN
      aL=xx-1
      IF(yp>=0) aL=xx+1
     ELSE
      aL=xx-1
     ENDIF
     chk=(dd>=0.AND.xx<=aL)
     IF(j==1.OR.chk) EXIT
     k=k+1
    ENDDO
    IF(j==0.AND.chk) THEN
!$OMP ATOMIC
     j=j+1
!$OMP CRITICAL
     ynew=yy
!$OMP END CRITICAL
     k=100
    ENDIF
   ENDDO
!$OMP END PARALLEL DO
  ELSEIF(knot==1) THEN
   xn=xchain(knot+1)
   yn=ychain(knot+1)
!$OMP PARALLEL DO PRIVATE(i,k,dd,aR,rnd,yy,chk)
 DO i=1,320
   IF(j==0) k=0
   DO WHILE(k<100)
    CALL RANDOM_NUMBER(rnd)
    yy=yold+djump*(1-2*rnd)
    dd=yy-LOG(do+EXP(xx))
    aR=(yy-yn)/(xx-xn)
    IF(aR/=1.AND.aR>0) THEN
     aR=(LOG(aR)+yy)/(1-aR)
    ELSEIF(aR==1) THEN
     aR=xn-1
     IF(yy>=0) aR=xn+1
    ELSE
     aR=xn-1
    ENDIF
    chk=(dd>=0.AND.xn<=aR)
    IF(j==1.OR.chk) EXIT
    k=k+1
   ENDDO
   IF(j==0.AND.chk) THEN
!$OMP ATOMIC
    j=j+1
!$OMP CRITICAL
    ynew=yy
!$OMP END CRITICAL
    k=100
   ENDIF
 ENDDO
!$OMP END PARALLEL DO
  ELSEIF(knot==No) THEN
!   WRITE(*,*) mpi_myproc,': start jump 2'
   xn=xchain(knot+1)
   yn=ychain(knot+1)
!$OMP PARALLEL DO PRIVATE(chk,i,k,dd,aL,aR,rnd,yy)
 DO i=1,320
   IF(j==0) k=0
   DO WHILE(k<100)
    CALL RANDOM_NUMBER(rnd)
    yy=yp+(yn-yp)*rnd
    dd=yy-LOG(do+EXP(xx))
    aL=(yy-yp)/(xx-xp)
    IF(aL/=1.AND.aL>0) THEN
     aL=(LOG(aL)+yp)/(1-aL)
    ELSEIF(aL==1) THEN
     aL=xx-1
     IF(yp>=0) aL=xx+1
    ELSE
     aL=xx-1
    ENDIF
    aR=(yy-yn)/(xx-xn)
    IF(aR/=1.AND.aR>0) THEN
     aR=(LOG(aR)+yy)/(1-aR)
    ELSEIF(aR==1) THEN
     aR=xn-1
     IF(yy>=0) aR=xn+1
    ELSE
     aR=xn-1
    ENDIF
    chk=(dd>=0.AND.xx<=aL.AND.xn<=aR)
    IF(j==1.OR.chk) EXIT
    k=k+1
   ENDDO
   IF(j==0.AND.chk) THEN
!$OMP ATOMIC
    j=j+1
!$OMP CRITICAL
    ynew=yy
!$OMP END CRITICAL
    k=100
   ENDIF
 ENDDO
!$OMP END PARALLEL DO
  ELSE
!   WRITE(*,*) mpi_myproc,': start jump 3'
   xn=xchain(knot+1)
   yn=ychain(knot+1)
!$OMP PARALLEL DO PRIVATE(chk,i,k,dd,xx,aL,aR,rnd,yy)
  DO i=1,320
   IF(j==0) k=0
   DO WHILE(k<100)
    CALL RANDOM_NUMBER(rnd)
    yy=yp+(yn-yp)*rnd
    CALL RANDOM_NUMBER(rnd)
    xx=xp+(xn-xp)*rnd
    dd=yy-LOG(do+EXP(xx))
    aL=(yy-yp)/(xx-xp)
    IF(aL/=1.AND.aL>0) THEN
     aL=(LOG(aL)+yp)/(1-aL)
    ELSEIF(aL==1) THEN
     aL=xx-1
     IF(yp>=0) aL=xx+1
    ELSE
     aL=xx-1
    ENDIF
    aR=(yy-yn)/(xx-xn)
    IF(aR/=1.AND.aR>0) THEN
     aR=(LOG(aR)+yy)/(1-aR)
    ELSEIF(aR==1) THEN
     aR=xn-1
     IF(yy>=0) aR=xn+1
    ELSE
     aR=xn-1
    ENDIF
    chk=(dd>0.AND.xx>xchain(knot-1).AND.xx<xchain(knot+1).AND.xx<=aL.AND.xn<=aR)
    IF(j>0.OR.chk) EXIT
    k=k+1
   ENDDO
   IF(j==0.AND.chk) THEN
!$OMP ATOMIC
    j=j+1
!$OMP CRITICAL
    ynew=yy
    xnew=xx
!$OMP END CRITICAL
    k=100
   ENDIF
 ENDDO
!$OMP END PARALLEL DO
  ENDIF
!  WRITE(*,*) mpi_myproc,': end jump'
  IF(j==0) THEN
   ok=0
   WRITE(*,*) mpi_myproc,': cannot move knot=',knot
  ELSE
   ok=1
   xchain(knot)=xnew
   ychain(knot)=ynew
   CALL fix_links
  ENDIF
 END SUBROUTINE modify_chain
 
 SUBROUTINE fix_links
   INTEGER :: i,j
   IF(knot==1) NSchain(1,4)=ychain(1)
!$OMP PARALLEL DO PRIVATE(i,j)
   DO i=1,2*linkpoints
    j=MAX(1,1+linkpoints*(knot-2)+i)
    IF(knot<Nknots.AND.i>linkpoints) THEN
     NSchain(j,3)=xchain(knot)+(i-linkpoints)*(xchain(knot+1)-xchain(knot))/linkpoints
     NSchain(j,4)=ychain(knot)+(i-linkpoints)*(ychain(knot+1)-ychain(knot))/linkpoints
     fchain(j)=1/(1+EXP(NSchain(j,4)-NSchain(j,3)))
    ENDIF
    IF(i<=linkpoints.AND.knot>1) THEN
     NSchain(j,3)=xchain(knot-1)+i*(xchain(knot)-xchain(knot-1))/linkpoints
     NSchain(j,4)=ychain(knot-1)+i*(ychain(knot)-ychain(knot-1))/linkpoints
     fchain(j)=1/(1+EXP(NSchain(j,4)-NSchain(j,3)))
    ENDIF
   ENDDO
!$OMP END PARALLEL DO
 END SUBROUTINE fix_links

 SUBROUTINE test_chain
  CALL RANDOM_NUMBER(rnd)
  IF(rnd<EXP(lnPmax-lnPold).AND.(lnPall-cut**2*0.5<lnPmax)) RETURN
  tst=.FALSE.
  WRITE(*,*) mpi_myproc,': rejected'
  xchain(knot)=xold
  ychain(knot)=yold
  CALL fix_links
  lnPmax=lnPold
 END SUBROUTINE test_chain

 SUBROUTINE startLOG
    LOGICAL :: exists
    INQUIRE(FILE='log.dat',EXIST=exists)
    IF(exists.AND.step>0) RETURN
    OPEN(UNIT=scratch,FILE='log.dat',FORM='FORMATTED',STATUS='REPLACE')
     WRITE(scratch,'(A10,3A16,A10,A16)') 'Step','0:Max. ln(P)','Max. M','ok','knot','Tot. Max.'
    CLOSE(UNIT=scratch)
 END SUBROUTINE startLOG

 SUBROUTINE saveLOG
    OPEN(unit=scratch2,file='log.dat',POSITION='APPEND')
     WRITE(scratch2,'(I10,2F16.6,I16,I10,F16.6)') step,lnPmax,mMaxAll,ok,knot,lnPall
    CLOSE(scratch2)
 END SUBROUTINE saveLOG

END MODULE Markov