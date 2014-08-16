PROGRAM TOOL
  USE Params
  USE Parallel
  USE Solver
  USE Markov
  USE Data
  USE Chain
  USE User
  IMPLICIT NONE
  INTEGER :: i,rndseed(1:20),rndsize,j,jj
  REAL(db) :: mseed, r2seed, lnPm,lnPb, lnPn
  CHARACTER(13) :: names
  LOGICAL :: exists
  NAMELIST /init/ refden, xseed, xMax, lnPmin, linkdots, withMPI, &
		    Love, MomI, MomQ, nDen, BE
  NAMELIST /Mchain/ CRUSTknots,COREknots,linkpoints, restart, step, cut, djump, step0, stepMax, MINknot, lnPall0
  NAMELIST /maps/ Nx,Ny,NM,NR,xL,xR,yD,yU,ML,MR,RD,RU,NLove,LoveL,LoveR,Ntrn,trnL,trnR,Nfrc,frcD,frcU, &
		    NusrX,NusrY,usrXL,usrXR,usrYU,usrYD
  NAMELIST /input/ fileEOS, dotsEOS, units, addBAND,addMAP,addLOVE,addMOM, addUSR
  NAMELIST /output/ writeevery, writeNS
  NAMELIST /Odata/  MmaxNS, &
		    M1maxNS,sM1maxNS,M2maxNS,sM2maxNS,M3maxNS,sM3maxNS, &
		    RsNS,sRsNS,MmnNS, &
		    M1,sM1,R1,sR1,&
		    R1inf,sR1inf,R2inf,sR2inf,z2,sz2, &
		    LoveLimit, Mlimit, Rlimit, &
		    glitchAmax, glitchRmax, glitchPmx
  NAMELIST /Ndata/ NucPh,fileEDF
  
  
  CALL init_all_mpi
  
  WRITE(*,*) mpi_myproc,': reading knobs'

   OPEN(unit=11,file='knobs',status='old',form='formatted')
    READ(11,init,END=11)
    READ(11,Mchain,END=11)
    READ(11,maps,END=11)
    READ(11,input,END=11)
    READ(11,output,END=11)
    READ(11,Odata,END=11)
    READ(11,Ndata,END=11)
11 CLOSE(11)

  IF(wflag.AND.MomQ) STOP 'Quadrupole moment is not suported yet. Set MomQ=F.'

  IF(MomI.AND.xseed==0.AND.Ntrn/=0) STOP 'No bare SQM star without crust, failed to get crust moment of inertia'

  IF(BE) nDen=.TRUE.
  
  IF(wflag.AND.Love) WRITE(*,*) 'Love number lambda will be calculated, output units km^5'
  IF(wflag.AND.MomI) WRITE(*,*) 'Moment of inertia for slow spin will be calculated, output units M_Sun*km^2'
  IF(wflag.AND.BE) WRITE(*,*) 	'Binding Energy will be calculated, output units M_Sun*c^2'

  IF(restart) THEN
   IF(wflag) WRITE(*,*) 'Restarting from step=',step
   addMAP=.TRUE.
   addBAND=.TRUE.
   addLOVE=.TRUE.
   addMOM=.TRUE.
   addUSR=.TRUE.
  ENDIF

  WRITE(names,'(A5,I6.6)') 'pmap.',mpi_myproc
  INQUIRE(FILE=TRIM('MAP/'//names),EXIST=exists)
  IF(.NOT.exists.AND.addMAP) addMAP=.FALSE.
  WRITE(names,'(A5,I6.6)') 'band.',mpi_myproc
  INQUIRE(FILE=TRIM('BND/'//names),EXIST=exists)
  IF(.NOT.exists.AND.addBAND) addBAND=.FALSE.
  WRITE(names,'(A5,I6.6)') 'lmap.',mpi_myproc
  INQUIRE(FILE=TRIM('LOV/'//names),EXIST=exists)
  IF(.NOT.exists.AND.addLOVE) addLOVE=.FALSE.
  WRITE(names,'(A5,I6.6)') 'ivst.',mpi_myproc
  INQUIRE(FILE=TRIM('MOM/'//names),EXIST=exists)
  IF(.NOT.exists.AND.addMOM) addMOM=.FALSE.
  WRITE(names,'(A5,I6.6)') 'umap.',mpi_myproc
  INQUIRE(FILE=TRIM('USR/'//names),EXIST=exists)
  IF(.NOT.exists.AND.addUSR) addUSR=.FALSE.

  IF(wflag) WRITE(*,*) 'Chain: crust knots = ',CRUSTknots,', core knots = ',COREknots,', points per link =',linkpoints
  IF(wflag) WRITE(*,*) 'Movable knots >= ',MINknot

  CALL alloc_chain

  rhos=EDtoKM*refden*(2.0D0/5)*(hbarc)**2*(3*pi**2*refden)**(2.0D0/3)/(2*mnc2)
  rss=rhos*Rkm*Rkm
  IF(wflag) WRITE(*,*) 'Reference energy density from neutron Fermi gas pressure ',rhos,' km^-2 (',rhos/EDtoKM,' MeV/fm^3)'
  
  IF (dotsEOS>0) THEN
   CALL loadEOS
  ELSE
   STOP 'no seed EOS'
  ENDIF

  IF (NucPh) THEN
   CALL loadEDF
  ELSE
   IF(wflag) WRITE(*,*) 'no EDF is used'
  ENDIF

  CALL init_MR(xseed,mseed,r2seed)
  IF(wflag) WRITE(*,*) mpi_myproc,': Seed EOS for n<=',nseed,' ln(rhoseed/rhos)<=',xseed,' yseed=',yseed,' iseed=',iseed
  IF(wflag) WRITE(*,*) mpi_myproc,': Seed star M=',mseed,'Msol R=',SQRT(r2seed),'km ',' I=',NSchain(1,6),' Msol*km^2'

  IF(wflag) WRITE(*,*) 'Init. random seed'
  call random_seed(size=rndsize)
  do i=1,rndsize
   call system_clock(count=rndseed(i))
   rndseed(i)=rndseed(i)+mpi_myproc
  enddo
  call random_seed(put=rndseed(1:rndsize))
  
  CALL init_chain(lnPm,lnPb)
  lnPn=MAX(lnPm,lnPb)
  WRITE(*,*) mpi_myproc,': lnPn=',lnPn

  IF(NM>0.OR.Nfrc>0) CALL init_grid
  IF(Nx>0) CALL init_gridEOS

  WRITE(*,*) mpi_myproc,': xchain=',xchain
  WRITE(*,*) mpi_myproc,': ychain=',ychain

  CALL Markov_start

  CALL kill_data
  CALL kill_chain
  CALL kill_grids

  IF(MomI) CALL kill_profile
  IF(NusrX>0) CALL kill_user
  
  CALL finish_mpi

END PROGRAM TOOL