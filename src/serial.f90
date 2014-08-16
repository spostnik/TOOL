MODULE Parallel
  USE Params, ONLY : wflag, withMPI,db
  IMPLICIT NONE
  SAVE
  INTEGER :: mpi_nprocs=1,mpi_myproc=0,MPI_MAX,MPI_COMM_WORLD,MPI_IERROR,MPI_DOUBLE_PRECISION,MPI_IN_PLACE
CONTAINS
!***********************************************************************
  SUBROUTINE init_all_mpi
    wflag=.TRUE.
    withMPI=.FALSE.
    mpi_myproc=0
    mpi_nprocs=1
    IF(wflag) WRITE(*,'(A,i5)') 'No MPI, single process: ',mpi_nprocs
  END SUBROUTINE init_all_mpi
!***********************************************************************
  SUBROUTINE finish_mpi
    IF(wflag) WRITE(*,'(A)') 'Finishing serial run...'
  END SUBROUTINE finish_mpi
!***********************************************************************
  SUBROUTINE mpi_allreduce(MPI_IN_PLAC,lnPll,ii1,mpi_double_PRECISIO,MPI_MA,mpi_comm_worl,mpi_ierro)
   INTEGER, INTENT(IN) :: MPI_IN_PLAC,ii1,mpi_double_PRECISIO,MPI_MA,mpi_comm_worl,mpi_ierro
   REAL(db), INTENT(IN) :: lnPll
   IF(wflag) WRITE(*,'(A)') 'How come you are reading this?! Use no MPI option and version!'
   STOP
  END SUBROUTINE mpi_allreduce
END MODULE Parallel