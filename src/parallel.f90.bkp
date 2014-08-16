MODULE Parallel
  USE Params, ONLY : wflag, withMPI
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  INTEGER :: mpi_nprocs=1,mpi_ierror,mpi_myproc=0,stts(MPI_STATUS_SIZE)
CONTAINS
!***********************************************************************
  SUBROUTINE init_all_mpi
    CALL mpi_init(mpi_ierror)
    CALL mpi_comm_size(mpi_comm_world,mpi_nprocs,mpi_ierror)
    CALL mpi_comm_rank(mpi_comm_world,mpi_myproc,mpi_ierror)
    IF(mpi_ierror/=0) STOP 'MPI init. error!'
    wflag=(mpi_myproc==0)
    IF(wflag) WRITE(*,'(A,i5)') 'Number of MPI proc. ',mpi_nprocs
    IF(mpi_nprocs==1) withMPI=.FALSE.
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  END SUBROUTINE init_all_mpi
!***********************************************************************
  SUBROUTINE finish_mpi
    INTEGER :: ierr
    CALL mpi_finalize(ierr)
  END SUBROUTINE finish_mpi
!***********************************************************************
END MODULE Parallel