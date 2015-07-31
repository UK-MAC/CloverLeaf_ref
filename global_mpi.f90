MODULE global_mpi_module

  USE definitions_module
  !USE MPI

  IMPLICIT NONE

  include "mpif.h"

CONTAINS

SUBROUTINE clover_sum(value)

  ! Only sums to the master

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  INTEGER :: err

  total=value

  CALL MPI_REDUCE(value,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_cart_comm,err)

  value=total

END SUBROUTINE clover_sum

SUBROUTINE clover_allsum(value)

  ! Global reduction for CG solver

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total, timer

  INTEGER :: err

  total=value

  CALL MPI_ALLREDUCE(value,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_cart_comm,err)

  value=total

END SUBROUTINE clover_allsum

SUBROUTINE clover_min(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: minimum

  INTEGER :: err

  minimum=value

  CALL MPI_ALLREDUCE(value,minimum,1,MPI_DOUBLE_PRECISION,MPI_MIN,mpi_cart_comm,err)

  value=minimum

END SUBROUTINE clover_min

SUBROUTINE clover_max(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: maximum

  INTEGER :: err

  maximum=value

  CALL MPI_ALLREDUCE(value,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,mpi_cart_comm,err)

  value=maximum

END SUBROUTINE clover_max

SUBROUTINE clover_allgather(value,values)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: values(parallel%max_task)

  INTEGER :: err

  values(1)=value ! Just to ensure it will work in serial

  CALL MPI_ALLGATHER(value,1,MPI_DOUBLE_PRECISION,values,1,MPI_DOUBLE_PRECISION,mpi_cart_comm,err)

END SUBROUTINE clover_allgather

SUBROUTINE clover_check_error(error)

  IMPLICIT NONE

  INTEGER :: error

  INTEGER :: maximum

  INTEGER :: err

  maximum=error

  CALL MPI_ALLREDUCE(error,maximum,1,MPI_INTEGER,MPI_MAX,mpi_cart_comm,err)

  error=maximum

END SUBROUTINE clover_check_error

SUBROUTINE clover_barrier

  INTEGER :: err

  CALL MPI_BARRIER(mpi_cart_comm,err)

END SUBROUTINE clover_barrier

SUBROUTINE clover_abort

  INTEGER :: ierr,err

  CALL MPI_ABORT(mpi_cart_comm,ierr,err)

END SUBROUTINE clover_abort

END MODULE global_mpi_module

