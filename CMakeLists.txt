cmake_minimum_required(VERSION 3.10)

project( 
       CLOVERLEAF
       VERSION 1.3.1
       LANGUAGES C Fortran)

FIND_PACKAGE(OpenMP COMPONENTS C Fortran)

IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  list(APPEND CMAKE_EXE_LINKER_FLAGS "${OpenMP_EXE_LINKER_FLAGS}")
ELSE()
  MESSAGE(WARNING "OpenMP not found")
  SET(ENABLE_OPENMP off)
ENDIF()

set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)

set(CLOVERLEAF_SRC
    ${PROJECT_SOURCE_DIR}/src/parse.f90
    ${PROJECT_SOURCE_DIR}/src/update_halo_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/hydro.f90
    ${PROJECT_SOURCE_DIR}/src/ideal_gas.f90
    ${PROJECT_SOURCE_DIR}/src/ideal_gas_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/initialise_chunk.f90
    ${PROJECT_SOURCE_DIR}/src/initialise_chunk_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/initialise.f90
    ${PROJECT_SOURCE_DIR}/src/pack_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/read_input.f90
    ${PROJECT_SOURCE_DIR}/src/report.f90
    ${PROJECT_SOURCE_DIR}/src/reset_field.f90
    ${PROJECT_SOURCE_DIR}/src/reset_field_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/revert.f90
    ${PROJECT_SOURCE_DIR}/src/revert_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/start.f90
    ${PROJECT_SOURCE_DIR}/src/timer.f90
    ${PROJECT_SOURCE_DIR}/src/timestep.f90
    ${PROJECT_SOURCE_DIR}/src/update_halo.f90
    ${PROJECT_SOURCE_DIR}/src/update_tile_halo.f90
    ${PROJECT_SOURCE_DIR}/src/update_tile_halo_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/viscosity.f90
    ${PROJECT_SOURCE_DIR}/src/viscosity_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/visit.f90
    ${PROJECT_SOURCE_DIR}/src/calc_dt_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/clover.f90
    ${PROJECT_SOURCE_DIR}/src/clover_leaf.f90
    ${PROJECT_SOURCE_DIR}/src/data.f90
    ${PROJECT_SOURCE_DIR}/src/definitions.f90
    ${PROJECT_SOURCE_DIR}/src/field_summary.f90
    ${PROJECT_SOURCE_DIR}/src/field_summary_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/flux_calc.f90
    ${PROJECT_SOURCE_DIR}/src/flux_calc_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/generate_chunk.f90
    ${PROJECT_SOURCE_DIR}/src/generate_chunk_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/accelerate.f90
    ${PROJECT_SOURCE_DIR}/src/accelerate_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/advec_cell_driver.f90
    ${PROJECT_SOURCE_DIR}/src/advec_cell_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/advec_mom_driver.f90
    ${PROJECT_SOURCE_DIR}/src/advec_mom_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/advection.f90
    ${PROJECT_SOURCE_DIR}/src/build_field.f90
    ${PROJECT_SOURCE_DIR}/src/calc_dt.f90
    ${PROJECT_SOURCE_DIR}/src/PdV_kernel.f90
    ${PROJECT_SOURCE_DIR}/src/PdV.f90
    ${PROJECT_SOURCE_DIR}/src/ideal_gas_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/initialise_chunk_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/pack_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/reset_field_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/revert_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/timer_c.c
    ${PROJECT_SOURCE_DIR}/src/update_halo_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/viscosity_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/calc_dt_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/field_summary_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/flux_calc_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/generate_chunk_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/accelerate_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/advec_cell_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/advec_mom_kernel_c.c
    ${PROJECT_SOURCE_DIR}/src/PdV_kernel_c.c)

add_executable(clover_leaf ${CLOVERLEAF_SRC})

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   if(CMAKE_BUILD_TYPE MATCHES "Release")
      set( CMAKE_C_FLAGS_RELEASE        "-O3 -march=native -funroll-loops" )
      set( CMAKE_Fortran_FLAGS_RELEASE  "-O3 -march=native -funroll-loops -cpp -ffree-line-length-none" )
   elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
      set( CMAKE_C_FLAGS_DEBUG          "-O0 -g -O -Wall -Wextra -fbounds-check")
      set( CMAKE_Fortran_FLAGS_DEBUG    "${CMAKE_C_FLAGS_DEBUG}")
   endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   if(CMAKE_BUILD_TYPE MATCHES "Release")
      set( CMAKE_C_FLAGS_RELEASE        "-O3 -no-prec-div -restrict -fno-alias" )
      set( CMAKE_Fortran_FLAGS_RELEASE  "-O3 -no-prec-div" )
   elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
      set( CMAKE_C_FLAGS_DEBUG          "-O0 -g -debug all -traceback")
      set( CMAKE_Fortran_FLAGS_DEBUG    "-O0 -g -debug all -check all -traceback -check noarg_temp_created")
   endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang") 
   if(CMAKE_BUILD_TYPE MATCHES "Release")
      set( CMAKE_C_FLAGS_RELEASE        "-O3 -mcpu=native" )
      set( CMAKE_Fortran_FLAGS_RELEASE  "-cpp -O3 -ffp-contract=fast -mcpu=native" )
   elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
      set( CMAKE_C_FLAGS_DEBUG          "-O0 -g")
      set( CMAKE_Fortran_FLAGS_DEBUG    "-cpp -O0 -g")
   endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
   if(CMAKE_BUILD_TYPE MATCHES "Release")
      set( CMAKE_C_FLAGS_RELEASE        "-em -h list=a" )
      set( CMAKE_Fortran_FLAGS_RELEASE  "-em -ra -h acc_model=fast_addr:no_deep_copy:auto_async_all" )
   elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
      set( CMAKE_C_FLAGS_DEBUG          "-O0 -g -em -eD")
      set( CMAKE_Fortran_FLAGS_DEBUG    "${CMAKE_C_FLAGS_DEBUG}")
   endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
   if(CMAKE_BUILD_TYPE MATCHES "Release")
      set( CMAKE_C_FLAGS_RELEASE        "-fastsse -Mipa=fast -Mlist" )
      set( CMAKE_Fortran_FLAGS_RELEASE  "${CMAKE_C_FLAGS_RELEASE}" )
   elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
      set( CMAKE_C_FLAGS_DEBUG          "-O0 -g -C -Mchkstk -Ktrap=fp -Mchkfpstk -Mchkptr")
      set( CMAKE_Fortran_FLAGS_DEBUG    "-O0 -g -C -Mchkstk -Ktrap=fp -Mchkfpstk")
   endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "XL")
   if(CMAKE_BUILD_TYPE MATCHES "Release")
      set( CMAKE_C_FLAGS_RELEASE        "-O3 -qipa=partition=large -g -qfullpath -Q -qlistopt -qattr=full -qlist -qreport -qxref=full -qsource -qsuppress=1506-224:1500-036 -qsrcmsg" )
      set( CMAKE_Fortran_FLAGS_RELEASE  "-O3 -qipa=partition=large -g -qfullpath -Q -qsigtrap -qlistopt -qattr=full -qlist -qreport -qxref=full -qsource -qsuppress=1506-224:1500-036" )
   elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
      set( CMAKE_C_FLAGS_DEBUG          "-O0 -g -qfullpath -qcheck -qflttrap=ov:zero:invalid:en -qsource -qinitauto=FF -qmaxmem=-1 -qsrcmsg")
      set( CMAKE_Fortran_FLAGS_DEBUG    "-O0 -g -qfullpath -qcheck -qflttrap=ov:zero:invalid:en -qsource -qinitauto=FF -qmaxmem=-1 -qinit=f90ptr -qsigtrap")
   endif()
endif()

#set(CMAKE_C_FLAGS         "${CMAKE_C_FLAGS}"         CACHE STRING "compiler flags" FORCE)
set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}"   CACHE STRING "compiler flags" FORCE)
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "compiler flags" FORCE)

#set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS}"         CACHE STRING "compiler flags" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG}"   CACHE STRING "compiler flags" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" CACHE STRING "compiler flags" FORCE)