    F90 := ifort
    FC  := ifort
    CC  := icc

    ifdef MPI
      MPIF90 := mpiifort
      F90    := mpiifort
      FC     := mpiifort
      CC     := mpiicc
    endif

    #FFLAGS   = -module $(mdir) -I $(mdir)
    #F90FLAGS = -module $(mdir) -I $(mdir)
    CFLAGS   = -std=c99

    ifdef OMP
      FFLAGS   += -openmp -openmp-report2 -g -debug inline-debug-info -parallel-source-info=2 -align array64byte -opt-assume-safe-padding -opt-streaming-cache-evict=0
      F90FLAGS += -openmp -openmp-report2 -g -debug inline-debug-info -parallel-source-info=2  -align array64byte -opt-assume-safe-padding -opt-streaming-cache-evict=0
      CFLAGS   += -openmp -openmp-report2 -g -debug inline-debug-info -parallel-source-info=2  -opt-assume-safe-padding -opt-streaming-cache-evict=0
      #FFLAGS   +=  -qopenmp -qopt-report=5 -qopt-report-file=myopts.txt -g -debug inline-debug-info -parallel-source-info=2 
      #F90FLAGS +=   -qopenmp -qopt-report=5 -qopt-report-file=myopts.txt -g -debug inline-debug-info -parallel-source-info=2 
      #CFLAGS   +=  -qopenmp -qopt-report=5 -qopt-report-file=myopts.txt -g -debug inline-debug-info -parallel-source-info=2 
    endif

    ifdef MIC
      FFLAGS   += -mmic -no-prec-sqrt -no-prec-div -fno-alias -fimf-precision=low -fimf-domain-exclusion=15
      F90FLAGS += -mmic -no-prec-sqrt -no-prec-div -fno-alias -fimf-precision=low -fimf-domain-exclusion=15
      CFLAGS   += -mmic -no-prec-sqrt -no-prec-div -fno-alias
    endif

    ifdef NDEBUG
      F90FLAGS += -O3
      FFLAGS   += -O3
      CFLAGS   += -O3
    else
      F90FLAGS += -g -traceback -O0 #-check all -warn all -u 
      FFLAGS   += -g -traceback -O0 #-check all -warn all -u 
      CFLAGS   += -g -Wcheck
    endif

    ifdef GPROF
      FFLAGS   += -p
      F90FLAGS += -p
      CFLAGS   += -p
    endif
