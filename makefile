include make-local.inc

CFLAGS     = -Wall \
             -L${PETSC_DIR}/lib/ \
             -L${PETSC_DIR}/${PETSC_ARCH}/lib/ -lpetsc \
             -I${PETSC_DIR}/include \
             -I${PETSC_DIR}/${PETSC_ARCH}/include \
             -I${MPI_INCLUDE}

# PETSC_LIB=-L${PETSC_DIR}/lib -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include \
#   -I${MPI_INCLUDE} #-lpetsc 
#             -lpetsc \
#             -I/usr/local/pkgs-modules/openmpi_1.6.5_intel_2013/include

FILE_MAT   = "system_data_binary/two_hole_nocrack_mat.bin"
FILE_VEC   = "system_data_binary/two_hole_nocrack_vec.bin"
FILE_SOL   = "system_data_binary/two_hole_nocrack_sol.bin"# LU solution
FILE_REF   = "system_data_binary/two_hole_nocrack_ref.bin"# LU +5cg iters
KSP_OPTS   =  -ksp_converged_reason -ksp_norm_type preconditioned -pc_type jacobi

# FFLAGS     = -I${PETSC_DIR}/include/petsc/finclude
# SOURCESC   = solve-fem-system.c
# SOURCESF   =
# OBJ        = $(SOURCESC:.c=.o)
# CLEANFILES = ${OBJ} test

HOSTNAME  := $(shell hostname)

all: exe solve-lu solve-ref time-all

test: exe test-solve test-rtol test-time

exe: ksp-solve.c ksp-solve.h
	# module load petsc_3.6.3
	$(CC) $(CFLAGS) ksp-solve.c -o ksp-solve
	mkdir -p $(HOSTNAME)

test-solve: ksp-solve
	./ksp-solve $(KSP_OPTS) -v 3 \
  -rmat $(FILE_MAT) -rrhs $(FILE_VEC) -rref $(FILE_SOL) \
  -ksp_type cg \
  -ksp_max_it 100 \
  -log_view :$(HOSTNAME)/test-100-cg.log \
  -wcsv $(HOSTNAME)/test-100-cg.log

test-rtol: ksp-solve
	for RTOL in "1e-5" "1e-4" "1e-3" ; do \
   for KSP in \
   cg cr ; \
   do rm -f $(HOSTNAME)/test-ksp-$$KSP-r$$RTOL.out ; \
   ./ksp-solve $(KSP_OPTS) \
  -rmat $(FILE_MAT) -rrhs $(FILE_VEC) -rref $(FILE_SOL) \
  -v 3 \
  -ksp_type $$KSP  -ksp_rtol $$RTOL \
  -ksp_max_it 100 \
  -log_view :$(HOSTNAME)/test-ksp-$$KSP-r$$RTOL.log \
  -wcsv $(HOSTNAME)/test-ksp-$$KSP-r$$RTOL.log ; \
   done ; done ;
	./log-parse.sh test

RTOL=1e-3
KSP=cr
test-time: ksp-solve
	./ksp-solve $(KSP_OPTS) -v 3 \
  -rmat $(FILE_MAT) -rrhs $(FILE_VEC) -rref $(FILE_SOL) \
  -ksp_type $(KSP) -ksp_rtol $(RTOL) \
  -log_view :$(HOSTNAME)/test-$(KSP)-r$(RTOL).log \
  -wcsv $(HOSTNAME)/test-ksp-$(KSP)-r$(RTOL).log

time-all: ksp-solve
	for RTOL in "1e-3" "1e-5" "1e-4" ; do \
   for KSP in \
   symmlq cr    cgs    lcd    cg     fcg \
   bicg   bcgs  ibcgs  bcgsl             \
   minres gmres lgmres dgmres tfqmr ;    \
   do rm -f $(HOSTNAME)/test-ksp-$$KSP-r$$RTOL.out ; \
   ./ksp-solve $(KSP_OPTS) \
  -rmat $(FILE_MAT) -rrhs $(FILE_VEC) \
  -v 1 \
  -ksp_type $$KSP -ksp_rtol $$RTOL \
  -log_view :$(HOSTNAME)/time-ksp-$$KSP-r$$RTOL.log \
  -wcsv $(HOSTNAME)/time-ksp-$$KSP-r$$RTOL.log ; \
   done ; done ;
	./log-parse.sh time

solve-lu: ksp-solve
	./ksp-solve $(KSP_OPTS) -v 3 \
  -rmat $(FILE_MAT) -rrhs $(FILE_VEC) -wsol $(FILE_SOL) \
  -ksp_type preonly -pc_type lu \
  -log_view :$(HOSTNAME)/solve-ksp-lu.log

solve-ref: ksp-solve
	./ksp-solve $(KSP_OPTS) -v 3 \
  -rmat $(FILE_MAT) -rrhs $(FILE_VEC) -ruin $(FILE_SOL) -wsol $(FILE_REF) \
  -ksp_type cg \
  -ksp_max_it 5 -ksp_rtol 1e-12 \
  -log_view :$(HOSTNAME)/solve-ksp-ref.log
# CG iterations don't improve the LU solution at all
