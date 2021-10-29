include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

TopOpt: TopOpt.o
	-${CLINKER} -o TopOpt TopOpt.o ${PETSC_LIB}
	${RM} TopOpt.o
hex888: hex888.o
	-${CLINKER} -o hex888 hex888.o ${PETSC_LIB}
	${RM} hex888.o
