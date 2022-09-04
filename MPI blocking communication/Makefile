CC       =  gcc

CFLAGS   = -O3 -Wall -pg

OBJS = boundaries.o field.o gauss_seidel.o print.o

sor_prof: laplace.c ${OBJS} 
	${CC} ${CFLAGS} -o sor_prof laplace.c ${OBJS} -lm


.PHONY: clean
clean:
	-rm *.o *.d$(EXECS)
