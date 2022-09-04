CC       =  gcc

CFLAGS   = -O3 -Wall -pg

CPPFLAGS = -MMD

SRC		 = field.c gauss_seidel.c print.c boundaries.c

OBJS     =	$(SRC:.c=.o)

DEP      =	$(SRC:.c=.d)

EXECS    =laplace

laplace: laplace.o	$(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(CPPFLAGS) -lm

-include	${DEP}

.PHONY: clean
clean:
	-rm *.o *.d$(EXECS)
