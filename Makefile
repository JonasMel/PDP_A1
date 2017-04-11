CC         =  mpicc
CCFLAGS    =  -O3 -march=native
CCGFLAGS   =  -g
LIBS       =  -lmpi -lm

BINS= assign1

assign1: Assign1_1.c
	$(CC) $(CCFLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) $(BINS)



