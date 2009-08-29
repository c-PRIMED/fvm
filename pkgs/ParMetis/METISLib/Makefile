include ../Makefile.in

.SUFFIXES: .c .lo .o

all: ../libmetis.a ../libmetis.so

CFLAGS = $(COPTIONS) $(OPTFLAGS) -I. $(INCDIR)


OBJS = coarsen.o fm.o initpart.o match.o ccgraph.o \
       pmetis.o pqueue.o refine.o util.o timing.o debug.o \
       bucketsort.o graph.o stat.o kmetis.o kwayrefine.o \
       kwayfm.o balance.o ometis.o srefine.o sfm.o separator.o \
       mincover.o mmd.o mesh.o meshpart.o frename.o fortran.o \
       myqsort.o compress.o parmetis.o estmem.o \
       mpmetis.o mcoarsen.o mmatch.o minitpart.o mbalance.o \
       mutil.o mkmetis.o mkwayrefine.o mkwayfmh.o \
       mrefine2.o minitpart2.o mbalance2.o mfm2.o \
       kvmetis.o kwayvolrefine.o kwayvolfm.o subdomains.o \
       mfm.o memory.o mrefine.o checkgraph.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c

SHLOBJS = $(OBJS:.o=.lo)

.c.lo:
	$(CC) $(CFLAGS) -fPIC -DPIC -c $< -o $@

../libmetis.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

../libmetis.so: $(SHLOBJS)
	$(CC) -shared $^ -Wl,-soname,libmetis.so.3.1.1 -o $@.3.1.1 
	ln -s libmetis.so.3.1.1 $@

clean:
	rm -f *.o *.lo

realclean:
	rm -f *.o *.lo; rm -f ../libmetis.*


checkin:
	@for file in *.[c,h]; \
	do \
	ci -u -m'Maintance' $$file;\
	done 

checkin2:
	@for file in *.[c,h]; \
	do \
	ci $$file;\
	rcs -U $$file;\
	co $$file;\
	done
