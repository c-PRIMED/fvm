include ../Makefile.in

.SUFFIXES: .c .lo .o

all: ../libparmetis.a ../libparmetis.so

CFLAGS = $(COPTIONS) $(OPTFLAGS) -I. $(INCDIR)


OBJS = comm.o util.o debug.o setup.o grsetup.o timer.o \
       node_refine.o initmsection.o order.o \
       xyzpart.o pspases.o frename.o \
       iintsort.o iidxsort.o ikeysort.o ikeyvalsort.o \
       kmetis.o gkmetis.o ometis.o \
       initpart.o match.o \
       kwayfm.o kwayrefine.o kwaybalance.o \
       remap.o stat.o fpqueue.o \
       ametis.o rmetis.o lmatch.o initbalance.o \
       mdiffusion.o diffutil.o wave.o \
       csrmatch.o redomylink.o balancemylink.o \
       selectq.o akwayfm.o serial.o move.o \
       mmetis.o mesh.o memory.o weird.o backcompat.o 

.c.o:
	$(CC) $(CFLAGS) -c $*.c

SHLOBJS = $(OBJS:.o=.lo)

.c.lo:
	$(CC) $(CFLAGS) -fPIC -DPIC -c $< -o $@


../libparmetis.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

../libparmetis.so: $(SHLOBJS)
	$(CC) -shared $^ -Wl,-soname,libparmetis.so.3.1.1 -o $@.3.1.1 \
	  $(LIBDIR) -lmpi
	ln -s libparmetis.so.3.1.1 $@

clean:
	rm -f *.o *.lo

realclean:
	rm -f *.o *.lo; rm -f ../libparmetis.*

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
