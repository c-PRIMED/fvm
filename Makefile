.PHONY	:	clean all test

all	:
	./tools/build/make.py

clean	:
	./tools/build/make.py clean

test	:
	./tools/build/make.py test
