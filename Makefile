.PHONY	:	clean all test mmh steele

all	:
	./tools/build/make.py

mmh	:
	./tools/build/make.py mmh

steele	:
	./tools/build/make.py steele

clean	:
	./tools/build/make.py clean

test	:
	./tools/build/make.py test
