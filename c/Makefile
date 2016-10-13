all: librgt.so

librgt.so: librgt.c
	$(CC) -Wall -g -fPIC -shared -o $@ $? -lc -O0

test: test_main.c librgt.o
	$(CC) -o $@ $?

testmodule: testmodule.c
	python setup.py build

clean:
	rm -f test librgt.o librgt.so *.pyc
	rm -rf build
