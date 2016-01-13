include ../inc/common.mk
CC = gcc
PY = ${HOME}/local/include/python2.7/
O = wWigIO.so
wWigIO:
	${CC}  -w -shared -fPIC -DMACHTYPE_x86_64 -p -I ${PY} -I ../inc wWigIO.c ../lib/jkweb.a -lm -lz -o ${O}

clean:
	if [ -f ${O} ]; then rm ${O}; fi
