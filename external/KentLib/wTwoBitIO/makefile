CC = gcc
PY = /usr/include/python2.6/
O = wTwoBitIO.so
wTwoBitIO:
	${CC}  -w -shared -fPIC -DMACHTYPE_x86_64 -p -I ${PY} -I ../inc wTwoBitIO.c ../lib/jkweb.a -lm -lz -o ${O}

clean:
	if [ -f ${O} ]; then rm ${O}; fi
