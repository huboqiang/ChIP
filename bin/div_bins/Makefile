cc=g++
exe=bed_read
obj=main.o bed_read.o gzstream.o 

$(exe):$(obj)
	$(cc) -o $(exe) $(obj) -lz
main.o:main.cpp bed_read.h
	$(cc) -c main.cpp -lz
bed_read.o:bed_read.cpp bed_read.h
	$(cc) -c bed_read.cpp -lz
gzstream.o:gzstream.cpp gzstream.h
	$(cc) -c gzstream.cpp -lz
clean:
	rm -rf *.o $(exe)
