cc=g++
exe=get_bin
obj=main.o get_bin.o

$(exe):$(obj)
		   $(cc) -o $(exe) $(obj)
main.o:main.cpp get_bin.h
		   $(cc) -c main.cpp
get_bin.o:get_bin.cpp get_bin.h
		   $(cc) -c get_bin.cpp
clean:
		   rm -rf *.o $(exe)
