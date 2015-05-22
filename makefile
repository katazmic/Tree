
#GSL_L	=	/usr/local/lib/
#GSL_I	=	/usr/local/include/


GSL_L	=	/usr/lib/
GSL_I	=	/usr/include/



clean: rm *.o

hrchl_tree: work.o hrchl_tree.o
	g++  -I$(GSL_I) -L$(GSL_L) -fopenmp -lgsl -lgslcblas -lm   work.o hrchl_tree.o -o hrchl_tree

work.o: work.cpp hrchl_tree.h
	g++ -c work.cpp

hrchl_tree.o: hrchl_tree.cpp hrchl_tree.h
	g++ -c hrchl_tree.cpp 
