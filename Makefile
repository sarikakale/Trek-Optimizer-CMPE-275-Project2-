# location of the Python header files
 

INCFLAG = /usr/include/python2.7
 
# location of the Boost Python include files and library
 
BOOSTFLAG = /usr/local/gcc/boost_1_60_0/lib
 
TARGET = TrekOptimizer
 
$(TARGET).so: $(TARGET).o
	g++ -shared -lgomp -fopenmp -Wl,--export-dynamic $(TARGET).o -L$(BOOSTFLAG) -lboost_python -lboost_chrono -L/usr/lib/python2.7/config -lpython2.7 -o $(TARGET).so -std=c++1y
 
$(TARGET).o: $(TARGET).cpp
	g++ -fopenmp -I $(INCFLAG) -L$(BOOSTFLAG) -lboost_python -lboost_chrono -L/usr/lib/python2.7/config -lpython2.7 -fPIC -c $(TARGET).cpp -std=c++1y
 