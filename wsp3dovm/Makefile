# OpenVolumeMesh was installed in 
# /usr/local/include/ resp. lib and will be found
# forexecution, adjust LD_LIBRARY_PATH properly

# boost 1.54 was installed in 
# /usr/include/boost/
# boost does not support pkg-config

CPPFLAGS=-std=c++11 -O3 -DNDEBUG -DBOOST_ALL_NO_LIB
LDFLAGS=-std=c++11 

LIBS=-lboost_system -lboost_filesystem -lboost_graph -lboost_chrono -L/usr/local/lib/OpenVolumeMesh/ -lOpenVolumeMesh

SRC=create_steinerpoints.cpp read_tet.cpp write_vtk.cpp statistics.cpp main.cpp
OBJ=$(SRC:.cpp=.o)

PROG=wsp3dovm

.PHONY: clean

$(PROG): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX)  -c $(CPPFLAGS) -o $@ $<
	
clean:
	- rm $(OBJ)
	- rm $(PROG)

