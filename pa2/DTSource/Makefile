# Makefile to create the library.  DataTank will call this script,
# so be careful when modifying.  This can be used as a blueprint for
# porting the library to other unix platforms.

include makefile.sources
OBJS = $(C_SRCS:.c=.o) $(F_SRCS:.f=.o) $(CXX_SRCS:.cpp=.o)

CC = cc
CXX =  c++

.SUFFIXES: .cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

CXXFLAGS = -O3 -g

libDT.a : $(OBJS)
	ar -r -u libDT.a $(OBJS)
	ranlib libDT.a

missing : $(OBJS)
	make

clean :
	rm -f $(OBJS)
	rm libDT.a

# Re-create the file of dependencies
depend : 
	rm -f makefile.depend
	gcc $(HEADERS) -MM $(CXX_SRCS) > makefile.depend

# Read dependencies
include makefile.depend
