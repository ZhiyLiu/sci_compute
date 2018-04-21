# Makefile to create the library.  DataTank will call this script,
# so be careful when modifying.  This can be used as a blueprint for
# porting the library to other unix platforms.

CXX_SRCS = main.cpp

OBJS = $(CXX_SRCS:.cpp=.o)

CC = cc
CXX =  c++

.SUFFIXES: .cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

HEADERS = -IDTSource

CXXFLAGS = -O3 $(HEADERS)


runme : $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o runme DTSource/libDT.a

missing : $(OBJS)
	make

clean :
	rm -f $(OBJS)
	rm -f runme
