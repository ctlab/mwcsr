LIBFORMAT  = static_pic
BOOST_LDFLAGS= -lboost_program_options -lboost_filesystem -lboost_system  -Xlinker -rpath -Xlinker
CCLNFLAGS = -lm -pthread -lrt -lpthread

CFLAGS = -I. -fPIC -pedantic -Wall -Wextra -Werror -Wdisabled-optimization -std=c++11 

LDFLAGS = $(BOOST_LDFLAGS)
EXEC=mwcs

all: CFLAGS+= -O3 -DNDEBUG
all: $(EXEC)

debug: CFLAGS+= -g3 
debug: $(EXEC)

.SUFFIXES: .o .cpp .d
.cpp.o:
	$(CXX) -c $(CFLAGS) $< -o $@
.cpp.d:
	$(CXX) -MM $(CFLAGS) -MF $@ $<

SRCS=$(wildcard *.cpp) $(wildcard */*.cpp)

OBJS=$(SRCS:.cpp=.o)

DEPS=$(SRCS:.cpp=.d)

$(EXEC): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(CCLNFLAGS)

clean:
	rm -f $(OBJS) $(DEPS) $(TARGETS) $(EXEC)

sinclude $(DEPS)
