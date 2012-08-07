CPP=g++ 
INCLUDE=-I ./include
LIBS=
CPPFLAGS=-g -O3 -w $(INCLUDE) $(LIBS)
PROGRAM=Kavosh

SRC=src/main.cpp\
    src/graph.cpp\
    src/ZeroOneTree.cpp\
    src/randomGenerator.cpp

OBJ=$(SRC:.cpp=.o)

NAUTY = nauty/nauty.o\
    	nauty/naugraph.o\
        nauty/nautil.o

all : $(PROGRAM)

%.o: src/%.cpp include/%.h
	$(CPP) $(CPPFLAGS) $< -o $@

nauty/%.o: nauty/%.c
	$(CPP) $(CPPFLAGS) -c $< -o $@

$(PROGRAM): $(OBJ) $(NAUTY)
	$(CPP) $(CPPFLAGS) -o $@ $(OBJ) $(NAUTY)
clean:
	@echo "Removing objects..."
	rm -f $(PROGRAM) *.log src/*.o *~ core
	rm -f nauty/*.o 
	rm -f result/*

install: $(PROGRAM)
	cp $(PROGRAM) ../bin

ar:
	make clean
	tar -czvf ../$(PROGRAM)"(`date`)".tar.gz *
