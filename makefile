
all:
	g++ -g -lm -Wall -o main molecular_system.cpp molecule.cpp parameter.cpp main.cpp 


molecule.o: molecule.cpp molecule.h
	g++ -c molecule.cpp

molecular_system.o: molecular_system.cpp molecular_system.h
	g++ -c molecular_system.cpp

parameter.o: parameter.cpp
	g++ -c parameter.cpp 
	

testunit.o: testunit.cpp
	g++ -c testunit.cpp 
	
testmain.o: testmain.cpp
	g++ -c testmain.cpp 

clean:
	rm -f *.o main testmain
