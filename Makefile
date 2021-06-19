udshp-dynamic-v2: udshp-dynamic.o GraphLoader.o DynamicGraph.o DynOpManager.o EdgeManager.o
	g++ -std=c++17 -O3 -g udshp-dynamic.o GraphLoader.o DynamicGraph.o DynOpManager.o EdgeManager.o -o udshp-dynamic-v2

udshp-dynamic.o: udshp-dynamic.cpp
	g++ -std=c++17 -O3 -g -c udshp-dynamic.cpp

DynOpManager.o: ./include/DynOpManager.cpp ./include/DynOpManager.h ./include/namespace.h
	g++ -std=c++17 -O3 -g -c ./include/DynOpManager.cpp

DynamicGraph.o: ./include/DynamicGraph.cpp ./include/DynamicGraph.h ./include/GraphLoader.cpp ./include/GraphLoader.h ./include/namespace.h
	g++ -std=c++17 -O3 -g -c ./include/DynamicGraph.cpp

EdgeManager.o: ./include/EdgeManager.h ./include/EdgeManager.cpp ./include/namespace.h
	g++ -std=c++17 -O3 -g -c ./include/EdgeManager.cpp

GraphLoader.o: ./include/GraphLoader.cpp ./include/GraphLoader.h ./include/namespace.h
	g++ -std=c++17 -O3 -g -c ./include/GraphLoader.cpp

clean:
	rm *.o udshp-dynamic-v2
