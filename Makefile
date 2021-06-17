udshp-dynamic-v2: udshp-dynamic.o GraphLoader.o DynamicGraph.o DynOpManager.o EdgeManager.o
	g++ -std=c++17 -O3 udshp-dynamic.o GraphLoader.o DynamicGraph.o DynOpManager.o EdgeManager.o -o udshp-dynamic-v2

udshp-dynamic.o: udshp-dynamic.cpp
	g++ -std=c++17 -O3 -c udshp-dynamic.cpp

DynOpManager.o: ./include/DynamicGraph.cpp ./include/DynamicGraph.h ./include/namespace.h
	g++ -std=c++17 -O3 -c ./include/DynOpManager.cpp

DynamicGraph.o: ./include/DynamicGraph.cpp ./include/DynamicGraph.h ./include/GraphLoader.cpp ./include/GraphLoader.h ./include/namespace.h
	g++ -std=c++17 -O3 -c ./include/DynamicGraph.cpp

EdgeManager.o: ./include/EdgeManager.h ./include/EdgeManager.cpp ./include/namespace.h
	g++ -std=c++17 -O3 -c ./include/EdgeManager.cpp

GraphLoader.o: ./include/GraphLoader.cpp ./include/GraphLoader.h ./include/namespace.h
	g++ -std=c++17 -O3 -c ./include/GraphLoader.cpp

clean:
	rm *.o udshp-dynamic-v2