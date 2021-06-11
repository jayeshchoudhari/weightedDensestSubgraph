udshp-dynamic-v2: udshp-dynamic.o GraphLoader.o DynamicGraph.o
	g++ -O3 udshp-dynamic.o GraphLoader.o DynamicGraph.o -o udshp-dynamic-v2

udshp-dynamic.o: udshp-dynamic.cpp ./include/namespace.h
	g++ -O3 -c udshp-dynamic.cpp

DynamicGraph.o: ./include/DynamicGraph.cpp ./include/DynamicGraph.h ./include/GraphLoader.cpp ./include/GraphLoader.h ./include/namespace.h
	g++ -O3 -c ./include/DynamicGraph.cpp

GraphLoader.o: ./include/GraphLoader.cpp ./include/GraphLoader.h ./include/namespace.h
	g++ -O3 -c ./include/GraphLoader.cpp