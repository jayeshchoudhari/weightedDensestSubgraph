# Densest Subgraph for Hypergraphs:

'densestHyperSubGraph-v2-MultiplyEdges-MultiInstance-newDelete-clearpending.cpp' is the latest file.

Compile: g++ -std=c++17 -Wall -g -o dhmi-nd-cp ./densestHyperSubGraph-v2-MultiplyEdges-MultiInstance-newDelete-clearpending.cpp


## Input file format
```
n max_cardinality
+ v1 v2 v3
+ v2 v3 v4
- v1 v2 v3
= timeval
+ v1 v3 v4
- v2 v3 v4
= timeval
```

Run: ./dhmi-nd-cp input-file output-file epsilon-value
