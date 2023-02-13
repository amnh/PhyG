# PhyGraph
Progam to perform phylogenetic searches on general graphs with diverse data types


## Installation

Precompiled binaries for linux (Intel), OSX (Intel and M1) are in the "bin" directory.

Documentation is in the "doc" directory.

Test files are in the "testData" directory.

## Make from Source

Compilation is verified up to ghc-9.4.4 with cabal v. 3.8.

```
cabal build PhyGraph:phyg --flags=super-optimization
```
