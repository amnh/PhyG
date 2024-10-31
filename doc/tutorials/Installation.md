Manual installation and building

-- Get Haskell ghc compiler

curl --proto '=https' --tlsv1.2 -sSf https://get-ghcup.haskell.org | sh

-- invoke ghcup

ghcup tui

install ghc-9.10.1

install cabal 3.12

-- update cabal libraries

cabal update

install LLVM via your package manager (e.g. homebrew, apt). The process for installation 
of LLVM varies from machine to machine. The user should consult their package manager
for information. 

-- build phyg (after you have installed LLVM, or delete Use-LLVM from flags option)

cabal build PhyG:phyg --with-compiler ghc-9.10.1 --flags="Super-Optimization Use-LLVM"

-- install manually

note the very very long path for phyg binary and copy the desired location

