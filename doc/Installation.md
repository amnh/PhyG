Manual installation and building

-- Get Haskell ghc compiler

curl --proto '=https' --tlsv1.2 -sSf https://get-ghcup.haskell.org | sh 

-- invoke ghcup 

ghcup tui

install ghc-9.10.1
install cabal 3.12

-- update cabal libraries 

cabal update

-- build phyg

cabal build PhyG:phyg --with-compiler ghc-9.10.1 --flags="Super-Optimization Use-LLVM"

-- install manually
-- note the very very long path for phyg binary and copy to a convenient location
