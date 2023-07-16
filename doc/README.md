# Phylogenetic Graph (PhyG)

Performs heuristic search of phylogenetic graph space via scoring abstract input data.

[![CI    Status        ][GitHub-actions-img]][GitHub-actions-ref]
[![Maintained          ][GitHub-support-img]][GitHub-support-ref]

[![Release             ][GitHub-release-img]][GitHub-release-ref]
[![Release Date        ][GitHub-tagdate-img]][GitHub-release-ref]

[![New Commits         ][GitHub-commits-img]][GitHub-commits-ref]
[![Code Size           ][GitHub-codelen-img]][GitHub-codelen-ref]

[![Author              ][GitHub-authors-img]][GitHub-authors-ref]
[![BSD3 License        ][GitHub-license-img]][GitHub-license-ref]
[![Contributor Covenant][GitHub-conduct-img]][GitHub-conduct-ref]

The `PhyG` package exposes the `phyg` executable and the `PhyG` sub-libraries.

## `phyg`

Phylogenetic Graph (PhyG) is a multi-platform programram designed to produce phylogenetic graphs from input data and graphs via heuristic searching of general phylogenetic graph space.
The bio-informatics framework libraries of the broader [Phylogenetic Haskell Analytic Network Engine (PHANE) project][GitHub-PHANE-Readme] are the foundation upon which PhyG is constructed.
PhyG offers vast functionality, including the optimization of unaligned sequences, and the ability to implement search strategies such as random addition sequence, swapping, and tree fusing.
Furthermore, PhyG can generate outputs in the form of implied alignments and graphical representations of cladograms and graphs.
What sets PhyG apart from other phylogenetic analysis programs, is the extension to broader classes of input data and phylogenetic graphs.
The phylogenetic graph inputs and outputs of PhyG include trees, as well as softwired and hardwired networks.


### [Funding provided by][GitHub-Funding]:

  * [American Museum of Natural History][Funding-0]

  * [DARPA SIMPLEX][Funding-1]

  * [Kleberg Foundation][Funding-2]


### [Installation instructions][GitHub-Install]

```
which ghcup || curl --proto '=https' --tlsv1.2 -sSf https://get-ghcup.haskell.org | sh
which ghc   || ghcup install ghc   latest
which cabal || ghcup install cabal latest
cabal update
cabal install
```

### [Publications][GitHub-PHANE-Papers]

[Funding-0]: https://www.amnh.org/our-research/computational-sciences
[Funding-1]: https://www.darpa.mil/program/simplifying-complexity-in-scientific-discovery
[Funding-2]: http://www.klebergfoundation.org/

[GitHub-actions-img]: https://github.com/amnh/PhyG/actions/workflows/integration-test-suite.yaml/badge.svg?branch=master
[GitHub-actions-ref]: https://github.com/AMNH/PhyG/actions
[GitHub-authors-img]: https://img.shields.io/badge/author-Ward%20Wheeler-blue.svg?color=134EA2
[GitHub-authors-ref]: https://github.com/AMNH/PhyG/tree/master/doc/AUTHORS.md
[GitHub-codelen-img]: https://img.shields.io/github/languages/code-size/AMNH/PhyG.svg?style=popout&color=yellowgreen
[GitHub-codelen-ref]: https://github.com/AMNH/PhyG/archive/master.zip
[GitHub-commits-img]: https://img.shields.io/github/commits-since/AMNH/PhyG/v0.1.2.svg?style=popout&color=yellowgreen
[GitHub-commits-ref]: https://github.com/AMNH/PhyG/commits/master
[GitHub-conduct-img]: https://img.shields.io/badge/Contributor%20Covenant-2.0-4baaaa.svg
[GitHub-conduct-ref]: https://github.com/AMNH/PhyG/blob/master/doc/Code_Of_Conduct.md
[GitHub-license-img]: https://img.shields.io/badge/license-BSD3-blue.svg?color=134EA2
[GitHub-license-ref]: https://github.com/AMNH/PhyG/blob/master/doc/LICENSE
[GitHub-release-img]: https://img.shields.io/github/release-pre/AMNH/PhyG.svg?style=popout&color=orange
[GitHub-release-ref]: https://github.com/AMNH/PhyG/releases/latest
[GitHub-tagdate-img]: https://img.shields.io/github/release-date-pre/AMNH/PhyG.svg?style=popout&color=orange
[GitHub-support-img]: https://img.shields.io/maintenance/yes/2023.svg?style=popout
[GitHub-support-ref]: https://github.com/AMNH/PhyG/graphs/contributors

[GitHub-Funding]: https://github.com/AMNH/PhyG/blob/master/doc/Funding.md
[GitHub-Install]: https://github.com/AMNH/PhyG/blob/master/doc/tutorials/Installation.md

[GitHub-PHANE-Readme]: https://github.com/AMNH/PHANE#readme
[GitHub-PHANE-Papers]: https://github.com/AMNH/PHANE/blob/master/doc/Publications.md
