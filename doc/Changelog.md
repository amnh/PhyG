Phylogenetic Haskell Analytic Network Engine (PHANE)
====================================================

# Changelog

PHANE uses [Semantic Versioning (v2.0.0)][1].
The changelog is available [on GitHub][2].


## Unreleased (`v0.2.0`)

### PhyGraph

  * Major codebase layout rearchitecting

  * Modularized and migrated sub-libraries to [PHANE project][GitHub-PHANE]

    - `PHANE-alphabet`

    - `PHANE-dynamic-character`

    - `PHANE-dynamic-character-element`

    - `PHANE-measure-units`

    - `PHANE-PhyloLib`

  * Added command line options

    - `--credits` To list financial and technical contributors

    - `--license` Emits the license under which the software is distributed

    - `--splash` A nice ASCII art title "splash screen"

    - `--version` Outputs the software version, commit hash, and compilation timestamp

  * Corrected defect in median extraction from preliminary contexts


### PhyloLib


## `v0.1.0`

  * Initial "alpha" state of PHANE domain


[1]: https://semver.org/spec/v2.0.0.html
[2]: https://github.com/wardwheeler/PhyGraph/blob/main/doc/Changelog.md
[GitHub-PHANE]: https://github.com/AMNH/PHANE#readme
