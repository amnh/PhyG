
## Changes from `measurability`

	* Major refacotring of SCMs and TCMs

	* Pre-computed specialized measures


## Changes from `dynamic-character-refactor`

	* Corrected computational defect in `Data.TCM.Overlap`

	* Updated `CompactMeasure` to correctly store, computute, and memoize TCMs

	* Fixed numerous defects in Ukkonen string alignment

	* Changed memory layout of dynamic characters

	* Changed gap state index from `n-1` to `0` (least significant bit)

	* Updated `Data.Alphabet` to handle gaps as first symbol

	* Updated C FFI code to handle gap state index as east significant bit

	* Corrected `insertGaps` defect

	* Refactored `DirectOptimization.Pairwise` sub-modules

	* Changed  `Input.*` modules to use `Data.Alphabet` for consistency and correctness

	* Added `test-dynamic-character` as a pairwise alignment testsuite

	* Improved efficiency and memory usage of pairwise string alignment tracback
