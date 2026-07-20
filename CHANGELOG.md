# Change Log
All notable changes to this project will be documented in this file.

## [1.1.0] - 2026-05-20

### Changed
In this version, the code has been completely restructured. Unnecessary logic has been removed, and all functions have been encapsulated to minimize memory errors (e.g., array out-of-bounds and unreleased memory) and improve the running speed.

### Added
- Compilation support for 32-bit systems.
- Output the complete GenBank file containing the sequence.
- Specify the topology for each sequence in one assembly.
- Complete RNA genome processing capability.

## [1.0.6] - 2026-04-22

### Changed
- Disable the TriTISA trainer in the metagenomic mode.
- Force the output of partial-5'-end ORFs with scores >= threshold for linear genomes. 

## [1.0.5] - 2026-04-03

### Added
- Added RNA sequence handling support for TriTISA+ and protein translation.
- Added translation table #15 (Blepharisma Nuclear Code).
- Incorporating the RBS score into the determination of the final score.
- Added version information display in the help menu.

### Changed
- Optimized the processing procedure for the phage genomes.
- Fixed some errors present in the Python scripts.
- Optimized the overlapping gene reporter.
- Fixed the calculation error of the Bayesian factor in TriTISA+.
- For the protein output, replace non-AGCT chars with "X"s (unknown amino acids) instead of reporting an error.
- Fixed the issue where the external SVM model could not be loaded in certain cases.

## [1.0.4] - 2026-03-23

### Added
- Added auto-detection of translation table before training the model.
- Added more translation tables, and now there are 5 tables in total (1, 4, 11, 16, 25).

## [1.0.3] - 2026-03-20

### Added
- Added EMBL format support for input files.
- Added third-party terms for TriTISA to the LICENSE file.

### Changed
- Fixed a possible array out-of-bounds issue in the function `mm_train`.
- TIS revision has now been moved to take place before the CDS training.

## [1.0.2] - 2026-03-10

### Added
- Added the parameter `minolen` to filter overprinted genes with a minimum overprinted length between two ORFs;
- Added the length/topology of each sequence and date information to the output (GFF3, GBK and MED).

### Changed
- Optimized the exception handling of all the modules to improve performance;
- Fix the issue where GBK output cannot be parsed by Biopython.

## [1.0.1] - 2026-03-09

### Added
- Pretrained markov model for TriTISA+ to complete TIS revision on sequences where the samples are too few to initiate self-training.

### Changed
- Fixed the bug that overprinted genes cannot be detected at the boundaries of the circular genomes;
- Fixed the incorrect example format in the README file;
- Fixed the issue where svm.h could not be compiled on machines without AVX support.
- Fixed compilation issues on machines without OpenMP support.