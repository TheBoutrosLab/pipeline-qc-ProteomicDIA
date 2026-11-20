# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.3.2] - 2026-11-20

### Changed

- Migrate the pipeline license from GPL to BSD 3-Clause and update license documentation; shared code and third-party tools retain their respective licenses.

## [1.3.1] - 2026-11-18

### Fixed

- Let Nextflow cache Apptainer images before launching parallel QC tasks ([#19](https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/issues/19)).

## [1.3.0] - 2026-11-13

### Added

- Add a DIA QC workflow schematic showing input files, Python modules, metrics, MS1 visualization, and cohort aggregation.

### Changed

- Align DIA acquisition, identification, PSM-feature, inference, and quantification metrics with the DDA QC definitions and shared 12-column output contract.
- Count all retained filtered idXML hits and use sample-level protein-group quantity metrics from the final aggregated DIA-NN report.
- Refactor acquisition, identification, protein inference, library, quantification, plotting, and metric collection scripts into focused classes.
- Format Nextflow processes consistently with the DDA QC modules.

### Removed

- Remove DIA isolation-pattern and window-cycle metrics, rank-one validated-hit metrics, and DIA-NN report-row, run-count, and cohort-union missingness metrics.

## [1.2.1] - 2026-11-12

### Added

- Support multiple DIA acquisitions per sample with arbitrary fraction labels, fraction-specific acquisition and MS1 outputs, and sample-level DIA-NN metrics from the final aggregated report.
- Update the ProteinProphet integration case to use published multi-fraction E. coli metapipeline outputs.

## [1.2.0] - 2026-11-03

### Added

- Read original Bruker `.d` acquisitions for DIA scan and MS1 QC through pyOpenMS 3.6 while retaining reconstructed mzML paths for search and DIA-NN sample mapping.
- Add a HeLa dia-PASEF ProteinProphet integration case with verified published inputs and map DIA-NN raw-acquisition run names to their sample names.

## [1.1.1] - 2026-10-28

### Added

- Count pre-validation MSBooster PSMs by missed cleavages, enzyme cleavage specificity, and precursor charge for all, target, and decoy groups.

## [1.1.0] - 2026-10-27

### Added

- Summarize published MSBooster PIN and PSMFeatureExtractor idXML features from an optional manifest, including target/decoy mass-error and retention-time-deviation distributions when available.
- Add stable `layer`, `source_process`, search-engine, fraction, and availability columns to stage and merged QC metrics tables.

### Changed

- Leave Nextflow log files uncompressed by removing the workflow finalizer call.

## [1.0.0] - 2026-10-26

### Added

- Initial Nextflow DSL2 DIA quality-control pipeline for completed metapipeline outputs, with versioned Boutros Lab configuration, container profiles, resource allocation, and a combined long-form QC metrics TSV.
- Per-sample mzML acquisition summaries and an offline interactive 3D MS1 peak map generated through separate pyOpenMS extraction and Plotly rendering processes.
- Unfiltered MSFragger search and final filtered identification summaries from idXML using pyOpenMS.
- Cohort protein-group summaries for ProteinProphet, Fido, Epifany, and ProtLoom; EasyPQP library and export summaries; and per-sample metrics from the DIA-NN main report.
- UCLA nf-test integration cases for all four protein-inference routes, with case-specific metapipeline inputs and shared E. coli mzML files.
- Run, input, output, and testing documentation.

[1.0.0]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/releases/tag/v1.0.0
[1.1.0]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/compare/v1.0.0...v1.1.0
[1.1.1]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/compare/v1.1.0...v1.1.1
[1.2.0]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/compare/v1.1.1...v1.2.0
[1.2.1]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/compare/v1.2.0...v1.2.1
[1.3.0]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/compare/v1.2.1...v1.3.0
[1.3.1]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/compare/v1.3.0...v1.3.1
[1.3.2]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/compare/v1.3.1...v1.3.2
[unreleased]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/compare/v1.3.2...HEAD
