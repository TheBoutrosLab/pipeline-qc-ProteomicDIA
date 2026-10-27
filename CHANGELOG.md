# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2026-10-27

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

[unreleased]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDIA/compare/v1.0.0...HEAD
