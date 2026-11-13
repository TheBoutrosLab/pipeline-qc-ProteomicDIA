# DIA QC integration tests

The UCLA nf-test suite is defined in `nftest.yml`. The ProteinProphet case uses six fractionated E. coli mzML files.
Fido, Epifany, and ProtLoom use two unfractionated E. coli mzML files. Each case reads
search idXML, filtered idXML, protein groups, the EasyPQP library and
export summary, and the DIA-NN report from its own metapipeline run directory
under `TEST_DATA_DIR`. Their case-local CSVs expand `$TEST_DATA_DIR` before
Nextflow reads them. ProtLoom uses a `groups.csv` manifest listing the per-sample
`*_groups.tsv` files from its cohort inference run. The mzML manifests use `sample_name,fraction_id,path`;
the ProteinProphet case has three rows per sample.

Load the lab's `nftest` and Nextflow environment, then run inside a Slurm job:

```bash
bash test/run-nftest.sh test-integration-Fido
bash test/run-nftest.sh test-integration-ProteinProphet
bash test/run-nftest.sh test-integration-Epifany
bash test/run-nftest.sh test-integration-ProtLoom
```

The runner reads a Git-ignored `.env` file with `TEST_DATA_DIR` for these cases. `.env-template` pins
Nextflow 24.10.6 and the v1 syntax parser. It writes work and output to
`/scratch/$SLURM_JOB_ID/<test-case>` and retains them for inspection.
These tests use the pipeline's configured container profile on the cluster.

The `test-integration-HeLa-ProteinProphet` case uses the original HeLa dia-PASEF `.d`
directory for acquisition and MS1 QC, with the diaTracer reconstructed mzML for sample
mapping. Its manifests reference the published HeLa search, validation, and library
outputs under `$TEST_DATA_DIR/hela_diapasef/output`, including the published DIA-NN
report. The report's `Run` value is the original Bruker `.d` stem, which the QC reader
maps to `hela_50ng_5_6min` through `acquisition_path`.

The `test-integration-ProteinProphet` case uses the E. coli DIA 1.5.0 output,
with three acquisitions per sample (`F1`–`F3`).
It uses the published combined DIA-NN report, which yields sample-level
quantification metrics.
