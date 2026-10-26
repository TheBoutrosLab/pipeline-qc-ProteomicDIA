# DIA QC integration tests

The UCLA nf-test suite is defined in `nftest.yml`. The ProteinProphet,
Fido, Epifany, and ProtLoom cases use the same two E. coli mzML files. Each case reads
search idXML, filtered idXML, protein groups, the EasyPQP library and
export summary, and the DIA-NN report from its own metapipeline run directory
under `TEST_DATA_DIR`. Their case-local CSVs expand `$TEST_DATA_DIR` before
Nextflow reads them. ProtLoom uses a `groups.csv` manifest listing the per-sample
`*_groups.tsv` files from its cohort inference run. The mzML manifests use `sample_name,fraction_id,path`;
only one row per sample is supported at present.

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
