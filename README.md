# pipeline-qc-ProteomicDIA

This Nextflow DSL2 pipeline collects quality-control metrics from a completed `metapipeline-ProteomicDIA` cohort. It summarizes each sample's mzML, unfiltered MSFragger search idXML, and final filtered idXML; the selected protein-inference route; the cohort EasyPQP library; and the DIA-NN report. It also writes an interactive 3D MS1 peak map for each mzML.

## How to run

1. Initialize Git submodules with `git submodule update --init --recursive`.
2. Copy [config/template.config](config/template.config) to a run-specific config and set all input paths, `inference_route`, `dataset_id`, and `output_dir`.
3. Run with Nextflow 24.10.6 and the v1 syntax parser:

   ```bash
   export NXF_VER=24.10.6
   export NXF_SYNTAX_PARSER=v1
   nextflow run . -c /path/to/qc.config -profile apptainer
   ```

The `standard` and `docker` profiles use Docker; `apptainer` uses Apptainer. The OpenMS processes use `ghcr.io/theboutroslab/pyopenms:3.5.0`; the HTML plot uses `multiqc/multiqc:v1.35`, which provides Plotly. Set `work_dir` to scratch storage when processing large mzML files.

## Inputs

Three CSV manifests identify the same set of samples. Paths in these CSVs may be absolute or relative to the CSV directory. The current pipeline supports one mzML per sample; `fraction_id` is retained for future fraction support.

**mzML CSV** (`mzml_csv`):

```csv
sample_name,fraction_id,path
sample_01,1,/data/sample_01.mzML
sample_02,1,/data/sample_02.mzML
```

The two-column `sample_name,mzML` form is also accepted.

**Unfiltered search idXML CSV** (`search_idxml_csv`):

```csv
sample_name,search_engine,idXML
sample_01,MSFragger,/data/sample_01_MSFragger.idXML
sample_02,MSFragger,/data/sample_02_MSFragger.idXML
```

**Final filtered idXML CSV** (`filtered_idxml_csv`):

```csv
sample_name,idXML
sample_01,/data/sample_01_MSFragger_filter.idXML
sample_02,/data/sample_02_MSFragger_filter.idXML
```

The search CSV can be the metapipeline's `validation.csv`. Create the filtered CSV from the final validation outputs; the metapipeline does not currently publish it. Sample names must be unique and identical across all three CSVs.

| Parameter | Description |
| --- | --- |
| `dataset_id` | Cohort identifier used in the output path and metric filenames. |
| `mzml_csv` | Per-sample mzML manifest. |
| `search_idxml_csv` | Unfiltered MSFragger search idXML manifest. |
| `filtered_idxml_csv` | Final filtered validation idXML manifest. |
| `inference_route` | `ProteinProphet`, `Fido`, `Epifany`, or `ProtLoom`. |
| `protein_groups` | Cohort ProteinProphet or Epifany idXML; cohort Fido groups JSON; or a ProtLoom CSV described below. |
| `library_tsv` | Published cohort EasyPQP library TSV. |
| `export_summary_tsv` | Validation EasyPQP export summary TSV with `metric,count` columns. |
| `diann_report` | DIA-NN main `report.tsv`. Its `Run` values must match a sample name, mzML filename, or mzML stem. |
| `output_dir` | Root directory for QC results. |
| `work_dir` | Optional Nextflow work directory. |
| `pyopenms_version` | pyOpenMS image tag; default `3.5.0`. |
| `docker_image_plotly` | Plotly image; default `multiqc/multiqc:v1.35`. |
| `base_resource_update` | Optional process CPU and memory overrides. |

ProteinProphet and Epifany use their cohort idXML outputs. Fido uses its cohort `*-Fido-groups.json`. For ProtLoom, set `protein_groups` to a CSV that points to the per-sample `*_groups.tsv` files produced by the cohort inference run:

```csv
sample_name,groups_tsv
sample_01,/data/validation/sample_01/ProtLoom-0.2.3/output/sample_01_groups.tsv
sample_02,/data/validation/sample_02/ProtLoom-0.2.3/output/sample_02_groups.tsv
```

## Pipeline steps

| Process | Input | Result |
| --- | --- | --- |
| `summarize_mzML` | One mzML per sample | Acquisition duration, scan counts and rates, cycle timing, isolation windows, and TIC. |
| `extract_MS1_matrix` | One mzML per sample | Bounded RT, m/z, and intensity matrix plus selection metadata. |
| `plot_MS1_landscape` | MS1 matrix and metadata | Offline interactive 3D stick spectrum in HTML. |
| `summarize_search_idXML` | One unfiltered search idXML per sample | Candidate and best-hit counts, search-score percentiles, and charge counts. |
| `summarize_filtered_idXML` | One final filtered idXML per sample | Accepted PSM and distinct peptide counts and score percentiles. |
| `summarize_ProteinGroups_*` | One selected inference route | Cohort protein-group counts, sizes, and probabilities when available. |
| `summarize_EasyPQP_library` | Library and export summary | Cohort transition, precursor, peptide, and export counts. |
| `summarize_DIANN_report` | DIA-NN main report | Per-sample quantified counts, intensity percentiles, and precursor missingness. |
| `collect_QCMetrics` | All stage metric TSVs | Combined QC metric table. |

The MS1 matrix contains up to 120 strongest positive peaks from each of 350 evenly spaced MS1 scans within m/z 300–1200. The HTML embeds Plotly and works offline. The plot colors encode within-sample log-intensity percentiles: yellow (0–25), orange (25–50), red (50–75), purple (75–95), and blue (95–100). Colors therefore do not represent the same absolute intensity across samples. Continuous-value metric summaries include a nonmissing count and the 5th, 25th, 50th, 75th, and 95th percentiles.

## Outputs

Results are published under `<output_dir>/qc-ProteomicDIA-1.0.0/<dataset_id>/`:

| Output | Description |
| --- | --- |
| `output/<dataset_id>_QC_metrics.tsv` | Combined long-form QC table. |
| `output/mzML/<sample_name>_mzML_metrics.tsv` | Per-sample mzML metrics. |
| `output/mzML/<sample_name>_MS1_matrix.tsv` and `.json` | Selected MS1 peaks and selection metadata. |
| `output/mzML/<sample_name>_MS1_landscape.html` | Interactive MS1 figure. |
| `output/search_idXML/<sample_name>_search_idXML_metrics.tsv` | Unfiltered search metrics. |
| `output/filtered_idXML/<sample_name>_filtered_idXML_metrics.tsv` | Final filtered identification metrics. |
| `output/ProteinGroups/<dataset_id>_ProteinGroups_metrics.tsv` | Cohort inference metrics. |
| `output/EasyPQP/<dataset_id>_EasyPQP_metrics.tsv` | Library and export metrics. |
| `output/DIA-NN/<dataset_id>_DIANN_metrics.tsv` | DIA-NN quantitation metrics. |
| `log-*/nextflow-log/` | Nextflow trace, timeline, report, and effective parameters. |

Each metric row records `dataset_id`, `sample_name` (empty for cohort metrics), `metric_name`, `group`, `value`, `unit`, and `source_file`. Programmed LC gradient duration and retention-time error are not derived from these inputs.

## Testing

The UCLA nf-test suite has integration cases for ProteinProphet, Fido, Epifany, and ProtLoom. All four use the same two *E. coli* mzML files and their route-specific metapipeline outputs. Set `TEST_DATA_DIR` in a local `.env` file using [`.env-template`](.env-template), then run a case inside a Slurm job:

```bash
bash test/run-nftest.sh test-integration-ProteinProphet
```

See [test/README.md](test/README.md) for the other case names and test-data layout. These tests require the lab's external metapipeline outputs; they are not self-contained fixtures.

## License

Licensed under GNU General Public License version 2. See [LICENSE](LICENSE).
