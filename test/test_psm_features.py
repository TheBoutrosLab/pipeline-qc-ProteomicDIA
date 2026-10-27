"""Check feature QC against representative published PIN and idXML inputs."""

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pyopenms as oms


ROOT = Path(__file__).resolve().parents[1]


def metric_rows(path):
    """Read a stage or cohort metric TSV."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class FeatureQCTests(unittest.TestCase):
    """Exercise both published validation feature formats."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def summarize(self, input_path, route, output_name):
        """Run the same entry point used by Nextflow."""
        output = self.folder / output_name
        subprocess.run([
            sys.executable, str(ROOT / "src/SummarizePSMFeatures.py"),
            "--input", str(input_path), "--output", str(output),
            "--dataset-id", "cohort", "--sample-name", "sample1",
            "--search-engine", "MSFragger", "--route", route,
            "--fraction-id", "1" if route == "MSBooster" else "",
        ], check=True)
        return output, metric_rows(output)

    def test_msbooster_pin_keeps_target_decoy_and_fraction(self):
        """Signed mass and RT differences retain their source and groups."""
        pin = self.folder / "run_01-fraction_1_edited.pin"
        pin.write_text(
            "SpecId\tLabel\tScanNr\tdeltamass\tCalcMass\tdeltaRT\tspectralAngle\n"
            "a\t1\t1\t0.01\t500\t1.5\t0.8\n"
            "b\t-1\t2\t-0.02\t500\t-2.5\t0.2\n"
            "c\t1\t3\t\t500\t\t0.6\n", encoding="utf-8")
        _, rows = self.summarize(pin, "MSBooster", "pin-metrics.tsv")
        lookup = {(row["metric_name"], row["group"]): row for row in rows}
        self.assertEqual(lookup[("feature_psm_count", "all")]["value"], "3")
        self.assertEqual(lookup[("feature_psm_count", "target")]["value"], "2")
        self.assertEqual(lookup[("feature_deltamass_missing_count", "target")]["value"], "1")
        mass_median = lookup[("feature_mass_error_ppm_signed_p50", "all")]["value"]
        self.assertAlmostEqual(float(mass_median), -10)
        self.assertEqual(lookup[("feature_deltart_p50", "target")]["value"], "1.5")
        self.assertEqual(lookup[("feature_spectralangle_p50", "decoy")]["value"], "0.2")
        self.assertTrue(all(row["layer"] == "psm_features" and
                            row["source_process"] == "MSBooster" and
                            row["fraction_id"] == "1" for row in rows))

    def test_openms_idxml_reads_declared_features(self):
        """The idXML reader uses OpenMS extra_features and hit metadata."""
        protein = oms.ProteinIdentification()
        protein.setIdentifier("run")
        search = oms.SearchParameters()
        search.setMetaValue("extra_features", "deltamass,CalcMass,deltaRT,spectralAngle")
        protein.setSearchParameters(search)
        identification = oms.PeptideIdentification()
        identification.setIdentifier("run")
        hits = []
        for label, delta in (("target", 0.01), ("decoy", -0.02)):
            hit = oms.PeptideHit()
            hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
            hit.setCharge(2)
            hit.setScore(1.0)
            hit.setMetaValue("target_decoy", label)
            hit.setMetaValue("deltamass", delta)
            hit.setMetaValue("CalcMass", 500.0)
            hit.setMetaValue("deltaRT", 1.0 if label == "target" else -2.0)
            hit.setMetaValue("spectralAngle", 0.8 if label == "target" else 0.2)
            hits.append(hit)
        identification.setHits(hits)
        path = self.folder / "features.idXML"
        oms.IdXMLFile().store(str(path), [protein], oms.PeptideIdentificationList([identification]))
        _, rows = self.summarize(path, "PSMFeatureExtractor", "idxml-metrics.tsv")
        lookup = {(row["metric_name"], row["group"]): row for row in rows}
        self.assertEqual(lookup[("feature_psm_count", "all")]["value"], "2")
        self.assertEqual(lookup[("feature_deltart_p50", "decoy")]["value"], "-2.0")
        self.assertEqual(lookup[("feature_mass_error_ppm_absolute_p50", "target")]["value"], "20.0")
        self.assertTrue(all(row["source_process"] == "PSMFeatureExtractor" and
                            row["fraction_id"] == "" for row in rows))

    def test_missing_features_are_unavailable_in_merged_table(self):
        """The merger preserves layer, source, and unavailable values."""
        pin = self.folder / "without-rt.pin"
        pin.write_text("SpecId\tLabel\tScanNr\nA\t1\t1\n", encoding="utf-8")
        stage, rows = self.summarize(pin, "MSBooster", "stage.tsv")
        missing = next(row for row in rows if row["metric_name"] == "feature_rt_delta_p50" and
                       row["group"] == "all")
        self.assertEqual((missing["status"], missing["value"]), ("unavailable", ""))
        merged = self.folder / "merged.tsv"
        subprocess.run([sys.executable, str(ROOT / "src/CollectQCMetrics.py"),
                        "--inputs", str(stage), "--output", str(merged)], check=True)
        self.assertEqual(metric_rows(merged), sorted(rows, key=lambda row: (
            row["dataset_id"], row["layer"], row["sample_name"],
            row["search_engine"], row["fraction_id"], row["metric_name"],
            row["group"], row["source_file"])))

    def test_merged_table_uses_pipeline_layer_order(self):
        """Layer ordering follows the pipeline rather than alphabetic names."""
        pin = self.folder / "one.pin"
        pin.write_text("SpecId\tLabel\tScanNr\nA\t1\t1\n", encoding="utf-8")
        stage, rows = self.summarize(pin, "MSBooster", "stage.tsv")
        layer_order = (
            "acquisition", "raw_identification", "psm_features",
            "validated_identification", "protein_inference", "spectral_library",
            "quantification",
        )
        with stage.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t")
            writer.writeheader()
            for layer in reversed(layer_order):
                writer.writerow({**rows[0], "layer": layer})
        merged = self.folder / "ordered.tsv"
        subprocess.run([sys.executable, str(ROOT / "src/CollectQCMetrics.py"),
                        "--inputs", str(stage), "--output", str(merged)], check=True)
        self.assertEqual([row["layer"] for row in metric_rows(merged)], list(layer_order))


if __name__ == "__main__":
    unittest.main()
