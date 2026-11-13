"""Helpers for cohort protein-group QC from OpenMS idXML."""

import pyopenms as oms

from QCMetrics import Metrics, as_text, finite


class ProteinGroupSummarizer:
    """Summarize cohort protein groups from an OpenMS idXML file."""

    def __init__(self, input_path, dataset_id, route):
        self.input_path = input_path
        self.route = route
        self.metrics = Metrics(dataset_id, "", input_path, "protein_inference")
        self.groups = []

    def load(self):
        """Read protein groups and validate their membership and probability."""
        protein_ids = []
        peptide_ids = oms.PeptideIdentificationList()
        oms.IdXMLFile().load(str(self.input_path), protein_ids, peptide_ids)
        for identification in protein_ids:
            selected = (identification.getProteinGroups() if self.route == "ProteinProphet"
                        else identification.getIndistinguishableProteins())
            for group in selected:
                probability = finite(group.probability)
                members = [as_text(member) for member in group.accessions]
                if probability is None or not members:
                    raise ValueError(f"Invalid {self.route} protein group in {self.input_path}")
                self.groups.append((members, probability))

    def summarize(self):
        """Count groups and summarize member counts and probabilities."""
        self.metrics.add("protein_group_count", len(self.groups))
        self.metrics.distribution(
            name="proteins_per_group",
            values=[len(members) for members, _ in self.groups],
        )
        self.metrics.distribution(
            name="protein_group_probability",
            values=[probability for _, probability in self.groups],
        )

    def run(self, output_path):
        """Read, summarize, and write cohort protein-group QC."""
        self.load()
        self.summarize()
        self.metrics.write(output_path)
