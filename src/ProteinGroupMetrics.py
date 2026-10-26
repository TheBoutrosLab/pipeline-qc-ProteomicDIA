"""Helpers for cohort protein-group QC from OpenMS idXML."""

import pyopenms as oms

from QCMetrics import Metrics, as_text, finite


def summarize_idxml(path, dataset_id, output, route):
    protein_ids = []
    peptide_ids = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(str(path), protein_ids, peptide_ids)
    groups = []
    for identification in protein_ids:
        selected = (identification.getProteinGroups() if route == "ProteinProphet"
                    else identification.getIndistinguishableProteins())
        for group in selected:
            probability = finite(group.probability)
            members = [as_text(member) for member in group.accessions]
            if probability is None or not members:
                raise ValueError(f"Invalid {route} protein group in {path}")
            groups.append((members, probability))
    metrics = Metrics(dataset_id, "", path)
    metrics.add("protein_group_count", len(groups), group=route)
    metrics.distribution("proteins_per_group", [len(members) for members, _ in groups], group=route)
    metrics.distribution("protein_group_probability", [probability for _, probability in groups], group=route)
    metrics.write(output)
