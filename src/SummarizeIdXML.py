#!/usr/bin/env python3
"""Summarize raw or filtered DIA peptide identifications from idXML."""
import argparse
from collections import Counter, defaultdict, namedtuple
import pyopenms as oms
from QCMetrics import Metrics, as_text as text, finite

HitSummary = namedtuple(
    'HitSummary',
    (
        'scores',
        'errors',
        'absolute_errors',
        'charges',
        'missing_charges',
        'missing_errors',
        'score_types',
        'directions',
        'total_hits',
        'spectra',
        'sequences',
    ),
)


def label_of(hit):
    """Return the target/decoy class of a peptide hit."""
    if not hit.metaValueExists('target_decoy'):
        return 'unknown'
    value = text(hit.getMetaValue('target_decoy')).lower()
    return value if value in ('target', 'decoy') else 'unknown'


def selected_hits(hits, higher_is_better, filtered):
    """Keep all validated hits or the best finite-score raw hit."""
    if filtered:
        return hits
    valid = [(finite(hit.getScore()), hit) for hit in hits]
    valid = [(score, hit) for score, hit in valid if score is not None]
    if not valid:
        return []
    rank = max if higher_is_better else min
    best = rank(valid, key=lambda pair: pair[0])
    return [best[1]]


def precursor_mass_error(identification, hit, charge):
    """Calculate precursor mass error in ppm when m/z and charge are valid."""
    observed_mz = finite(identification.getMZ())
    theoretical_mz = (finite(hit.getSequence().getMZ(charge)) if charge > 0 else None)
    if observed_mz is None or theoretical_mz is None or theoretical_mz <= 0:
        return None
    return (observed_mz - theoretical_mz) / theoretical_mz * 1_000_000


def score_rows(identifications, filtered):
    """Collect selected hit scores and attributes by target/decoy class."""
    scores = defaultdict(list)
    signed_errors = defaultdict(list)
    absolute_errors = defaultdict(list)
    charges = defaultdict(Counter)
    missing_charges = Counter()
    missing_errors = Counter()
    score_types, directions = set(), set()
    total_hits, represented_spectra = 0, set()
    sequences = set()

    for identification in identifications:
        hits = list(identification.getHits())
        total_hits += len(hits)
        if not hits:
            continue
        score_types.add(text(identification.getScoreType()))
        higher = identification.isHigherScoreBetter()
        directions.add(higher)
        spectrum_ref = (text(identification.getMetaValue('spectrum_reference'))
                        if identification.metaValueExists('spectrum_reference') else '')
        if spectrum_ref:
            merge_index = (text(identification.getMetaValue('id_merge_index'))
                           if identification.metaValueExists('id_merge_index') else '')
            represented_spectra.add((merge_index, spectrum_ref))
        for hit in selected_hits(hits, higher, filtered):
            label = label_of(hit)
            score = finite(hit.getScore())
            sequence = text(hit.getSequence().toString())
            if sequence:
                sequences.add(sequence)
            for group in ('all', label):
                scores[group].append(score)
            charge = hit.getCharge()
            if charge > 0:
                for group in ('all', label):
                    charges[group][charge] += 1
            else:
                for group in ('all', label):
                    missing_charges[group] += 1
            error = precursor_mass_error(identification, hit, charge)
            if error is not None:
                for group in ('all', label):
                    signed_errors[group].append(error)
                    absolute_errors[group].append(abs(error))
            else:
                for group in ('all', label):
                    missing_errors[group] += 1
    return HitSummary(
        scores=scores,
        errors=signed_errors,
        absolute_errors=absolute_errors,
        charges=charges,
        missing_charges=missing_charges,
        missing_errors=missing_errors,
        score_types=score_types,
        directions=directions,
        total_hits=total_hits,
        spectra=represented_spectra,
        sequences=sequences,
    )


class IdXMLSummarizer:
    """Summarize raw search or final validated peptide identifications."""

    def __init__(self, input_path, dataset_id, sample_name, search_engine, stage):
        self.input_path = input_path
        self.filtered = stage == 'filtered'
        layer = 'validated_identification' if self.filtered else 'raw_identification'
        self.metrics = Metrics(
            dataset_id=dataset_id,
            layer=layer,
            source=input_path,
            sample_name=sample_name,
            search_engine=search_engine,
        )
        self.identifications = None
        self.hits = None

    def load(self):
        """Read idXML and collect selected hit values."""
        protein_ids = []
        self.identifications = oms.PeptideIdentificationList()
        oms.IdXMLFile().load(self.input_path, protein_ids, self.identifications)
        self.hits = score_rows(self.identifications, self.filtered)

    def summarize_counts(self):
        """Count identification records, hits, peptides, and unknown labels."""
        self.metrics.add('peptide_identification_count', len(self.identifications))
        self.metrics.add(
            name='accepted_psm_count' if self.filtered else 'candidate_psm_count',
            value=self.hits.total_hits,
        )
        self.metrics.add(
            name='distinct_peptide_count' if self.filtered else 'distinct_top_peptide_count',
            value=len(self.hits.sequences),
        )
        self.metrics.add('target_decoy_unknown_count', len(self.hits.scores['unknown']))
        if not self.filtered:
            self.metrics.add('top_hit_count', len(self.hits.scores['all']))
            if self.hits.spectra:
                self.metrics.add('spectrum_with_candidate_count', len(self.hits.spectra))
            else:
                self.metrics.unavailable('spectrum_with_candidate_count')

    def summarize_score_metadata(self):
        """Report the single score type and direction used in this idXML."""
        prefix = 'validated_score' if self.filtered else 'search_score'
        if len(self.hits.score_types) == 1:
            self.metrics.add(
                name=prefix + '_type',
                value=next(iter(self.hits.score_types)),
                unit='',
            )
        elif self.hits.score_types:
            raise ValueError('Mixed score types in one idXML; separate them before QC')
        if len(self.hits.directions) == 1:
            self.metrics.add(
                name=prefix + '_higher_is_better',
                value=int(next(iter(self.hits.directions))),
                unit='boolean',
            )
        elif self.hits.directions:
            raise ValueError('Mixed score directions in one idXML; separate them before QC')

    def summarize_hit_distributions(self):
        """Report scores and, for raw search, mass error and precursor charge."""
        prefix = 'validated_score' if self.filtered else 'search_score'
        for group in ('all', 'target', 'decoy', 'unknown'):
            self.metrics.distribution(
                name=prefix,
                values=self.hits.scores[group],
                group=group,
            )
            if not self.filtered:
                self.summarize_raw_group(group)

    def summarize_raw_group(self, group):
        """Report mass error and charge categories for one target/decoy group."""
        self.metrics.distribution(
            name='precursor_mass_error_signed',
            values=self.hits.errors[group],
            unit='ppm',
            group=group,
        )
        self.metrics.distribution(
            name='precursor_mass_error_absolute',
            values=self.hits.absolute_errors[group],
            unit='ppm',
            group=group,
        )
        self.metrics.add(
            name='precursor_mass_error_missing_count',
            value=self.hits.missing_errors[group],
            group=group,
        )
        self.metrics.add(
            name='precursor_charge_missing_count',
            value=self.hits.missing_charges[group],
            group=group,
        )
        total = sum(self.hits.charges[group].values())
        for charge, count in sorted(self.hits.charges[group].items()):
            self.metrics.add(
                name='precursor_charge_count',
                value=count,
                group=f'{group}:z{charge}',
            )
            self.metrics.add(
                name='precursor_charge_fraction',
                value=count / total,
                unit='fraction',
                group=f'{group}:z{charge}',
            )

    def run(self, output_path):
        """Read identifications and write their QC metrics."""
        self.load()
        self.summarize_counts()
        self.summarize_score_metadata()
        self.summarize_hit_distributions()
        self.metrics.write(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input', 'output', 'dataset-id', 'sample-name', 'search-engine'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--stage', choices=('search', 'filtered'), required=True)
    args = parser.parse_args()
    summarizer = IdXMLSummarizer(
        input_path=args.input,
        dataset_id=args.dataset_id,
        sample_name=args.sample_name,
        search_engine=args.search_engine,
        stage=args.stage,
    )
    summarizer.run(args.output)
