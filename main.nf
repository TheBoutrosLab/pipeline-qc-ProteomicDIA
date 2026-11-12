#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { summarize_mzML } from './module/summarize_mzML'
include { extract_MS1_matrix } from './module/extract_MS1_matrix'
include { plot_MS1_landscape } from './module/plot_MS1_landscape'
include { summarize_search_idXML } from './module/summarize_search_idXML'
include { summarize_filtered_idXML } from './module/summarize_filtered_idXML'
include { summarize_PSMFeatures } from './module/summarize_PSMFeatures'
include { summarize_ProteinGroups_ProteinProphet } from './module/summarize_ProteinGroups_ProteinProphet'
include { summarize_ProteinGroups_Fido } from './module/summarize_ProteinGroups_Fido'
include { summarize_ProteinGroups_Epifany } from './module/summarize_ProteinGroups_Epifany'
include { summarize_ProteinGroups_ProtLoom } from './module/summarize_ProteinGroups_ProtLoom'
include { summarize_EasyPQP_library } from './module/summarize_EasyPQP_library'
include { summarize_DIANN_report } from './module/summarize_DIANN_report'
include { collect_QCMetrics } from './module/collect_QCMetrics'

def resolveInput(def csv, String value) {
    if (!value?.trim()) error "Empty input path in ${csv}"
    def result = file(value.startsWith('/') ? value : "${csv.parent}/${value}")
    if (!result.isFile()) error "Missing input file: ${result}"
    return result
}
def resolveAcquisition(def csv, String value) {
    if (!value?.trim()) error "Empty acquisition path in ${csv}"
    def result = file(value.startsWith('/') ? value : "${csv.parent}/${value}")
    def name = result.name.toLowerCase()
    if (!((result.isFile() && name.endsWith('.mzml')) ||
          (result.isDirectory() && name.endsWith('.d')))) {
        error "Acquisition input must be an mzML file or Bruker .d directory: ${result}"
    }
    return result
}
workflow {
    mzmlRows = Channel.fromPath(params.mzml_csv).splitCsv(header: true).toList().map { [rows: it] }
    searchRows = Channel.fromPath(params.search_idxml_csv).splitCsv(header: true).toList().map { [rows: it] }
    filteredRows = Channel.fromPath(params.filtered_idxml_csv).splitCsv(header: true).toList().map { [rows: it] }
    inputs = mzmlRows.combine(searchRows).combine(filteredRows).map { mzPayload, rawPayload, finalPayload ->
        def mzRows = mzPayload.rows
        def rawRows = rawPayload.rows
        def finalRows = finalPayload.rows
        if (!mzRows || !rawRows || !finalRows) error 'All three sample manifests must contain rows'
        if (!(mzRows[0].keySet().containsAll(['sample_name', 'fraction_id', 'path'])))
            error 'mzML manifest needs sample_name,fraction_id,path columns'
        def safe = { value -> value ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/ }
        def mz = [:]
        mzRows.each { row ->
            def sample = row.sample_name?.trim()
            def fraction = row.fraction_id?.trim()
            if (!fraction)
                error "Missing fraction_id in ${params.mzml_csv}: ${row}"
            def key = [sample, fraction]
            if (!safe(sample) || !safe(fraction) || mz.containsKey(key))
                error "Invalid or duplicate sample_name/fraction_id in ${params.mzml_csv}: ${key}"
            def path = row.path?.trim()
            def reconstructed = resolveInput(file(params.mzml_csv), path)
            def acquisition = resolveAcquisition(file(params.mzml_csv),
                row.acquisition_path?.trim() ?: path)
            mz[key] = [reconstructed, acquisition]
        }
        def asMap = { rows, source ->
            def values = [:]
            rows.each { row ->
                def sample = row.sample_name?.trim()
                if (!safe(sample) || values.containsKey(sample))
                    error "Invalid or duplicate sample_name in ${source}: ${sample}"
                values[sample] = resolveInput(file(source), row.idXML?.trim())
            }
            values
        }
        def raw = asMap(rawRows, params.search_idxml_csv)
        def filtered = asMap(finalRows, params.filtered_idxml_csv)
        def names = mz.keySet().collect { it[0] } as Set
        if (names != raw.keySet() || names != filtered.keySet())
            error 'mzML, search idXML, and filtered idXML manifests must contain the same sample names'
        rawRows.each { row ->
            if (row.search_engine && row.search_engine != 'MSFragger')
                error "DIA QC expects MSFragger search idXML: ${row.search_engine}"
        }
        [acquisitions: mz.keySet().sort().collect { key ->
            tuple(key[0], key[1], mz[key][1])
        }, samples: names.sort().collect { sample ->
            tuple(sample, raw[sample], filtered[sample])
        }, pairs: mz.keySet()]
    }
    acquisitions = inputs.flatMap { it.acquisitions }
    sampleIds = inputs.flatMap { it.samples }
    summarize_mzML(acquisitions)
    extract_MS1_matrix(acquisitions)
    plot_MS1_landscape(extract_MS1_matrix.out.matrix)
    summarize_search_idXML(sampleIds.map { sample, raw, filtered -> tuple(sample, raw) })
    summarize_filtered_idXML(sampleIds.map { sample, raw, filtered -> tuple(sample, filtered) })
    if (params.psm_features_csv) {
        featureRows = Channel.fromPath(params.psm_features_csv).splitCsv(header: true)
            .toList().map { [rows: it] }
        featureSamples = inputs.map { [names: it.samples.collect { sample, raw, filtered -> sample },
                                          pairs: it.pairs] }
        featureEntries = featureRows.combine(featureSamples).flatMap { rowPayload, samplePayload ->
            def rows = rowPayload.rows
            def sampleNames = samplePayload.names
            def pairs = samplePayload.pairs
            def safe = { value -> value ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/ }
            def seen = [] as Set
            rows.collect { row ->
                def sample = row.sample_name?.trim()
                def engine = row.search_engine?.trim()
                def route = row.route?.trim()
                def fraction = row.fraction_id?.trim() ?: ''
                if (!safe(sample) || !(sample in sampleNames) || engine != 'MSFragger' ||
                    !(route in ['MSBooster', 'PSMFeatureExtractor']) ||
                    (route == 'MSBooster' && (!fraction || !pairs.contains([sample, fraction]))) ||
                    (fraction && (!safe(fraction) || !pairs.contains([sample, fraction]))) ||
                    !seen.add([sample, engine, route, fraction]))
                    error "Invalid or duplicate PSM feature manifest row: ${row}"
                def input = resolveInput(file(params.psm_features_csv), row.path?.trim())
                if ((route == 'MSBooster' && !input.name.endsWith('.pin')) ||
                    (route == 'PSMFeatureExtractor' && !input.name.endsWith('.idXML')))
                    error "Wrong PSM feature file type for ${route}: ${input}"
                tuple(sample, engine, route, fraction, input)
            }
        }
        summarize_PSMFeatures(featureEntries)
        featureMetrics = summarize_PSMFeatures.out.metrics
    } else {
        featureMetrics = Channel.empty()
    }
    if (params.inference_route == 'ProteinProphet') {
        summarize_ProteinGroups_ProteinProphet(file(params.protein_groups))
        groupMetrics = summarize_ProteinGroups_ProteinProphet.out.metrics
    } else if (params.inference_route == 'Fido') {
        summarize_ProteinGroups_Fido(file(params.protein_groups))
        groupMetrics = summarize_ProteinGroups_Fido.out.metrics
    } else if (params.inference_route == 'Epifany') {
        summarize_ProteinGroups_Epifany(file(params.protein_groups))
        groupMetrics = summarize_ProteinGroups_Epifany.out.metrics
    } else {
        protloomGroups = Channel.fromPath(params.protein_groups).splitCsv(header: true).toList().map { rows ->
            if (!rows || rows.any { !it.sample_name || !it.groups_tsv })
                error 'ProtLoom groups CSV needs sample_name,groups_tsv rows'
            if ((rows*.sample_name as Set).size() != rows.size())
                error 'ProtLoom groups CSV has duplicate sample names'
            rows.collect { resolveInput(file(params.protein_groups), it.groups_tsv) }
        }
        summarize_ProteinGroups_ProtLoom(protloomGroups)
        groupMetrics = summarize_ProteinGroups_ProtLoom.out.metrics
    }
    summarize_EasyPQP_library(file(params.library_tsv), file(params.export_summary_tsv))
    summarize_DIANN_report(file(params.diann_report), file(params.mzml_csv))
    allMetrics = summarize_mzML.out.metrics
        .mix(summarize_search_idXML.out.metrics)
        .mix(summarize_filtered_idXML.out.metrics)
        .mix(featureMetrics)
        .mix(groupMetrics)
        .mix(summarize_EasyPQP_library.out.metrics)
        .mix(summarize_DIANN_report.out.metrics)
        .collect()
    collect_QCMetrics(allMetrics)
}
