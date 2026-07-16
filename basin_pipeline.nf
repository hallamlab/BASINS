#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def projectRootDir = new File(projectDir.toString())
def defaultConfigPath = "${projectDir}/basin_pipeline_nextflow.yml"
def paramsMap = [:]

try {
    paramsMap = workflow.params ? new LinkedHashMap(workflow.params) : [:]
} catch( Throwable ignored ) {
    paramsMap = [:]
}

def inlineKeys = ['paths', 'resources', 'biochem', 'biochem_pre_asv', 'environments', 'config_root', 'pipeline_config']
def hasInlineConfig = inlineKeys.any { paramsMap.containsKey(it) }
def config
File configFile = null
File configRoot = projectRootDir

def providedConfigPath = System.getenv('BASIN_PIPELINE_CONFIG') ?: paramsMap.pipeline_config ?: paramsMap.config
if( providedConfigPath ) {
    def configPath = file(providedConfigPath)
    configFile = configPath.toFile()
    if( !configFile.exists() ) {
        exit 1, "Config file not found: ${configFile}"
    }
    config = new groovy.yaml.YamlSlurper().parse(configFile)
    configRoot = configFile.parentFile ?: projectRootDir
    log.info "Loaded config from ${configFile}"
} else if( hasInlineConfig ) {
    config = paramsMap
    if( paramsMap.containsKey('config_root') ) {
        configRoot = file(paramsMap.config_root).toFile()
    }
    log.info "Using inline Nextflow parameters as configuration."
} else {
    configFile = file(defaultConfigPath).toFile()
    if( !configFile.exists() ) {
        exit 1, "Config file not found: ${defaultConfigPath}"
    }
    config = new groovy.yaml.YamlSlurper().parse(configFile)
    configRoot = configFile.parentFile ?: projectRootDir
    log.info "Loaded default config from ${configFile}"
}

def resolvePath(String pathValue) {
    if( !pathValue ) {
        return null
    }
    def candidate = new File(pathValue)
    if( candidate.isAbsolute() ) {
        return candidate.canonicalPath
    }
    return new File(configRoot, pathValue).canonicalPath
}

def resolveOptionalPath = { value, File root ->
    if( value == null ) {
        return null
    }
    def text = value.toString().trim()
    if( !text ) {
        return null
    }
    def candidate = new File(text)
    if( candidate.isAbsolute() ) {
        return candidate.canonicalPath
    }
    return new File(root, text).canonicalPath
}

def resolveOutputRelative = { String pathValue, String outputRoot ->
    if( !pathValue ) {
        return null
    }
    def candidate = new File(pathValue)
    if( candidate.isAbsolute() ) {
        return candidate.canonicalPath
    }
    return new File(outputRoot, pathValue).canonicalPath
}

def listFromConfig = { raw ->
    if( raw instanceof List ) {
        return raw.collect { it.toString().trim() }.findAll { it }
    }
    if( raw ) {
        return raw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    return []
}

def outputDir = resolvePath(config.paths?.output_dir)
assert outputDir : "paths.output_dir must be provided in the YAML config"

def workDirOverride = config.paths?.work_dir ? resolvePath(config.paths.work_dir) : null
def condaCacheOverride = config.paths?.conda_cache_dir ? resolvePath(config.paths.conda_cache_dir) : null
if( workDirOverride ) {
    def workDirFile = new File(workDirOverride)
    workDirFile.mkdirs()
    workflow.workDir = java.nio.file.Paths.get(workDirFile.canonicalPath)
    log.info "Using custom Nextflow work directory: ${workflow.workDir}"
}
def resolvedCondaCacheDir = condaCacheOverride ?: new File(outputDir, ".conda_cache").canonicalPath
def condaCacheDirFile = new File(resolvedCondaCacheDir)
condaCacheDirFile.mkdirs()
System.setProperty('NXF_CONDA_CACHEDIR', condaCacheDirFile.canonicalPath)
log.info "Using Conda cache directory: ${condaCacheDirFile.canonicalPath}"

int hostThreads = Runtime.runtime.availableProcessors()
int pipelineThreads = config.resources?.threads ? (config.resources.threads as int) : hostThreads

def biochemConfig = config.biochem ?: (config.biochem_pre_asv ?: [:])
boolean biochemEnabled = biochemConfig.containsKey('enabled') ? (biochemConfig.enabled as boolean) : true
def biochemTableAPath = biochemConfig.table_a ? resolveOptionalPath(biochemConfig.table_a, configRoot) : null
def biochemTableBPath = biochemConfig.table_b ? resolveOptionalPath(biochemConfig.table_b, configRoot) : null
if( biochemEnabled && (!biochemTableAPath || !new File(biochemTableAPath).exists()) ) {
    exit 1, "biochem.table_a not found: ${biochemTableAPath}"
}
if( biochemEnabled && (!biochemTableBPath || !new File(biochemTableBPath).exists()) ) {
    exit 1, "biochem.table_b not found: ${biochemTableBPath}"
}

def biochemOutputRoot = biochemConfig.output_root ? resolveOutputRelative(biochemConfig.output_root.toString(), outputDir) : outputDir
def biochemProcessingDirAbs = new File(biochemOutputRoot, 'biochem_processing').canonicalPath
def biochemStratMetricsDirAbs = new File(biochemProcessingDirAbs, 'stratification_metrics').canonicalPath
def biochemPcaDirAbs = new File(biochemOutputRoot, 'env_pca').canonicalPath
def biochemSelectkDirAbs = new File(biochemOutputRoot, 'env_compartments_selectk').canonicalPath
def biochemGmmDirAbs = new File(biochemOutputRoot, 'env_compartments_gmm').canonicalPath
def biochemO2DirAbs = new File(biochemOutputRoot, 'env_o2_soft_compartments').canonicalPath
def biochemHybridDirAbs = new File(biochemOutputRoot, 'env_hybrid_soft_compartments').canonicalPath
def biochemCompareDirAbs = new File(biochemOutputRoot, 'env_compare_compartments').canonicalPath
def biochemSplitDirAbs = new File(biochemOutputRoot, 'env_o2_split_by_gmm').canonicalPath
def biochemStratIndexDirAbs = new File(biochemOutputRoot, 'env_stratification_index').canonicalPath
def biochemStateTransitionsDirAbs = new File(biochemOutputRoot, 'env_state_transitions').canonicalPath
def biochemSuccessionDirAbs = new File(biochemOutputRoot, 'env_succession_graphs').canonicalPath
def biochemFeatureAssocDirAbs = new File(biochemOutputRoot, 'env_compartment_feature_assoc').canonicalPath
def biochemEofPcaDirAbs = new File(biochemOutputRoot, 'eof_pca').canonicalPath
def biochemEofStatesDirAbs = new File(biochemOutputRoot, 'eof_states').canonicalPath
def biochemEofPlotsDirAbs = new File(biochemOutputRoot, 'eof_plots').canonicalPath
def biochemWithinGmmDirAbs = new File(biochemGmmDirAbs, 'within_gmm_hdbscan').canonicalPath
def biochemMergedOxygenPath = new File(biochemProcessingDirAbs, '02_oxygen_best_available.tsv').canonicalPath
def biochemDensityPath = new File(biochemProcessingDirAbs, '02_oxygen_best_available_density.tsv').canonicalPath
def biochemDensityCleanedFile = biochemConfig.cleaned_density_filename ?: '02_oxygen_best_available_density_RJM.tsv'
def biochemDensityCleanedPath = new File(biochemProcessingDirAbs, biochemDensityCleanedFile.toString()).canonicalPath
def biochemFeatureCols = biochemConfig.feature_cols ?: 'Oxygen,Nitrate,Nitrite,Nitrous Oxide,Ammonium,Hydrogen Sulfide,Methane,Phosphate,Silicate,Temperature,Salinity,Density,Fe,Dimethyl Sulfide'
def biochemGmmKRaw = biochemConfig.gmm_k
boolean biochemGmmKAuto = (biochemGmmKRaw == null) || (biochemGmmKRaw.toString().trim().equalsIgnoreCase('auto'))
def biochemGmmK = biochemGmmKAuto ? 5 : (biochemGmmKRaw as int)
def biochemEofPcs = biochemConfig.eof_pcs ?: '1,2,4'
List<String> biochemCleanKeepCols = listFromConfig(biochemConfig.clean_keep_cols)
List<String> biochemCleanDropCols = listFromConfig(biochemConfig.clean_drop_cols)
Map<String, String> biochemCleanRenameMap = [:]
if( biochemConfig.clean_rename_map instanceof Map ) {
    biochemConfig.clean_rename_map.each { k, v ->
        if( k != null && v != null ) {
            def kk = k.toString().trim()
            def vv = v.toString().trim()
            if( kk && vv ) {
                biochemCleanRenameMap[kk] = vv
            }
        }
    }
}

def defaultBiochemEnvPath = new File("${projectDir}/processes/shared_envs/biochem.yml").canonicalPath
def resolveBiochemStepEnv = { String envKey ->
    def stepEnvCfg = config.environments ? (config.environments[envKey] ?: config.environments.biochem) : null
    def stepEnvPath = stepEnvCfg ? resolveOptionalPath(stepEnvCfg, configRoot) : defaultBiochemEnvPath
    def stepEnvFile = file(stepEnvPath)
    if( !stepEnvFile.exists() ) {
        exit 1, "Biochem step conda environment YAML not found for ${envKey}: ${stepEnvPath}"
    }
    log.info "Using biochem Conda/Mamba env for ${envKey}: ${stepEnvPath}"
    return stepEnvPath
}

def biochemMergeCondaEnvPath = resolveBiochemStepEnv('biochem_merge')
def biochemDensityCondaEnvPath = resolveBiochemStepEnv('biochem_density')
def biochemStratMetricsCondaEnvPath = resolveBiochemStepEnv('biochem_strat_metrics')
def biochemCustomCleanCondaEnvPath = resolveBiochemStepEnv('biochem_custom_clean')
def biochemEigenvectorsCondaEnvPath = resolveBiochemStepEnv('biochem_eigenvectors')
def biochemSelectkCondaEnvPath = resolveBiochemStepEnv('biochem_selectk')
def biochemGmmCondaEnvPath = resolveBiochemStepEnv('biochem_gmm')
def biochemO2SoftCondaEnvPath = resolveBiochemStepEnv('biochem_o2_soft')
def biochemHybridCondaEnvPath = resolveBiochemStepEnv('biochem_hybrid')
def biochemCompareCondaEnvPath = resolveBiochemStepEnv('biochem_compare')
def biochemSplitCondaEnvPath = resolveBiochemStepEnv('biochem_split_o2_by_gmm')
def biochemStratAnomalyCondaEnvPath = resolveBiochemStepEnv('biochem_strat_anomaly')
def biochemStateTransitionsCondaEnvPath = resolveBiochemStepEnv('biochem_state_transitions')
def biochemSuccessionCondaEnvPath = resolveBiochemStepEnv('biochem_succession')
def biochemFeatureAssocCondaEnvPath = resolveBiochemStepEnv('biochem_feature_assoc')
def biochemEofPipelineCondaEnvPath = resolveBiochemStepEnv('biochem_eof_pipeline')
def biochemEofStateCondaEnvPath = resolveBiochemStepEnv('biochem_eof_state_cluster')
def biochemEofModeCondaEnvPath = resolveBiochemStepEnv('biochem_eof_mode_plots')
def biochemWithinGmmCondaEnvPath = resolveBiochemStepEnv('biochem_within_gmm')

def scriptPath = { String relativePath ->
    def scriptFile = new File("${projectDir}/${relativePath}")
    if( !scriptFile.exists() ) {
        exit 1, "${relativePath} not found in project directory"
    }
    return scriptFile.canonicalPath
}

def biochemMergeTablesScriptPath = scriptPath('processes/merge_tables/merge_tables_ctd_nearest_depth.py')
def biochemCalcDensityScriptPath = scriptPath('processes/density/env_calc_density.py')
def biochemStratMetricsScriptPath = scriptPath('processes/stratification_metrics/env_stratification_metrics.py')
def biochemCustomCleanerScriptPath = scriptPath('processes/custom_clean/custom_density_cleaner.py')
def biochemEigenvectorsScriptPath = scriptPath('processes/eigenvectors/env_eigenvectors.py')
def biochemSelectkScriptPath = scriptPath('processes/selectk/env_compartments_selectk.py')
def biochemGmmScriptPath = scriptPath('processes/gmm/env_compartments_gmm.py')
def biochemO2SoftScriptPath = scriptPath('processes/o2_soft/env_compartments_o2_soft.py')
def biochemHybridScriptPath = scriptPath('processes/hybrid/env_hybrid_compartment_builder.py')
def biochemCompareScriptPath = scriptPath('processes/compare_compartments/env_compare_compartments.py')
def biochemSplitScriptPath = scriptPath('processes/split_o2_by_gmm/env_split_o2_by_gmm.py')
def biochemStratAnomalyScriptPath = scriptPath('processes/stratification_anomaly/env_stratification_anomaly_detection.py')
def biochemStateTransitionScriptPath = scriptPath('processes/state_transitions/env_state_transition_analysis.py')
def biochemSuccessionScriptPath = scriptPath('processes/succession_graph/env_succession_graph.py')
def biochemFeatureAssocScriptPath = scriptPath('processes/feature_assoc/env_compartment_feature_assoc.py')
def biochemEofPipelineScriptPath = scriptPath('processes/eof_pipeline/env_eof_pipeline.py')
def biochemEofStateScriptPath = scriptPath('processes/eof_state_cluster/eof_state_clustering.py')
def biochemEofModePlotScriptPath = scriptPath('processes/eof_mode_plots/eof_mode_plots.py')
def biochemWithinGmmScriptPath = scriptPath('processes/within_gmm_hdbscan/env_within_gmm_hdbscan.py')

workflow {
    if( !biochemEnabled ) {
        log.info "biochem.enabled is false; no BASIN stages will run."
        return
    }

    b01 = BIOCHEM_MERGE()
    b02 = BIOCHEM_DENSITY(b01.done)
    b03 = BIOCHEM_STRAT_METRICS(b02.done)
    b04 = BIOCHEM_CUSTOM_CLEAN(b03.done)
    b05 = BIOCHEM_EIGENVECTORS(b04.done)
    b06 = BIOCHEM_SELECTK(b05.done)
    b07 = BIOCHEM_GMM(b06.done, b06.selected_k)
    b08 = BIOCHEM_O2_SOFT(b07.done)
    b09 = BIOCHEM_HYBRID(b08.done)
    b10 = BIOCHEM_COMPARE(b09.done)
    b11 = BIOCHEM_SPLIT_O2_BY_GMM(b10.done)
    b12 = BIOCHEM_STRAT_ANOMALY(b11.done)
    b13 = BIOCHEM_STATE_TRANSITIONS(b12.done)
    b14 = BIOCHEM_SUCCESSION_GRAPH(b13.done)
    b15 = BIOCHEM_FEATURE_ASSOC(b14.done)
    b16 = BIOCHEM_EOF_PIPELINE(b15.done)
    b17 = BIOCHEM_EOF_STATE_CLUSTER(b16.done)
    b18 = BIOCHEM_EOF_MODE_PLOTS(b17.done)
    BIOCHEM_WITHIN_GMM_HDBSCAN(b18.done)
}

process BIOCHEM_MERGE {
    cpus pipelineThreads
    conda "${biochemMergeCondaEnvPath}"

    output:
    path("biochem_merge.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemProcessingDirAbs}"
python "${biochemMergeTablesScriptPath}" \\
  --table-a "${biochemTableAPath}" \\
  --table-b "${biochemTableBPath}" \\
  --outdir "${biochemProcessingDirAbs}"
[[ -f "${biochemMergedOxygenPath}" ]] || { echo "Missing ${biochemMergedOxygenPath}" >&2; exit 1; }
touch biochem_merge.done
"""
}

process BIOCHEM_DENSITY {
    cpus pipelineThreads
    conda "${biochemDensityCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_density.done"), emit: done

    script:
    """
set -euo pipefail
python "${biochemCalcDensityScriptPath}" \\
  --input "${biochemMergedOxygenPath}" \\
  --salinity-col Salinity \\
  --temperature-col Temperature \\
  --depth-col Depth \\
  --latitude-col Latitude \\
  --longitude-col Longitude \\
  --sigma0
[[ -f "${biochemDensityPath}" ]] || { echo "Missing ${biochemDensityPath}" >&2; exit 1; }
touch biochem_density.done
"""
}

process BIOCHEM_STRAT_METRICS {
    cpus pipelineThreads
    conda "${biochemStratMetricsCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_strat_metrics.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemStratMetricsDirAbs}"
python "${biochemStratMetricsScriptPath}" \\
  --input "${biochemDensityPath}" \\
  --output-dir "${biochemStratMetricsDirAbs}" \\
  --salinity-col Salinity \\
  --temperature-col Temperature \\
  --depth-col Depth \\
  --latitude-col Latitude \\
  --longitude-col Longitude \\
  --profile-cols Cruise \\
  --date-col Date \\
  --layer-split-mode mld125
[[ -f "${biochemStratMetricsDirAbs}/stratification_summary.tsv" ]] || { echo "Missing ${biochemStratMetricsDirAbs}/stratification_summary.tsv" >&2; exit 1; }
touch biochem_strat_metrics.done
"""
}

process BIOCHEM_CUSTOM_CLEAN {
    cpus pipelineThreads
    conda "${biochemCustomCleanCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_custom_clean.done"), emit: done
    path("biochem_density_cleaned.tsv"), emit: cleaned_density

    script:
    def biochemRenamePairs = biochemCleanRenameMap ? biochemCleanRenameMap.collect { k, v -> "${k}:${v}" }.join(',') : ''
    def biochemKeepArg = biochemCleanKeepCols && !biochemCleanKeepCols.isEmpty() ? """ --keep-cols "${biochemCleanKeepCols.join(',')}" """ : ''
    def biochemDropArg = biochemCleanDropCols && !biochemCleanDropCols.isEmpty() ? """ --drop-cols "${biochemCleanDropCols.join(',')}" """ : ''
    def biochemRenameArg = biochemRenamePairs ? """ --rename-map "${biochemRenamePairs}" """ : ''
    """
set -euo pipefail
python "${biochemCustomCleanerScriptPath}" --input "${biochemDensityPath}" --output "${biochemDensityCleanedPath}"${biochemKeepArg}${biochemDropArg}${biochemRenameArg}
[[ -f "${biochemDensityCleanedPath}" ]] || { echo "Missing ${biochemDensityCleanedPath}" >&2; exit 1; }
ln -sf "${biochemDensityCleanedPath}" biochem_density_cleaned.tsv
touch biochem_custom_clean.done
"""
}

process BIOCHEM_EIGENVECTORS {
    cpus pipelineThreads
    conda "${biochemEigenvectorsCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_eigenvectors.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemPcaDirAbs}"
python "${biochemEigenvectorsScriptPath}" \\
  --input "${biochemDensityCleanedPath}" \\
  --outdir "${biochemPcaDirAbs}" \\
  --feature-cols "${biochemFeatureCols}" \\
  --pc-selection \\
  --anchor-depths
[[ -f "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" ]] || { echo "Missing ${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" >&2; exit 1; }
[[ -f "${biochemPcaDirAbs}/tables/pc_keep_decision.csv" ]] || { echo "Missing ${biochemPcaDirAbs}/tables/pc_keep_decision.csv" >&2; exit 1; }
touch biochem_eigenvectors.done
"""
}

process BIOCHEM_SELECTK {
    cpus pipelineThreads
    conda "${biochemSelectkCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_selectk.done"), emit: done
    path("selected_k.txt"), emit: selected_k

    script:
    """
set -euo pipefail
mkdir -p "${biochemSelectkDirAbs}"
python "${biochemSelectkScriptPath}" \\
  --eigenvectors "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" \\
  --pc-keep "${biochemPcaDirAbs}/tables/pc_keep_decision.csv" \\
  --outdir "${biochemSelectkDirAbs}" \\
  --sep "," \\
  --stability-block-col Cruise \\
  --min-cluster-frac 0.02
[[ -s "${biochemSelectkDirAbs}/SELECTED_K.txt" ]] || { echo "Missing ${biochemSelectkDirAbs}/SELECTED_K.txt" >&2; exit 1; }
ln -sf "${biochemSelectkDirAbs}/SELECTED_K.txt" selected_k.txt
touch biochem_selectk.done
"""
}

process BIOCHEM_GMM {
    cpus pipelineThreads
    conda "${biochemGmmCondaEnvPath}"

    input:
    path(prev_done)
    path(selected_k_file)

    output:
    path("biochem_gmm.done"), emit: done

    script:
    def biochemGmmKAutoFlag = biochemGmmKAuto ? '1' : '0'
    """
set -euo pipefail
mkdir -p "${biochemGmmDirAbs}"
gmm_k="${biochemGmmK}"
if [[ "${biochemGmmKAutoFlag}" == "1" ]]; then
  gmm_k="\$(tr -d '[:space:]' < "${selected_k_file}")"
  case "\$gmm_k" in
    ''|*[!0-9]*) echo "Invalid selected K: \$gmm_k" >&2; exit 1 ;;
  esac
fi
python "${biochemGmmScriptPath}" \\
  --eigenvectors "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" \\
  --pc-keep "${biochemPcaDirAbs}/tables/pc_keep_decision.csv" \\
  --outdir "${biochemGmmDirAbs}" \\
  --sep "," \\
  --pc-use-mode keep \\
  --standardize-pc-space \\
  --episodic-smoothing \\
  --random-state 42 \\
  --K "\$gmm_k" \\
  --matrix-cleaned "${biochemPcaDirAbs}/tables/matrix_cleaned_with_sparse.csv"
[[ -f "${biochemGmmDirAbs}/tables/compartments_assignments_smoothed.csv" ]] || { echo "Missing ${biochemGmmDirAbs}/tables/compartments_assignments_smoothed.csv" >&2; exit 1; }
touch biochem_gmm.done
"""
}

process BIOCHEM_O2_SOFT {
    cpus pipelineThreads
    conda "${biochemO2SoftCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_o2_soft.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemO2DirAbs}"
python "${biochemO2SoftScriptPath}" \\
  --input "${biochemPcaDirAbs}/tables/matrix_cleaned.csv" \\
  --outdir "${biochemO2DirAbs}" \\
  --o2-col Oxygen \\
  --T-oxic-dyso 90 \\
  --T-dyso-sub 20 \\
  --T-sub-anox 1 \\
  --episodic-smoothing \\
  --episodic-block-col Cruise \\
  --episodic-sort-cols Depth_anchored \\
  --episodic-sticky-prob 0.85 \\
  --episodic-apply-to all
[[ -f "${biochemO2DirAbs}/tables/o2_compartments_assignments_smoothed.csv" ]] || { echo "Missing ${biochemO2DirAbs}/tables/o2_compartments_assignments_smoothed.csv" >&2; exit 1; }
touch biochem_o2_soft.done
"""
}

process BIOCHEM_HYBRID {
    cpus pipelineThreads
    conda "${biochemHybridCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_hybrid.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemHybridDirAbs}"
python "${biochemHybridScriptPath}" \\
  --gmm-assign "${biochemGmmDirAbs}/tables/compartments_assignments_smoothed.csv" \\
  --o2-assign "${biochemO2DirAbs}/tables/o2_compartments_assignments_smoothed.csv" \\
  --outdir "${biochemHybridDirAbs}" \\
  --join-key cruise_year_month_depth \\
  --make-plot \\
  --cmap rainbow \\
  --make-umap \\
  --umap-n-neighbors 15 \\
  --umap-min-dist 0.1 \\
  --umap-seed 42
[[ -f "${biochemHybridDirAbs}/tables/cruise_composition_hybrid.csv" ]] || { echo "Missing ${biochemHybridDirAbs}/tables/cruise_composition_hybrid.csv" >&2; exit 1; }
touch biochem_hybrid.done
"""
}

process BIOCHEM_COMPARE {
    cpus pipelineThreads
    conda "${biochemCompareCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_compare.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemCompareDirAbs}"
python "${biochemCompareScriptPath}" \\
  --matrix-cleaned "${biochemPcaDirAbs}/tables/matrix_cleaned_with_sparse.csv" \\
  --eigenvectors "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" \\
  --assignments "${biochemGmmDirAbs}/tables/compartments_assignments_smoothed.csv" \\
  --o2-assignments "${biochemO2DirAbs}/tables/o2_compartments_assignments_smoothed.csv" \\
  --o2-compartment-col compartment_name \\
  --outdir "${biochemCompareDirAbs}" \\
  --sep-matrix "," \\
  --sep-eig "," \\
  --sep-assign "," \\
  --sep-o2-assign "," \\
  --pca-tables-dir "${biochemPcaDirAbs}/tables" \\
  --key-mode composite \\
  --key-cols "Cruise,Year,Month,Day,Depth" \\
  --pc-cols "PC1,PC2"
[[ -f "${biochemCompareDirAbs}/tables/umap_embedding.csv" ]] || { echo "Missing ${biochemCompareDirAbs}/tables/umap_embedding.csv" >&2; exit 1; }
touch biochem_compare.done
"""
}

process BIOCHEM_SPLIT_O2_BY_GMM {
    cpus pipelineThreads
    conda "${biochemSplitCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_split.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemSplitDirAbs}"
python "${biochemSplitScriptPath}" \\
  --matrix-cleaned "${biochemPcaDirAbs}/tables/matrix_cleaned_with_sparse.csv" \\
  --eigenvectors "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" \\
  --assignments "${biochemGmmDirAbs}/tables/compartments_assignments_smoothed.csv" \\
  --o2-assignments "${biochemO2DirAbs}/tables/o2_compartments_assignments_smoothed.csv" \\
  --o2-compartment-col compartment_name \\
  --outdir "${biochemSplitDirAbs}" \\
  --sep-matrix "," \\
  --sep-eig "," \\
  --sep-assign "," \\
  --sep-o2-assign "," \\
  --key-mode composite \\
  --key-cols "Cruise,Year,Month,Day,Depth" \\
  --pc-cols "PC1,PC2" \\
  --plots \\
  --plot-formats "pdf,png,svg" \\
  --umap-embedding "${biochemCompareDirAbs}/tables/umap_embedding.csv" \\
  --reassign \\
  --borderline-mode other_or_low_conf \\
  --borderline-max-prob 0.70 \\
  --core-min-prob 0.90 \\
  --reassign-radius-quantile 0.95 \\
  --reassign-min-core-n 30 \\
  --min-subcluster-size 20
[[ -f "${biochemSplitDirAbs}/tables/merged_o2_split_by_gmm.csv" ]] || { echo "Missing ${biochemSplitDirAbs}/tables/merged_o2_split_by_gmm.csv" >&2; exit 1; }
touch biochem_split.done
"""
}

process BIOCHEM_STRAT_ANOMALY {
    cpus pipelineThreads
    conda "${biochemStratAnomalyCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_strat_anomaly.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemStratIndexDirAbs}"
python "${biochemStratAnomalyScriptPath}" \\
  --input "${biochemPcaDirAbs}/tables/matrix_cleaned.csv" \\
  --sample-id-col cruise_year_month_depth \\
  --date-col date \\
  --month-col Month \\
  --year-col Year \\
  --depth-col Depth \\
  --output-dir "${biochemStratIndexDirAbs}" \\
  --pea-metrics "${biochemStratMetricsDirAbs}/stratification_summary.tsv"
[[ -f "${biochemStratIndexDirAbs}/stratification_timeseries.tsv" ]] || { echo "Missing ${biochemStratIndexDirAbs}/stratification_timeseries.tsv" >&2; exit 1; }
touch biochem_strat_anomaly.done
"""
}

process BIOCHEM_STATE_TRANSITIONS {
    cpus pipelineThreads
    conda "${biochemStateTransitionsCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_state_transitions.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemStateTransitionsDirAbs}"
python "${biochemStateTransitionScriptPath}" \\
  --o2 "${biochemHybridDirAbs}/tables/cruise_composition_o2.csv" \\
  --gmm "${biochemHybridDirAbs}/tables/cruise_composition_gmm.csv" \\
  --hybrid "${biochemHybridDirAbs}/tables/cruise_composition_hybrid.csv" \\
  --outdir "${biochemStateTransitionsDirAbs}" \\
  --changepoint-metric braycurtis \\
  --changepoint-threshold 0.35 \\
  --coupling \\
  --coupling-method spearman \\
  --coupling-cluster-threshold 0.30 \\
  --coupling-edge-threshold 0.60 \\
  --strat-timeseries "${biochemStratIndexDirAbs}/stratification_timeseries.tsv" \\
  --eof-states "${biochemStratIndexDirAbs}/stratification_timeseries.tsv" \\
  --eof-state-col anomaly_type
touch biochem_state_transitions.done
"""
}

process BIOCHEM_SUCCESSION_GRAPH {
    cpus pipelineThreads
    conda "${biochemSuccessionCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_succession.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemSuccessionDirAbs}"
python "${biochemSuccessionScriptPath}" \\
  --o2 "${biochemHybridDirAbs}/tables/cruise_composition_o2.csv" \\
  --gmm "${biochemHybridDirAbs}/tables/cruise_composition_gmm.csv" \\
  --hybrid "${biochemHybridDirAbs}/tables/cruise_composition_hybrid.csv" \\
  --outdir "${biochemSuccessionDirAbs}" \\
  --make-plots \\
  --top-n 1 \\
  --no-keep-self \\
  --min-prob 0.1
touch biochem_succession.done
"""
}

process BIOCHEM_FEATURE_ASSOC {
    cpus pipelineThreads
    conda "${biochemFeatureAssocCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_feature_assoc.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemFeatureAssocDirAbs}"
python "${biochemFeatureAssocScriptPath}" \\
  --matrix-cleaned "${biochemPcaDirAbs}/tables/matrix_cleaned_with_sparse.csv" \\
  --assignments-gmm "${biochemGmmDirAbs}/tables/compartments_assignments_smoothed.csv" \\
  --assignments-o2 "${biochemO2DirAbs}/tables/o2_compartments_assignments_smoothed.csv" \\
  --assignments-hybrid "${biochemHybridDirAbs}/tables/compartments_assignments_hybrid.csv" \\
  --outdir "${biochemFeatureAssocDirAbs}" \\
  --sep-matrix "," \\
  --sep-assign "," \\
  --bootstrap-B 500 \\
  --top-n-each-side 8 \\
  --min-n-comp 20 \\
  --min-n-rest 50 \\
  --depth-adjust \\
  --hybrid-split-table "${biochemSplitDirAbs}/tables/merged_o2_split_by_gmm.csv"
touch biochem_feature_assoc.done
"""
}

process BIOCHEM_EOF_PIPELINE {
    cpus pipelineThreads
    conda "${biochemEofPipelineCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_eof_pipeline.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemEofPcaDirAbs}"
python "${biochemEofPipelineScriptPath}" \\
  --matrix-cleaned "${biochemPcaDirAbs}/tables/matrix_cleaned_with_sparse.csv" \\
  --core-loadings "${biochemPcaDirAbs}/tables/pca_loadings.csv" \\
  --outdir "${biochemEofPcaDirAbs}" \\
  --sep "," \\
  --pc-selection
[[ -f "${biochemEofPcaDirAbs}/tables/eof_eigenvectors_scores_by_cruise.csv" ]] || { echo "Missing ${biochemEofPcaDirAbs}/tables/eof_eigenvectors_scores_by_cruise.csv" >&2; exit 1; }
touch biochem_eof_pipeline.done
"""
}

process BIOCHEM_EOF_STATE_CLUSTER {
    cpus pipelineThreads
    conda "${biochemEofStateCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_eof_state_cluster.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemEofStatesDirAbs}"
python "${biochemEofStateScriptPath}" \\
  --scores "${biochemEofPcaDirAbs}/tables/eof_eigenvectors_scores_by_cruise.csv" \\
  --pcs "${biochemEofPcs}" \\
  --k auto \\
  --k-min 2 \\
  --k-max 10 \\
  --covariance-type full \\
  --standardize-pc-space \\
  --n-init 20 \\
  --max-iter 500 \\
  --cv-folds 5 \\
  --stability-R 200 \\
  --stability-block-col Cruise \\
  --stability-oob-min 10 \\
  --stability-min-ari 0.25 \\
  --min-cluster-frac 0.01 \\
  --select-by icl \\
  --select-delta 5 \\
  --sep "," \\
  --outdir "${biochemEofStatesDirAbs}" \\
  --sticky-smoothing \\
  --time-col date \\
  --sticky-prob 0.85 \\
  --apply-to all
touch biochem_eof_state_cluster.done
"""
}

process BIOCHEM_EOF_MODE_PLOTS {
    cpus pipelineThreads
    conda "${biochemEofModeCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_eof_mode_plots.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemEofPlotsDirAbs}"
python "${biochemEofModePlotScriptPath}" \\
  --loadings "${biochemEofPcaDirAbs}/tables/eof_pca_loadings.csv" \\
  --explained "${biochemEofPcaDirAbs}/tables/eof_pca_explained_variance.csv" \\
  --outdir "${biochemEofPlotsDirAbs}" \\
  --eofs "${biochemEofPcs}" \\
  --top-n 100 \\
  --sep ","
touch biochem_eof_mode_plots.done
"""
}

process BIOCHEM_WITHIN_GMM_HDBSCAN {
    cpus pipelineThreads
    conda "${biochemWithinGmmCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_within_gmm.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemWithinGmmDirAbs}"
python "${biochemWithinGmmScriptPath}" \\
  --eigenvectors "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" \\
  --assignments "${biochemGmmDirAbs}/tables/compartments_assignments_smoothed.csv" \\
  --outdir "${biochemWithinGmmDirAbs}" \\
  --sep "," \\
  --pc-cols "PC1,PC2" \\
  --standardize-pc-space \\
  --hdbscan-min-cluster-size 10 \\
  --hdbscan-metric euclidean \\
  --min-rows-per-component 10 \\
  --high-conf-only \\
  --high-conf-maxprob 0.80 \\
  --strict-unique-ids
touch biochem_within_gmm.done
"""
}
