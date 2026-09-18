#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def projectRootDir = new File(projectDir.toString())
def defaultConfigPath = "${projectDir}/basin_pipeline_nextflow.yml"
def paramsMap = [:]

try {
    paramsMap = new LinkedHashMap(params)
} catch( Throwable ignored ) {
    paramsMap = [:]
}

def inlineKeys = ['paths', 'resources', 'biochem', 'biochem_pre_asv', 'continuous_sections', 'gapseq_media', 'genome_modeling', 'master_summary', 'environments', 'config_root', 'pipeline_config']
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

def readPathManifest = { String manifestPath, String manifestName ->
    File manifestFile = new File(manifestPath)
    List<String> lines = manifestFile.readLines('UTF-8')
        .findAll { line -> line.trim() && !line.trim().startsWith('#') }
    if( !lines ) {
        exit 1, "${manifestName} is empty: ${manifestFile}"
    }
    String header = lines[0].replaceFirst(/^\uFEFF/, '')
    String delimiter = header.contains('\t') ? '\t' : ','
    List<String> columns = header.split(java.util.regex.Pattern.quote(delimiter), -1)
        .collect { it.trim().toLowerCase() }
    if( columns != ['label', 'filepath'] ) {
        exit 1, "${manifestName} must have exactly two columns in this order: label, filepath"
    }
    List<Map> records = []
    Set<String> labels = [] as Set
    lines.drop(1).eachWithIndex { line, index ->
        List<String> fields = line.split(java.util.regex.Pattern.quote(delimiter), -1)
            .collect { it.trim() }
        int lineNumber = index + 2
        if( fields.size() != 2 || !fields[0] || !fields[1] ) {
            exit 1, "${manifestName} line ${lineNumber} must contain a non-empty label and filepath"
        }
        String label = fields[0]
        if( !(label ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/) ) {
            exit 1, "${manifestName} line ${lineNumber} has unsafe label '${label}'; use letters, numbers, period, underscore, or hyphen"
        }
        String pathText = fields[1]
        boolean isGlob = pathText.find(/[*?\[\]{}]/) != null
        File candidate = new File(pathText)
        String resolvedText = (candidate.isAbsolute() ? candidate : new File(manifestFile.parentFile, pathText))
            .absoluteFile.toPath().normalize().toString()
        List<File> matches
        if( isGlob ) {
            int wildcardAt = resolvedText.findIndexOf { character -> '*?[]{}'.contains(character.toString()) }
            int separatorAt = resolvedText.substring(0, wildcardAt).lastIndexOf(File.separator)
            File searchRoot = new File(separatorAt > 0 ? resolvedText.substring(0, separatorAt) : File.separator)
            if( !searchRoot.isDirectory() ) {
                exit 1, "${manifestName} glob search root not found on line ${lineNumber}: ${searchRoot}"
            }
            def matcher = java.nio.file.FileSystems.default.getPathMatcher("glob:${resolvedText}")
            def stream = java.nio.file.Files.walk(searchRoot.toPath())
            matches = stream
                .filter { path -> java.nio.file.Files.isRegularFile(path) && matcher.matches(path.toAbsolutePath().normalize()) }
                .map { path -> path.toFile() }
                .sorted()
                .collect(java.util.stream.Collectors.toList())
            stream.close()
            if( !matches ) {
                exit 1, "${manifestName} glob matched no files for '${label}' on line ${lineNumber}: ${pathText}"
            }
        } else {
            File resolved = new File(resolvedText)
            if( !resolved.isFile() ) {
                exit 1, "${manifestName} filepath not found for '${label}': ${resolved}"
            }
            matches = [resolved]
        }
        matches.each { matched ->
            String recordLabel = label
            if( label.equalsIgnoreCase('auto') || (isGlob && matches.size() > 1) ) {
                String basename = matched.name
                    .replaceFirst(/(?i)\.gz$/, '')
                    .replaceFirst(/(?i)\.(fna|fa|fasta|faa)$/, '')
                recordLabel = label.equalsIgnoreCase('auto') ? basename : "${label}_${basename}"
            }
            if( !labels.add(recordLabel) ) {
                exit 1, "${manifestName} produces duplicate genome label '${recordLabel}'; use a distinct label prefix"
            }
            records << [label: recordLabel, filepath: matched.canonicalPath]
        }
    }
    if( !records ) {
        exit 1, "${manifestName} contains no data rows: ${manifestFile}"
    }
    return records
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

def publicOutputDir = resolvePath(config.paths?.output_dir)
assert publicOutputDir : "paths.output_dir must be provided in the YAML config"
def runtimeDir = config.paths?.runtime_dir ? resolvePath(config.paths.runtime_dir) : new File(publicOutputDir, '.basin').canonicalPath
def outputDir = new File(runtimeDir, 'publication_staging').canonicalPath
new File(outputDir).mkdirs()
log.info "Using BASINS publication staging directory: ${outputDir}"

def workDirOverride = config.paths?.work_dir ? resolvePath(config.paths.work_dir) : null
def condaCacheOverride = config.paths?.conda_cache_dir ? resolvePath(config.paths.conda_cache_dir) : null
if( workDirOverride ) {
    def workDirFile = new File(workDirOverride)
    workDirFile.mkdirs()
    workflow.workDir = java.nio.file.Paths.get(workDirFile.canonicalPath)
    log.info "Using custom Nextflow work directory: ${workflow.workDir}"
}
def resolvedCondaCacheDir = condaCacheOverride ?: new File(runtimeDir, "conda_cache").canonicalPath
def condaCacheDirFile = new File(resolvedCondaCacheDir)
condaCacheDirFile.mkdirs()
System.setProperty('NXF_CONDA_CACHEDIR', condaCacheDirFile.canonicalPath)
log.info "Using Conda cache directory: ${condaCacheDirFile.canonicalPath}"

int hostThreads = Runtime.runtime.availableProcessors()
int pipelineThreads = config.resources?.threads ? (config.resources.threads as int) : hostThreads

def biochemConfig = config.biochem ?: (config.biochem_pre_asv ?: [:])
def continuousSectionsConfig = config.continuous_sections ?: [:]
boolean biochemEnabled = biochemConfig.containsKey('enabled') ? (biochemConfig.enabled as boolean) : true
def gapseqMediaConfig = config.gapseq_media ?: [:]
boolean gapseqMediaEnabled = gapseqMediaConfig.containsKey('enabled') ? (gapseqMediaConfig.enabled as boolean) : false
def genomeModelingConfig = config.genome_modeling ?: [:]
boolean genomeModelingEnabled = genomeModelingConfig.containsKey('enabled') ? (genomeModelingConfig.enabled as boolean) : false
def masterSummaryConfig = config.master_summary ?: [:]
boolean masterSummaryEnabled = masterSummaryConfig.containsKey('enabled') ? (masterSummaryConfig.enabled as boolean) : true
def biochemTableAPath = biochemConfig.table_a ? resolveOptionalPath(biochemConfig.table_a, configRoot) : null
def biochemTableBPath = biochemConfig.table_b ? resolveOptionalPath(biochemConfig.table_b, configRoot) : null
double biochemCtdMaxDepthDiff = biochemConfig.ctd_max_depth_diff_m != null ? (biochemConfig.ctd_max_depth_diff_m as double) : 10.0
if( biochemEnabled && (!biochemTableAPath || !new File(biochemTableAPath).exists()) ) {
    exit 1, "biochem.table_a not found: ${biochemTableAPath}"
}
if( biochemEnabled && (!biochemTableBPath || !new File(biochemTableBPath).exists()) ) {
    exit 1, "biochem.table_b not found: ${biochemTableBPath}"
}
if( biochemEnabled && biochemCtdMaxDepthDiff < 0 ) {
    exit 1, "biochem.ctd_max_depth_diff_m must be nonnegative"
}

def configuredBiochemOutput = biochemConfig.output_root?.toString()?.trim()
def biochemOutputRoot = outputDir
if( configuredBiochemOutput ) {
    def configuredFile = new File(configuredBiochemOutput)
    if( configuredFile.isAbsolute() ) {
        def publicRoot = new File(publicOutputDir).canonicalFile
        def requested = configuredFile.canonicalFile
        def publicPrefix = publicRoot.path + File.separator
        if( requested.path == publicRoot.path ) {
            biochemOutputRoot = outputDir
        } else if( requested.path.startsWith(publicPrefix) ) {
            def relative = publicRoot.toPath().relativize(requested.toPath()).toString()
            biochemOutputRoot = new File(outputDir, relative).canonicalPath
        } else {
            exit 1, "biochem.output_root must be empty, relative, or located under paths.output_dir for atomic publication: ${requested}"
        }
    } else {
        biochemOutputRoot = new File(outputDir, configuredBiochemOutput).canonicalPath
    }
}
def biochemProcessingDirAbs = new File(biochemOutputRoot, 'biochem_processing').canonicalPath
def biochemStratMetricsDirAbs = new File(biochemProcessingDirAbs, 'stratification_metrics').canonicalPath
def biochemPcaDirAbs = new File(biochemOutputRoot, 'env_pca').canonicalPath
def biochemContinuousSectionsDirAbs = new File(biochemOutputRoot, continuousSectionsConfig.output_dir ?: 'env_continuous_sections').canonicalPath
def biochemMissingnessSensitivityDirAbs = new File(biochemOutputRoot, 'env_missingness_sensitivity').canonicalPath
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
def biochemGapseqMediaDirAbs = new File(biochemOutputRoot, 'gapseq_media').canonicalPath
def biochemGenomeModelingDirAbs = new File(biochemOutputRoot, 'genome_modeling').canonicalPath
def biochemWithinGmmDirAbs = new File(biochemGmmDirAbs, 'within_gmm_hdbscan').canonicalPath
def biochemMasterSummaryDirAbs = new File(biochemOutputRoot, 'master_summary').canonicalPath
def biochemMergedOxygenPath = new File(biochemProcessingDirAbs, '02_oxygen_best_available.tsv').canonicalPath
def biochemDensityPath = new File(biochemProcessingDirAbs, '02_oxygen_best_available_density.tsv').canonicalPath
def biochemDensityCleanedFile = biochemConfig.cleaned_density_filename ?: '02_oxygen_best_available_density_RJM.tsv'
def biochemDensityCleanedPath = new File(biochemProcessingDirAbs, biochemDensityCleanedFile.toString()).canonicalPath
def gapseqNutrientsPath = gapseqMediaConfig.nutrients_tsv ? resolveOptionalPath(gapseqMediaConfig.nutrients_tsv, configRoot) : null
def gapseqTemplatePath = gapseqMediaConfig.template ? resolveOptionalPath(gapseqMediaConfig.template, configRoot) : null
def gapseqTemplateArg = gapseqTemplatePath ? "--template \"${gapseqTemplatePath}\"" : ""
def gapseqMappingPath = gapseqMediaConfig.compound_mapping ? resolveOptionalPath(gapseqMediaConfig.compound_mapping, configRoot) : new File("${projectDir}/processes/gapseq_media/compound_mapping.tsv").canonicalPath
def gapseqTransformMethod = gapseqMediaConfig.transform_method ?: 'ordinal'
double gapseqMinEffectiveN = gapseqMediaConfig.min_effective_n != null ? (gapseqMediaConfig.min_effective_n as double) : 1.0
double gapseqDetectionLimit = gapseqMediaConfig.detection_limit != null ? (gapseqMediaConfig.detection_limit as double) : 0.0
def gapseqDepthBaselineRaw = gapseqMediaConfig.depth_baseline_m ?: []
List<Double> gapseqDepthBaselineM = (
    gapseqDepthBaselineRaw instanceof List
        ? gapseqDepthBaselineRaw
        : gapseqDepthBaselineRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
).collect { it as double }.unique().sort()
def gapseqDepthBaselineArg = gapseqDepthBaselineM.join(',')
def genomeManifestPath = genomeModelingConfig.genomes_manifest ? resolveOptionalPath(genomeModelingConfig.genomes_manifest, configRoot) : null
def genomeMetadataPath = genomeModelingConfig.genome_metadata_tsv ? resolveOptionalPath(genomeModelingConfig.genome_metadata_tsv, configRoot) : null
def genomeMetadataArg = genomeMetadataPath ? "--genome-metadata \"${genomeMetadataPath}\"" : ""
def genomeAbundancePath = genomeModelingConfig.abundance_tsv ? resolveOptionalPath(genomeModelingConfig.abundance_tsv, configRoot) : null
def genomeAbundanceSampleCol = genomeModelingConfig.abundance_sample_col ?: 'sample'
def genomeAbundanceGenomeCol = genomeModelingConfig.abundance_genome_col ?: 'genome'
def genomeAbundanceValueCol = genomeModelingConfig.abundance_value_col ?: 'read_count'
def genomeAbundanceNormalization = genomeModelingConfig.abundance_normalization ?: 'auto'
def genomeAbundanceSeqkitPath = genomeModelingConfig.abundance_seqkit_tsv ? resolveOptionalPath(genomeModelingConfig.abundance_seqkit_tsv, configRoot) : null
def genomeAbundanceSeqkitFileCol = genomeModelingConfig.abundance_seqkit_file_col ?: 'file'
def genomeAbundanceSeqkitCountCol = genomeModelingConfig.abundance_seqkit_count_col ?: 'num_seqs'
def genomeAbundanceSeqkitArg = genomeAbundanceSeqkitPath ? "--abundance-seqkit \"${genomeAbundanceSeqkitPath}\" --abundance-seqkit-file-col \"${genomeAbundanceSeqkitFileCol}\" --abundance-seqkit-count-col \"${genomeAbundanceSeqkitCountCol}\"" : ""
def genomeAbundanceArg = genomeAbundancePath ? "--abundance \"${genomeAbundancePath}\" --abundance-sample-col \"${genomeAbundanceSampleCol}\" --abundance-genome-col \"${genomeAbundanceGenomeCol}\" --abundance-value-col \"${genomeAbundanceValueCol}\" --abundance-normalization \"${genomeAbundanceNormalization}\" ${genomeAbundanceSeqkitArg}" : ""
def genomeTranscriptAbundancePath = genomeModelingConfig.transcript_abundance_tsv ? resolveOptionalPath(genomeModelingConfig.transcript_abundance_tsv, configRoot) : null
def genomeTranscriptAbundanceSampleCol = genomeModelingConfig.transcript_abundance_sample_col ?: 'sample'
def genomeTranscriptAbundanceGenomeCol = genomeModelingConfig.transcript_abundance_genome_col ?: 'genome'
def genomeTranscriptAbundanceValueCol = genomeModelingConfig.transcript_abundance_value_col ?: 'read_count'
def genomeTranscriptAbundanceNormalization = genomeModelingConfig.transcript_abundance_normalization ?: 'auto'
def genomeTranscriptAbundanceSeqkitPath = genomeModelingConfig.transcript_abundance_seqkit_tsv ? resolveOptionalPath(genomeModelingConfig.transcript_abundance_seqkit_tsv, configRoot) : null
def genomeTranscriptAbundanceSeqkitFileCol = genomeModelingConfig.transcript_abundance_seqkit_file_col ?: 'file'
def genomeTranscriptAbundanceSeqkitCountCol = genomeModelingConfig.transcript_abundance_seqkit_count_col ?: 'num_seqs'
def genomeTranscriptAbundanceSeqkitArg = genomeTranscriptAbundanceSeqkitPath ? "--transcript-abundance-seqkit \"${genomeTranscriptAbundanceSeqkitPath}\" --transcript-abundance-seqkit-file-col \"${genomeTranscriptAbundanceSeqkitFileCol}\" --transcript-abundance-seqkit-count-col \"${genomeTranscriptAbundanceSeqkitCountCol}\"" : ""
def genomeTranscriptAbundanceArg = genomeTranscriptAbundancePath ? "--transcript-abundance \"${genomeTranscriptAbundancePath}\" --transcript-abundance-sample-col \"${genomeTranscriptAbundanceSampleCol}\" --transcript-abundance-genome-col \"${genomeTranscriptAbundanceGenomeCol}\" --transcript-abundance-value-col \"${genomeTranscriptAbundanceValueCol}\" --transcript-abundance-normalization \"${genomeTranscriptAbundanceNormalization}\" ${genomeTranscriptAbundanceSeqkitArg}" : ""
int genomeAgreementPermutations = genomeModelingConfig.agreement_permutations != null ? (genomeModelingConfig.agreement_permutations as int) : 10000
int genomeAgreementBootstrapIterations = genomeModelingConfig.agreement_bootstrap_iterations != null ? (genomeModelingConfig.agreement_bootstrap_iterations as int) : 5000
int genomeAgreementRandomState = genomeModelingConfig.agreement_random_state != null ? (genomeModelingConfig.agreement_random_state as int) : 42
double genomeExpectedCompartmentTopFraction = genomeModelingConfig.expected_compartment_top_fraction != null ? (genomeModelingConfig.expected_compartment_top_fraction as double) : 0.25
if( genomeAgreementPermutations < 1 || genomeAgreementBootstrapIterations < 1 ) {
    exit 1, "genome_modeling agreement iteration counts must be positive"
}
if( genomeExpectedCompartmentTopFraction <= 0 || genomeExpectedCompartmentTopFraction >= 1 ) {
    exit 1, "genome_modeling.expected_compartment_top_fraction must be between zero and one"
}
def genomeManifestRows = []
def genomeTaxonomy = genomeModelingConfig.taxonomy ?: 'auto'
def genomeAligner = genomeModelingConfig.aligner ?: 'diamond'
def genomeMmseqs2Bin = genomeModelingConfig.mmseqs2_bin ? resolveOptionalPath(genomeModelingConfig.mmseqs2_bin, configRoot) : ''
int genomeCpusPerTask = genomeModelingConfig.cpus_per_genome != null ? (genomeModelingConfig.cpus_per_genome as int) : Math.min(8, pipelineThreads)
int genomeMaxForks = genomeModelingConfig.max_parallel_genomes != null ? (genomeModelingConfig.max_parallel_genomes as int) : Math.max(1, (int)(pipelineThreads / Math.max(1, genomeCpusPerTask)))
def genomeSolver = genomeModelingConfig.solver ?: 'glpk'
double genomeFractionOptimum = genomeModelingConfig.fraction_of_optimum != null ? (genomeModelingConfig.fraction_of_optimum as double) : 1.0
double genomeFluxThreshold = genomeModelingConfig.flux_threshold != null ? (genomeModelingConfig.flux_threshold as double) : 1e-9
double genomeSupplementMaxFlux = genomeModelingConfig.supplement_max_flux != null ? (genomeModelingConfig.supplement_max_flux as double) : 10.0
boolean genomeRunFva = genomeModelingConfig.containsKey('run_fva') ? (genomeModelingConfig.run_fva as boolean) : true
boolean genomeRunCounterfactuals = genomeModelingConfig.containsKey('run_counterfactuals') ? (genomeModelingConfig.run_counterfactuals as boolean) : true
def genomeCounterfactualCompounds = genomeModelingConfig.counterfactual_compounds ?: 'all'
if( gapseqMediaEnabled ) {
    if( !gapseqNutrientsPath || !new File(gapseqNutrientsPath).isFile() ) {
        exit 1, "gapseq_media.nutrients_tsv not found: ${gapseqNutrientsPath}"
    }
    if( gapseqTemplatePath && !new File(gapseqTemplatePath).isFile() ) {
        exit 1, "gapseq_media.template not found: ${gapseqTemplatePath}"
    }
    if( !new File(gapseqMappingPath).isFile() ) {
        exit 1, "gapseq_media.compound_mapping not found: ${gapseqMappingPath}"
    }
    if( !(gapseqTransformMethod in ['ordinal', 'binary']) ) {
        exit 1, "gapseq_media.transform_method must be ordinal or binary"
    }
    if( gapseqDepthBaselineM.any { it.isNaN() || it.isInfinite() || it < 0 } ) {
        exit 1, "gapseq_media.depth_baseline_m must contain finite nonnegative depths"
    }
}
if( genomeModelingEnabled ) {
    if( !gapseqMediaEnabled ) {
        exit 1, "genome_modeling.enabled requires gapseq_media.enabled"
    }
    if( !genomeManifestPath || !new File(genomeManifestPath).isFile() ) {
        exit 1, "genome_modeling.genomes_manifest not found: ${genomeManifestPath}"
    }
    if( genomeMetadataPath && !new File(genomeMetadataPath).isFile() ) {
        exit 1, "genome_modeling.genome_metadata_tsv not found: ${genomeMetadataPath}"
    }
    if( genomeAbundancePath && !new File(genomeAbundancePath).isFile() ) {
        exit 1, "genome_modeling.abundance_tsv not found: ${genomeAbundancePath}"
    }
    if( genomeAbundanceSeqkitPath && !new File(genomeAbundanceSeqkitPath).isFile() ) {
        exit 1, "genome_modeling.abundance_seqkit_tsv not found: ${genomeAbundanceSeqkitPath}"
    }
    if( genomeTranscriptAbundancePath && !new File(genomeTranscriptAbundancePath).isFile() ) {
        exit 1, "genome_modeling.transcript_abundance_tsv not found: ${genomeTranscriptAbundancePath}"
    }
    if( genomeTranscriptAbundanceSeqkitPath && !new File(genomeTranscriptAbundanceSeqkitPath).isFile() ) {
        exit 1, "genome_modeling.transcript_abundance_seqkit_tsv not found: ${genomeTranscriptAbundanceSeqkitPath}"
    }
    if( !(genomeAbundanceNormalization in ['auto', 'input_fragment_fpm', 'provided_fpkm', 'provided_tpm', 'median_ratio', 'raw', 'relative']) ) {
        exit 1, "genome_modeling.abundance_normalization must be auto, input_fragment_fpm, provided_fpkm, provided_tpm, median_ratio, raw, or relative"
    }
    if( !(genomeTranscriptAbundanceNormalization in ['auto', 'input_fragment_fpm', 'provided_fpkm', 'provided_tpm', 'median_ratio', 'raw', 'relative']) ) {
        exit 1, "genome_modeling.transcript_abundance_normalization must be auto, input_fragment_fpm, provided_fpkm, provided_tpm, median_ratio, raw, or relative"
    }
    if( genomeAbundancePath && (genomeAbundanceNormalization in ['auto', 'input_fragment_fpm']) && !genomeAbundanceSeqkitPath ) {
        exit 1, "genome_modeling.abundance_normalization=${genomeAbundanceNormalization} requires abundance_seqkit_tsv; select median_ratio explicitly to use the legacy normalization"
    }
    if( genomeTranscriptAbundancePath && (genomeTranscriptAbundanceNormalization in ['auto', 'input_fragment_fpm']) && !genomeTranscriptAbundanceSeqkitPath ) {
        exit 1, "genome_modeling.transcript_abundance_normalization=${genomeTranscriptAbundanceNormalization} requires transcript_abundance_seqkit_tsv; select median_ratio explicitly to use the legacy normalization"
    }
    if( genomeAbundanceSeqkitPath && !genomeAbundancePath ) {
        exit 1, "genome_modeling.abundance_seqkit_tsv requires abundance_tsv"
    }
    if( genomeTranscriptAbundanceSeqkitPath && !genomeTranscriptAbundancePath ) {
        exit 1, "genome_modeling.transcript_abundance_seqkit_tsv requires transcript_abundance_tsv"
    }
    genomeManifestRows = readPathManifest(genomeManifestPath, 'genome_modeling.genomes_manifest')
    genomeManifestRows.each { row ->
        if( !(row.filepath ==~ /(?i).+\.(fna|fa|fasta|faa)(\.gz)?$/) ) {
            exit 1, "Unsupported genome extension for '${row.label}': ${row.filepath}"
        }
    }
    if( !(genomeTaxonomy in ['auto', 'Bacteria', 'Archaea']) ) {
        exit 1, "genome_modeling.taxonomy must be auto, Bacteria, or Archaea"
    }
    if( !(genomeAligner in ['blast', 'diamond', 'mmseqs2']) ) {
        exit 1, "genome_modeling.aligner must be blast, diamond, or mmseqs2"
    }
    if( genomeAligner == 'mmseqs2' && (!genomeMmseqs2Bin || !new File(genomeMmseqs2Bin).isFile()) ) {
        exit 1, "genome_modeling.aligner mmseqs2 requires an executable genome_modeling.mmseqs2_bin"
    }
    if( genomeCpusPerTask < 1 || genomeMaxForks < 1 ) {
        exit 1, "genome_modeling CPU and parallel-genome settings must be positive"
    }
    if( genomeSupplementMaxFlux <= 0 ) {
        exit 1, "genome_modeling.supplement_max_flux must be positive"
    }
}
def biochemFeatureCols = biochemConfig.feature_cols ?: 'Oxygen,Nitrate,Nitrite,Nitrous Oxide,Ammonium,Hydrogen Sulfide,Methane,Phosphate,Silicate,Temperature,Salinity,Density,Fe,Dimethyl Sulfide'
boolean biochemContinuousSectionsEnabled = continuousSectionsConfig.containsKey('enabled') ? (continuousSectionsConfig.enabled as boolean) : true
def biochemContinuousVariables = continuousSectionsConfig.variables ?: biochemFeatureCols
if( biochemContinuousVariables instanceof List ) {
    biochemContinuousVariables = biochemContinuousVariables.join(',')
}
def biochemContinuousDateCol = continuousSectionsConfig.date_col ?: 'date'
def biochemContinuousDepthCol = continuousSectionsConfig.depth_col ?: 'Depth_anchored'
double biochemContinuousMaxDepth = continuousSectionsConfig.maximum_depth_m != null ? (continuousSectionsConfig.maximum_depth_m as double) : 210.0
double biochemContinuousDepthStep = continuousSectionsConfig.depth_step_m != null ? (continuousSectionsConfig.depth_step_m as double) : 1.0
int biochemContinuousTimeStepDays = continuousSectionsConfig.time_step_days != null ? (continuousSectionsConfig.time_step_days as int) : 7
double biochemContinuousMaxTimeSupport = continuousSectionsConfig.max_time_support_days != null ? (continuousSectionsConfig.max_time_support_days as double) : 90.0
double biochemContinuousMaxDepthSupport = continuousSectionsConfig.max_depth_support_m != null ? (continuousSectionsConfig.max_depth_support_m as double) : 30.0
int biochemContinuousMinimumSamples = continuousSectionsConfig.minimum_samples != null ? (continuousSectionsConfig.minimum_samples as int) : 4
int biochemContinuousLevels = continuousSectionsConfig.levels != null ? (continuousSectionsConfig.levels as int) : 64
def biochemContinuousFormatsRaw = continuousSectionsConfig.formats ?: ['pdf', 'png', 'svg']
def biochemContinuousFormats = biochemContinuousFormatsRaw instanceof List ? biochemContinuousFormatsRaw.join(',') : biochemContinuousFormatsRaw.toString()
if( biochemContinuousMaxDepth <= 0 || biochemContinuousDepthStep <= 0 || biochemContinuousTimeStepDays < 1 || biochemContinuousMaxTimeSupport <= 0 || biochemContinuousMaxDepthSupport <= 0 || biochemContinuousMinimumSamples < 3 || biochemContinuousLevels < 8 ) {
    exit 1, "continuous_sections contains invalid grid, support, or sample settings"
}
double biochemMissingnessCutoff = biochemConfig.missingness_cutoff != null ? (biochemConfig.missingness_cutoff as double) : 0.20
boolean biochemMissingnessSensitivityEnabled = biochemConfig.containsKey('missingness_sensitivity_enabled') ? (biochemConfig.missingness_sensitivity_enabled as boolean) : false
def biochemMissingnessThresholdsRaw = biochemConfig.missingness_sensitivity_thresholds ?: [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
List<Double> biochemMissingnessThresholds = (
    biochemMissingnessThresholdsRaw instanceof List
        ? biochemMissingnessThresholdsRaw
        : biochemMissingnessThresholdsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
).collect { it as double }.unique().sort()
int biochemMissingnessPcParallelReps = biochemConfig.missingness_sensitivity_pc_parallel_replicates != null ? (biochemConfig.missingness_sensitivity_pc_parallel_replicates as int) : 100
int biochemMissingnessPcStabilityReps = biochemConfig.missingness_sensitivity_pc_stability_replicates != null ? (biochemConfig.missingness_sensitivity_pc_stability_replicates as int) : 100
int biochemMissingnessGmmStabilityReps = biochemConfig.missingness_sensitivity_gmm_stability_replicates != null ? (biochemConfig.missingness_sensitivity_gmm_stability_replicates as int) : 100
double biochemMissingnessMinGmmAri = biochemConfig.missingness_sensitivity_min_gmm_ari != null ? (biochemConfig.missingness_sensitivity_min_gmm_ari as double) : 0.70
if( biochemMissingnessCutoff.isNaN() || biochemMissingnessCutoff.isInfinite() || biochemMissingnessCutoff <= 0.0 || biochemMissingnessCutoff >= 1.0 ) {
    exit 1, "biochem.missingness_cutoff must be finite and strictly between 0 and 1"
}
if( !biochemMissingnessThresholds || biochemMissingnessThresholds.any { it.isNaN() || it.isInfinite() || it <= 0.0 || it >= 1.0 } ) {
    exit 1, "biochem.missingness_sensitivity_thresholds must contain finite values strictly between 0 and 1"
}
if( [biochemMissingnessPcParallelReps, biochemMissingnessPcStabilityReps, biochemMissingnessGmmStabilityReps].any { it < 1 } ) {
    exit 1, "biochem missingness-sensitivity replicate counts must be positive"
}
if( [biochemMissingnessMinGmmAri].any { it.isNaN() || it.isInfinite() || it < 0.0 || it > 1.0 } ) {
    exit 1, "biochem missingness-sensitivity feasibility thresholds must be finite and between 0 and 1"
}
def biochemMissingnessThresholdsCsv = biochemMissingnessThresholds.join(',')
def biochemGmmKRaw = biochemConfig.gmm_k
boolean biochemGmmKAuto = (biochemGmmKRaw == null) || (biochemGmmKRaw.toString().trim().equalsIgnoreCase('auto'))
def biochemGmmK = biochemGmmKAuto ? 5 : (biochemGmmKRaw as int)
def biochemEofPcs = biochemConfig.eof_pcs ?: '1,2'
int biochemEofKMin = biochemConfig.eof_k_min != null ? (biochemConfig.eof_k_min as int) : 2
int biochemEofKMax = biochemConfig.eof_k_max != null ? (biochemConfig.eof_k_max as int) : 4
def biochemEofCovarianceType = biochemConfig.eof_covariance_type ?: 'tied'
double biochemEofStabilityMinAri = biochemConfig.eof_stability_min_ari != null ? (biochemConfig.eof_stability_min_ari as double) : 0.50
double biochemEofMinClusterFrac = biochemConfig.eof_min_cluster_frac != null ? (biochemConfig.eof_min_cluster_frac as double) : 0.10
int biochemEofMinClusterN = biochemConfig.eof_min_cluster_n != null ? (biochemConfig.eof_min_cluster_n as int) : 5
double biochemEofAssignmentProbThreshold = biochemConfig.eof_assignment_prob_threshold != null ? (biochemConfig.eof_assignment_prob_threshold as double) : 0.80
if( biochemEofAssignmentProbThreshold < 0.0 || biochemEofAssignmentProbThreshold > 1.0 ) {
    exit 1, "biochemical_analysis.eof_assignment_prob_threshold must be between 0 and 1"
}
int biochemEofBaselineMonths = biochemConfig.eof_baseline_months != null ? (biochemConfig.eof_baseline_months as int) : 24
int biochemEofBaselineMinCruises = biochemConfig.eof_baseline_min_cruises != null ? (biochemConfig.eof_baseline_min_cruises as int) : 5
int biochemPeaBootstrapIterations = biochemConfig.pea_bootstrap_iterations != null ? (biochemConfig.pea_bootstrap_iterations as int) : 500
int biochemPeaRandomState = biochemConfig.pea_random_state != null ? (biochemConfig.pea_random_state as int) : 42
int biochemPhysicalRegimeKMax = biochemConfig.physical_regime_k_max != null ? (biochemConfig.physical_regime_k_max as int) : 6
double biochemDeepIntrusionQuantile = biochemConfig.deep_intrusion_quantile != null ? (biochemConfig.deep_intrusion_quantile as double) : 0.90
double biochemOxygenLowCompartmentMax = biochemConfig.oxygen_low_compartment_max != null ? (biochemConfig.oxygen_low_compartment_max as double) : 90.0
int biochemOxygenIntrusionBottomN = biochemConfig.oxygen_intrusion_bottom_n != null ? (biochemConfig.oxygen_intrusion_bottom_n as int) : 3
double biochemOxygenIntrusionOnsetThreshold = biochemConfig.oxygen_intrusion_onset_threshold != null ? (biochemConfig.oxygen_intrusion_onset_threshold as double) : 4.5
double biochemOxygenIntrusionPersistenceThreshold = biochemConfig.oxygen_intrusion_persistence_threshold != null ? (biochemConfig.oxygen_intrusion_persistence_threshold as double) : 4.5
int biochemOxygenIntrusionEndConsecutive = biochemConfig.oxygen_intrusion_end_consecutive != null ? (biochemConfig.oxygen_intrusion_end_consecutive as int) : 2
def biochemRenewalNitrateRequestedCol = biochemConfig.renewal_nitrate_col ?: 'Nitrate'
int biochemRenewalNitrateMinDepths = biochemConfig.renewal_nitrate_min_depths != null ? (biochemConfig.renewal_nitrate_min_depths as int) : 2
int biochemRenewalNitrateBottomN = biochemConfig.renewal_nitrate_bottom_n != null ? (biochemConfig.renewal_nitrate_bottom_n as int) : 3
double biochemRenewalNitrateDetectionLimit = biochemConfig.renewal_nitrate_detection_limit != null ? (biochemConfig.renewal_nitrate_detection_limit as double) : 0.0
boolean biochemRenewalBridgeEnabled = biochemConfig.renewal_bridge_enabled != null ? (biochemConfig.renewal_bridge_enabled as boolean) : true
int biochemRenewalBridgeMaxCruises = biochemConfig.renewal_bridge_max_cruises != null ? (biochemConfig.renewal_bridge_max_cruises as int) : 2
double biochemRenewalBridgeMaxDays = biochemConfig.renewal_bridge_max_days != null ? (biochemConfig.renewal_bridge_max_days as double) : 150.0
def biochemRenewalBridgeFlag = biochemRenewalBridgeEnabled ? '--renewal-bridge-enabled' : ''
if( biochemRenewalNitrateBottomN < 1 ) {
    exit 1, "biochem.renewal_nitrate_bottom_n must be positive"
}
if( biochemRenewalNitrateMinDepths < 1 || biochemRenewalNitrateMinDepths > biochemRenewalNitrateBottomN ) {
    exit 1, "biochem.renewal_nitrate_min_depths must be between 1 and renewal_nitrate_bottom_n"
}
if( biochemRenewalNitrateDetectionLimit < 0.0 ) {
    exit 1, "biochem.renewal_nitrate_detection_limit must be nonnegative"
}
if( biochemRenewalBridgeMaxCruises < 1 || biochemRenewalBridgeMaxDays <= 0.0 ) {
    exit 1, "biochem renewal bridge limits must be positive"
}
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
// Stratification precedes the canonical clean/rename stage. When the renewal
// setting uses a canonical name, resolve it back to the corresponding source
// column so the stratification process reads the actual upstream table.
def biochemRenewalNitrateCol = biochemCleanRenameMap.find {
    sourceName, canonicalName -> canonicalName == biochemRenewalNitrateRequestedCol
}?.key ?: biochemRenewalNitrateRequestedCol

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
def biochemMissingnessSensitivityCondaEnvPath = resolveBiochemStepEnv('biochem_missingness_sensitivity')
def biochemEigenvectorsCondaEnvPath = resolveBiochemStepEnv('biochem_eigenvectors')
def biochemContinuousSectionsCondaEnvPath = resolveBiochemStepEnv('biochem_continuous_sections')
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
def biochemGapseqMediaCondaEnvPath = resolveBiochemStepEnv('biochem_gapseq_media')
def defaultGenomeModelingEnvPath = new File("${projectDir}/processes/shared_envs/gapseq_modeling.yml").canonicalPath
def defaultGenomeSimulationEnvPath = new File("${projectDir}/processes/shared_envs/cobrapy_simulation.yml").canonicalPath
def genomeModelingCondaEnvPath = config.environments?.genome_modeling ? resolveOptionalPath(config.environments.genome_modeling, configRoot) : defaultGenomeModelingEnvPath
def genomeSimulationCondaEnvPath = config.environments?.genome_simulation ? resolveOptionalPath(config.environments.genome_simulation, configRoot) : defaultGenomeSimulationEnvPath
if( !new File(genomeModelingCondaEnvPath).isFile() || !new File(genomeSimulationCondaEnvPath).isFile() ) {
    exit 1, "Genome modeling Conda environment YAML is missing"
}
def biochemMasterSummaryCondaEnvPath = resolveBiochemStepEnv('master_summary')

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
def biochemMissingnessSensitivityScriptPath = scriptPath('processes/missingness_sensitivity/run_missingness_sensitivity.py')
def biochemEigenvectorsScriptPath = scriptPath('processes/eigenvectors/env_eigenvectors.py')
def biochemContinuousSectionsScriptPath = scriptPath('processes/continuous_sections/continuous_time_depth_sections.py')
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
def biochemLocalCruiseEofScriptPath = scriptPath('processes/eof_state_cluster/local_cruise_eof.py')
def biochemEofStateInterpretationScriptPath = scriptPath('processes/eof_state_cluster/eof_state_interpretation.py')
def biochemCruiseGroupSeasonBenchmarkScriptPath = scriptPath('processes/eof_state_cluster/cruise_group_season_benchmark.py')
def biochemEofModePlotScriptPath = scriptPath('processes/eof_mode_plots/eof_mode_plots.py')
def biochemWithinGmmScriptPath = scriptPath('processes/within_gmm_hdbscan/env_within_gmm_hdbscan.py')
def biochemGapseqMediaScriptPath = scriptPath('processes/gapseq_media/build_gapseq_media.py')
def genomeReconstructionScriptPath = scriptPath('processes/gapseq_modeling/run_gapseq_reconstruction.sh')
def genomeCompareScriptPath = scriptPath('processes/gapseq_modeling/compare_gapseq_media.py')
def genomeCombineScriptPath = scriptPath('processes/gapseq_modeling/combine_gapseq_results.py')
def genomeNicheBiplotScriptPath = scriptPath('processes/gapseq_modeling/plot_genome_niche_biplot.py')
def biochemMasterSummaryScriptPath = scriptPath('processes/master_summary/build_basin_summary.py')

workflow {
    if( !biochemEnabled ) {
        log.info "biochem.enabled is false; no BASINS stages will run."
        return
    }

    b01 = BIOCHEM_MERGE()
    b02 = BIOCHEM_DENSITY(b01.done)
    b03 = BIOCHEM_STRAT_METRICS(b02.done)
    b04 = BIOCHEM_CUSTOM_CLEAN(b03.done)
    if( biochemMissingnessSensitivityEnabled ) {
        b04s = BIOCHEM_MISSINGNESS_SENSITIVITY(b04.done)
        selected_missingness_cutoff = b04s.selected_cutoff.map { cutoff_file -> cutoff_file.text.trim() }
        b05 = BIOCHEM_EIGENVECTORS(b04s.done, selected_missingness_cutoff)
    } else {
        b05 = BIOCHEM_EIGENVECTORS(b04.done, biochemMissingnessCutoff.toString())
    }
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
    if( gapseqMediaEnabled ) {
        b19 = BIOCHEM_GAPSEQ_MEDIA(b18.done)
    } else {
        b19 = b18
    }
    b20 = BIOCHEM_WITHIN_GMM_HDBSCAN(b19.done)
    if( biochemContinuousSectionsEnabled ) {
        b21 = BIOCHEM_CONTINUOUS_SECTIONS(b20.done)
        continuous_done = b21.done
    } else {
        continuous_done = b20.done
    }
    if( genomeModelingEnabled ) {
        genome_inputs = Channel.fromList(genomeManifestRows)
            .map { row -> tuple(row.label, file(row.filepath, checkIfExists: true)) }
        reconstructed = GAPSEQ_RECONSTRUCT(genome_inputs)
        media_ready = b20.done.collect()
        compared = GAPSEQ_COMPARE_MEDIA(reconstructed.models, media_ready)
        result_dirs = compared.results.map { genome_id, result_dir -> result_dir }.collect()
        combined = GAPSEQ_COMBINE_RESULTS(result_dirs)
        final_done = combined.done
    } else {
        final_done = b20.done
    }
    if( masterSummaryEnabled ) {
        MASTER_SUMMARY(final_done, continuous_done)
    }
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
"\${CONDA_PREFIX}/bin/python" "${biochemMergeTablesScriptPath}" \\
  --table-a "${biochemTableAPath}" \\
  --table-b "${biochemTableBPath}" \\
  --outdir "${biochemProcessingDirAbs}" \\
  --max-depth-diff ${biochemCtdMaxDepthDiff}
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
"\${CONDA_PREFIX}/bin/python" "${biochemCalcDensityScriptPath}" \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemStratMetricsScriptPath}" \\
  --input "${biochemDensityPath}" \\
  --output-dir "${biochemStratMetricsDirAbs}" \\
  --salinity-col Salinity \\
  --temperature-col Temperature \\
  --depth-col Depth \\
  --latitude-col Latitude \\
  --longitude-col Longitude \\
  --profile-cols Cruise \\
  --date-col Date \\
  --layer-split-mode mld125 \\
  --pea-bootstrap-iterations ${biochemPeaBootstrapIterations} \\
  --pea-random-state ${biochemPeaRandomState} \\
  --physical-regime-k-max ${biochemPhysicalRegimeKMax} \\
  --deep-intrusion-quantile ${biochemDeepIntrusionQuantile} \\
  --deep-intrusion-oxygen-col Oxygen \\
  --oxygen-low-compartment-max ${biochemOxygenLowCompartmentMax} \\
  --oxygen-intrusion-bottom-n ${biochemOxygenIntrusionBottomN} \\
  --oxygen-intrusion-onset-threshold ${biochemOxygenIntrusionOnsetThreshold} \\
  --oxygen-intrusion-persistence-threshold ${biochemOxygenIntrusionPersistenceThreshold} \\
  --oxygen-intrusion-end-consecutive ${biochemOxygenIntrusionEndConsecutive} \\
  --renewal-nitrate-col "${biochemRenewalNitrateCol}" \\
  --renewal-nitrate-bottom-n ${biochemRenewalNitrateBottomN} \\
  --renewal-nitrate-min-depths ${biochemRenewalNitrateMinDepths} \\
  --renewal-nitrate-detection-limit ${biochemRenewalNitrateDetectionLimit} \\
  ${biochemRenewalBridgeFlag} \\
  --renewal-bridge-max-cruises ${biochemRenewalBridgeMaxCruises} \\
  --renewal-bridge-max-days ${biochemRenewalBridgeMaxDays}
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
"\${CONDA_PREFIX}/bin/python" "${biochemCustomCleanerScriptPath}" --input "${biochemDensityPath}" --output "${biochemDensityCleanedPath}"${biochemKeepArg}${biochemDropArg}${biochemRenameArg}
[[ -f "${biochemDensityCleanedPath}" ]] || { echo "Missing ${biochemDensityCleanedPath}" >&2; exit 1; }
ln -sf "${biochemDensityCleanedPath}" biochem_density_cleaned.tsv
touch biochem_custom_clean.done
"""
}

process BIOCHEM_MISSINGNESS_SENSITIVITY {
    cpus pipelineThreads
    conda "${biochemMissingnessSensitivityCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_missingness_sensitivity.done"), emit: done
    path("selected_missingness_cutoff.txt"), emit: selected_cutoff

    script:
    """
set -euo pipefail
export XDG_CACHE_HOME="\$PWD/.cache"
export MPLCONFIGDIR="\$PWD/.cache/matplotlib"
mkdir -p "\$XDG_CACHE_HOME" "\$MPLCONFIGDIR" "${biochemMissingnessSensitivityDirAbs}"
"\${CONDA_PREFIX}/bin/python" "${biochemMissingnessSensitivityScriptPath}" \\
  --input "${biochemDensityCleanedPath}" \\
  --outdir "${biochemMissingnessSensitivityDirAbs}" \\
  --feature-cols "${biochemFeatureCols}" \\
  --thresholds "${biochemMissingnessThresholdsCsv}" \\
  --eigenvectors-script "${biochemEigenvectorsScriptPath}" \\
  --selectk-script "${biochemSelectkScriptPath}" \\
  --gmm-script "${biochemGmmScriptPath}" \\
  --pc-parallel-replicates ${biochemMissingnessPcParallelReps} \\
  --pc-stability-replicates ${biochemMissingnessPcStabilityReps} \\
  --gmm-stability-replicates ${biochemMissingnessGmmStabilityReps} \\
  --min-cluster-frac 0.02 \\
  --min-gmm-stability-ari ${biochemMissingnessMinGmmAri}
[[ -s "${biochemMissingnessSensitivityDirAbs}/SELECTED_MISSINGNESS_CUTOFF.txt" ]] || { echo "Missing selected missingness cutoff" >&2; exit 1; }
ln -sf "${biochemMissingnessSensitivityDirAbs}/SELECTED_MISSINGNESS_CUTOFF.txt" selected_missingness_cutoff.txt
touch biochem_missingness_sensitivity.done
"""
}

process BIOCHEM_EIGENVECTORS {
    cpus pipelineThreads
    conda "${biochemEigenvectorsCondaEnvPath}"

    input:
    path(prev_done)
    val(missingness_cutoff)

    output:
    path("biochem_eigenvectors.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemPcaDirAbs}"
"\${CONDA_PREFIX}/bin/python" "${biochemEigenvectorsScriptPath}" \\
  --input "${biochemDensityCleanedPath}" \\
  --outdir "${biochemPcaDirAbs}" \\
  --feature-cols "${biochemFeatureCols}" \\
  --dropna-col-thresh "${missingness_cutoff}" \\
  --pc-selection \\
  --anchor-depths
[[ -f "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" ]] || { echo "Missing ${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" >&2; exit 1; }
[[ -f "${biochemPcaDirAbs}/tables/pc_keep_decision.csv" ]] || { echo "Missing ${biochemPcaDirAbs}/tables/pc_keep_decision.csv" >&2; exit 1; }
touch biochem_eigenvectors.done
"""
}

process BIOCHEM_CONTINUOUS_SECTIONS {
    cpus 1
    conda "${biochemContinuousSectionsCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_continuous_sections.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemContinuousSectionsDirAbs}"
export MPLCONFIGDIR="\$PWD/.matplotlib"
mkdir -p "\$MPLCONFIGDIR"
"\${CONDA_PREFIX}/bin/python" "${biochemContinuousSectionsScriptPath}" \
  --input "${biochemPcaDirAbs}/tables/matrix_cleaned_with_sparse.csv" \
  --outdir "${biochemContinuousSectionsDirAbs}" \
  --variables "${biochemContinuousVariables}" \
  --date-col "${biochemContinuousDateCol}" \
  --depth-col "${biochemContinuousDepthCol}" \
  --renewal-events "${biochemStratMetricsDirAbs}/stratification_nitrate_renewal_events.tsv" \
  --renewal-date-col start_date \
  --o2-assignments "${biochemO2DirAbs}/tables/o2_compartments_assignments_smoothed.csv" \
  --gmm-assignments "${biochemGmmDirAbs}/tables/compartments_assignments_smoothed.csv" \
  --hybrid-assignments "${biochemHybridDirAbs}/tables/compartments_assignments_hybrid.csv" \
  --maximum-depth ${biochemContinuousMaxDepth} \
  --depth-step ${biochemContinuousDepthStep} \
  --time-step-days ${biochemContinuousTimeStepDays} \
  --max-time-support-days ${biochemContinuousMaxTimeSupport} \
  --max-depth-support-m ${biochemContinuousMaxDepthSupport} \
  --minimum-samples ${biochemContinuousMinimumSamples} \
  --levels ${biochemContinuousLevels} \
  --formats "${biochemContinuousFormats}"
[[ -f "${biochemContinuousSectionsDirAbs}/tables/continuous_section_variable_audit.tsv" ]] || { echo "Missing continuous-section audit table" >&2; exit 1; }
[[ -f "${biochemContinuousSectionsDirAbs}/tables/continuous_time_depth_grids.tsv.gz" ]] || { echo "Missing continuous-section interpolation grid" >&2; exit 1; }
[[ -f "${biochemContinuousSectionsDirAbs}/tables/continuous_compartment_manifest.tsv" ]] || { echo "Missing continuous-compartment manifest" >&2; exit 1; }
touch biochem_continuous_sections.done
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
"\${CONDA_PREFIX}/bin/python" "${biochemSelectkScriptPath}" \\
  --eigenvectors "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" \\
  --pc-keep "${biochemPcaDirAbs}/tables/pc_keep_decision.csv" \\
  --outdir "${biochemSelectkDirAbs}" \\
  --sep "," \\
  --covariance-type tied \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemGmmScriptPath}" \\
  --eigenvectors "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" \\
  --pc-keep "${biochemPcaDirAbs}/tables/pc_keep_decision.csv" \\
  --outdir "${biochemGmmDirAbs}" \\
  --sep "," \\
  --pc-use-mode keep \\
  --standardize-pc-space \\
  --covariance-type tied \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemO2SoftScriptPath}" \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemHybridScriptPath}" \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemCompareScriptPath}" \\
  --matrix-cleaned "${biochemPcaDirAbs}/tables/matrix_cleaned_with_sparse.csv" \\
  --eigenvectors "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" \\
  --assignments "${biochemGmmDirAbs}/tables/compartments_assignments_smoothed.csv" \\
  --o2-assignments "${biochemO2DirAbs}/tables/o2_compartments_assignments_smoothed.csv" \\
  --hybrid-assignments "${biochemHybridDirAbs}/tables/compartments_assignments_hybrid.csv" \\
  --o2-compartment-col compartment_name \\
  --hybrid-compartment-col compartment_name \\
  --outdir "${biochemCompareDirAbs}" \\
  --sep-matrix "," \\
  --sep-eig "," \\
  --sep-assign "," \\
  --sep-o2-assign "," \\
  --sep-hybrid-assign "," \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemSplitScriptPath}" \\
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
  --min-subcluster-size 20 \\
  --curtain-hide-other \\
  --curtain-time-subdivisions-per-month 4 \\
  --curtain-time-sigma-months 0.75 \\
  --curtain-depth-step-m 1.0 \\
  --curtain-maximum-depth-m 210.0 \\
  --curtain-contour-visual-depth-sigma-m 5.0 \\
  --curtain-renewal-events "${biochemStratMetricsDirAbs}/stratification_nitrate_renewal_events.tsv" \\
  --curtain-renewal-date-col start_date
[[ -f "${biochemSplitDirAbs}/tables/merged_o2_split_by_gmm.csv" ]] || { echo "Missing ${biochemSplitDirAbs}/tables/merged_o2_split_by_gmm.csv" >&2; exit 1; }
[[ -f "${biochemSplitDirAbs}/tables/hybrid_compartment_time_depth_curtain_cells.csv" ]] || { echo "Missing hybrid time-depth curtain audit table" >&2; exit 1; }
[[ -f "${biochemSplitDirAbs}/tables/hybrid_compartment_time_depth_curtain_grid.csv" ]] || { echo "Missing hybrid time-depth curtain display grid" >&2; exit 1; }
[[ -f "${biochemSplitDirAbs}/tables/hybrid_compartment_time_depth_curtain_renewals.csv" ]] || { echo "Missing hybrid time-depth curtain renewal audit table" >&2; exit 1; }
[[ -f "${biochemSplitDirAbs}/plots/hybrid_compartment_time_depth_curtain.pdf" ]] || { echo "Missing hybrid time-depth curtain PDF" >&2; exit 1; }
[[ -f "${biochemSplitDirAbs}/plots/hybrid_compartment_time_depth_curtain.png" ]] || { echo "Missing hybrid time-depth curtain PNG" >&2; exit 1; }
[[ -f "${biochemSplitDirAbs}/plots/hybrid_compartment_time_depth_curtain.svg" ]] || { echo "Missing hybrid time-depth curtain SVG" >&2; exit 1; }
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
"\${CONDA_PREFIX}/bin/python" "${biochemStratAnomalyScriptPath}" \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemStateTransitionScriptPath}" \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemSuccessionScriptPath}" \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemFeatureAssocScriptPath}" \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemEofPipelineScriptPath}" \\
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
"\${CONDA_PREFIX}/bin/python" "${biochemLocalCruiseEofScriptPath}" \\
  --profiles "${biochemEofPcaDirAbs}/tables/eof_cruise_feature_matrix.tsv" \\
  --metadata "${biochemEofPcaDirAbs}/tables/eof_eigenvectors_scores_by_cruise.csv" \\
  --outdir "${biochemEofStatesDirAbs}/tables" \\
  --baseline-months ${biochemEofBaselineMonths} \\
  --baseline-min-cruises ${biochemEofBaselineMinCruises}
"\${CONDA_PREFIX}/bin/python" "${biochemEofStateScriptPath}" \\
  --scores "${biochemEofStatesDirAbs}/tables/local_eof_scores_by_cruise.csv" \\
  --pcs "${biochemEofPcs}" \\
  --k auto \\
  --k-min ${biochemEofKMin} \\
  --k-max ${biochemEofKMax} \\
  --covariance-type "${biochemEofCovarianceType}" \\
  --standardize-pc-space \\
  --n-init 20 \\
  --max-iter 500 \\
  --cv-folds 5 \\
  --stability-R 200 \\
  --stability-block-col Cruise \\
  --stability-oob-min 10 \\
  --stability-min-ari ${biochemEofStabilityMinAri} \\
  --min-cluster-frac ${biochemEofMinClusterFrac} \\
  --min-cluster-n ${biochemEofMinClusterN} \\
  --lowconf-maxprob ${biochemEofAssignmentProbThreshold} \\
  --select-by icl \\
  --select-delta 5 \\
  --sep "," \\
  --outdir "${biochemEofStatesDirAbs}" \\
  --time-col date
"\${CONDA_PREFIX}/bin/python" "${biochemEofStateInterpretationScriptPath}" \\
  --assignments "${biochemEofStatesDirAbs}/tables/gmm_selected_assignments.tsv" \\
  --cruise-profiles "${biochemEofPcaDirAbs}/tables/eof_cruise_feature_matrix.tsv" \\
  --stratification "${biochemStratIndexDirAbs}/stratification_physical_biochem_timeseries.tsv" \\
  --outdir "${biochemEofStatesDirAbs}" \\
  --formats pdf,png,svg
"\${CONDA_PREFIX}/bin/python" "${biochemCruiseGroupSeasonBenchmarkScriptPath}" \\
  --assignments "${biochemEofStatesDirAbs}/tables/cruise_group_assignments_interpreted.tsv" \\
  --matrix "${biochemPcaDirAbs}/tables/matrix_cleaned_with_sparse.csv" \\
  --outdir "${biochemEofStatesDirAbs}"
[[ -f "${biochemEofStatesDirAbs}/tables/cruise_group_season_redundancy.tsv" ]] || { echo "Missing cruise-group season benchmark output" >&2; exit 1; }
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
"\${CONDA_PREFIX}/bin/python" "${biochemEofModePlotScriptPath}" \\
  --loadings "${biochemEofPcaDirAbs}/tables/eof_pca_loadings.csv" \\
  --explained "${biochemEofPcaDirAbs}/tables/eof_pca_explained_variance.csv" \\
  --outdir "${biochemEofPlotsDirAbs}" \\
  --eofs "${biochemEofPcs}" \\
  --top-n 100 \\
  --sep ","
touch biochem_eof_mode_plots.done
"""
}

process GAPSEQ_RECONSTRUCT {
    tag "${genome_id}"
    cpus genomeCpusPerTask
    maxForks genomeMaxForks
    conda "${genomeModelingCondaEnvPath}"

    input:
    tuple val(genome_id), path(genome_file)

    output:
    tuple val(genome_id), path("${genome_id}_reconstruction"), emit: models

    script:
    """
set -euo pipefail
bash "${genomeReconstructionScriptPath}" \\
  "${genome_file}" \\
  "${genome_id}" \\
  "${genomeTaxonomy}" \\
  "${genomeAligner}" \\
  "${genomeMmseqs2Bin}" \\
  "${task.cpus}" \\
  "${genome_id}_reconstruction"
mkdir -p "${biochemGenomeModelingDirAbs}/reconstructions/${genome_id}"
cp -a "${genome_id}_reconstruction/." "${biochemGenomeModelingDirAbs}/reconstructions/${genome_id}/"
"""
}

process GAPSEQ_COMPARE_MEDIA {
    tag "${genome_id}"
    cpus 1
    maxForks genomeMaxForks
    conda "${genomeSimulationCondaEnvPath}"

    input:
    tuple val(genome_id), path(reconstruction_dir)
    path(media_ready)

    output:
    tuple val(genome_id), path("${genome_id}_environment_comparison"), emit: results

    script:
    def fvaFlag = genomeRunFva ? '--run-fva' : '--no-run-fva'
    def counterfactualFlag = genomeRunCounterfactuals ? '--run-counterfactuals' : '--no-run-counterfactuals'
    """
set -euo pipefail
export XDG_CACHE_HOME="\$PWD/.cache"
export MPLCONFIGDIR="\$PWD/.cache/matplotlib"
mkdir -p "\$XDG_CACHE_HOME" "\$MPLCONFIGDIR"
mapfile -t predicted_media < <(
  find -L "${reconstruction_dir}" -maxdepth 1 -type f -name '*-medium.csv' -print
)
if (( \${#predicted_media[@]} != 1 )); then
  printf 'Expected exactly one gapseq-predicted medium, found %d in %s\\n' \\
    "\${#predicted_media[@]}" "${reconstruction_dir}" >&2
  exit 1
fi
"\${CONDA_PREFIX}/bin/python" "${genomeCompareScriptPath}" \\
  --model "${reconstruction_dir}/${genome_id}.xml" \\
  --supplement-medium "\${predicted_media[0]}" \\
  --media-manifest "${biochemGapseqMediaDirAbs}/tables/gapseq_media_manifest.tsv" \\
  --media-root "${biochemGapseqMediaDirAbs}" \\
  --outdir "${genome_id}_environment_comparison" \\
  --solver "${genomeSolver}" \\
  --fraction-of-optimum ${genomeFractionOptimum} \\
  --flux-threshold ${genomeFluxThreshold} \\
  --supplement-max-flux ${genomeSupplementMaxFlux} \\
  ${fvaFlag} \\
  ${counterfactualFlag} \\
  --counterfactual-compounds "${genomeCounterfactualCompounds}"
rm -rf "${biochemGenomeModelingDirAbs}/comparisons/${genome_id}"
mkdir -p "${biochemGenomeModelingDirAbs}/comparisons/${genome_id}"
cp -a "${genome_id}_environment_comparison/." "${biochemGenomeModelingDirAbs}/comparisons/${genome_id}/"
"""
}

process GAPSEQ_COMBINE_RESULTS {
    cpus 1
    conda "${genomeSimulationCondaEnvPath}"

    input:
    path(result_dirs)

    output:
    path("gapseq_modeling.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemGenomeModelingDirAbs}/combined"
mkdir -p "${biochemGenomeModelingDirAbs}/plots"
"\${CONDA_PREFIX}/bin/python" "${genomeCombineScriptPath}" \\
  --results-root "${biochemGenomeModelingDirAbs}/comparisons" \\
  --reconstructions-root "${biochemGenomeModelingDirAbs}/reconstructions" \\
  --outdir "${biochemGenomeModelingDirAbs}/combined" \\
  ${genomeMetadataArg} \\
  --genomes-manifest "${genomeManifestPath}"
"\${CONDA_PREFIX}/bin/python" "${genomeNicheBiplotScriptPath}" \\
  --scores "${biochemPcaDirAbs}/tables/eigenvectors_scores.csv" \\
  --explained "${biochemPcaDirAbs}/tables/pca_explained_variance.csv" \\
  --hybrid-assignments "${biochemHybridDirAbs}/tables/hybrid_per_sample_joined.csv" \\
  --media-manifest "${biochemGapseqMediaDirAbs}/tables/gapseq_media_manifest.tsv" \\
  --growth "${biochemGenomeModelingDirAbs}/combined/combined_growth_comparison.tsv" \\
  --nutrient-importance "${biochemGenomeModelingDirAbs}/combined/combined_nutrient_importance_summary.tsv" \\
  --genome-metadata "${genomeMetadataPath}" \\
  ${genomeAbundanceArg} \\
  ${genomeTranscriptAbundanceArg} \\
  --agreement-permutations ${genomeAgreementPermutations} \\
  --agreement-bootstrap-iterations ${genomeAgreementBootstrapIterations} \\
  --agreement-random-state ${genomeAgreementRandomState} \\
  --expected-compartment-top-fraction ${genomeExpectedCompartmentTopFraction} \\
  --table-outdir "${biochemGenomeModelingDirAbs}/combined" \\
  --plot-outdir "${biochemGenomeModelingDirAbs}/plots"
[[ -s "${biochemGenomeModelingDirAbs}/combined/combined_results_inventory.tsv" ]]
[[ -s "${biochemGenomeModelingDirAbs}/combined/combined_reconstruction_qc.tsv" ]]
[[ -s "${biochemGenomeModelingDirAbs}/combined/combined_reconstruction_summary.tsv" ]]
if [[ -n "${genomeMetadataPath ?: ''}" && "${genomeRunCounterfactuals}" == "true" ]]; then
  [[ -s "${biochemGenomeModelingDirAbs}/combined/combined_counterfactual_taxonomic_summary.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/combined_counterfactual_lineage_summary.tsv" ]]
fi
[[ -s "${biochemGenomeModelingDirAbs}/combined/genome_predicted_niche_positions.tsv" ]]
[[ -s "${biochemGenomeModelingDirAbs}/combined/species_predicted_niche_positions.tsv" ]]
[[ -s "${biochemGenomeModelingDirAbs}/combined/species_representative_selection_audit.tsv" ]]
[[ -s "${biochemGenomeModelingDirAbs}/combined/species_representative_position_sensitivity.tsv" ]]
[[ -s "${biochemGenomeModelingDirAbs}/combined/genome_niche_recipe_weights.tsv" ]]
[[ -s "${biochemGenomeModelingDirAbs}/plots/genome_hybrid_compartment_niche_biplot.pdf" ]]
[[ -s "${biochemGenomeModelingDirAbs}/plots/genome_hybrid_compartment_niche_biplot.png" ]]
[[ -s "${biochemGenomeModelingDirAbs}/plots/genome_hybrid_compartment_niche_biplot.svg" ]]
if [[ -n "${genomeAbundancePath ?: ''}" ]]; then
  [[ -s "${biochemGenomeModelingDirAbs}/combined/genome_observed_abundance_positions.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/species_observed_abundance_positions.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/species_predicted_observed_positions.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/species_predicted_observed_agreement.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/genome_observed_hybrid_affinity.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/species_observed_hybrid_affinity.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/genome_expected_compartment_performance.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/species_expected_compartment_performance.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/genome_hybrid_recipe_growth_ranks.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/species_hybrid_recipe_growth_ranks.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/expected_compartment_performance_summary.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/lineage_expected_compartment_enrichment.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/plots/genome_hybrid_compartment_predicted_observed_biplot.pdf" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/plots/genome_hybrid_compartment_predicted_observed_biplot.png" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/plots/genome_hybrid_compartment_predicted_observed_biplot.svg" ]]
fi
if [[ -n "${genomeTranscriptAbundancePath ?: ''}" ]]; then
  [[ -s "${biochemGenomeModelingDirAbs}/combined/genome_observed_expression_positions.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/species_observed_expression_positions.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/species_predicted_observed_expression_positions.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/combined/species_predicted_observed_expression_agreement.tsv" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/plots/genome_hybrid_compartment_predicted_observed_expression_biplot.pdf" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/plots/genome_hybrid_compartment_predicted_observed_expression_biplot.png" ]]
  [[ -s "${biochemGenomeModelingDirAbs}/plots/genome_hybrid_compartment_predicted_observed_expression_biplot.svg" ]]
fi
touch gapseq_modeling.done
"""
}

process BIOCHEM_GAPSEQ_MEDIA {
    cpus 1
    conda "${biochemGapseqMediaCondaEnvPath}"

    input:
    path(prev_done)

    output:
    path("biochem_gapseq_media.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemGapseqMediaDirAbs}"
"\${CONDA_PREFIX}/bin/python" "${biochemGapseqMediaScriptPath}" \\
  --matrix "${biochemPcaDirAbs}/tables/matrix_cleaned_with_sparse.csv" \\
  --cruise-groups "${biochemEofStatesDirAbs}/tables/gmm_selected_assignments.tsv" \\
  --o2-assignments "${biochemO2DirAbs}/tables/o2_compartments_assignments_smoothed.csv" \\
  --gmm-assignments "${biochemGmmDirAbs}/tables/compartments_assignments_smoothed.csv" \\
  --hybrid-assignments "${biochemHybridDirAbs}/tables/compartments_assignments_hybrid.csv" \\
  --nutrients "${gapseqNutrientsPath}" \\
  ${gapseqTemplateArg} \\
  --mapping "${gapseqMappingPath}" \\
  --outdir "${biochemGapseqMediaDirAbs}" \\
  --method "${gapseqTransformMethod}" \\
  --min-effective-n ${gapseqMinEffectiveN} \\
  --detection-limit ${gapseqDetectionLimit} \\
  --depth-baseline-m "${gapseqDepthBaselineArg}"
[[ -s "${biochemGapseqMediaDirAbs}/tables/gapseq_media_manifest.tsv" ]] || { echo "Missing gapseq media manifest" >&2; exit 1; }
touch biochem_gapseq_media.done
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
"\${CONDA_PREFIX}/bin/python" "${biochemWithinGmmScriptPath}" \\
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

process MASTER_SUMMARY {
    cpus 1
    conda "${biochemMasterSummaryCondaEnvPath}"

    input:
    path(prev_done)
    path(continuous_done)

    output:
    path("master_summary.done"), emit: done

    script:
    """
set -euo pipefail
mkdir -p "${biochemMasterSummaryDirAbs}"
"\${CONDA_PREFIX}/bin/python" "${biochemMasterSummaryScriptPath}" \
  --input-root "${biochemOutputRoot}" \
  --output-dir "${biochemMasterSummaryDirAbs}"
[[ -f "${biochemMasterSummaryDirAbs}/basin_run_overview.tsv" ]] || { echo "Missing BASINS master summary" >&2; exit 1; }
touch master_summary.done
"""
}
