#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                            R N A S E Q   P I P E L I N E
========================================================================================
 RNA-Seq pipeline for in-house use to run on Puhti supercomputing facility by csc.fi
 
 #### Contact Person
 Bhagwan Yadav @Bhagwan <bhagwan.yadav@helsinki.fi>
 
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     RNAseq : RNA-Seq pipeline
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    RNAseqPipeline [help / run / resume]

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference

    Options:
      --singleEnd                   Specifies that the input is single end reads
    Strandedness:
      --forward_stranded            The library is forward stranded
      --reverse_stranded            The library is reverse stranded
      --unstranded                  The default behaviour

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --star_index                  Path to STAR index
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF file
      --bed12                       Path to bed12 file
      --downloadFasta               If no STAR / Fasta reference is supplied, a URL can be supplied to download a Fasta file at the start of the pipeline.
      --downloadGTF                 If no GTF reference is supplied, a URL can be supplied to download a Fasta file at the start of the pipeline.
      --saveReference               Save the generated reference files the the Results directory.
      --saveAlignedIntermediates    Save the BAM files from the Aligment step  - not done by default

    Trimming options
      --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2 [int]   Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed

    Presets:
      --pico                        Sets trimming and standedness settings for the SMARTer Stranded Total RNA-Seq Kit - Pico Input kit. Equivalent to: --forward_stranded --clip_r1 3 --three_prime_clip_r2 3

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --sampleLevel                 Used to turn of the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples
      --rlocation                   Location to save R-libraries used in the pipeline. Default value is ~/R/nxtflow_libs/
      --clusterOptions              Extra SLURM options, used in conjunction with base.config
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '1.0'

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

// Configurable variables
params.name = false
params.project = false
params.ercc = false    //true or false
params.genomedir = false
if (params.genomedir){
  if (!file("${params.genomedir}/Genome").isEmpty()){
     params.genome = "${params.genomedir}/Genome"
  } else {
     params.genome = false
  }
  params.fasta = "${params.genomedir}/*.fa"
  if (!file("${params.genomedir}/SAindex").isEmpty()){
     params.star_index = "${params.genomedir}/SAindex"
  } else {
     params.star_index = false
  }

  params.gtf = "${params.genomedir}/*.gtf"
  if (!file("${params.genomedir}/*.bed").isEmpty()){
     params.bed12 = "${params.genomedir}/*.bed"
  } else {
     params.bed12 = false
	}
}
else {
  params.fasta = false
  params.gtf = false
  params.bed12 = false
  params.star_index = false
  params.genome = false
}
params.sjdbOverhang = 149
if (params.strandedness == 'unstranded'){
	params.unstranded = true
	params.forward_stranded = false
	params.reverse_stranded = false
}
if (params.strandedness == 'forward_stranded'){
        params.unstranded = false
        params.forward_stranded = true
        params.reverse_stranded = false
}
if (params.strandedness == 'reverse_stranded'){
	params.unstranded = false
        params.forward_stranded = false
        params.reverse_stranded = true
}
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.download_fasta = false
params.download_gtf = false
if (params.star_index){
params.saveReference = false}
else {
params.saveReference = true}
params.saveTrimmed = true
params.saveAlignedIntermediates = true
params.reads = false
params.outdir = './results'
params.workdir = '/.work'
params.email = true
params.plaintext_email = true

// R library locations
params.rlocation = false
if (params.rlocation){
    nxtflow_libs = file(params.rlocation)
    nxtflow_libs.mkdirs()
}

multiqc_config = file(params.multiqc_config)
params.sampleLevel = false

// Custom trimming options
params.trim = true
if (params.trim){
  params.clip_r1 = 0
  params.clip_r2 = 0
  params.three_prime_clip_r1 = 0
  params.three_prime_clip_r2 = 0

// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
}

forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded

// Preset trimming options
params.pico = false
if (params.pico){
  clip_r1 = 3
  clip_r2 = 0
  three_prime_clip_r1 = 0
  three_prime_clip_r2 = 3
  forward_stranded = true
  reverse_stranded = false
  unstranded = false
}

if(unstranded){
	Strandedness = 0
} else if(forward_stranded){
	Strandedness = 1
}else if(reverse_stranded){
	Strandedness = 2
}

// Choose aligner
params.aligner = 'star'
if ( params.aligner != 'star' ){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star' "
}

// Validate inputs
if ( params.genomedir ){
    genomeDir = Channel
        .fromPath(params.genomedir)
        .ifEmpty { exit 1, "GenomeDir not found: ${params.genomedir}" }
}

if( params.star_index && params.aligner == 'star' ){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}

if ( params.genome ){
    genome = Channel
        .fromPath(params.genome)
        .ifEmpty { exit 1, "Genome not found: ${params.genome}" }
}

if ( params.fasta ){
    fasta = file(params.fasta)
    if( fasta.isEmpty() ) exit 1, "Fasta file not found: ${params.fasta}"

}

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_makeBED12; gtf_star; gtf_dupradar; gtf_featureCounts; gtf_stringtieFPKM }
}
else if ( !params.download_gtf ){
    exit 1, "No GTF annotation specified!"
}
if( params.bed12 ){
    bed12 = Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
        .into {bed_rseqc; bed_genebody_coverage}
}
if( workflow.profile == 'standard' && !params.project ) exit 1, "No project ID found! Use --project"

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input read files
 */
params.singleEnd = false
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2)
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_trimming; untrimmed_reads}






// Header log info
log.info "========================================="
log.info " RNAseq Pipeline"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']       = params.genome
summary['GenomeDir']       = params.genomedir
if( params.pico ) summary['Library Prep'] = "SMARTer Stranded Total RNA-Seq Kit - Pico Input"
summary['Strandedness'] = ( unstranded ? 'None' : forward_stranded ? 'Forward' : reverse_stranded ? 'Reverse' : 'None' )
if (params.trim){
	summary['Trim R1'] = clip_r1
	summary['Trim R2'] = clip_r2
	summary["Trim 3' R1"] = three_prime_clip_r1
	summary["Trim 3' R2"] = three_prime_clip_r2
}
if(params.aligner == 'star'){
    summary['Aligner'] = "STAR"
    if(params.star_index)          summary['STAR Index']   = params.star_index
    else if(params.fasta)          summary['Fasta Ref']    = params.fasta
    else if(params.download_fasta) summary['Fasta URL']    = params.download_fasta
}
if(params.gtf)                 summary['GTF Annotation']  = params.gtf
else if(params.download_gtf)   summary['GTF URL']         = params.download_gtf
if(params.bed12)               summary['BED Annotation']  = params.bed12
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Output dir']     = params.outdir
summary['Working dir']    = workDir
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['R libraries']    = params.rlocation
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = (workflow.profile == 'standard')
if(params.project) summary['Project'] = params.project
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
 * PREPROCESSING - Download Fasta
 */
if((params.aligner == 'star' && !params.star_index) && !params.fasta && params.download_fasta){
    process downloadFASTA {
        tag "${params.download_fasta}"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        output:
        file "*.{fa,fasta}" into fasta

        script:
        """
        curl -O -L ${params.download_fasta}
        if [ -f *.tar.gz ]; then
            tar xzf *.tar.gz
        elif [ -f *.gz ]; then
            gzip -d *.gz
        fi
        """
    }
}
/*
 * PREPROCESSING - Download GTF
 */
if(!params.gtf && params.download_gtf){
    process downloadGTF {
        tag "${params.download_gtf}"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        output:
        file "*.gtf" into gtf_makeSTARindex, gtf_makeHisatSplicesites, gtf_makeHISATindex, gtf_makeBED12, gtf_star, gtf_dupradar, gtf_featureCounts, gtf_stringtieFPKM

        script:
        """
        curl -O -L ${params.download_gtf}
        if [ -f *.tar.gz ]; then
            tar xzf *.tar.gz
        elif [ -f *.gz ]; then
            gzip -d *.gz
        fi
        """
    }
}



/*
 * PREPROCESSING - Build STAR index
 */
if(params.aligner == 'star' && !params.star_index){
    process makeSTARindex {
        tag fasta
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta
        file gtf from gtf_makeSTARindex

        output:
        file "star" into genomeDir

        script:
        """
        ulimit -c unlimited
        ulimit -s unlimited
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN 16 \\
            --sjdbGTFfile $gtf \\
            --sjdbOverhang "${params.sjdbOverhang}" \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta
        """
    }

}




/*
 * PREPROCESSING - Build BED12 file
 */
if(!params.bed12){
    process makeBED12 {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeBED12

        output:
        file "${gtf.baseName}.bed" into bed_rseqc, bed_genebody_coverage

        script: // This script is bundled with the pipeline, in RNAseqPipeline/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results
    file '.command.out' into fastqc_stdout

    script:
    """
    fastqc -t 20 -q $reads
    fastqc --version
    """
}


/*
 * STEP 2 - Trim Galore!
 */
if (params.trim){
 process trim_galore {
     tag "$prefix"
     publishDir "${params.outdir}/trim_galore", mode: 'copy',
         saveAs: {filename ->
             if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
             else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
             else params.saveTrimmed ? filename : null
         }

     input:
     set val(prefix), file(reads) from read_files_trimming

     output:
     file "*.{fq,fastq}.gz" into trimmed_reads
     file "*trimming_report.txt" into trimgalore_results, trimgalore_logs
     file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports


     script:
     c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
     c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
     tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
     tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
     if (params.singleEnd) {
         """
         trim_galore --cores 16 --fastqc --gzip $c_r1 $tpc_r1 $reads
         """
     } else {
         """
         trim_galore --cores 16 --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
         """
     }
 }
}



/*
 * STEP 3 - align with STAR
 */
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    if(percent_aligned.toFloat() <= '5'.toFloat() ){
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}




if(params.aligner == 'star' && params.trim){
    process star_v1 {
        tag "$prefix"
        publishDir "${params.outdir}/STAR", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bam") == -1) "logs/$filename"
                else params.saveAlignedIntermediates ? filename : null
            }
        
        input:
        file reads from trimmed_reads
        file genome from genomeDir.collect()
        file gtf from gtf_star.collect()

        output:
        set file("*Log.final.out"), file ("*.{bam,bai}") into star_aligned
        file "*.out" into alignment_logs
        file "*SJ.out.tab"
        file "*Log.out" into star_log

        script:
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_001)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        ulimit -c unlimited
        ulimit -s unlimited
        STAR --genomeDir $genome \\
            --genomeLoad NoSharedMemory \\
            --sjdbGTFfile $gtf \\
            --readFilesIn $reads  \\
            --runThreadN 20 \\
            --twopassMode Basic \\
            --outWigType bedGraph \\
            --outSAMtype BAM SortedByCoordinate \\
            --readFilesCommand zcat \\
            --runDirPerm All_RWX \\
            --outFileNamePrefix "${prefix}_" \\
            --outReadsUnmapped Fastx \\
            --outSAMattributes All
        """
   }

// Filter removes all 'aligned' channels that fail the check
    star_aligned
        .filter { logs, bams -> check_log(logs) }
        .flatMap {  logs, bams -> bams }
    .into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM; bam_geneBodyCoverage }
}






if(params.aligner == 'star' && !params.trim){
    process star_v2 {
        tag "$prefix"
        publishDir "${params.outdir}/STAR", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bam") == -1) "logs/$filename"
                else params.saveAlignedIntermediates ? filename : null
            }

        input:
        set val(prefix), file(reads) from untrimmed_reads
        file genome from genomeDir.collect()
        file gtf from gtf_star.collect()
        
        output:
        set file("*Log.final.out"), file("*.{bam,bai}") into star_aligned
        file "*.out" into alignment_logs
        file "*SJ.out.tab"
        file "*Log.out" into star_log
       
        script:
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_001)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        ulimit -c unlimited
        ulimit -s unlimited
        STAR --genomeDir $genome \\
            --genomeLoad NoSharedMemory \\
            --sjdbGTFfile $gtf \\
            --readFilesIn $reads  \\
            --runThreadN 20 \\
            --twopassMode Basic \\
            --outWigType bedGraph \\
            --outSAMtype BAM SortedByCoordinate \\
            --readFilesCommand zcat \\
            --runDirPerm All_RWX \\
            --outFileNamePrefix "${prefix}_" \\
            --outReadsUnmapped Fastx \\
            --outSAMattributes All
        """

   }

// Filter removes all 'aligned' channels that fail the check
    star_aligned
        .filter { logs, bams -> check_log(logs) }
        .flatMap {  logs, bams -> bams }
    .into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM; bam_geneBodyCoverage }
}




/*
 * STEP 4 - RSeQC analysis
 */
process rseqc {
    tag "${bam_rseqc.baseName - '.sorted'}"
    publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else if (filename.indexOf("log.txt") > -1) false
            else "$filename"
        }

    input:
    file bam_rseqc
    file bed12 from bed_rseqc.collect()

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    def strandRule = ''
    if (forward_stranded && !unstranded){
        strandRule = params.singleEnd ? '-d ++,--' : '-d 1++,1--,2+-,2-+'
    } else if (reverse_stranded && !unstranded){
        strandRule = params.singleEnd ? '-d +-,-+' : '-d 1+-,1-+,2++,2--'
    }
    """
    samtools index $bam_rseqc
    infer_experiment.py -i $bam_rseqc -s 20000 -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
    junction_annotation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    bam_stat.py -i $bam_rseqc > ${bam_rseqc.baseName}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12 > ${bam_rseqc.baseName}.junction_annotation_log.txt
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc.baseName}.read_duplication
    echo "Filename $bam_rseqc RseQC version: "\$(read_duplication.py --version)
    """
}

/*
 * Step 4.1 Rseqc genebody_coverage
 */
process genebody_coverage {
    tag "${bam_geneBodyCoverage.baseName - '.sorted'}" 
       publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
            else if (filename.indexOf("geneBodyCoverage.r") > 0)                "geneBodyCoverage/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.txt") > 0)              "geneBodyCoverage/data/$filename"
            else "$filename"
        }

    input:
    file bam_geneBodyCoverage
    file bed12 from bed_genebody_coverage.collect()
    
    output:
    file "*.{txt,pdf,r,xls}" into genebody_coverage_results
    
    script:
    """
    cat <(samtools view -H ${bam_geneBodyCoverage}) \\
        <(samtools view ${bam_geneBodyCoverage} | shuf -n 10000) \\
        | samtools sort - -o ${bam_geneBodyCoverage.baseName}_subsamp_sorted.bam 
    samtools index ${bam_geneBodyCoverage.baseName}_subsamp_sorted.bam
    geneBody_coverage.py -i ${bam_geneBodyCoverage.baseName}_subsamp_sorted.bam -o ${bam_geneBodyCoverage.baseName}.rseqc -r $bed12
    """
}

/*
 * STEP 5 - preseq analysis
 */
process preseq {
    tag "${bam_preseq.baseName - '.sorted'}"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    input:
    file bam_preseq

    output:
    file "${bam_preseq.baseName}.ccurve.*" into preseq_results
    file '.command.log' into preseq_stdout

    script:
    """
    preseq lc_extrap -v -B $bam_preseq -o ${bam_preseq.baseName}.ccurve.txt
    ##plot_complexity_curves.py ${bam_preseq.baseName}.ccurve.txt -o ${bam_preseq.baseName}.ccurve.pdf
    echo "File name: $bam_preseq  preseq version: "\$(preseq)
    """
}


/*
 * STEP 6 Mark duplicates
 */
process markDuplicates {
    tag "${bam_markduplicates.baseName - '.sorted'}"
    publishDir "${params.outdir}/markDuplicates", mode: 'copy',
        saveAs: {filename -> filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"}

    input:
    file bam_markduplicates

    output:
    file "${bam_markduplicates.baseName}.markDups.bam" into bam_md
    file "${bam_markduplicates.baseName}.markDups_metrics.txt" into picard_results
    file "${bam_markduplicates.baseName}.bam.bai"
    file '.command.log' into markDuplicates_stdout

    script:
    """
    ulimit -c unlimited
    java -Xmx60g -Djava.io.tmpdir='/scratch/project_2001209/.work' -jar \${PICARD_HOME}/picard.jar MarkDuplicates \\
        INPUT=$bam_markduplicates \\
        OUTPUT=${bam_markduplicates.baseName}.markDups.bam \\
        METRICS_FILE=${bam_markduplicates.baseName}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT

    # Print version number to standard out
      echo "File name: $bam_markduplicates Picard version "\$(java -Xmx60g -jar \${PICARD_HOME}/picard.jar  MarkDuplicates --version 2>&1)
      samtools index $bam_markduplicates
         """
}


/*
 * STEP 7 - dupRadar
 */
process dupradar {
    tag "${bam_md.baseName - '.sorted.markDups'}"
    publishDir "${params.outdir}/dupradar", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
            else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
            else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
            else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
            else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
            else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
            else "$filename"
        }

    input:
    file bam_md
    file gtf from gtf_dupradar.collect()

    output:
    file "*.{pdf,txt}" into dupradar_results
    file '.command.log' into dupradar_stdout

    script: // This script is bundled with the pipeline, in RNAseqPipeline/bin/
    def paired = params.singleEnd ? 'FALSE' :  'TRUE'
    def rlocation = params.rlocation ?: ''
    def strandedVal = Strandedness ?: 0
    """
    module load r-env
    export R_LIBS_USER='/projappl/project_2001209/Softwares/Rlib'
    dupRadar.r $bam_md $gtf $strandedVal $paired $rlocation
    """
}


/*
 * STEP 8 Feature counts
 */

if (params.ercc){
process featureCountsV1 {
    tag "${bam_featurecounts.baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    input:
    file bam_featurecounts
    file gtf from gtf_featureCounts.collect()

    output:
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
    file '.command.log' into featurecounts_stdout

    script:
    def featureCounts_direction = 0
    if (forward_stranded && !unstranded) {
        featureCounts_direction = 1
    } else if (reverse_stranded && !unstranded){
        featureCounts_direction = 2
    }
    """
    ulimit -c unlimited
    featureCounts -T 20 -a $gtf -g gene_id -o ${bam_featurecounts.baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
    """
}
}
else {
process featureCountsV2 {
    tag "${bam_featurecounts.baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_biotype_counts.txt") > 0) "biotype_counts/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    input:
    file bam_featurecounts
    file gtf from gtf_featureCounts.collect()

    output:
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
    file "${bam_featurecounts.baseName}_biotype_counts.txt" into featureCounts_biotype
    file '.command.log' into featurecounts_stdout

    script:
    def featureCounts_direction = 0
    if (forward_stranded && !unstranded) {
        featureCounts_direction = 1
    } else if (reverse_stranded && !unstranded){
        featureCounts_direction = 2
    }
    """
    ulimit -c unlimited
    featureCounts -T 20 -a $gtf -g gene_id -o ${bam_featurecounts.baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
    featureCounts -T 20 -a $gtf -g gene_biotype -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
    cut -f 1,7 ${bam_featurecounts.baseName}_biotype.featureCounts.txt > ${bam_featurecounts.baseName}_biotype_counts.txt
    """
}
}

/*
 * STEP 9 - Merge featurecounts
 */
process merge_featureCounts {
    tag "${input_files[0].baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy'

    input:
    file input_files from featureCounts_to_merge.collect()

    output:
    file 'merged_gene_counts.txt'

    script:
    """
    merge_featurecounts.py -o merged_gene_counts.txt -i $input_files
    """
}


/*
 * STEP 10 - stringtie FPKM
 */
process stringtieFPKM {
    tag "${bam_stringtieFPKM.baseName - '.sorted'}"
    publishDir "${params.outdir}/stringtieFPKM", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
            else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
            else "$filename"
        }

    input:
    file bam_stringtieFPKM
    file gtf from gtf_stringtieFPKM.collect()

    output:
    file "${bam_stringtieFPKM.baseName}_transcripts.gtf"
    file "${bam_stringtieFPKM.baseName}.gene_abund.txt"
    file "${bam_stringtieFPKM}.cov_refs.gtf"
    file ".command.log" into stringtie_log, stringtie_stdout

    script:
    def st_direction = ''
    if (forward_stranded && !unstranded){
        st_direction = "--fr"
    } else if (reverse_stranded && !unstranded){
        st_direction = "--rf"
    }
    """
    stringtie $bam_stringtieFPKM \\
        $st_direction \\
        -o ${bam_stringtieFPKM.baseName}_transcripts.gtf \\
        -v \\
        -G $gtf \\
        -A ${bam_stringtieFPKM.baseName}.gene_abund.txt \\
        -C ${bam_stringtieFPKM}.cov_refs.gtf \\
        -e \\
        -b ${bam_stringtieFPKM.baseName}_ballgown

    echo "File name: $bam_stringtieFPKM Stringtie version "\$(stringtie --version)
    """
}
def num_bams
bam_count.count().subscribe{ num_bams = it }

/*
 * STEP 11 - edgeR MDS and heatmap
 */
process sample_correlation {
    tag "${input_files[0].toString() - '.sorted_gene.featureCounts.txt' - 'Aligned'}"
    publishDir "${params.outdir}/sample_correlation", mode: 'copy'

    input:
    file input_files from geneCounts.collect()
    bam_count

    output:
    file "*.{txt,pdf}" into sample_correlation_results

    when:
    num_bams > 2 && (!params.sampleLevel)

    script: // This script is bundled with the pipeline, in RNAseqPipeline/bin/
    def rlocation = params.rlocation ?: ''
    """
    edgeR_heatmap_MDS.r "rlocation=$rlocation" $input_files
    """
}

/*
 * Parse software version numbers
 */

if(params.trim){
software_versions = [
  'FastQC': null, 'Trim Galore!': null, 'Star': null, 'StringTie': null,
  'Preseq': null, 'featureCounts': null, 'dupRadar': null, 'Picard MarkDuplicates': null,
  'Nextflow': "v$workflow.nextflow.version"
]
process get_software_V1 {
    cache false
    executor 'local'

    input:
    val fastqc from fastqc_stdout.collect()
    val trim_galore from trimgalore_logs.collect()
    val star from star_log.collect()
    val stringtie from stringtie_stdout.collect()
    val preseq from preseq_stdout.collect()
    val featurecounts from featurecounts_stdout.collect()
    val dupradar from dupradar_stdout.collect()
    val markDuplicates from markDuplicates_stdout.collect()

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    exec:
    software_versions['FastQC'] = fastqc[0].getText().find(/FastQC v(\S+)/) { match, version -> "v$version" }
    software_versions['Trim Galore!'] = trim_galore[0].getText().find(/Trim Galore version: (\S+)/) {match, version -> "v$version"}
    if(params.aligner == 'star') software_versions['Star'] = star[0].getText().find(/STAR v.+(\d+\.\d+\.\d+)/) { match, version -> "v$version" }
    software_versions['StringTie'] = stringtie[0].getText().find(/StringTie (\S+)/) { match, version -> "v"+version.replaceAll(/\.$/, "") }
    software_versions['Preseq'] = preseq[0].getText().find(/Version: (\S+)/) { match, version -> "v$version" }
    software_versions['featureCounts'] = featurecounts[0].getText().find(/\s+v([\.\d]+)/) {match, version -> "v$version"}
    software_versions['dupRadar'] = dupradar[0].getText().find(/dupRadar\_(\S+)/) {match, version -> "v$version"}
    software_versions['Picard MarkDuplicates'] = markDuplicates[0].getText().find(/Picard version ([\d\.]+)/) {match, version -> "v$version"}

    def sw_yaml_file = task.workDir.resolve('software_versions_mqc.yaml')
    sw_yaml_file.text  = """
    id: 'rnaseqpipeline'
    section_name: 'RNAseq Pipeline'
    plot_type: 'html'
    description: 'are collected at run time from the software output.'
    data: |
        <dl class=\"dl-horizontal\">
${software_versions.collect{ k,v -> "            <dt>$k</dt><dd>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</dd>" }.join("\n")}
        </dl>
    """.stripIndent()
  }
}

if(!params.trim){
software_versions = [
  'FastQC': null, 'Star': null, 'StringTie': null,
  'Preseq': null, 'featureCounts': null, 'dupRadar': null, 'Picard MarkDuplicates': null,
  'Nextflow': "v$workflow.nextflow.version"
]

process get_software_V2 {
    cache false
    executor 'local'

    input:
    val fastqc from fastqc_stdout.collect()
    val star from star_log.collect()
    val stringtie from stringtie_stdout.collect()
    val preseq from preseq_stdout.collect()
    val featurecounts from featurecounts_stdout.collect()
    val dupradar from dupradar_stdout.collect()
    val markDuplicates from markDuplicates_stdout.collect()

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    exec:
    software_versions['FastQC'] = fastqc[0].getText().find(/FastQC v(\S+)/) { match, version -> "v$version" }
    if(params.aligner == 'star') software_versions['Star'] = star[0].getText().find(/STAR v.+(\d+\.\d+\.\d+?)/) { match, version -> "v$version" }
    software_versions['StringTie'] = stringtie[0].getText().find(/StringTie (\S+)/) { match, version -> "v"+version.replaceAll(/\.$/, "") }
    software_versions['Preseq'] = preseq[0].getText().find(/Version: (\S+)/) { match, version -> "v$version" }
    software_versions['featureCounts'] = featurecounts[0].getText().find(/\s+v([\.\d]+)/) {match, version -> "v$version"}
    software_versions['dupRadar'] = dupradar[0].getText().find(/dupRadar\_(\S+)/) {match, version -> "v$version"}
    software_versions['Picard MarkDuplicates'] = markDuplicates[0].getText().find(/Picard version ([\d\.]+)/) {match, version -> "v$version"}

    def sw_yaml_file = task.workDir.resolve('software_versions_mqc.yaml')
    sw_yaml_file.text  = """
    id: 'rnaseqpipeline'
    section_name: 'RNAseq Pipeline'
    plot_type: 'html'
    description: 'are collected at run time from the software output.'
    data: |
        <dl class=\"dl-horizontal\">
${software_versions.collect{ k,v -> "            <dt>$k</dt><dd>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</dd>" }.join("\n")}
        </dl>
    """.stripIndent()
  }
}


/*
 *STEP 12 MultiQC
 */

process multiqc {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    file ('software_versions/*') from software_versions_yaml 

    output:
    file "*_report.*" into multiqc_report
    file "*_data"
    file '.command.err' into multiqc_stderr
    val prefix into multiqc_prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename ${params.outdir}
    """
}


/*
process multiqc {
    executor 'local'

    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('trimgalore/*') from trimgalore_results.collect().ifEmpty([])
    file ('alignment/*') from alignment_logs.collect().ifEmpty([])
    file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    file ('rseqc/*') from genebody_coverage_results.collect().ifEmpty([])
    file ('preseq/*') from preseq_results.collect().ifEmpty([])
    file ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    file ('featureCounts/*') from featureCounts_logs.collect().ifEmpty([])
    file ('stringtie/*') from stringtie_log.collect().ifEmpty([])
    file ('sample_correlation_results/*') from sample_correlation_results.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    file '.command.err' into multiqc_stderr
    val prefix into multiqc_prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}
*/
multiqc_stderr.subscribe { stderr ->
  software_versions['MultiQC'] = stderr.getText().find(/This is MultiQC v(\S+)/) { match, version -> "v$version" }
}






/*
 * STEP 13 - Output Description HTML
 */

process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    val prefix from multiqc_prefix

    output:
    file "results_description.html"

    script:
    def rlocation = params.rlocation ?: ''
    """
    markdown_to_html.r $baseDir/docs/output.md results_description.html $rlocation
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[RNAseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[RNAseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['software_versions'] = software_versions
    email_fields['software_versions']['Nextflow Build'] = workflow.nextflow.build
    email_fields['software_versions']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp
    
    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    if (params.email) {
        try {
          // Try to send HTML e-mail using sendmail
            def html_email = "To: $params.email\nSubject: $subject\nMime-Version: 1.0\nContent-Type: text/html\n\n$email_html";
            
            [ 'sendmail', '-t' ].execute() << html_email
            log.debug "Sent summary e-mail using sendmail"
        } catch (all) {
         // Catch failures and try with plaintext
            [ 'mail', '-s', subject, params.email ].execute() << email_txt
            log.debug "Sendmail failed, failing back to sending summary e-mail using mail"
        }
    log.info "Sent summary e-mail to $params.email"
   }


/*
    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "RNAseq summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "RNAseq summary e-mail to $params.email (mail)"
        }
    }
*/
    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "RNAseq Pipeline Complete"

}
