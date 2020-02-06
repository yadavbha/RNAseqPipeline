#!/bin/bash -l

##Set parameters below

#set project name
projectName='project_2001209'

##set directories
dataDir='/scratch/project_2001209/RawNGSdata/RNAseqFastqData/'
genomeDir='/scratch/project_2001209/GenomeDir'
refGeneDir='/scratch/project_2001209/RefHuGenome'
outDir='/scratch/project_2001209/RNAseqResults'
workDir='/scratch/project_2001209/.work'

## Path to different files
fastaFile="$refGeneDir/GRCh38_13/*.fa"
gtfFile="$refGeneDir/GRCh38_13/*.gtf"
fastaERCC92="$refGeneDir/ERCC/*.fa"
gtfERCC="$refGeneDir/ERCC/*.gtf"

##set other parameters
trimReads=true  #Trim Adaptors and filter phred 20
sjdbOverhang=49 #default value is 149. Because our nucleotide length is 50i
strandedness='unstranded'  ## unstranded, forward_stranded, reverse_stranded
addERCC92=true

##set if
downloadFasta=false
downloadGTF=false

email='bhagwan.yadav@gmail.com'  #Don't work because of server blocking
##End Parameters
saveTrimmed=false
if [ $trimReads ]
then
   saveTrimmed=$trimReads
fi

##clean workdir

################################## Don't Change anything down unless it is required #######################

ulimit -c unlimited
ulimit -s unlimited

export HOME='/projappl/project_2001209'
export NXF_HOME='/projappl/project_2001209/Softwares/RNAseqPipeline'
export PICARD_HOME='/projappl/project_2001209/Softwares/conda/condawockstrom/share/picard-2.21.4-0/'
export NXF_OPTS='-Xms512M -Xmx10G'
export LANGUAGE='en_US.UTF-8'
export LC_ALL='en_US.UTF-8'
export LANG='en_US.UTF-8'
export LC_CTYPE='en_US.UTF-8'
export R_LIBS_USER='/projappl/project_2001209/Softwares/Rlib'
module load gcc/9.1.0
module load r-env

script_path='/projappl/project_2001209/Softwares/RNAseqPipeline/main.nf'

if [ $addERCC92 ]
then
       fastaFile="$refGeneDir/GRCh38_13_with_ERCC92_spikes/*.fa"
       gtfFile="$refGeneDir/GRCh38_13_with_ERCC92_spikes/*.gtf"
fi


if [[ $1 == 'help' || $1 == '--help' || $1 == '-h' ]]
then
	cmd="nextflow run $script_path -name $name --help"
	eval $cmd

elif [ $1 == 'resume' ]
then
     if [ -f "$genomeDir/SAindex" ]
     then
        cmd="nextflow run $script_path -profile standard \
        --project $projectName \
        --genomedir $genomeDir \
        --reads '${dataDir}*_R{1,2}_*.fastq.gz' \
	--strandedness $strandedness \
	--trim $trimReads \
        --outdir $outDir \
        --workdir $workDir \
        --ercc $addERCC92 \
	--download_fasta $downloadFasta \
	--download_gtf $downloadGTF \
	--saveTrimmed $saveTrimmed \
        --email $email \
        -resume"
        echo "Starting nextflow... Command:"
        echo $cmd
        echo "-----"
        eval $cmd
     else
	cmd="nextflow run $script_path -profile standard \
	--project $projectName \
	--fasta $fastaFile \
        --gtf $gtfFile \
        --sjdbOverhang $sjdbOverhang \
	--reads '${dataDir}*_R{1,2}_*.fastq.gz' \
        --strandedness $strandedness \
        --trim $trimReads \
	--outdir $outDir \
	--workdir $workDir \
        --ercc $addERCC92 \
        --download_fasta $downloadFasta \
        --download_gtf $downloadGTF \
        --saveTrimmed $saveTrimmed \
	--email $email \
	-resume"
	echo "Starting nextflow... Command:"
	echo $cmd
	echo "-----"
	eval $cmd
     fi

elif [ $1 == 'run' ]

#rm -r $workDir/*  #clean workdir

then
     if [ -f "$genomeDir/SAindex" ]
     then
        cmd="nextflow run $script_path -profile standard \
        --project $projectName \
        --genomedir $genomeDir \
        --reads '${dataDir}*_R{1,2}_*.fastq.gz' \
        --strandedness $strandedness \
        --trim $trimReads \
        --outdir $outDir \
	--workdir $workDir \
        --ercc $addERCC92 \
        --download_fasta $downloadFasta \
        --download_gtf $downloadGTF \
        --saveTrimmed $saveTrimmed \
        --email $email"
        echo "Starting nextflow... Command:"
        echo $cmd
        echo "-----"
        eval $cmd
     else
        cmd="nextflow run $script_path -profile standard \
        --project $projectName \
        --fasta $fastaFile \
        --gtf $gtfFile \
        --sjdbOverhang $sjdbOverhang \
        --reads '${dataDir}*_R{1,2}_*.fastq.gz' \
        --strandedness $strandedness \
        --trim $trimReads \
        --outdir $outDir \
        --workdir $workDir \
        --ercc $addERCC92 \
        --download_fasta $downloadFasta \
        --download_gtf $downloadGTF \
        --saveTrimmed $saveTrimmed \
        --email $email"
        echo "Starting nextflow... Command:"
        echo $cmd
        echo "-----"
        eval $cmd
     fi 
else
	ech "Please enter correct parameter [help/run/resume]"

fi

rm -r $outDir/RNAseq_timeline.html.*
rm -r $outDir/RNAseq_trace.txt.*

#st1=$(tail -1 $NXF_HOME/.nextflow.log|grep -cim1 'Execution complete')

#if [ $st1 -eq 1 ]
#then
#	if [ -d "$outDir/multiQC" ]; then
#		rm -r $outDir/multiQC
#		mkdir $outDir/multiQC
#	fi
#	cd $outDir
#	multiqc $outDir -f --filename "multiqc_report" -o "$outDir/multiQC" -c "$NXF_HOME/conf/" --cl_config "multiqc_config.yaml"
#	$NXF_HOME/bin/markdown_to_html.r "$NXF_HOME/docs/output.md" "$outDir/Documentation/results_description.html" $R_LIBS_USER

#fi



