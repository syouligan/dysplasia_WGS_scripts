#!/bin/bash

###Filter VCF files to high quality calls
ncore='15'

#define project directory
scratchDir='/share/ScratchGeneral/scoyou/projects/andreas_OImin/'

#Real dataset file paths
sampleDir='input_vcfs/'
resultsDir='project_results/'
QCDir='logs/'
#tool='VCF_annotations'

#Set input paths
samplePath=$scratchDir$sampleDir
outPath=$scratchDir$sampleDir
#outPath=$scratchDir$sampleDir$tool

logPath=$scratchDir$sampleDir
#logPath=$scratchDir$sampleDir$tool

inFileExt='_decomp_norm.vcf.gz'

        ##make an array containing names of directories containing samples
        sample_arr=( $(ls $samplePath) )
        echo "# in samples array ${#sample_arr[@]}"
        echo "names in samples array ${sample_arr[@]}"

        #within each folder
        #for sample in ${sample_arr[@]}; do

        #Runs loop for only the first sample
        for sample in ${sample_arr[0]}; do

        ##define input directory, define and make output and log directories
                inPath=$samplePath$sample
                echo "inPath $inPath"

                logDir=$logPath$sample/$QCDir
                mkdir -p $logDir
                echo $logDir

                ##make an array containing names of vcf files within each sample directory
                fastq_arr=( $(ls $inPath/"F"??$inFileExt) )
                echo "# in fastq array ${#fastq_arr[@]}"
                echo "names in fastq array ${fastq_arr[@]}"

                ##define each input file within array
                inFile1=${fastq_arr[0]}
                echo "inFile $inFile1"

                #submit jobs to the cluster                
command1="vep -i $inFile1 --cache --offline --keep_csq --fasta /share/ClusterShare/biodata/contrib/scoyou/GRCh37/vep/homo_sapiens_dep/96_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz  --verbose --port 3337 --dir_cache /share/ClusterShare/biodata/contrib/scoyou/GRCh37/vep/ --fork 6 --o $outPath/$sample/$sample"_decomp_norm_VEP.vcf.gz" --force_overwrite --stats_text --vcf --variant_class --nearest gene --regulatory --exclude_null_alleles --no_check_alleles --terms SO --symbol --biotype --domains --compress_output gzip --gencode_basic --failed 0 --pick"

                echo "command $command1"
                qsub -P OsteoporosisandTranslationalResearch -N "VEP_"$sample -b y -wd $logDir -j y -R y -pe smp $ncore -V $command1

        done;

