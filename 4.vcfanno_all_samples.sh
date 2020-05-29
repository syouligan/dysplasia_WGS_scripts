#!/bin/bash

###Filter VCF files to high quality calls
ncore='20'

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

inFileExt='_decomp_norm_VEP_bgz.vcf.gz'

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
                fastq_arr=( $(ls $inPath/*$inFileExt) )
                echo "# in fastq array ${#fastq_arr[@]}"
                echo "names in fastq array ${fastq_arr[@]}"

                ##define each input file within array
                inFile1=${fastq_arr[0]}
                echo "inFile $inFile1"

                #submit jobs to the cluster                
                command1="vcfanno -p 15 /home/scoyou/projects/andreas_OImin/scripts/vcfanno_config_1.conf $inFile1 | bgzip > $outPath$sample/$sample'_decomp_norm_VEP_vcfanno.vcf.gz'"

                echo "command $command1"
                qsub -P OsteoporosisandTranslationalResearch -N "vcfanno_"$sample -b y -wd $logDir -j y -R y -pe smp $ncore -V $command1

        done;

