# Pipeline_SV
Pipeline for structural variants detection using whole genome long reads


nota: mancano i seguenti file 
- repeatmasker.bed
- segdups.bed
- campione.bam e relativo bai (da specificare nel samplesheet)
- reference.fa


comando per lanciare: 
nextflow run main.nf -profile local -resume



file da scaricare: 
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed 
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
