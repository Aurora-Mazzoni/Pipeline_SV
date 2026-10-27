nextflow.enable.dsl=2

//dati, percorsi locali specificati nel config
params.samplesheet = ""
params.reference = ""
params.cleaning_script = ""
params.stats_script = ""
params.speedxx_script = ""   
params.merge_script = ""
params.convert_script = ""
params.repeat_masker = ""    
params.seg_dups = ""
params.speedxx_checkpoint = "./speedxx_checkpoint" // cartella dove speedxx salva i blocchi già fatti, per non perdere progressi se si blocca
params.jasmine_vcf = ""      // FACOLTATIVO: vcf di jasmine già pronto (non compresso), salta i caller

// SNIFFLES (solo esecuzione)
process RUN_SNIFFLES2 {
    input:
    tuple val(sample_id), path(bam_files)
    path reference //file path perché è uno solo per tutti i sample

    output:
    //file zippato .vcf.gz e il suo indice .tbi di tabix
    tuple val(sample_id), path("${sample_id}_sniffles.vcf")
    path "${sample_id}_log_risorse_sniffles.txt" //salva anche il log

    script:
    """
    /usr/bin/time -v -o ${sample_id}_log_risorse_sniffles.txt sniffles \\
        -i ${bam_files[0]} \\
        -v ${sample_id}_sniffles.vcf \\
        --reference ${reference} \\
        --threads ${task.cpus}
    """
}

// CUTESV (solo esecuzione)
process RUN_CUTESV {
    input:
    tuple val(sample_id), path(bam_files)
    path reference //file path perché è uno solo per tutti i sample

    output:
    //file zippato .vcf.gz e il suo indice .tbi di tabix
    tuple val(sample_id), path("${sample_id}_cutesv.vcf")
    path "${sample_id}_log_risorse_cutesv.txt" //salva anche il log

    script:
    """
    mkdir -p cutesv_tmp

    /usr/bin/time -v -o ${sample_id}_log_risorse_cutesv.txt cuteSV \\
        ${bam_files[0]} \\
        ${reference} \\
        ${sample_id}_cutesv.vcf \\
        cutesv_tmp \\
        --sample ${sample_id} \\
        --genotype \\
        --threads ${task.cpus}
    """
}

// SVIM (solo esecuzione e senza filtro)
process RUN_SVIM {
    input:
    tuple val(sample_id), path(bam_files)
    path reference //file path perché è uno solo per tutti i sample

    output:
    // Passa al processo successivo direttamente il variants.vcf estratto
    tuple val(sample_id), path("svim_outdir_${sample_id}/variants.vcf")
    path "${sample_id}_log_risorse_svim.txt" //salva anche il log

    script:
    """
    mkdir -p svim_outdir_${sample_id} 

    /usr/bin/time -v -o ${sample_id}_log_risorse_svim.txt svim alignment \\
        svim_outdir_${sample_id} \\
        ${bam_files[0]} \\
        ${reference} \\
        --skip_consensus
    """
}

// COMPRESSIONI E FILTRO SVIM
process FILTER_AND_ZIP {    
    publishDir path: { "results/tool_singoli/${sample_id}" }, mode: 'copy', pattern: '*.vcf.gz*'

    input:
    tuple val(sample_id), path(sniffles_vcf), path(cutesv_vcf), path(svim_vcf)

    output:
    //emit manda direttamente a jasmine
    tuple val(sample_id), path(sniffles_vcf), path(cutesv_vcf), path("svim_filtered.vcf"), emit: jasmine_input
    // salvataggio file zippati 
    tuple val(sample_id), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: zipped_files

    script:
    """
    # Filtra svim (QUAL >= 15), zippa in bgzip (-O z) e salva tutto 
    bcftools view -i 'QUAL >= 15' ${svim_vcf} -o svim_filtered.vcf

    #zippa e indicizza 
    bgzip -c ${sniffles_vcf} > ${sample_id}_sniffles.vcf.gz
    bgzip -c ${cutesv_vcf} > ${sample_id}_cutesv.vcf.gz
    bgzip -c svim_filtered.vcf > ${sample_id}_svim.vcf.gz

    tabix -p vcf ${sample_id}_sniffles.vcf.gz
    tabix -p vcf ${sample_id}_cutesv.vcf.gz
    tabix -p vcf ${sample_id}_svim.vcf.gz    
    """
}

// UNION WITH JASMINE
process UNION_JASMINE {
    publishDir path: {"results/jasmine/${sample_id}"}, mode: 'copy'

    input:
    // riceve i 3 file da emit del processo precedente
    tuple val(sample_id), path(sniffles_raw), path(cutesv_raw), path(svim_raw)

    output:
    tuple val(sample_id), path("jasmine_union_3tools.vcf")

    script:
    """
    # file con i percorsi (richiesto per l'esecuzione di jasmine)
    echo "${sniffles_raw}" > lista_jasmine.txt
    echo "${cutesv_raw}" >> lista_jasmine.txt
    echo "${svim_raw}" >> lista_jasmine.txt

    #Jasmine
    jasmine file_list=lista_jasmine.txt out_file=jasmine_union_3tools.vcf max_dist=1000
    """
}

// SOLO PULIZIA DEL FILE JASMINE 
process CLEAN_JASMINE {
    input: 
    tuple val(sample_id), path("jasmine_union_3tools.vcf")
    path py_script

    output:
    tuple val(sample_id), path("jasmine_union_3tools_clean.vcf") 

    script: 
    """
    # pulizia con script + sorting e zip
    python3 ${py_script}
    """
}

//ORDINAMENTO E COMPRESSIONE DEL FILE JASMINE
process SORT_JASMINE {
    publishDir "results/jasmine", mode: 'copy'

    input: 
    tuple val(sample_id), path(clean_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_jasmine_clean_sorted.vcf.gz"), path("${sample_id}_jasmine_clean_sorted.vcf.gz.tbi")

    script: 
    """
    bcftools sort ${clean_vcf} -O z -o ${sample_id}_jasmine_clean_sorted.vcf.gz
    tabix -p vcf ${sample_id}_jasmine_clean_sorted.vcf.gz
    """
}

// STATISTICHE VCF JASMINE 
process STAT_VCF {
    publishDir "results/stats", mode: 'copy' 

    input:
    tuple val(sample_id), path(clean_vcf), path(clean_tbi)
    path stats_script

    output:
    tuple val(sample_id), path("${sample_id}_stats.txt")

    script:
    """
    python3 ${stats_script} ${clean_vcf} ${sample_id}_stats.txt
    """
}

//EXTRACTING FEATURES FOR VARTRUSTML

//ESTRAI LE REGIONI DELLE SV  
process SV_REGIONS_EXTRACT {                                                        
    publishDir path: {"results/features/${sample_id}"}, mode:"copy"                               
    input:
    tuple val(sample_id), path(vcf), path(tbi)
    output:
    tuple val(sample_id), path("${sample_id}_sv_regions.bed")                                 
    script:
    """
    bcftools query  -f '%CHROM\t%POS0\t%END\t%SVTYPE\t%SVLEN\n' ${vcf} > ${sample_id}_sv_regions.bed
    """ 
 }

 //ESTRAI GC CONTENT
 process GC_CONTENT {
    publishDir path: {"results/features/${sample_id}"}, mode: 'copy'

    input:
    tuple val(sample_id), path(sv_bed)
    path reference 

    output:
    tuple val(sample_id), path("${sample_id}_GC.txt")

    script:
    """
    bedtools nuc -fi ${reference} -bed ${sv_bed} > ${sample_id}_GC.txt
    """
}

//ESTRAI COVERAGE 
process SV_COVERAGE {
    publishDir path: {"results/features/${sample_id}"}, mode: 'copy'

    input:
    tuple val(sample_id), path(bam_files), path(sv_bed) // BAM + BED stesso sample
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.sv.regions.bed.gz") // mosdepth zippa l'output di default

    script:
    """
    mosdepth --threads ${task.cpus} --by ${sv_bed} --use-median --mapq 30 --fasta ${reference} ${sample_id}.sv ${bam_files[0]}
    """
}

// ESTRAZIONE DELLE FEATURE PER VARTRUSTML (3 script della tutor, in ordine)
process EXTRACT_FEATURES {
    publishDir path: {"results/features/${sample_id}"}, mode: 'copy'

    input:
    // per ogni campione: bam+bai, vcf pulito e ordinato + indice, GC (da GC_CONTENT), coverage (da SV_COVERAGE)
    tuple val(sample_id), path(bam_files), path(vcf), path(tbi), path(gc_txt), path(cov_bed)
    path speedxx_script
    path merge_script
    path convert_script
    path repeat_masker
    path seg_dups
    path ckpt_dir   // cartella dei checkpoint di speedxx (sta fuori da work/, così sopravvive alle interruzioni)

    output:
    // 4 file finali per VarTrustML, uno per tipo di SV
    tuple val(sample_id), path("${sample_id}_features_DEL.csv"), path("${sample_id}_features_DUP.csv"), path("${sample_id}_features_INS.csv"), path("${sample_id}_features_INV.csv"), emit: per_tipo
    // tabella completa con tutte le colonne (utile per controllare)
    tuple val(sample_id), path("${sample_id}_features.tsv"), emit: tabella

    script:
    """
    # 1) speedxx: misure sulle read attorno a ogni SV (gli si passa il BAM anche se l'opzione si chiama --cram_file)
    python3 ${speedxx_script} --cram_file ${bam_files[0]} --vcf_file ${vcf} --out_file ${sample_id}_extra_bam.tsv --threads ${task.cpus} --checkpoint_dir ${ckpt_dir}

    # 2) merge_fixed: unisce VCF, GC, coverage, ripetuti e la tabella di speedxx
    python3 ${merge_script} --repeat_masker ${repeat_masker} --segmental_duplications ${seg_dups} --vcf_file ${vcf} --gc_content ${gc_txt} --coverage_mosdepth ${cov_bed} --sv_feature_matrix ${sample_id}_extra_bam.tsv --output_path ${sample_id}_features.tsv

    # 3) convert_sv: tiene le colonne del modello e divide per tipo di SV
    python3 ${convert_script} ${sample_id}_features.tsv -o ${sample_id}_features
    """
}


//------------------------------------------------------------------------------------

workflow {
    bam_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> 
            // Ricrea la tuple: sample_id e la lista dei due file (bam, bai)
            tuple(row.sample_id, [ file(row.bam), file(row.bai) ]) 
        }
    ref_ch = file(params.reference, checkIfExists: true) 
    py_ch  = file(params.cleaning_script, checkIfExists: true)

    if (params.jasmine_vcf) {
        // SCORCIATOIA: usa un vcf di jasmine già pronto e salta i 3 caller e jasmine
        // (pensato per UN solo campione: il vcf viene associato al sample_id della samplesheet)
        jasmine_ch = bam_ch.map { sample_id, bam_files ->
            tuple(sample_id, file(params.jasmine_vcf, checkIfExists: true))
        }
    } else {
        RUN_SNIFFLES2(bam_ch, ref_ch) 
        RUN_CUTESV(bam_ch, ref_ch) 
        RUN_SVIM(bam_ch, ref_ch) 

        //  join() guarda il 'sample_id' di ogni output e aspetta di avere tutti
        // e 3 i pezzi dello stesso sample prima di passarli al channel di Jasmine (unisce i file grezzi!)
        raw_vcf_ch = RUN_SNIFFLES2.out[0]
                        .join(RUN_CUTESV.out[0])
                        .join(RUN_SVIM.out[0])
        filtered_and_zipped_ch = FILTER_AND_ZIP(raw_vcf_ch)
        //usa l'outut di emit del processo di filtro e zip
        jasmine_ch = UNION_JASMINE(filtered_and_zipped_ch.jasmine_input)
    }

    clean_raw_ch = CLEAN_JASMINE(jasmine_ch, py_ch)
    clean_sorted_ch = SORT_JASMINE(clean_raw_ch)    
    stats_ch = file(params.stats_script, checkIfExists: true)
    STAT_VCF(clean_sorted_ch,stats_ch)

    sv_regions_ch = SV_REGIONS_EXTRACT(clean_sorted_ch)
    gc_content_ch = GC_CONTENT(sv_regions_ch, ref_ch)
    coverage_ch = SV_COVERAGE(bam_ch.join(sv_regions_ch), ref_ch)

    // file per EXTRACT_FEATURES (uguali per tutti i campioni)
    speedxx_ch = file(params.speedxx_script, checkIfExists: true)
    merge_ch   = file(params.merge_script, checkIfExists: true)
    convert_ch = file(params.convert_script, checkIfExists: true)
    rmsk_ch    = file(params.repeat_masker, checkIfExists: true)
    segdup_ch  = file(params.seg_dups, checkIfExists: true)

    // join per sample_id: (sample_id, [bam, bai], vcf, tbi, gc, coverage)
    features_input_ch = bam_ch
        .join(clean_sorted_ch)
        .join(gc_content_ch)
        .join(coverage_ch)

    // cartella dei checkpoint: viene creata se non esiste
    ckpt_ch = file(params.speedxx_checkpoint)
    ckpt_ch.mkdirs()

    EXTRACT_FEATURES(features_input_ch, speedxx_ch, merge_ch, convert_ch, rmsk_ch, segdup_ch, ckpt_ch)
}