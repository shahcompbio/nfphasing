/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// set default params
params.skip_longphase = false
params.skip_whatshap = false
params.skip_qc = false

// Define inputs
println "Running with the following parameters:"
println "Input samplesheet: ${params.input}"
println "Output directory: ${params.outdir}"
println "Skip longphase: ${params.skip_longphase}"
println "Skip whatshap: ${params.skip_whatshap}"
println "Skip QC: ${params.skip_qc}"

no_file_name = file(params.no_file).name
qc_script = file("${workflow.projectDir}/subworkflows/local/phasingqc/phasingqc.py")

workflow NFPHASING {

    take:
    samplesheet // channel: samplesheet read in from --input
    main:

    inputs = Channel
         .fromPath(samplesheet)
         .splitCsv(header: true)
         .map { 
            row -> [
                file(row.ref),
                file(row.ref_index),
                row.sample,
                file(row.bam),
                file(row.bam_index), 
                file(row.snv), 
                row.sv ? file(row.sv) : file(params.no_file) // use no_file as sentinel if sv is not provided
            ]
            }
    ch_whatshap_phase = inputs
        .filter { 
            ref, ref_index, sample, bam, bam_index, snv, sv -> 
            sv.name == no_file_name
            } // filter out rows where sv is provided
        .map {
            ref, ref_index, sample, bam, bam_index, snv, sv -> 
            [ ref, ref_index, sample, bam, bam_index, snv ] 
        }
    ch_longphase_phase = inputs
        .map {
            ref, ref_index, sample, bam, bam_index, snv, sv -> 
            [ ref, ref_index, sample, bam, bam_index, snv, sv ] 
        }

    longphase = longphase_phase(ch_longphase_phase)
    whatshap = whatshap_phase(ch_whatshap_phase)

    ch_phased_vcf = longphase.concat( whatshap)
    gzipped_vcf = bgzip(ch_phased_vcf)
    tabix_indexed_vcf = tabix(gzipped_vcf)

}


process bgzip {
    publishDir params.outdir, mode: 'copy'
    container "quay.io/biocontainers/samtools:1.22--h96c455f_0"

    input:
    tuple val (sample), path (vcf)

    output:
    tuple val (sample), path (gz_vcf)

    script:
    gz_vcf = "${vcf.baseName}.vcf.gz"
    """
    bgzip -c "${vcf}" > "${gz_vcf}"
    """
}

process tabix {
    publishDir params.outdir, mode: 'copy'
    container "quay.io/biocontainers/tabix:0.2.5--0"

    input:
    tuple val (sample), path (vcf)

    output:
    tuple val (sample), path (tbi)

    script:
    tbi = "${vcf}.tbi"
    """
    tabix -p vcf "$vcf"
    """
}


process longphase_phase {
    container "quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0"

    input:
    tuple path (ref), path (ref_index), val (sample), path (bam), path (bam_index), path (snv), path (sv)

    output:
    tuple val (sample), path (phased_vcf)

    script:
    phased_prefix = "longphase_${sample}"
    phased_vcf = "${phased_prefix}.vcf"
    sv = sv.name == no_file_name ? null : sv
    """
    longphase phase \\
        -b "$bam" \\
        -r "$ref" \\
        -s "$snv" \\
        ${sv ? "--sv-file $sv" : ""} \\
        -t ${task.cpus} \\
        -o "$phased_prefix" \\
        --ont
    """
}

process whatshap_phase {
    container "quay.io/biocontainers/whatshap:2.8--py312hf731ba3_0"

    input:
    tuple path (ref), path (ref_index), val (sample) , path (bam) , path (bam_index) , path (snv)

    output:
    tuple val (sample), path (phased_vcf)

    script:
    phased_vcf = "whatshap_${sample}.vcf"
    """
    whatshap phase \\
        --reference "$ref" \\
        --output "$phased_vcf" \\
        --ignore-read-groups \\
        "$snv" \\
        "$bam"
    """
}

process phasing_qc {
    container "quay.io/shahlab_singularity/haplotagqc:latest"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val (sample), path (snv), path (script)

    output:
    path "*.txt", emit: phasing_qc_txt
    path "*.png", emit: phasing_qc_png


    script:
    """
    python phasingqc.py \\
        --vcf "$snv" \\
        --sample "${sample}}"
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/