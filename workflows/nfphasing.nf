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

// Define inputs
println "Running with the following parameters:"
println "Input samplesheet: ${params.input}"

no_file_name = file(params.no_file).name

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
            [ ref, sample, bam, snv, sv ] 
        }
    longphase = longphase_phase(ch_longphase_phase)
    whatshap = whatshap_phase(ch_whatshap_phase)

    ch_phased_vcf = longphase.phased_vcf.concat( whatshap.phased_vcf)
    gzipped_vcf = bgzip(ch_phased_vcf)
    tabix_indexed_vcf = tabix(gzipped_vcf.gz_vcf)


}

process longphase_phase {
    container "quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0"

    input:
    tuple path(ref), val(sample), path(bam), path(snv), path(sv)

    output:
    path phased_vcf, emit: phased_vcf

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

process bgzip {
    publishDir "results"
    container "quay.io/biocontainers/samtools:1.22--h96c455f_0"

    input:
    path vcf

    output:
    path gz_vcf, emit: gz_vcf

    script:
    gz_vcf = "${vcf.baseName}.vcf.gz"
    """
    bgzip -c "${vcf}" > "${gz_vcf}"
    """
}

process tabix {
    publishDir "results"
    container "quay.io/biocontainers/tabix:0.2.5--0"

    input:
    path vcf

    output:
    path tbi, emit: tbi

    script:
    tbi = "${vcf}.tbi"
    """
    tabix -p vcf "$vcf"
    """
}

process index_bam {
    publishDir "results"
    container "quay.io/biocontainers/samtools:1.22--h96c455f_0"

    input:
    path bam

    output:
    path indexed_bam, emit: indexed_bam

    script:
    indexed_bam = "${bam}.bai"
    """
    samtools index "$bam"
    """
}

process whatshap_phase {
    publishDir "results"
    container "quay.io/biocontainers/whatshap:2.8--py312hf731ba3_0"

    input:
    tuple path(ref), path(ref_index), val(sample), path(bam), path(bam_index), path(snv)

    output:
    path phased_vcf, emit: phased_vcf

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/