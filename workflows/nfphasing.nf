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

params.run_whatshap_phase = false
params.run_longphase_phase = false

println "Running with the following parameters:"
println "Run longphase phase: ${params.run_longphase_phase}"
println "Run whatshap phase: ${params.run_whatshap_phase}"


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
                row.sample,
                file(row.bam), 
                file(row.snv), 
                row.sv ? file(row.sv) : file(params.no_file)
            ]
            }
    ch_whatshap_phase = inputs.map { 
        ref, sample, bam, snv, sv -> 
            [ ref, sample, bam, snv ] 
        }
    ch_longphase_phase = inputs.map {
        ref, sample, bam, snv, sv -> 
            [ ref, sample, bam, snv, sv ] 
        }
    longphase = longphase_phase(ch_longphase_phase)
    gzipped_vcf = bgzip(longphase.phased_vcf)
    tabix_indexed_vcf = tabix(gzipped_vcf.gz_vcf)


}

process longphase_phase {
    publishDir "results"
    container "quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0"

    input:
    tuple path(ref), val(sample), path(bam), path(snv), path(sv)

    output:
    path phased_vcf, emit: phased_vcf

    script:
    phased_prefix = "longphase_${sample}"
    phased_vcf = "${phased_prefix}.vcf"
    sv = sv.name == params.no_file ? null : sv
    """
    longphase phase \
        -b "$bam" \
        -r "$ref" \
        -s "$snv" \
        ${sv ? "--sv-file $sv" : ""} \
        -t ${task.cpus} \
        -o "$phased_prefix" \
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/