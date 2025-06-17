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
            .view()
    ch_whatshap_phase = inputs.map { 
        ref, sample, bam, snv, sv -> 
            [ ref, sample, bam, snv ] 
        }
    ch_longphase_phase = inputs.map {
        ref, sample, bam, snv, sv -> 
            [ ref, sample, bam, snv, sv ] 
        }
    longphase_results = longphase_phase(ch_longphase_phase)
    


}

process longphase_phase {
    publishDir "results"
    container "quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0"

    input:
    tuple path(ref), val(sample), path(bam), path(snv), path(sv)

    output:
    path phased_vcf, emit: phased_vcf
    val(sample), emit: sample

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

process tabix {
    publishDir "results"

    input:
    path phased_vcf, emit: phased_vcf
    val(sample), emit: sample

    output:
    path "${phased_vcf}.tbi", emit: tbi

    script:
    """
    tabix -p vcf "$phased_vcf"
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
