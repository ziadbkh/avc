/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAvc.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { ALIGNANDSORT                  } from '../modules/local/alignandsort'
include { MERGEBAMS                     } from '../modules/local/mergebams.nf'
include { HAPLOTYPECALLER               } from '../modules/local/haplotypecaller.nf'
include { MERGEVCFS                     } from '../modules/local/mergevcfs.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SPLITFASTQ                     } from '../modules/local/splitfastq.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow AVC {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    
    //INPUT_CHECK (    ch_input )
    //ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    


    /*Channel
        .fromFilePairs(params.input + '*_R{1,2}*.fastq.gz', checkIfExists:true)
        .map{[it[0], it[1][0].simpleName, it[1]]}
        .set{ch_input_fastq}
    */
    
    Channel
        .fromPath( params.fasta, checkIfExists:true).first().set{ch_reference_genome}

    Channel
        .fromPath( params.fasta + '.{bwt,sa,ann,amb,pac}', checkIfExists:true).toSortedList().set{ch_reference_genome_extra_bwa}
    
    
    ch_reference_genome.map{"${it.getParent()}/${it.getBaseName()}.dict"}
    .combine(Channel
            .fromPath( params.fasta + '.fai', checkIfExists:true)
    )
    .set{ch_reference_genome_extra_gatk}


    Channel
        .fromPath( params.intervals + '*scattered.interval_list', checkIfExists:true).set{ch_intervals}

    intervals_files = file(params.intervals + '*scattered.interval_list', checkIfExists:true)

    Channel
        .fromPath( params.input, checkIfExists:true )
        .splitCsv(header: true)
        .map { 
            row -> {
            if (row.fastq2)
                {    
                    [row.sample_id, 
                    row.type.toLowerCase(), 
                    true, 
                    [file(row.fastq1, checkIfExists:true), file(row.fastq2, checkIfExists:true)] ] 
                }else{
                    [row.sample_id, row.type.toLowerCase(), false, [file(row.fastq1, checkIfExists:true)] ]
                }
            }
           
        }
        .set{ch_meta_complete}
    

    ch_meta_complete
    .map{[it[0], it[2], it[3]]}
    .branch {
        small: it[2][0].size() <= params.fastq_max_size
        large: it[2][0].size() > params.fastq_max_size
    }
    .set {
        ch_meta_all
    }
    
    ch_meta_complete
    .map{[it[0], it[3], Math.ceil(it[3][0].size() / params.fastq_max_size).intValue()]}
    .groupTuple()
    .map{[it[0], it[2].sum()]}
    .set{ch_sample_file_cnt}

    //ch_sample_file_cnt.view()
    //ch_meta_all.large.view()
    //ch_meta_all.small.view()

    ch_meta_complete
    .map{[it[0], it[1]]}
    .branch {
        tumor: it[1] == "t" || it[1] == "tumor" || it[1] == "0"
        normal: true 
    }
    .set{ch_sample_type}


    SPLITFASTQ(
        ch_meta_all.large,
        Channel.value(params.fastq_max_size)
    )

    
    SPLITFASTQ.out.splited_fastq.
        transpose()
        .map{[it[0], it[2].getName().replace("R1", "").replace("R2", ""), it[1], it[2]]}
        .groupTuple(by:[0, 1, 2], size:2, remainder:true)
        .map{[it[0], it[2], it[3]]}
        .mix(ch_meta_all.small)
        .set{ch_all_fastq}
    
    //ch_all_fastq.view()
    
    ALIGNANDSORT(
        ch_all_fastq,
        ch_reference_genome.combine(ch_reference_genome_extra_bwa).flatten().toSortedList()
    )
    
    ch_sample_file_cnt
    .cross(ALIGNANDSORT.out.aligned_and_sorted_bam)
    .map{ tuple(groupKey(it[0][0], it[0][1]), it[1][1], it[1][2] )}
    .groupTuple()
    .branch {
        singelton: it[1].size() <= 1
        multiple: it[1].size() > 1
    }
    .set {
        ch_sample_bams
    }
    
    //ch_sample_bams.multiple.view()
    //ch_sample_bams.singelton.view()
    
    MERGEBAMS(
        ch_sample_bams.multiple
    )
    
    
    MERGEBAMS.out.merged_bam
    .mix(
        ch_sample_bams
        .singelton
        .map{[it[0], it[1][0], it[2][0]]}
    ).set{
        ch_all_bams
    }
    
    ch_all_bams
    .join(ch_sample_type.normal)
    .map{[it[0], it[1], it[2]]}
    .set{ch_normal_bams}

    
    HAPLOTYPECALLER(
        ch_normal_bams,
        ch_reference_genome
            .combine(ch_reference_genome_extra_gatk)
            .flatten()
            .toSortedList(),
        ch_intervals

    )
    
    HAPLOTYPECALLER.out.gvcf
    .groupTuple(size: intervals_files.size())
    .set{ch_sample_gvcf}
    
    MERGEVCFS(
        ch_sample_gvcf
    )
    
    if (params.publish_dir_mode == "move"){
        ch_sample_gvcf
        .map{it[0]}
        .join(ch_normal_bams)
        .mix(
            ch_sample_type.tumor
            .map{it[0]}
            .join(ch_all_bams)
        ).set{ch_output_files}
    }
    else{
        ch_all_bams.set{ch_output_files}
    }

    ch_output_files.subscribe{
        file("${params.outdir}/bams/").mkdirs()
        
        println("Migrating bam fils for the sample ${it[0]}")
        if (params.publish_dir_mode == "move"){
            it[1].moveTo("${params.outdir}/bams/${it[1].getName()}")
            it[2].moveTo("${params.outdir}/bams/${it[2].getName()}")
        }else{
            it[1].copyTo("${params.outdir}/bams/${it[1].getName()}")
            it[2].copyTo("${params.outdir}/bams/${it[2].getName()}")
        }
    }

    /*
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    */
    
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
