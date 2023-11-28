//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow AVC_SW {
    take:
        

    main:
    
    

    emit:
        reads                                     // channel: [ val(meta), [ reads ] ]
        versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
    }
