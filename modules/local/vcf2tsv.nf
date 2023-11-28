//Process
process VCF2TSV {    
    input:
        path (vcfanno)
        path (pythonconfig)
        path (preftsv)

    output:
        path "*.tsv" , emit:  vcf2tsv_out 
    
    script:
        """
        set -x
        id
        NO_COLOR=1
        vcf2tsv --nprocs 8 \
                --config $pythonconfig \
                --debug --keep-header \ 
                --prefer $preftsv \ 
                --sort --chunksize 250000 \
                $vcfanno
       
        """
}