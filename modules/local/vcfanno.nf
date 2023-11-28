
process vcfanno {
    tag "$sample_id"
    label 'process_low'

    conda "bioconda::vcfanno=0.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcfanno:0.3.3--h9ee0642_0':
        'biocontainers/vcfanno:0.3.3--h9ee0642_0' }"

    input:
    tuple val(sample_id), path(vcf), path(tbi), path(specific_resources)
    path toml
    path lua
    path resources

    output:
    tuple val(sample_id), path("*.vcf.gz") , emit: vcf
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    def lua_cmd = lua ? "--lua ${lua}" : ""
    """
    vcfanno \\
        -p ${task.cpus} \\
        ${args} \\
        ${lua_cmd} \\
        ${toml} \\
        ${vcf} \\
        > ${prefix}.vcf

    bgzip ${args2} --threads ${task.cpus} ${prefix}.vcf
    tabix ${prefix}.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(echo \$(vcfanno 2>&1 | grep version | cut -f3 -d' ' ))
    END_VERSIONS
    """


  input:
    path vcffile
    path vcftomldir 
    path vcftomlfilename 
    path luafile
    path hg38  
    path hg38liftover
    path tbi


  output:
   path "vcfano.*.gz" , emit: vcffile
   path "vcfano.*.gz.tbi" , emit: tbi

  script:
      """
      set -x
      
      find_sources.sh $hg38
      ls -al $hg38
      export ANN_VCF=vcfano.${vcffile}
      echo "OUTPUT VCF: \$ANN_VCF"
      export IRELATE_MAX_CHUNK=12000
      export IRELATE_MAX_GAP=1000
      vcfanno -p 4 -lua $luafile $vcftomlfilename $vcffile > \$ANN_VCF
      bgzip \$ANN_VCF
      tabix \$ANN_VCF.gz
      
      """

        // # find_sources.sh Assets/vcfanno/hg38/sources/MGRB.hg38_liftover.vcf.gz 1>&2 && \
        // # ANN_VCF=$input.vep.vcfanno.vcf && \
        // # echo "OUTPUT VCF: $ANN_VCF" 1>&2 && \
        // # IRELATE_MAX_CHUNK=12000 IRELATE_MAX_GAP=1000 vcfanno -p 4 -lua vcfanno/vcfanno_functions.lua vcfanno_MGRB.toml pedcbioportal_strand_3.vcf.gz > $ANN_VCF && \
        // # bgzip 1>&2 $ANN_VCF && \
        // # tabix 1>&2 $ANN_VCF.gz && \
        // # bcftools view -h $ANN_VCF.gz | tail -n1 | egrep '^#CHROM' --REMOVED

}