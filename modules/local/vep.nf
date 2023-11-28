
// NOTE: there is an nfcore: https://nf-co.re/modules/ENSEMBLVEP
process vep{
    tag "$sample_id"
    label 'process_medium'

    conda "bioconda::ensembl-vep=110.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:110.0--pl5321h2a3209d_0' :
        'biocontainers/ensembl-vep:110.0--pl5321h2a3209d_0' }"

    input:
    tuple val(sample_id), path(vcf), path(custom_extra_files)
    val   genome
    val   species
    val   cache_version
    path  cache
    tuple val(meta2), path(fasta)
    path  extra_files

    output:
    tuple val(meta), path("*.vcf.gz")  , optional:true, emit: vcf
    tuple val(meta), path("*.tab.gz")  , optional:true, emit: tab
    tuple val(meta), path("*.json.gz") , optional:true, emit: json
    path "*.summary.html"              , emit: report
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta $fasta" : ""
    """
    vep \\
        -i $vcf \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
        $reference \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --stats_file ${prefix}.summary.html \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
    
    
  /*  
        input:
        path refgnome
        path inputvcf
        path vepcache

    output:
        path "*.vep.vcf.gz" , emit: vepvcf
        path "*.vep.html" , emit: vephtml

    script:
        """
        # Set echo comands debug options
        set -x 
        ls -al $refgnome
        ls -al $inputvcf
        # Full vcf input file name , including extension
        FILENAME="$inputvcf"
        # Strip the extension using first . as delimiter (As shell parameter)
        PREFIX="\${FILENAME%%.*}"
        echo "Filename prefix : \$PREFIX"
        mkdir -p ./plugins
        #Put a file in here as it does a copy and fails if empty!
        touch ./plugins/filetest
        ls ./plugins
        #--vep_plugins_dir=./plugins
        run_vep.sh  --vep_plugins_dir=./plugins --options=--merged --prefix=\$PREFIX --input_vcf=$inputvcf --ref_genome=$refgnome --vep_cache=$vepcache
        ls -al .


        """
    */

    //   echo \${inputvcf}
    //     echo \${inputvcf%.*}
    //     PREFIX=${$inputvcf%:*}
//           --output_file MYPREFIX.vep.vcf \
//   --stats_file MYPREFIX.vep.html

        // --input_vcf=vt_norm/$input.vcf.gz 
        // --ref_genome=Assets/genome/hg38/GRCh38_no_alt.fa 
        // --vep_cache=Assets/vep/cache/homo_sapiens_merged_vep_100_GRCh38.tar.gz
// 
// run_vep.sh 1>&2 --vep_plugins_dir=./plugins/ --options=--merged --prefix=`echo $input.vcf | sed 's/\.vcf$//'` --input_vcf=vt_norm/$input.vcf.gz --ref_genome=Assets/genome/hg38/GRCh38_no_alt.fa --vep_cache=Assets/vep/cache/homo_sapiens_merged_vep_100_GRCh38.tar.gz

//| sed 's/\.vcf$//'
        //  run_vep.sh  --vep_plugins_dir=./plugins/ --options=--merged --prefix=`echo $inputvcf | sed 's/\.vcf$//'` --input_vcf=vt_norm/$input.vcf.gz --ref_genome=Assets/genome/hg38/GRCh38_no_alt.fa --vep_cache=Assets/vep/cache/homo_sapiens_merged_vep_100_GRCh38.tar.gz


//  # mkdir "./plugins" && \
//         # ls ./plugins 1>&2 && \
//         # run_vep.sh 1>&2 --vep_plugins_dir=./plugins/ --options=--merged --prefix=`echo $input.vcf | sed 's/\.vcf$//'` --input_vcf=vt_norm/$input.vcf.gz --ref_genome=Assets/genome/hg38/GRCh38_no_alt.fa --vep_cache=Assets/vep/cache/homo_sapiens_merged_vep_100_GRCh38.tar.gz


}