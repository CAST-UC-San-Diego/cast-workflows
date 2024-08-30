
version 1.0

workflow advntr_single_sample {
    input {
        String bam_file
        String region_file
        String google_project
        String vntr_id
        Int sleep_seconds
        Int mem
    }
    call download_input {
        input :
        bam_file = bam_file,
        region_file = region_file,
        google_project = google_project,
    }
    call genotype {
        input :
        region_file = region_file,
        vntr_id = vntr_id,
        target_bam_file = download_input.target_bam_file,
        target_bam_index_file = download_input.target_bam_index_file,
        sleep_seconds = sleep_seconds,
        mem = mem
    }
    
    call sort_index {
        input :
        vcf = genotype.genotype_output,
        mem = mem
    }
    output {
        File out_vcf = sort_index.out_vcf
        File out_vcf_index = sort_index.out_vcf_index
    }
    meta {
        description: "This workflow calls adVNTR to genotype VNTRs for a single sample"
    }
}



task sort_index {
  input {
    File vcf
    Int mem
  }

  String basename = basename(vcf, ".vcf")

  command <<<
    echo "Update the contig name in the header"
    bcftools view -h ~{vcf} | grep "^##" > header.txt
    cat header.txt | grep "contig" |sed 's/##contig=<ID=/##contig=<ID=chr/g' >> header.txt
    bcftools view -h ~{vcf} | grep -v "^##" >> header.txt
    cat header.txt
    echo "Header created. Now running reheader"
    bcftools reheader -h header.txt ~{vcf} > ~{basename}_rh.vcf
    echo "Add IDs to entries"
    bcftools annotate --set-id +'%CHROM\_%POS' ~{basename}_rh.vcf > ~{basename}_rh_id.vcf
    # Sort and index
    bcftools sort -Oz ~{basename}_rh_id.vcf  > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
  >>>

  runtime {
    docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    memory: mem + "GB"
    disks: "local-disk ${mem} SSD"
  }

  output {
    File out_vcf = "${basename}.sorted.vcf.gz"
    File out_vcf_index = "${basename}.sorted.vcf.gz.tbi"
  }
}

task download_input {
    input {
        String bam_file
        File region_file
        String google_project
    }


    # Names of all the intermediate files generated to get the target bam file.
    String unsorted_target_bam = "target_region_~{sample_id}_unsorted.bam"
    String sorted_target_bam = "~{sample_id}.bam"
    String sorted_target_bam_index = "~{sample_id}.bam.bai"
    String sample_id = sub(basename(bam_file), ".bam", "")

    command <<<
        ls -lh .
        export HTSLIB_CONFIGURE_OPTIONS="--enable-gcs"
        echo "pwd $(pwd)"
        /google-cloud-sdk/bin/gcloud --version
        /google-cloud-sdk/bin/gcloud config list --format='text(core.project)'
        export gcloud_token=$(/google-cloud-sdk/bin/gcloud auth application-default print-access-token --project ~{google_project})
        export GCS_OAUTH_TOKEN=${gcloud_token}
        export GCS_REQUESTER_PAYS_PROJECT="~{google_project}"
        samtools view -hb -o ~{unsorted_target_bam} --use-index -L ~{region_file} ~{bam_file}
        samtools sort -o ~{sorted_target_bam} ~{unsorted_target_bam}
        samtools index ~{sorted_target_bam}
    >>>

    runtime {
        docker:"sarajava/samtools:1.13_gcli"
        maxRetries: 3
        preemptible: 3
    }
    output {
        File target_bam_file = "~{sorted_target_bam}"
        File target_bam_index_file = "~{sorted_target_bam_index}"
    }
}

task genotype {
    input {
        File region_file
        String vntr_id
        File target_bam_file
        File target_bam_index_file
        Int sleep_seconds
        Int mem
    }

    # Provide the names of all the output files being generated by AdVNTR including intermediate files.
    String logging = "./log_~{sample_id}.bam.log"
    String filtering_out = "./filtering_out_~{sample_id}.unmapped.fasta.txt"
    String keywords = "./keywords_~{sample_id}.unmapped.fasta.txt"
    String unmapped = "./~{sample_id}.unmapped.fasta"
    String genotype_output = "./~{sample_id}.vcf"
    String sample_id = sub(basename(target_bam_file), ".bam", "")

    # VNTR_db is placed in the docker file. So the path is within the docker image.
    String vntr_db = "/adVNTR/vntr_db/p_vntrs_g_vntrs.db"
    #String vntr_db = "/adVNTR/vntr_db/p_vntrs_g_vntrs_lt_500bp.db"

    # To get all p-vntr ids
    #            -vid $(cat /adVNTR/vntr_db/phenotype_associated_vntrs_comma.txt | tr -d \\r\\n ) \

    command <<<
        sleep ~{sleep_seconds}
        echo "Num processors on this device is"
        nproc
        samtools --version
        # To suppress warnings on pysam about old index file
        touch ~{target_bam_index_file}
        echo "Num reads in input bam file $(samtools view -c ~{target_bam_file})"
        if [[ "~{vntr_id}" == "ALL" ]] ; then
                /usr/bin/time -v advntr genotype  \
                --alignment_file ~{target_bam_file} \
                --models ~{vntr_db}  \
                --working_directory . \
                --outfmt vcf \
                --threads 24 \
                --pacbio > ~{genotype_output}
        else
                /usr/bin/time -v advntr genotype  \
                --alignment_file ~{target_bam_file} \
                --models ~{vntr_db}  \
                --working_directory . \
                -vid ~{vntr_id} \
                --outfmt vcf \
                --threads 24 \
                --pacbio > ~{genotype_output}
        fi
    >>>

    runtime {
        docker:"sarajava/advntr:1.5.0_v16"
        memory: mem + "GB"
        cpu: "16"
    }

    output {
        File? log_file = "~{logging}"
        File? filtering_out = "~{filtering_out}"
        File? keywords = "~{keywords}"
        File genotype_output = "~{genotype_output}"
    }
}
