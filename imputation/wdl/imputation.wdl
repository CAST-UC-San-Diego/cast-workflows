version 1.0

workflow imputation {
    input {
        String vcf
        String vcf_index
        String ref_panel
        String ref_panel_bref
        String ref_panel_index
        String genetic_map
        String out_prefix
        String GOOGLE_PROJECT = ""
        String chrom
        Boolean skip_subset_vcf
        Int? mem
        Int? window_size
        Int? overlap
        File samples_file
	File regions_file
    }
    # skip_subset_vcf is deprecated and should be removed

    call subset_vcf {
        input:
            samples_file=samples_file,
            regions_file=regions_file,
            ref_panel=ref_panel,
            ref_panel_index=ref_panel_index,
            vcf=vcf,
            vcf_index=vcf_index,
            mem=mem,
            GOOGLE_PROJECT=GOOGLE_PROJECT,
            out_prefix=out_prefix
        }

    call beagle {
        input :
          vcf=subset_vcf.outfile,
          vcf_index=subset_vcf.outfile_index,
          ref_panel=ref_panel_bref,
          ref_panel_index=ref_panel_index,
          genetic_map=genetic_map,
          out_prefix=out_prefix,
          chrom=chrom,
          mem=mem,
          window_size=window_size,
          overlap=overlap
    }
    call sort_index_beagle {
        input :
            vcf=beagle.outfile,
            mem=mem
    }
    output {
        File outfile = sort_index_beagle.outvcf
        File outfile_index = sort_index_beagle.outvcf_index
    }
    meta {
      description: "Run Beagle on a subset of samples on a single chromesome with default parameters"
    }
}


task subset_vcf {
    input {
        String vcf
        String vcf_index
        String ref_panel
        String ref_panel_index
        File? samples_file
	File? regions_file
        Int? mem
        String GOOGLE_PROJECT = ""
        String out_prefix=out_prefix
    }

    command <<<
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        echo "memory on device"
        grep MemTotal /proc/meminfo
        echo "disk space on device"
        df -h
        #echo "Check if the bam file is corrupted or not"

        # To get a list of sample ids
        #bcftools query -l ~{vcf} > sample_ids.txt
        # The "bcftools head" command was to check the header for the labeling if contigs e.g. chr21 vs 21.
        # bcftools head ~{vcf} > header.txt
        # Subsetting region for each chromesome
        echo "Select reference samples from the actual samples (to exclude later)"
        bcftools query -l ~{ref_panel} > ref_sample_ids.txt
        echo "Exclude reference samples from the query. Otherwise beagle will give an error."
        grep -v -x -f ref_sample_ids.txt ~{samples_file} > samples_file_clean.txt
        echo "Subsetting the target region and samples from the vcf file"
        bcftools view -Oz -I -R ~{regions_file} -S samples_file_clean.txt ~{vcf} > ~{out_prefix}.vcf.gz && tabix -p vcf ~{out_prefix}.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
	memory: mem + "GB"
        #bootDiskSizeGb: mem
	disks: "local-disk " + mem + " SSD"
        maxRetries: 3
        preemptible: 3
    }

    output {
        File outfile = "${out_prefix}.vcf.gz"
        File outfile_index = "${out_prefix}.vcf.gz.tbi"
        File ref_sample_ids = "ref_sample_ids.txt"
    }
}

task index_vcf {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf")

    command <<<
        bgzip -c ~{vcf} > ~{basename}.vcf.gz && tabix -p vcf ~{basename}.vcf.gz
        tabix -p vcf ~{vcf}
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
        File outfile = "${basename}.vcf.gz"
        File outfile_index = "${basename}.vcf.gz.tbi"
    }
}

task beagle {
    input {
        File vcf
        File vcf_index
        File ref_panel
        File ref_panel_index
        File genetic_map
        String out_prefix
        String chrom
        Int? mem
        Int? window_size
        Int? overlap
    }

    command <<<
        echo "vcf: ~{vcf}"
        java -Xmx~{mem}g -jar /beagle.jar \
            gt=~{vcf} \
            ap=true \
            gp=true \
            ref=~{ref_panel} \
            window=~{window_size} \
            overlap=~{overlap} \
            chrom=~{chrom} \
            out=~{out_prefix}_output \
            map=~{genetic_map}
    >>>

    #file upto 300mb use mem=25
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/beagle:latest"
	disks: "local-disk " + mem + " SSD"
	memory: mem + "GB"
    }

    output {
       File outfile = "${out_prefix}_output.vcf.gz"
    }
}

task sort_index_beagle {
    input {
      File vcf
      Int? mem 
    }

    String basename = basename(vcf, ".vcf.gz")

    command <<<
        df -h /cromwell_root
        echo "Sorting and compressing"
        zcat ~{vcf} | vcf-sort | bgzip -c > ~{basename}.sorted.vcf.gz
        df -h /cromwell_root
        echo "Extracting TRs"
        bcftools view -Oz -i 'ID="."' ~{basename}.sorted.vcf.gz > ~{basename}.sorted_TR.vcf.gz
        tabix -p vcf ~{basename}.sorted_TR.vcf.gz
        echo "Number of TRs in the genotyped file"
        bcftools view -i 'ID="."' ~{basename}.sorted_TR.vcf.gz | grep -v "^#" | wc -l
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
	memory: mem + "GB"
	disks: "local-disk " + mem + " SSD"
    }

    output {
    File outvcf = "${basename}.sorted_TR.vcf.gz"
    File outvcf_index = "${basename}.sorted_TR.vcf.gz.tbi"
  }
}
