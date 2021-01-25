version 1.0

import "../wdl_structs.wdl"

# General tasks

task CompressVcf {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        File vcf
        String vcfCompressedPath = sub(vcf, "$", ".gz")
    }

    command {
        bgzip \
        -c \
        ~{vcf} \
        > ~{vcfCompressedPath} \
    }

    output {
        File vcfCompressed = "~{vcfCompressedPath}"
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task IndexVcf {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        File vcfCompressed
    }

    command {
        gatk \
        IndexFeatureFile \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=4" \
        --feature-file ~{vcfCompressed}
    }

    output {
        IndexedVcf vcfCompressedIndexed = object {
                vcf : "${vcfCompressed}", 
                vcfIndex : "${vcfCompressed}.tbi"
            }
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}


# Prep tasks

task RenameMetadata {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String pairName
        File callerVcf
        String renameMetaVcfPath = sub(callerVcf, "$", ".rename_metadata.vcf")
        String tool
        File callerVcf
    }

    command {
        python \
        rename_metadata.py \
        ~{callerVcf} \
        ~{renameMetaVcfPath} \
        ~{tool}
    }

    output {
        File renameMetaVcf = "~{renameMetaVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task MergePrepSupport {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath =  sub(renameMetaVcf, ".rename_metadata.vcf$", ".merge_prep.vcf") 
        String tool
    }

    command {
        python \
        merge_prep.py \
        --vcf ~{renameMetaVcf} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool} \
        --support \
    }

    output {
        File prepCallerVcf = "~{prepCallerVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task MergePrep{
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath = sub(renameMetaVcf, ".rename_metadata.vcf$", ".merge_prep.vcf")
        String tool
    }

    command {
        python \
        merge_prep.py \
        --vcf ~{renameMetaVcf} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool} \
    }

    output {
        File prepCallerVcf = "~{prepCallerVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task RenameVcf {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        File prepCallerVcf
        String pairName
        String renameVcfPath = sub(prepCallerVcf, ".merge_prep.vcf$", ".rename.vcf")
        String normal
        String tumor
        String tool
    }

    command {
        python \
        rename_vcf.py \
        ~{prepCallerVcf} \
        ~{renameVcfPath} \
        ~{normal} \
        ~{tumor} \
        ~{tool} \
    }

    output {
        File renameVcf = "~{renameVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task SplitMultiAllelic {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String pairName
        String splitVcfPath
        IndexedReference referenceFa
        IndexedVcf vcfCompressedIndexed
    }

    command {
        bcftools \
        norm \
        -m \
        -any \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{vcfCompressedIndexed} \
    }

    output {
        File splitVcf = "~{splitVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task SplitMnv {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        File splitVcf
        String mnvVcfPath = sub(splitVcf, ".split.vcf", ".split_mnvs.vcf")
        String tool
    }

    command {
        python \
        split_mnv.py \
        ~{splitVcf} \
        ~{mnvVcfPath} \
        ~{tool} \
    }

    output {
        File mnvVcf = "~{mnvVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task RemoveContig {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String mnvVcfPath
        String removeChromVcfPath = "~{mnvVcfPath}"
        File removeChromVcf
    }

    command {
        python2.7 \
        remove_contig.py \
        ~{removeChromVcf} \
        ~{removeChromVcfPath} \
    }

    output {
        File removeContigVcf = "~{removeChromVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task Gatk4MergeSortVcf {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String sortedVcfPath
        Array[File] tempVcfs
        IndexedReference referenceFa
    }

    command {
        gatk \
        SortVcf \
        --java-options "-Xmx8196m -XX:ParallelGCThreads=4" \
        -SD ~{referenceFa.dict} \
        -I ~{sep=" -I " tempVcfs} \
        -O ~{sortedVcfPath}
    }

    output {
        IndexedVcf sortedVcf = object {
                vcf : "~{sortedVcfPath}", 
                vcfIndex : "~{sortedVcfPath}.idx"
            }
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}






