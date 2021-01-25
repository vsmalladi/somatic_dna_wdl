version 1.0

import "../wdl_structs.wdl"

# General tasks

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


