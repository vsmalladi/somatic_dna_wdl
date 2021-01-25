version 1.0

import "merge_vcf.wdl" as merge_vcf
import "../wdl_structs.wdl"

workflow MergeVcf {
    input {
        String tumor
        String normal
        String tool
        File callerVcf
        String pairName
        String pysamDockerImage
        Int threads
        Int memory_gb
    }
    
    call merge_vcf.RenameMetadata {
        input:
            callerVcf = callerVcf,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    if (tool == 'manta') {
        call merge_vcf.MergePrepSupport {
            input:
                tool = tool,
                renameMetaVcf = RenameMetadata.renameMetaVcf,
                memory_gb = memory_gb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
    }
    if (tool != 'manta') {
        call merge_vcf.MergePrep {
            input:
                tool = tool,
                renameMetaVcf = RenameMetadata.renameMetaVcf,
                memory_gb = memory_gb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
        call merge_vcf.RenameVcf {
        input:
            pairName = pairName,
            tumor = tumor,
            normal = normal,
            tool = tool,
            prepCallerVcf = MergePrep.prepCallerVcf,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
        }
    }
    call merge_vcf.RenameVcf {
        input:
            pairName = pairName,
            tumor = tumor,
            normal = normal,
            tool = tool,
            prepCallerVcf = MergePrep.prepCallerVcf,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
}