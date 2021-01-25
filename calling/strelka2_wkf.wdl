version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Strelka2 {
    # command 
    #   run Strelka2 caller
    input {
        String tumor
        String normal
        String intHVmem
        IndexedReference referenceFa
        IndexedVcf candidateSmallIndels
        Bam normalFinalBam
        File callRegions
        Bam tumorFinalBam
        File jsonLog
        String pairName
        String strelka2SnvPath = "~{pairName}.snv.strelka2.v2.9.3.vcf"
        String strelka2IndelPath = "~{pairName}.indel.strelka2.v2.9.3.vcf"
        # resources
        Int memory_gb
        Int threads
        String strelka2DockerImage
        String pysamDockerImage
        String gatkDockerImage
    }
    
    call calling.Strelka2 {
        input:
            intHVmem = intHVmem,
            referenceFa = referenceFa,
            callRegions = callRegions,
            candidateSmallIndels = candidateSmallIndels,
            normalFinalBam = normalFinalBam,
            tumorFinalBam = tumorFinalBam,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = strelka2DockerImage
    }
    
    call calling.AddVcfCommand as strelka2SnvAddVcfCommand {
        input:
            inVcf = Strelka2.strelka2Snvs.vcf,
            jsonLog = jsonLog,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call calling.ReorderVcfColumns as strelka2SnvReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = strelka2SnvAddVcfCommand.outVcf,
            orderedVcfPath = strelka2SnvPath,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call calling.AddVcfCommand as strelka2IndelAddVcfCommand {
        input:
            inVcf = Strelka2.strelka2Indels.vcf,
            jsonLog = jsonLog,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call calling.ReorderVcfColumns as strelka2IndelReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = strelka2IndelAddVcfCommand.outVcf,
            orderedVcfPath = strelka2IndelPath,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    output {
        IndexedVcf strelka2Snvs = Strelka2.strelka2Snvs
        IndexedVcf strelka2Indels = Strelka2.strelka2Indels
        File strelka2Snv = strelka2SnvReorderVcfColumns.orderedVcf
        File strelka2Indel = strelka2IndelReorderVcfColumns.orderedVcf
    }
}