version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Strelka2 {
    # command 
    #   run Strelka2 caller
    input {
        String tumor
        String normal
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        IndexedVcf candidateSmallIndels
        Bam normalFinalBam
        IndexedTable callRegions
        Bam tumorFinalBam
        String pairName
        String strelka2SnvPath = "~{pairName}.snv.strelka2.v2.9.3.vcf"
        String strelka2IndelPath = "~{pairName}.indel.strelka2.v2.9.3.vcf"
        # resources
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 20
        Int memoryGb = 4
        Int threads = 8
        # remove definition after replacing the command step for gcp
        File jsonLog = "gs://nygc-comp-s-fd4e-input/internal/strelka_v2.9.3_Strelka2.json"
    }
    
    call calling.Strelka2 {
        input:
            intHVmem = intHVmem,
            referenceFa = referenceFa,
            callRegions = callRegions,
            pairName = pairName,
            candidateSmallIndels = candidateSmallIndels,
            normalFinalBam = normalFinalBam,
            tumorFinalBam = tumorFinalBam,
            memoryGb = memoryGb,
            diskSize = diskSize,
            threads = threads
    }
    
    call calling.AddVcfCommand as strelka2SnvAddVcfCommand {
        input:
            inVcf = Strelka2.strelka2Snvs.vcf,
            jsonLog = jsonLog,
            memoryGb = 2,
            diskSize = 1
    }
    
    call calling.ReorderVcfColumns as strelka2SnvReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = strelka2SnvAddVcfCommand.outVcf,
            orderedVcfPath = strelka2SnvPath,
            memoryGb = 2,
            diskSize = 1
    }
    
    call calling.AddVcfCommand as strelka2IndelAddVcfCommand {
        input:
            inVcf = Strelka2.strelka2Indels.vcf,
            jsonLog = jsonLog,
            memoryGb = 2,
            diskSize = 1
    }
    
    call calling.ReorderVcfColumns as strelka2IndelReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = strelka2IndelAddVcfCommand.outVcf,
            orderedVcfPath = strelka2IndelPath,
            memoryGb = 2,
            diskSize = 1
    }
    
    output {
        IndexedVcf strelka2Snvs = Strelka2.strelka2Snvs
        IndexedVcf strelka2Indels = Strelka2.strelka2Indels
        File strelka2Snv = strelka2SnvReorderVcfColumns.orderedVcf
        File strelka2Indel = strelka2IndelReorderVcfColumns.orderedVcf
    }
}