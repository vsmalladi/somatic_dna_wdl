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
        String strelka2SnvPath = "~{pairName}.snv.strelka2.vcf"
        String strelka2IndelPath = "~{pairName}.indel.strelka2.vcf"
        # resources
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 20
        Int memoryGb = 4
        Int threads = 8
        File configureStrelkaSomaticWorkflow
        File strelkaJsonLog
    }

    call calling.Strelka2 {
        input:
            configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
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
            jsonLog = strelkaJsonLog,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as strelka2SnvReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = strelka2SnvAddVcfCommand.outVcf,
            orderedVcfPath = strelka2SnvPath,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.AddVcfCommand as strelka2IndelAddVcfCommand {
        input:
            inVcf = Strelka2.strelka2Indels.vcf,
            jsonLog = strelkaJsonLog,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as strelka2IndelReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = strelka2IndelAddVcfCommand.outVcf,
            orderedVcfPath = strelka2IndelPath,
            memoryGb = 4,
            diskSize = 10
    }

    output {
        IndexedVcf strelka2Snvs = Strelka2.strelka2Snvs
        IndexedVcf strelka2Indels = Strelka2.strelka2Indels
        File strelka2Snv = strelka2SnvReorderVcfColumns.orderedVcf
        File strelka2Indel = strelka2IndelReorderVcfColumns.orderedVcf
    }
}
