version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow SvabaPon {
    # command
    #   run Svaba caller
    input {
        String tumor
        String sampleId
        BwaReference bwaReference
        IndexedTable callRegions
        IndexedReference referenceFa
        File dbsnpIndels
        Bam tumorFinalBam
        # resources
        Int threads = 8
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB")) + 40
        Int memoryGb = 32
        File svabaJsonLog
        
        Boolean highMem = false
    }
    
    Int lowCallMemoryGb = 16
    
    if (highMem) {
        Int highCallMemoryGb = 32
    }
    Int callMemoryGb = select_first([highCallMemoryGb, lowCallMemoryGb])

    call calling.SvabaWgsPon {
        input:
            sampleId=tumor,
            bwaReference = bwaReference,
            callRegions = callRegions,
            tumorFinalBam = tumorFinalBam,
            dbsnpIndels = dbsnpIndels,
            memoryGb = callMemoryGb,
            threads = threads,
            diskSize = diskSize

    }

    # indel
    call calling.ReheaderVcf as indelReheaderVcf {
        input:
            inVcf = SvabaWgsPon.svabaIndelGz,
            sampleIds = [tumor],
            referenceFa = referenceFa,
            outVcfPath = "~{tumor}.indel.svaba.v0.2.1.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    output {
        File svabaIndel = indelReheaderVcf.reheaderedVcf
    }
}
