version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Manta {
    # command 
    #   run Manta caller
    input {
        String tumor
        String normal
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        IndexedTable callRegions
        Bam normalFinalBam
        Bam tumorFinalBam
        # reorder columns
        String pairName 
        String mantaPath = "~{pairName}.manta.v1.4.0.vcf"
        String filteredMantafPath = "~{pairName}.manta.v1.4.0.filtered.vcf"
        # resources
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 20
        Int memoryGb = 64
        Int threads = 8
        # remove definition after replacing the command step for gcp
        File jsonLog = "gs://nygc-comp-s-fd4e-input/mutect2_4.0.5.1_COLO-829-NovaSeq_80--COLO-829BL-NovaSeq_40.json"

    }
    
    call calling.MantaWgs {
        input:
            intHVmem = intHVmem,
            referenceFa = referenceFa,
            callRegions = callRegions,
            normalFinalBam = normalFinalBam,
            tumorFinalBam = tumorFinalBam,
            memoryGb = memoryGb,
            diskSize = diskSize,
            threads = threads
    }
    
    call calling.AddVcfCommand as mantaAddVcfCommand {
        input:
            inVcf = MantaWgs.somaticSV.vcf,
            jsonLog = jsonLog,
            memoryGb = 2,
            diskSize = 1
    }
    
    call calling.ReorderVcfColumns as mantaReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = mantaAddVcfCommand.outVcf,
            orderedVcfPath = mantaPath,
            memoryGb = 2,
            diskSize = 1
    }
    
    call calling.FilterNonpass {
        input:
            referenceFa = referenceFa,
            pairName = pairName,
            vcf = mantaReorderVcfColumns.orderedVcf,
            memoryGb = 8,
            threads = 4,
            diskSize = 5
            
    }
    
    call calling.ReorderVcfColumns as mantaFilteredReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = FilterNonpass.outVcf,
            orderedVcfPath = filteredMantafPath,
            memoryGb = 2,
            diskSize = 1
    }
    
    output {
        IndexedVcf candidateSmallIndels = MantaWgs.candidateSmallIndels
        IndexedVcf diploidSV = MantaWgs.diploidSV
        IndexedVcf somaticSV = MantaWgs.somaticSV
        IndexedVcf candidateSV = MantaWgs.candidateSV
        File unfilteredMantaSV = mantaReorderVcfColumns.orderedVcf
        File filteredMantaSV = mantaFilteredReorderVcfColumns.orderedVcf
    }
    
}