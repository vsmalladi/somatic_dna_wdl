version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Manta {
    # command 
    #   run Manta caller
    input {
        String tumor
        String normal
        String intHVmem
        IndexedReference referenceFa
        File callRegions
        Bam normalFinalBam
        Bam tumorFinalBam
        File jsonLog
        # reorder columns
        String pairName 
        String mantaPath = "~{pairName}.manta.v1.4.0.vcf"
        String filteredMantafPath = "~{pairName}.manta.v1.4.0.filtered.vcf"
        # resources
        Int memory_gb
        Int threads
        String mantaDockerImage
        String pysamDockerImage
        String gatkDockerImage

    }
    
    call calling.MantaWgs {
        input:
            intHVmem = intHVmem,
            referenceFa = referenceFa,
            callRegions = callRegions,
            normalFinalBam = normalFinalBam,
            tumorFinalBam = tumorFinalBam,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = mantaDockerImage
    }
    
    call calling.AddVcfCommand as mantaAddVcfCommand {
        input:
            inVcf = MantaWgs.somaticSV.vcf,
            jsonLog = jsonLog,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call calling.ReorderVcfColumns as mantaReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = mantaAddVcfCommand.outVcf,
            orderedVcfPath = mantaPath,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call calling.FilterNonpass {
        input:
            referenceFa = referenceFa,
            pairName = pairName,
            vcf = mantaReorderVcfColumns.orderedVcf,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = gatkDockerImage
            
    }
    
    call calling.ReorderVcfColumns as mantaFilteredReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = FilterNonpass.outVcf,
            orderedVcfPath = FilterNonpass.outVcf,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
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