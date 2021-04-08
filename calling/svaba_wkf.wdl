version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Svaba {
    # command 
    #   run Svaba caller
    input {
        String tumor
        String normal
        String pairName
        BwaReference bwaReference
        File dbsnpIndels
        Bam normalFinalBam
        Bam tumorFinalBam
        # resources
        Int threads = 8
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 40
        Int memoryGb = 32
        # remove definition after replacing the command step for gcp
        File jsonLog = "gs://nygc-comp-s-fd4e-input/mutect2_4.0.5.1_COLO-829-NovaSeq_80--COLO-829BL-NovaSeq_40.json"
    }
    
    call calling.SvabaWgs {
        input:
            pairName=pairName,
            bwaReference = bwaReference,
            normalFinalBam = normalFinalBam,
            tumorFinalBam = tumorFinalBam,
            dbsnpIndels = dbsnpIndels,
            memoryGb = memoryGb,
            threads = threads,
            diskSize = diskSize
            
    }
    
    # sv
    call calling.AddVcfCommand as svAddVcfCommand {
        input:
            inVcf = SvabaWgs.svabaGz,
            jsonLog = jsonLog,
            memoryGb = 2,
            diskSize = 1
    }
    
    call calling.ReorderVcfColumns as svReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = svAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.sv.svaba.v0.2.1.vcf",
            memoryGb = 2,
            diskSize = 1
    }
    
    # indel
    call calling.AddVcfCommand as indelAddVcfCommand {
        input:
            inVcf = SvabaWgs.svabaIndelGz,
            jsonLog = jsonLog,
            memoryGb = 2,
            diskSize = 1
    }
    
    call calling.ReorderVcfColumns as indelReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = indelAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.indel.svaba.v0.2.1.vcf",
            memoryGb = 2,
            diskSize = 1
    }
    
    output {
        File svabaRawGermlineIndel = SvabaWgs.svabaRawGermlineIndel
        File svabaRawGermlineSv = SvabaWgs.svabaRawGermlineSv
        File svabaSv = svReorderVcfColumns.orderedVcf
        File svabaIndel = indelReorderVcfColumns.orderedVcf
    }
}