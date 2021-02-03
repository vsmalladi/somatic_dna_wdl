version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Mutect2 {
    # command 
    #   run Mutect2 caller
    input {
        String tumor
        String normal
        Int threads
        Int memoryGb
        String svabaDockerImage
        String gatkDockerImage
        String pysamDockerImage
        String pairName
        IndexedReference referenceFa
        File dbsnp
        Bam normalFinalBam
        Bam tumorFinalBam
        File jsonLog
    }
    
    call calling.SvabaWgs {
        input:
            pairName=pairName,
            referenceFa = referenceFa,
            normalFinalBam = normalFinalBam,
            tumorFinalBam = tumorFinalBam,
            dbsnp = dbsnp,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = svabaDockerImage
            
    }
    
    # sv
    call calling.AddVcfCommand as svAddVcfCommand {
        input:
            inVcf = SvabaWgs.svabaGz,
            jsonLog = jsonLog,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call calling.ReorderVcfColumns as svReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = svAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.sv.svaba.v0.2.1.vcf",
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    # indel
    call calling.AddVcfCommand as indelAddVcfCommand {
        input:
            inVcf = SvabaWgs.svabaIndelGz,
            jsonLog = jsonLog,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call calling.ReorderVcfColumns as indelReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = indelAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.indel.svaba.v0.2.1.vcf",
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
}