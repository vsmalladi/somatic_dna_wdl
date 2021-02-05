version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Lancet {
    # command 
    #   run Lancet caller
    input {
        String tumor
        String normal
        Array[String]+ listOfChroms
        Map[String, File] chromBeds
        Int threads
        Int memoryGb
        String lancetDockerImage
        String dockerImage
        String gatkDockerImage
        String pysamDockerImage
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
        File jsonLog
    }
    
    scatter(chrom in listOfChroms) {
        call calling.LancetExome {
            input:
                chrom = chrom,
                chromBed = chromBeds[chrom],
                referenceFa = referenceFa,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = lancetDockerImage
        }
    }
    
    call calling.Gatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{pairName}.lancet.v1.0.7.sorted.vcf",
            tempChromVcfs = LancetExome.lancetChromVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call calling.AddVcfCommand as lancetAddVcfCommand {
        input:
            inVcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            jsonLog = jsonLog,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call calling.ReorderVcfColumns as lancetReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = lancetAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.lancet.v1.0.7.vcf",
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    output {
        File lancet = lancetReorderVcfColumns.orderedVcf
    }
}