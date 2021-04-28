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
        Map[String, File] chromBedsWgs
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
        # resources
        Int threads = 8
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 20
        Int memoryGb = 40
        # remove definition after replacing the command step for gcp
        File jsonLog = "gs://nygc-comp-s-fd4e-input/internal/lancet_1.0.7_LancetWGSRegional.json"
    }
    
    scatter(chrom in listOfChroms) {
        call calling.LancetWGSRegional {
            input:
                chrom = chrom,
                chromBed = chromBedsWgs[chrom],
                referenceFa = referenceFa,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                pairName = pairName,
                threads = threads,
                memoryGb = memoryGb,
                diskSize = diskSize
        }
    }
    
    call calling.Gatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{pairName}.lancet.v1.0.7.sorted.vcf",
            tempChromVcfs = LancetWGSRegional.lancetChromVcf,
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }
    
    call calling.AddVcfCommand as lancetAddVcfCommand {
        input:
            inVcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            jsonLog = jsonLog,
            memoryGb = 2,
            diskSize = 1
    }
    
    call calling.ReorderVcfColumns as lancetReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = lancetAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.lancet.v1.0.7.vcf",
            memoryGb = 2,
            diskSize = 1
    }
    
    output {
        File lancet = lancetReorderVcfColumns.orderedVcf
    }
}