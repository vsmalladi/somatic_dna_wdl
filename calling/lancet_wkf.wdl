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
        File lancetJsonLog

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
            sortedVcfPath = "~{pairName}.lancet.sorted.vcf",
            tempChromVcfs = LancetWGSRegional.lancetChromVcf,
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }

    call calling.AddVcfCommand as lancetAddVcfCommand {
        input:
            inVcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            jsonLog = lancetJsonLog,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as lancetReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = lancetAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.lancet.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    output {
        File lancet = lancetReorderVcfColumns.orderedVcf
    }
}
