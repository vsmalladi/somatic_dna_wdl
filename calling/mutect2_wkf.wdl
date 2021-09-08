version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Mutect2 {
    # command
    #   run Mutect2 caller
    input {
        String tumor
        String normal
        Array[String]+ listOfChroms
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        # Add when Exome are added
        # Map[String, File] chromBeds
        Bam tumorFinalBam
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 20
        Int memoryGb = 8
        File mutectJsonLog
        File mutectJsonLogFilter
        
        Boolean highMem = false

    }
    
    Int callMemoryGb = 4
    Int filterMemoryGb = 4
    Int filterDiskSize = 5
    
    if (highMem) {
        Int callMemoryGb = 8
        Int filterMemoryGb = 4
        Int filterDiskSize = 10
    }

    scatter(chrom in listOfChroms) {
        call calling.Mutect2Wgs {
            input:
                chrom = chrom,
                tumor = tumor,
                normal = normal,
                pairName = pairName,
                referenceFa = referenceFa,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                memoryGb = callMemoryGb,
                diskSize = diskSize
        }

        call calling.Mutect2Filter {
            input:
                chrom = chrom,
                pairName = pairName,
                referenceFa = referenceFa,
                # mutect2ChromRawVcf = select_first([Mutect2Wgs.mutect2ChromRawVcf, Mutect2Exome.mutect2ChromRawVcf])
                mutect2ChromRawVcf = Mutect2Wgs.mutect2ChromRawVcf,
                memoryGb = callMemoryGb,
                diskSize = filterDiskSize
        }
    }

    # filtered
    call calling.Gatk4MergeSortVcf as filteredGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{pairName}.mutect2.v4.0.5.1.sorted.vcf",
            tempChromVcfs = Mutect2Filter.mutect2ChromVcf,
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }

    call calling.AddVcfCommand as filteredAddVcfCommand {
        input:
            inVcf = filteredGatk4MergeSortVcf.sortedVcf.vcf,
            jsonLog = mutectJsonLogFilter,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as filteredReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = filteredAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.mutect2.v4.0.5.1.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    # unfiltered
    call calling.Gatk4MergeSortVcf as unfilteredGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{pairName}.mutect2.v4.0.5.1.unfiltered.sorted.vcf",
            tempChromVcfs = Mutect2Wgs.mutect2ChromRawVcf,
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }

    call calling.AddVcfCommand as unfilteredAddVcfCommand {
        input:
            inVcf = unfilteredGatk4MergeSortVcf.sortedVcf.vcf,
            jsonLog = mutectJsonLog,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as unfilteredReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = unfilteredAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.mutect2.v4.0.5.1.unfiltered.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    output {
        File mutect2 = filteredReorderVcfColumns.orderedVcf
        File mutect2_unfiltered = unfilteredReorderVcfColumns.orderedVcf
    }
}
