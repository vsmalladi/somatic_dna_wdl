version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Mutect2 {
    # command
    #   run Mutect2 caller
    input {
        Boolean local = false
        String library
        String tumor
        String normal
        Array[String]+ callerIntervals
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        # Exome
        File invertedIntervalListBed
        Bam tumorFinalBam
        Int diskSize = 20
        File mutectJsonLog
        File mutectJsonLogFilter
        
        Boolean highMem = false

    }
    
    Int lowCallMemoryGb = 20
    Int lowFilterMemoryGb = 20
    Int lowFilterDiskSize = 10
    
    if (highMem) {
        Int highCallMemoryGb = 40
        Int highFilterMemoryGb = 20
        Int highFilterDiskSize = 20
    }
    
    if (local) {
        Int localCallMemoryGb = 48
    }
    Int callMemoryGb = select_first([localCallMemoryGb, highCallMemoryGb, lowCallMemoryGb])
    Int filterMemoryGb = select_first([highFilterMemoryGb, lowFilterMemoryGb])
    Int filterDiskSize = select_first([highFilterDiskSize, lowFilterDiskSize])

    scatter(callerInterval in callerIntervals) {
        if (library == 'WGS') {
            call calling.Mutect2Wgs {
                input:
                    chrom = callerInterval,
                    tumor = tumor,
                    normal = normal,
                    pairName = pairName,
                    referenceFa = referenceFa,
                    normalFinalBam = normalFinalBam,
                    tumorFinalBam = tumorFinalBam,
                    memoryGb = callMemoryGb,
                    diskSize = diskSize
            }
        }
        
        if (library == 'Exome') {
            call calling.Mutect2Exome {
                input:
                    chrom = callerInterval,
                    tumor = tumor,
                    normal = normal,
                    pairName = pairName,
                    referenceFa = referenceFa,
                    normalFinalBam = normalFinalBam,
                    tumorFinalBam = tumorFinalBam,
                    invertedIntervalListBed = invertedIntervalListBed,
                    memoryGb = callMemoryGb,
                    diskSize = diskSize
            }
        }
        
        File mutect2ChromRawVcfInput = select_first([Mutect2Wgs.mutect2ChromRawVcf, Mutect2Exome.mutect2ChromRawVcf])
        File mutect2ChromRawStats = select_first([Mutect2Wgs.mutect2ChromRawStats, Mutect2Exome.mutect2ChromRawStats])
        

        call calling.Mutect2Filter {
            input:
                chrom = callerInterval,
                pairName = pairName,
                referenceFa = referenceFa,
                mutect2ChromRawVcf = mutect2ChromRawVcfInput,
                #mutect2ChromRawVcf = mutect2ChromRawVcfInput,
                mutect2ChromRawStats = mutect2ChromRawStats,
                memoryGb = filterMemoryGb,
                diskSize = filterDiskSize
        }
    }

    # filtered
    call calling.Gatk4MergeSortVcf as filteredGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{pairName}.mutect2.sorted.vcf",
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
            orderedVcfPath = "~{pairName}.mutect2.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    # unfiltered
    call calling.Gatk4MergeSortVcf as unfilteredGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{pairName}.mutect2.unfiltered.sorted.vcf",
            tempChromVcfs = mutect2ChromRawVcfInput,
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
            orderedVcfPath = "~{pairName}.mutect2.unfiltered.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    output {
        File mutect2 = filteredReorderVcfColumns.orderedVcf
        File mutect2_unfiltered = unfilteredReorderVcfColumns.orderedVcf
    }
}
