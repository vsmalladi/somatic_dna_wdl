version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Mutect2Pon {
    # command
    #   run Mutect2 caller
    input {
        String library
        String tumor
        Array[String]+ callerIntervals
        File invertedIntervalListBed
        IndexedReference referenceFa
        # Files for when Exome PONs are added
        Bam tumorFinalBam
        Int diskSize = 10
        Int memoryGb = 8
        File mutectJsonLog
        File mutectJsonLogFilter
        
        Boolean highMem = false

    }
    
    Int lowCallMemoryGb = 4
    Int lowFilterMemoryGb = 4
    Int lowFilterDiskSize = 10
    
    if (highMem) {
        Int highCallMemoryGb = 8
        Int highFilterMemoryGb = 4
        Int highFilterDiskSize = 20
    }
    Int callMemoryGb = select_first([highCallMemoryGb, lowCallMemoryGb])
    Int filterMemoryGb = select_first([highFilterMemoryGb, lowFilterMemoryGb])
    Int filterDiskSize = select_first([highFilterDiskSize, lowFilterDiskSize])

    scatter(callerInterval in callerIntervals) {
        if (library == 'WGS') {
            call calling.Mutect2WgsPon {
                input:
                    chrom = callerInterval,
                    tumor = tumor,
                    referenceFa = referenceFa,
                    tumorFinalBam = tumorFinalBam,
                    memoryGb = callMemoryGb,
                    diskSize = diskSize
            }
        }
        
        if (library == 'Exome') {
            call calling.Mutect2ExomePon {
                input:
                    chrom = callerInterval,
                    invertedIntervalListBed = invertedIntervalListBed,
                    tumor = tumor,
                    referenceFa = referenceFa,
                    tumorFinalBam = tumorFinalBam,
                    memoryGb = callMemoryGb,
                    diskSize = diskSize
            }
        }
        
        File mutect2ChromRawVcfInput = select_first([Mutect2WgsPon.mutect2ChromRawVcf, Mutect2ExomePon.mutect2ChromRawVcf])
        File mutect2ChromRawStats = select_first([Mutect2WgsPon.mutect2ChromRawStats, Mutect2ExomePon.mutect2ChromRawStats])

        call calling.Mutect2Filter {
            input:
                chrom = callerInterval,
                pairName = tumor,
                referenceFa = referenceFa,
                mutect2ChromRawStats = mutect2ChromRawStats,
                mutect2ChromRawVcf = mutect2ChromRawVcfInput,
                memoryGb = callMemoryGb,
                diskSize = filterDiskSize
        }
    }

    # filtered
    call calling.Gatk4MergeSortVcf as filteredGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{tumor}.mutect2.sorted.vcf",
            tempChromVcfs = Mutect2Filter.mutect2ChromVcf,
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }

    output {
        File mutect2 = filteredGatk4MergeSortVcf.sortedVcf.vcf
    }
}
