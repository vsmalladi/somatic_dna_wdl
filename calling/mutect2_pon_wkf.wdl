version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Mutect2Pon {
    # command
    #   run Mutect2 caller
    input {
        String tumor
        Array[String]+ listOfChroms
        IndexedReference referenceFa
        # Add when Exome are added
        # Map[String, File] chromBeds
        Bam tumorFinalBam
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB")) + 20
        Int memoryGb = 8
        File mutectJsonLog
        File mutectJsonLogFilter
        
        Boolean highMem = false

    }
    
    Int lowCallMemoryGb = 4
    Int lowFilterMemoryGb = 4
    Int lowFilterDiskSize = 5
    
    if (highMem) {
        Int highCallMemoryGb = 8
        Int highFilterMemoryGb = 4
        Int highFilterDiskSize = 10
    }
    Int callMemoryGb = select_first([highCallMemoryGb, lowCallMemoryGb])
    Int filterMemoryGb = select_first([highFilterMemoryGb, lowFilterMemoryGb])
    Int filterDiskSize = select_first([highFilterDiskSize, lowFilterDiskSize])

    scatter(chrom in listOfChroms) {
        call calling.Mutect2WgsPon {
            input:
                chrom = chrom,
                tumor = tumor,
                referenceFa = referenceFa,
                tumorFinalBam = tumorFinalBam,
                memoryGb = callMemoryGb,
                diskSize = diskSize
        }

        call calling.Mutect2Filter {
            input:
                chrom = chrom,
                pairName = tumor,
                referenceFa = referenceFa,
                # mutect2ChromRawVcf = select_first([Mutect2WgsPon.mutect2ChromRawVcf])
                mutect2ChromRawVcf = Mutect2WgsPon.mutect2ChromRawVcf,
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
