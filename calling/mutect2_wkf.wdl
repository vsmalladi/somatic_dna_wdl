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
        Int threads
        Int memoryGb
        String mutect2GatkDockerImage
        String pysamDockerImage
        String gatkDockerImage
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        # Add when Exome are added
        # Map[String, File] chromBeds
        Bam tumorFinalBam
        File jsonLog
    }
    
    scatter(chrom in listOfChroms) {
        call calling.Mutect2Wgs {
            input:
                chrom = chrom,
                referenceFa = referenceFa,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = mutect2GatkDockerImage
        }
        
        call calling.Mutect2Filter {
            input:
                chrom = chrom,
                referenceFa = referenceFa,
                # mutect2ChromRawVcf = select_first([Mutect2Wgs.mutect2ChromRawVcf, Mutect2Exome.mutect2ChromRawVcf])
                mutect2ChromRawVcf = Mutect2Wgs.mutect2ChromRawVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = mutect2GatkDockerImage
        }
    }
    
    # filtered
    call calling.Gatk4MergeSortVcf as filteredGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{pairName}.mutect2.v4.0.5.1.sorted.vcf",
            tempChromVcfs = Mutect2Filter.mutect2ChromVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call calling.AddVcfCommand as filteredAddVcfCommand {
        input:
            inVcf = filteredGatk4MergeSortVcf.sortedVcf.vcf,
            jsonLog = jsonLog,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call calling.ReorderVcfColumns as filteredReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = filteredAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.mutect2.v4.0.5.1.vcf",
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    # unfiltered
    call calling.Gatk4MergeSortVcf as unfilteredGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{pairName}.mutect2.v4.0.5.1.unfiltered.sorted.vcf",
            tempChromVcfs = Mutect2Wgs.mutect2ChromRawVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call calling.AddVcfCommand as unfilteredAddVcfCommand {
        input:
            inVcf = unfilteredGatk4MergeSortVcf.sortedVcf.vcf,
            jsonLog = jsonLog,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call calling.ReorderVcfColumns as unfilteredReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = unfilteredAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.mutect2.v4.0.5.1.unfiltered.vcf",
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    output {
        File mutect2 = filteredReorderVcfColumns.orderedVcf
        File mutect2_unfiltered = unfilteredReorderVcfColumns.orderedVcf
    }
}

