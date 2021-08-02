version 1.0

import "../wdl_structs.wdl"
import "baf.wdl"
import "../test/tests.wdl" as tests


workflow Baf {
    # command 
    input {
        Array[String]+ listOfChroms
        String sampleId
        String pairName
        Bam normalFinalBam
        Bam tumorFinalBam
        File? germlineVcf
        IndexedReference referenceFa
    }
    
    if ( size(germlineVcf) > 0 ) {
        
        call baf.FilterForHetSnps {
            input:
                sampleId = sampleId,
                referenceFa = referenceFa,
                germlineVcf = germlineVcf
        }
        
        call baf.FilterBaf {
            input:
                sampleId = sampleId,
                hetVcf = FilterForHetSnps.hetVcf
        }
        
        scatter(chrom in listOfChroms) {
            call baf.AlleleCounts {
                input:
                    pairName = pairName,
                    referenceFa = referenceFa,
                    normalFinalBam = normalFinalBam,
                    tumorFinalBam = tumorFinalBam,
                    knownHetVcf = FilterBaf.knownHetVcf,
                    chrom = chrom
            }
            
            call baf.CalcBaf {
                input:
                    pairName = pairName,
                    alleleCountsTxt = AlleleCounts.alleleCountsTxt
            }
        }
        
        call tests.ConcateTables {
            input:
                tables = CalcBaf.bafTxt,
                outputTablePath = "~{pairName}.haplotypecaller.gatk.v4.1.8.0.alleles.txt"
        }
    }
    
    output {
        File? alleleCountsTxt = ConcateTables.outputTable
    }
}