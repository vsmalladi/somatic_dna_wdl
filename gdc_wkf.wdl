version 1.0

import "test/tests.wdl" as tests
import "wdl_structs.wdl"
import "test/gdc_wkf.wdl" as gdcComparePair

struct GdcRelationship {
    String pairId
    String tumor
    String normal
    File nygcVcf
    File gdcSnvVcf
    File gdcIndelVcf
}


workflow GdcCompare {
    input {
        String name = "gdc_comparison"
        File chromLengths
        File cosmicCensus
        Array[String] listOfChroms
        
        Array[GdcRelationship] gdcRelationships
        IndexedReference referenceFa
    }
 
    
    scatter (gdcRelationship in gdcRelationships) {
        
        call gdcComparePair.GdcComparePair {
            input:
                pairId = gdcRelationship.pairId,
                nygcVcf = gdcRelationship.nygcVcf,
                gdcIndelVcf = gdcRelationship.gdcIndelVcf,
                gdcSnvVcf = gdcRelationship.gdcSnvVcf,
                referenceFa = referenceFa
        }
    }
    
    
    call tests.ConcateTables as finalVcfConcateTables {
        input:
            tables = GdcComparePair.summaryOutputTable,
            outputTablePath = "~{name}.finalVcf.csv"
    }
    
    
    output {
        File finalVcfTable = finalVcfConcateTables.outputTable
        Array[File] detailedOutputVcfTable = GdcComparePair.detailedOutputTable

    }
}

