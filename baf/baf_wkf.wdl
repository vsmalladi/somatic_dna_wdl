version 1.0

import "../wdl_structs.wdl"
import "../annotate/germline_annotate_wkf.wdl" as germlineAnnotate
import "baf.wdl"


workflow Baf {
    # command 
    input {
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
        
        call baf.AlleleCounts {
            input:
                pairName = pairName,
                referenceFa = referenceFa,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                knownHetVcf = FilterBaf.knownHetVcf
        }
        
        call baf.CalcBaf {
            input:
                pairName = pairName,
                alleleCountsTxt = AlleleCounts.alleleCountsTxt
        }
    }
    
    output {
        File? alleleCountsTxt = CalcBaf.bafTxt
    }
}