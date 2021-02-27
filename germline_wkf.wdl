version 1.0

import "wdl_structs.wdl"
import "germline/germline_wkf.wdl" as germline

workflow GermlineAll {
    # command 
    input {
        Array[pairInfo]+ pairInfos
        
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        
        IndexedVcf MillsAnd1000G
        IndexedVcf omni
        IndexedVcf hapmap
        IndexedVcf onekG
        IndexedVcf dbsnp
    }
    
    scatter(pairInfo in pairInfos) {
        call germline.Germline {
            input:
                finalBam = pairInfo.normalFinalBam,
                normal = pairInfo.normal,
                referenceFa = referenceFa,
                listOfChroms = listOfChroms,
                MillsAnd1000G = MillsAnd1000G,
                omni = omni,
                hapmap = hapmap,
                onekG = onekG,
                dbsnp = dbsnp    
        }
    }
    
    output {
        Array[IndexedVcf] haplotypecallerNormVcf = Germline.haplotypecallerNormVcf
    }
    
}