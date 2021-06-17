version 1.0

import "wdl_structs.wdl"
import "germline/germline_wkf.wdl" as germline
import "annotate/germline_annotate_wkf.wdl" as germlineAnnotate

workflow GermlineAll {
    # command 
    input {
        Array[SampleBamInfo]+ normalSampleBamInfos
        
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        
        File excludeIntervalList
        Array[File] scatterIntervalsHcs
        
        IndexedVcf MillsAnd1000G
        IndexedVcf omni
        IndexedVcf hapmap
        IndexedVcf onekG
        IndexedVcf dbsnp
        
        IndexedVcf whitelist
        IndexedVcf nygcAf
        IndexedVcf pgx
        IndexedTable rwgsPgxBed
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf chdWhitelistVcf
        
        # annotation:
        String vepGenomeBuild
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding
        
        # Public
        File cancerResistanceMutations
        File vepCache
        File annotations
        File plugins
        String vepGenomeBuild
        IndexedReference vepFastaReference
        
        # NYGC-only
        IndexedVcf hgmdGene
        IndexedVcf hgmdUd10
        IndexedVcf hgmdPro
        IndexedVcf omimVcf
        
        # Public
        IndexedVcf chdGenesVcf
        IndexedVcf chdEvolvingGenesVcf
        IndexedVcf chdWhitelistVcf
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf masterMind
        
        # post annotation
        File cosmicCensus
        
        File ensemblEntrez
        String library
    }
    
    scatter (normalSampleBamInfo in normalSampleBamInfos) {
        call germline.Germline {
            input:
                finalBam = normalSampleBamInfo.finalBam,
                normal = normalSampleBamInfo.sampleId,
                referenceFa = referenceFa,
                listOfChroms = listOfChroms,
                MillsAnd1000G = MillsAnd1000G,
                omni = omni,
                hapmap = hapmap,
                onekG = onekG,
                dbsnp = dbsnp,
                nygcAf = nygcAf,
                excludeIntervalList=excludeIntervalList,
                scatterIntervalsHcs=scatterIntervalsHcs,
                pgx=pgx,
                rwgsPgxBed=rwgsPgxBed,
                whitelist=whitelist,
                chdWhitelistVcf=chdWhitelistVcf,
                deepIntronicsVcf=deepIntronicsVcf,
                clinvarIntronicsVcf=clinvarIntronicsVcf
        }
        
    }
    
    output {
        Array[IndexedVcf] haplotypecallerVcf = Germline.haplotypecallerVcf 
        Array[IndexedVcf] haplotypecallerFinalFiltered = Germline.haplotypecallerFinalFiltered 
    }
    
}