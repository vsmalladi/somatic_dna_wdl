version 1.0

import "annotate/annotate_wkf.wdl" as annotate
import "wdl_structs.wdl"

workflow Annotate {
    # command
    #   Annotate variants in VCFs
    input {
        Array[MergedPairVcfInfo]+ mergedPairVcfInfos
        IndexedReference referenceFa
        
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
    
        IndexedReference referenceFa

    }
    
    scatter(mergedPairVcfInfo in mergedPairVcfInfos) {
        call annotate.Annotate {
            input:
                unannotatedVcf = mergedPairVcfInfo.unannotatedVcf,
                referenceFa = referenceFa,
                tumor = mergedPairVcfInfo.tumor,
                normal = mergedPairVcfInfo.normal,
                pairName = mergedPairVcfInfo.pairId,
                vepGenomeBuild = vepGenomeBuild,
                cosmicCoding = cosmicCoding,
                cosmicNoncoding = cosmicNoncoding,
                # Public
                cancerResistanceMutations = cancerResistanceMutations,
                vepCache = vepCache,
                annotations = annotations,
                plugins = plugins,
                vepFastaReference = vepFastaReference,
                # NYGC-only
                hgmdGene = hgmdGene,
                hgmdUd10 = hgmdUd10,
                hgmdPro = hgmdPro,
                omimVcf = omimVcf,
                # Public
                chdGenesVcf = chdGenesVcf,
                chdEvolvingGenesVcf = chdEvolvingGenesVcf,
                chdWhitelistVcf = chdWhitelistVcf,
                deepIntronicsVcf = deepIntronicsVcf,
                clinvarIntronicsVcf = clinvarIntronicsVcf,
                masterMind = masterMind,
                # post annotation
                cosmicCensus = cosmicCensus,
                ensemblEntrez = ensemblEntrez,
                library = library
                
        }
    }
    
    output {   
        Array[PairVcfInfo] pairVcfInfos = Annotate.pairVcfInfo
    }
}