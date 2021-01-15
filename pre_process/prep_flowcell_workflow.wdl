version 1.0

import "prep_flowcell.wdl" as prepFlowcell
import "merge_flowcell.wdl" as mergeFlowcells
import "wdl_structs.wdl"

workflow RunFlowcell {
    # command 
    #     prep flowcell
    input {
        Array[Fastqs]+ listOfFastqPairs
        BwaReference bwaReference
        #    command merge flowcell
        String sampleId
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf DbSnp
        File chromFile
        IndexedReference indexedReference
        # resources
        #    prep flowcell
        Int mem
        Int threads
        String bwaDockerImage
        String shortAlignDockerImage
        String gatkDockerImage
        #    merge flowcell
        String novosortDockerImage
        String samtoolsDockerImage
    }
    
    scatter(fastqs in listOfFastqPairs) {
        call prepFlowcell.AlignBwaMem {
            input:
                fastqs = fastqs,
                bwaReference = bwaReference,
                mem = mem,
                threads = threads,
                dockerImage = bwaDockerImage
        }
        
        call prepFlowcell.ShortAlignMark {
            input:
                laneBam = AlignBwaMem.laneBam,
                bamBase = fastqs.bamBase,
                mem = mem,
                dockerImage = shortAlignDockerImage
        }
        
        call prepFlowcell.Fixmate {
            input:
                laneBamMark = ShortAlignMark.laneBamMark,
                bamBase = fastqs.bamBase,
                mem = mem,
                dockerImage = gatkDockerImage
        }
    }
    
    call mergeFlowcells.NovosortMarkDup {
        input:
            laneBams = Fixmate.laneFixmateBam,
            laneBamsLists = Fixmate.laneFixmateBamPath,
            sampleId = sampleId,
            dockerImage = novosortDockerImage
    }
    
    call mergeFlowcells.IndexBam as novosortMarkDupIndexed {
        input:
            bam = NovosortMarkDup.mergedDedupBamOnly,
            dockerImage = samtoolsDockerImage
    }
    
    call mergeFlowcells.Bqsr38 {
        input:
            mergedDedupBam = novosortMarkDupIndexed.indexedBam,
            MillsAnd1000G = MillsAnd1000G,
            Indels = Indels,
            DbSnp = DbSnp,
            chromFile = chromFile,
            sampleId = sampleId
            
    }
    call mergeFlowcells.PrintReads {
        input:
            indexedReference = indexedReference,
            mergedDedupBam = novosortMarkDupIndexed.indexedBam,
            recalGrp = Bqsr38.recalGrp,
            sampleId = sampleId
            
    }
    #output {
    #    Bam finalBam = PrintReads.finalBam
    #}
}
