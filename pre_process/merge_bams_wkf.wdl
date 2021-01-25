version 1.0

import "merge_flowcell.wdl" as mergeFlowcells
import "../wdl_structs.wdl"

workflow MergeBams {
    # command 
    #     merge lane level BAMs
    input {
        #    command merge flowcell
        Array[Bam] laneFixmateBams
        Array[String] laneFixmateBamPaths
        String sampleId
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf DbSnp
        File chromFile
        IndexedReference indexedReference
        # resources
        Int mem
        Int threads
        String gatkDockerImage
        String novosortDockerImage
        String samtoolsDockerImage
    }
    
    call mergeFlowcells.NovosortMarkDup {
        input:
            laneBams = laneFixmateBams,
            laneBamsLists = laneFixmateBamPaths,
            sampleId = sampleId,
            dockerImage = novosortDockerImage,
            mem = mem,
            threads = threads
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
            indexedReference = indexedReference,
            Indels = Indels,
            DbSnp = DbSnp,
            chromFile = chromFile,
            sampleId = sampleId,
            mem = mem,
            threads = threads,
            dockerImage = gatkDockerImage
            
    }
    call mergeFlowcells.PrintReads {
        input:
            indexedReference = indexedReference,
            mergedDedupBam = novosortMarkDupIndexed.indexedBam,
            recalGrp = Bqsr38.recalGrp,
            sampleId = sampleId,
            mem = mem,
            threads = threads,
            dockerImage = gatkDockerImage
            
    }
    
    output {
        Bam mergedDedupBam = novosortMarkDupIndexed.indexedBam
        Bam finalBam = PrintReads.finalBam
    }
}