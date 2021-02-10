version 1.0

import "pre_process/pre_process.wdl" as preProcess
import "wdl_structs.wdl"

workflow Preprocess {
    # command 
    #   align FASTQ files
    #   and 
    #   merge lane level BAMs
    input {
        Array[sampleInfo]+ sampleInfos
        BwaReference bwaReference
        #    command merge flowcell
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        IndexedTable callRegions
        IndexedReference referenceFa
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
    scatter(sampleInfo in sampleInfos) {
        call preProcess.Preprocess {
            input:
                sampleId = sampleInfo.sampleId,
                listOfFastqPairs = sampleInfo.listOfFastqPairs,
                bwaReference = bwaReference,
                mem = mem,
                threads = threads,
                bwaDockerImage = bwaDockerImage,
                shortAlignDockerImage = shortAlignDockerImage,
                gatkDockerImage = gatkDockerImage,
                novosortDockerImage = novosortDockerImage,
                samtoolsDockerImage = samtoolsDockerImage,
                MillsAnd1000G = MillsAnd1000G,
                Indels = Indels,
                dbsnp = dbsnp,
                callRegions = callRegions,
                referenceFa = referenceFa
        }
    }
    
    output {
        Array[Bam] mergedDedupBam = Preprocess.mergedDedupBam
        Array[Bam] finalBam = Preprocess.finalBam
    }
}