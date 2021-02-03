version 1.0

import "alignment_analysis.wdl" as alignment_analysis
import "../wdl_structs.wdl"

workflow Kourami {
    input {
        String sampleId
        String chr6Contigs
        File mergedHlaPanel
        File hlaGem
        Bam finalBam
        
        String memoryPerThread = "1G"
        
        String gemDockerImage
        String samtoolsDockerImage
        String pythonDockerImage
        String seqkitDockerImage
        String bwaDockerImage
        String kouramiDockerImage
        Int threads
        Int memoryGb
    }
    
    call alignment_analysis.GemSelect {
        input:
            sampleId = sampleId,
            chr6Contigs = chr6Contigs,
            finalBam = finalBam,
            hlaGem = hlaGem,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gemDockerImage
            
    }
    
    call alignment_analysis.LookUpMates {
        input:
            sampleId = sampleId,
            r2File = GemSelect.r2File,
            r2MappedFastq = GemSelect.r2MappedFastq,
            r1File = GemSelect.r1File,
            r1MappedFastq = GemSelect.r1MappedFastq,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pythonDockerImage
    }
    
    call alignment_analysis.GetMates {
        input:
            sampleId = sampleId,
            finalBam = finalBam,
            r1UnmappedFile = LookUpMates.r1UnmappedFile,
            r2UnmappedFile = LookUpMates.r2UnmappedFile,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = samtoolsDockerImage
    }
    
    call alignment_analysis.SortFastqs as r1SortFastqs {
        input:
            fastqPairId = "R1",
            sampleId = sampleId,
            chr6MappedFastq = GemSelect.r1MappedFastq,
            chr6MappedMatesFastq = GetMates.r1UnmappedFastq,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = seqkitDockerImage
    }
    
    call alignment_analysis.SortFastqs as r2SortFastqs {
        input:
            fastqPairId = "R2",
            sampleId = sampleId,
            chr6MappedFastq = GemSelect.r2MappedFastq,
            chr6MappedMatesFastq = GetMates.r2UnmappedFastq,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = seqkitDockerImage
    }
    
    call alignment_analysis.AlignToPanel {
        input:
            sampleId = sampleId,
            mergedHlaPanel = mergedHlaPanel,
            memoryPerThread = memoryPerThread,
            r1SortedFastq = r1SortFastqs.sortedFastq,
            r2SortedFastq = r2SortFastqs.sortedFastq,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bwaDockerImage
    }
    
    call alignment_analysis.Kourami {
        input:
            sampleId = sampleId,
            kouramiBam = AlignToPanel.kouramiBam,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = kouramiDockerImage
    }
    
    output {
        File result = Kourami.result
    }
}