version 1.0

import "alignment_analysis.wdl" as alignmentAnalysis
import "../wdl_structs.wdl"

workflow Kourami {
    input {
        String sampleId
        # mergedHlaPanel
        BwaReference kouramiReference
        File kouramiFastaGem1Index
        Bam finalBam
        Int diskSize = ceil( size(finalBam.bam, "GB")) + 20
    }
    
    call alignmentAnalysis.GetChr6Contigs {
        input:
            finalBam = finalBam,
            diskSize = diskSize
    }
    
    call alignmentAnalysis.GemSelect {
        input:
            sampleId = sampleId,
            finalBam = finalBam,
            chr6Contigs = GetChr6Contigs.chr6Contigs,
            kouramiFastaGem1Index = kouramiFastaGem1Index,
            diskSize = diskSize
    }
    
    call alignmentAnalysis.LookUpMates {
        input:
            sampleId = sampleId,
            r2File = GemSelect.r2File,
            r2MappedFastq = GemSelect.r2MappedFastq,
            r1File = GemSelect.r1File,
            r1MappedFastq = GemSelect.r1MappedFastq
    }
    
    call alignmentAnalysis.GetMates {
        input:
            sampleId = sampleId,
            finalBam = finalBam,
            r1UnmappedFile = LookUpMates.r1UnmappedFile,
            r2UnmappedFile = LookUpMates.r2UnmappedFile,
            diskSize = diskSize
    }
    
    call alignmentAnalysis.SortFastqs as r1SortFastqs {
        input:
            fastqPairId = "R1",
            sampleId = sampleId,
            chr6MappedFastq = GemSelect.r1MappedFastq,
            chr6MappedMatesFastq = GetMates.r1UnmappedFastq
    }
    
    call alignmentAnalysis.SortFastqs as r2SortFastqs {
        input:
            fastqPairId = "R2",
            sampleId = sampleId,
            chr6MappedFastq = GemSelect.r2MappedFastq,
            chr6MappedMatesFastq = GetMates.r2UnmappedFastq
    }
    
    call alignmentAnalysis.AlignToPanel {
        input:
            sampleId = sampleId,
            kouramiReference = kouramiReference,
            r1SortedFastq = r1SortFastqs.sortedFastq,
            r2SortedFastq = r2SortFastqs.sortedFastq
    }
    
    call alignmentAnalysis.Kourami {
        input:
            sampleId = sampleId,
            kouramiBam = AlignToPanel.kouramiBam
    }
    
    output {
        File result = Kourami.result
    }
}