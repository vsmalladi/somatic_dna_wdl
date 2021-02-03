version 1.0

import "merge_vcf.wdl" as merge_vcf
import "../wdl_structs.wdl"

workflow PrepMergeVcf {
    input {
        String tumor
        String normal
        String tool
        File callerVcf
        String pairName
        String pysamDockerImage
        String gatkDockerImage
        String bgzipDockerImage
        String bcftoolsDockerImage
        Int threads
        Int memoryGb
        IndexedReference referenceFa
    }
    
    call merge_vcf.RenameMetadata {
        input:
            callerVcf = callerVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    if (tool == 'manta') {
        call merge_vcf.MergePrepSupport {
            input:
                tool = tool,
                renameMetaVcf = RenameMetadata.renameMetaVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
    }
    if (tool != 'manta') {
        call merge_vcf.MergePrep {
            input:
                tool = tool,
                renameMetaVcf = RenameMetadata.renameMetaVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }

    }
    
    call merge_vcf.RenameVcf {
        input:
            pairName = pairName,
            tumor = tumor,
            normal = normal,
            tool = tool,
            prepCallerVcf = select_first([MergePrepSupport.prepCallerVcf, MergePrep.prepCallerVcf]),
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call merge_vcf.CompressVcf as renameCompressVcf {
        input:
            vcf = RenameVcf.renameVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bgzipDockerImage
    }
    
    call merge_vcf.IndexVcf as renameIndexVcf {
        input:
            vcfCompressed = renameCompressVcf.vcfCompressed,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call merge_vcf.SplitMultiAllelic as prepSplitMultiAllelic {
        input:
            vcfCompressedIndexed = renameIndexVcf.vcfCompressedIndexed,
            splitVcfPath = sub(renameIndexVcf.vcfCompressedIndexed.vcf, ".rename.vcf.gz$", ".split.vcf"),
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bcftoolsDockerImage
    }
    
    call merge_vcf.SplitMnv {
        input:
            tool = tool,
            mnvVcfPath = sub(renameIndexVcf.vcfCompressedIndexed.vcf, ".rename.vcf.gz$", ".split.vcf"),
            splitVcf = prepSplitMultiAllelic.splitVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
            
    }
    
    if (tool != 'svaba') {
        call merge_vcf.RemoveContig {
            input:
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage, 
                mnvVcfPath = sub(renameIndexVcf.vcfCompressedIndexed.vcf, ".rename.vcf.gz$", ".split.vcf"),
                removeChromVcf = SplitMnv.mnvVcf
        }
    }
    
    call merge_vcf.Gatk4MergeSortVcf {
        input:
            tempVcfs = [select_first([RemoveContig.removeContigVcf, SplitMnv.mnvVcf])],
            sortedVcfPath = sub(select_first([RemoveContig.removeContigVcf, SplitMnv.mnvVcf]), "$", ".gz"),
            referenceFa = referenceFa,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
            
    }
    
    output {
        IndexedVcf preppedVcf = Gatk4MergeSortVcf.sortedVcf
    }
}
