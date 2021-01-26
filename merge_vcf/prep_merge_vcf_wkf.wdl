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
        Int memory_gb
        IndexedReference referenceFa
    }
    
    call merge_vcf.RenameMetadata {
        input:
            callerVcf = callerVcf,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    if (tool == 'manta') {
        call merge_vcf.MergePrepSupport {
            input:
                tool = tool,
                renameMetaVcf = RenameMetadata.renameMetaVcf,
                memory_gb = memory_gb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
    }
    if (tool != 'manta') {
        call merge_vcf.MergePrep {
            input:
                tool = tool,
                renameMetaVcf = RenameMetadata.renameMetaVcf,
                memory_gb = memory_gb,
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
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call merge_vcf.CompressVcf as renameCompressVcf {
        input:
            vcf = RenameVcf.renameVcf,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = bgzipDockerImage
    }
    
    call merge_vcf.IndexVcf as renameIndexVcf {
        input:
            vcfCompressed = renameCompressVcf.vcfCompressed,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call merge_vcf.SplitMultiAllelic as prepSplitMultiAllelic {
        input:
            vcfCompressedIndexed = renameIndexVcf.vcfCompressedIndexed,
            splitVcfPath = sub(renameIndexVcf.vcfCompressedIndexed.vcf, ".rename.vcf.gz$", ".split.vcf"),
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = bcftoolsDockerImage
    }
    
    call merge_vcf.SplitMnv {
        input:
            tool = tool,
            mnvVcfPath = sub(renameIndexVcf.vcfCompressedIndexed.vcf, ".rename.vcf.gz$", ".split.vcf"),
            splitVcf = prepSplitMultiAllelic.splitVcf,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
            
    }
    
    if (tool != 'svaba') {
        call merge_vcf.RemoveContig {
            input:
                memory_gb = memory_gb,
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
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = gatkDockerImage
            
    }
    
    output {
        IndexedVcf preppedVcf = Gatk4MergeSortVcf.sortedVcf
    }
}
