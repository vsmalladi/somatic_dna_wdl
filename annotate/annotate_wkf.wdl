version 1.0

import "../merge_vcf/merge_vcf.wdl" as merge_vcf
import "annotate.wdl" as annotate
import "../wdl_structs.wdl"

workflow Annotate {
    input {
        String tumor
        String normal
        String pairName
        File unannotatedVcf
        
        String vepVersion
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding
        IndexedVcf vepClinvarXMLTrait
        IndexedTable vepCaddSnp
        IndexedVcf vepGnomadExomes
        IndexedVcf vcfCompressed
        IndexedTable vepCaddIndel
        IndexedVcf vepGnomadGenomes
        IndexedReference vepFastaReference
        
        File cosmicCensus
        
        File cancerResistanceMutations
        String genome
        
        File ensemblEntrez
        String library
        
        String vepFathmmDockerImage
        String pysamDockerImage
        String gatkDockerImage
        String bgzipDockerImage
        String bcftoolsDockerImage
        Int threads
        Int memoryGb
        IndexedReference referenceFa
    }
        
    call merge_vcf.CompressVcf as unannotatedCompressVcf {
        input:
            vcf = unannotatedVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bgzipDockerImage
    }

    call merge_vcf.IndexVcf as unannotatedIndexVcf {
        input:
            vcfCompressed = unannotatedCompressVcf.vcfCompressed,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call annotate.Vep {
        input:
            pairName = pairName,
            vcfCompressed = unannotatedIndexVcf.vcfCompressedIndexed,
            vepVersion = vepVersion,
            cosmicCoding = cosmicCoding,
            cosmicNoncoding = cosmicNoncoding,
            vepClinvarXMLTrait = vepClinvarXMLTrait,
            vepCaddSnp = vepCaddSnp,
            vepGnomadExomes = vepGnomadExomes,
            vepCaddIndel = vepCaddIndel,
            vepGnomadGenomes = vepGnomadGenomes,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = vepFathmmDockerImage
    }
    
    call annotate.AddCosmic {
        input:
            pairName = pairName,
            cosmicCensus = cosmicCensus,
            vcfAnnotatedVep = Vep.vcfAnnotatedVep,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call annotate.AddCancerResistanceMutations {
        input:
            pairName = pairName,
            cancerResistanceMutations = cancerResistanceMutations,
            genome = genome,
            vcfAnnotatedCancerGeneCensus = AddCosmic.vcfAnnotatedCancerGeneCensus,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call annotate.AnnotateId {
        input:
            pairName = pairName,
            vcfAnnotatedResistance = AddCancerResistanceMutations.vcfAnnotatedResistance,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call annotate.RenameCsqVcf {
        input:
            pairName = pairName,
            vcfAnnotatedId = AnnotateId.vcfAnnotatedId,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call annotate.MainVcf {
        input:
            pairName = pairName,
            vcfAnnotated = RenameCsqVcf.vcfCsqRenamed,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call annotate.TableVcf {
        input:
            tumor = tumor,
            normal = normal,
            pairName = pairName,
            mainVcf = MainVcf.mainVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call annotate.VcfToMaf {
        input:
            tumor = tumor,
            normal = normal,
            pairName = pairName,
            mainVcf = MainVcf.mainVcf,
            library = library,
            vepVersion = vepVersion,
            ensemblEntrez = ensemblEntrez,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    output {
        File supplementalVcf = RenameCsqVcf.vcfCsqRenamed
        File finalVcf = MainVcf.mainVcf
        File vcfAnnotatedTxt = TableVcf.vcfAnnotatedTxt
        File maf = VcfToMaf.maf
    }
}