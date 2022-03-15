version 1.0

import "./wdl_structs.wdl"
import "pon_wkf.wdl" as ponWkf
#import "widdleware/tasks/label_and_transfer.wdl" as dataTransfer

workflow PonWrapper {
    input {
        Boolean production = true
        #Array[SampleBamInfo]+ tumorInfos

        String project
        String bamPath
        String baiPath
        String sampleId
        

        # strelka2
        File strelkaJsonLog
        File configureStrelkaSomaticWorkflow
        #   mutect2
        File mutectJsonLog
        Array[String]+ listOfChroms
        IndexedReference referenceFa
        #   Manta
        IndexedTable callRegions
        File mantaJsonLog
        #   Svaba
        File dbsnpIndels
        File svabaJsonLog
        File mutectJsonLogFilter
        BwaReference svabaIndexedReference
        File refCache
        Boolean highMem = false
        
        # annotation:
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
    }


    call downloadFiles {
        input:
            project = project,
            bamPath = bamPath,
            baiPath = baiPath
    }

    SampleBamInfo tumorInfo = object {
       sampleId : sampleId,
       finalBam :  {
            "bam": downloadFiles.bam,
            "bamIndex": downloadFiles.bai
        }
    }

    call ponWkf.CallingPon as ponWkf {
        input:
            production = production,
            tumorInfos = [tumorInfo],
            # strelka2
            strelkaJsonLog = strelkaJsonLog,
            configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
            #   mutect2
            mutectJsonLog = mutectJsonLog,
            listOfChroms = listOfChroms,
            referenceFa = referenceFa,
            #   Manta
            callRegions = callRegions,
            mantaJsonLog = mantaJsonLog,
            #   Svaba
            dbsnpIndels = dbsnpIndels,
            svabaJsonLog = svabaJsonLog,
            mutectJsonLogFilter = mutectJsonLogFilter,
            svabaIndexedReference = svabaIndexedReference,
            refCache = refCache,
            highMem = highMem,
        
            # annotation:
            cosmicCoding = cosmicCoding,
            cosmicNoncoding = cosmicNoncoding,
        
            # Public
            cancerResistanceMutations = cancerResistanceMutations,
            vepCache = vepCache,
            annotations = annotations,
            plugins = plugins,
            vepGenomeBuild = vepGenomeBuild,
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
            masterMind = masterMind
    }


    output {
            # General
        Array[File] vcfAnnotatedVep = ponWkf.vcfAnnotatedVep
        
        # Mutect2
        Array[File] mutect2 = ponWkf.mutect2
        # Manta 
        Array[File] filteredMantaSV = ponWkf.filteredMantaSV
        # Strelka2

        # Svaba
        Array[File] svabaSvGz = ponWkf.svabaSvGz
        Array[File] svabaIndel = ponWkf.svabaIndel
    }
}


task downloadFiles {
    input {
        String bamPath
        String baiPath
        String project
        String baseName = basename(bamPath)
    }

    command {
        set -e
        mkdir output
        gsutil -u ~{project} -m cp ~{bamPath} output/
        gsutil -u ~{project} -m cp ~{baiPath} output/

    }

    output {
        File bam = "output/~{baseName}"
        File bai = "output/~{baseName}.crai"
    }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli"
    disks: "local-disk 100 HDD"
  }

}
