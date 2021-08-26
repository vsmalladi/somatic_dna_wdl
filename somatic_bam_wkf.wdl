version 1.0

import "./wdl_structs.wdl"
import "calling/calling_wkf.wdl" as calling
import "calling/calling.wdl" as callingTasks
import "merge_vcf/merge_vcf_wkf.wdl" as mergeVcf
import "alignment_analysis/kourami_wfk.wdl" as kourami
import "alignment_analysis/msi_wkf.wdl" as msi
import "pre_process/conpair_wkf.wdl" as conpair
import "pre_process/qc.wdl" as qc
import "annotate/annotate_wkf.wdl" as annotate
import "annotate/annotate_cnv_sv_wkf.wdl" as annotate_cnv_sv
import "germline/germline_wkf.wdl" as germline
import "annotate/germline_annotate_wkf.wdl" as germlineAnnotate
import "baf/baf_wkf.wdl" as baf

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Minita Shah (mshah@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================


# for wdl version 1.0

task GetIndex {
    input {
        String sampleId
        Array[String] sampleIds
    }

    command {
        python /get_index.py \
        --sample-id ~{sampleId} \
        --sample-ids ~{sep=' ' sampleIds}
    }

    output {
        Int index = read_int(stdout())
    }

    runtime {
        docker: "gcr.io/nygc-internal-tools/workflow_utils:2.0"
    }
}


workflow SomaticBamWorkflow {
    input {
        BwaReference bwaReference
        IndexedReference referenceFa
        
        Boolean production = true

        Array[pairInfo]+ pairInfos
        Array[SampleBamInfo]+ normalSampleBamInfos

        # For Tumor-Normal QC
        File markerBedFile
        File markerTxtFile

        # calling
        Array[String]+ listOfChromsFull
        Array[String]+ listOfChroms
        IndexedTable callRegions
        File dbsnpIndels
        Map[String, File] chromBedsWgs
        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
        File svabaJsonLog
        File mutectJsonLog
        File mutectJsonLogFilter
        File configureStrelkaSomaticWorkflow
        
        #   BicSeq2
        Int readLength
        Int coordReadLength
        Map[Int, Map[String, File]] uniqCoords
        File bicseq2ConfigFile
        File bicseq2SegConfigFile
        Map[String, File] chromFastas
        Int tumorMedianInsertSize = 400
        Int normalMedianInsertSize = 400
        Int lambda = 4
        
        # Gridss
        String bsGenome
        File ponTarGz
        Array[File] gridssAdditionalReference
        
        # merge callers
        File intervalListBed

        String library
        File ponWGSFile
        File ponExomeFile
        IndexedVcf gnomadBiallelic

        IndexedVcf germFile

        # kourami
        BwaReference kouramiReference
        File kouramiFastaGem1Index

        # mantis
        File mantisBed
        File intervalListBed
        IndexedReference referenceFa
        
        # annotation:
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
        
        # annotate cnv
        File cytoBand
        File dgv
        File thousandG
        File cosmicUniqueBed
        File cancerCensusBed
        File ensemblUniqueBed
        
        # annotate sv
        String vepGenomeBuild
        # gap,DGV,1000G,PON,COSMIC
        File gap
        File dgvBedpe
        File thousandGVcf
        File svPon
        File cosmicBedPe
        
        # post annotation
        File cosmicCensus
        
        File ensemblEntrez
        
        # germline
        
        File excludeIntervalList
        Array[File] scatterIntervalsHcs
        
        IndexedVcf MillsAnd1000G
        IndexedVcf omni
        IndexedVcf hapmap
        IndexedVcf onekG
        IndexedVcf dbsnp
        
        IndexedVcf whitelist
        IndexedVcf nygcAf
        IndexedVcf pgx
        IndexedTable rwgsPgxBed
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf chdWhitelistVcf
    }
    
    scatter (normalSampleBamInfo in normalSampleBamInfos) {
        String normalSampleIds = normalSampleBamInfo.sampleId

        call kourami.Kourami {
            input:
                sampleId = normalSampleBamInfo.sampleId,
                kouramiReference = kouramiReference,
                finalBam = normalSampleBamInfo.finalBam,
                kouramiFastaGem1Index = kouramiFastaGem1Index
        }
        
        call germline.Germline {
            input:
                finalBam = normalSampleBamInfo.finalBam,
                normal = normalSampleBamInfo.sampleId,
                referenceFa = referenceFa,
                listOfChroms = listOfChroms,
                MillsAnd1000G = MillsAnd1000G,
                hapmap = hapmap,
                omni = omni,
                onekG = onekG,
                dbsnp = dbsnp,
                nygcAf = nygcAf,
                excludeIntervalList = excludeIntervalList,
                scatterIntervalsHcs = scatterIntervalsHcs,
                pgx = pgx,
                rwgsPgxBed = rwgsPgxBed,
                whitelist = whitelist,
                chdWhitelistVcf = chdWhitelistVcf,
                deepIntronicsVcf = deepIntronicsVcf,
                clinvarIntronicsVcf = clinvarIntronicsVcf
        }
        
        call germlineAnnotate.GermlineAnnotate as filteredGermlineAnnotate {
            input:
                unannotatedVcf = Germline.haplotypecallerFinalFiltered,
                production = production,
                referenceFa = referenceFa,
                normal = normalSampleBamInfo.sampleId,
                listOfChroms = listOfChroms,
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
        
        call germlineAnnotate.GermlineAnnotate as unFilteredGermlineAnnotate {
            input:
                unannotatedVcf = Germline.haplotypecallerVcf,
                haplotypecallerAnnotatedVcfPath = "~{normalSampleBamInfo.sampleId}.haplotypecaller.gatk.v4.1.8.0.annotated.unfiltered.vcf",
                production = production,
                referenceFa = referenceFa,
                normal = normalSampleBamInfo.sampleId,
                listOfChroms = listOfChroms,
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
    
    scatter(pairInfo in pairInfos) {
        call GetIndex as germlineGetIndex {
            input:
                sampleIds = normalSampleIds,
                sampleId = pairInfo.normal
        }
        
        # tumor insert size
        Int tumorDiskSize = ceil(size(pairInfo.tumorFinalBam.bam, "GB")) + 30
                      
        call qc.MultipleMetrics as tumorMultipleMetrics {
            input:
                referenceFa = referenceFa,
                finalBam = pairInfo.tumorFinalBam,
                sampleId = pairInfo.tumor,
                diskSize = tumorDiskSize
        }
        
        call callingTasks.GetInsertSize as tumorGetInsertSize {
            input:
                insertSizeMetrics = tumorMultipleMetrics.insertSizeMetrics
        }
        
        # normal insert size
        Int normalDiskSize = ceil(size(pairInfo.normalFinalBam.bam, "GB")) + 30
                      
        call qc.MultipleMetrics as normalMultipleMetrics {
            input:
                referenceFa = referenceFa,
                finalBam = pairInfo.normalFinalBam,
                sampleId = pairInfo.normal,
                diskSize = normalDiskSize
        }
        
        call callingTasks.GetInsertSize as normalGetInsertSize {
            input:
                insertSizeMetrics = normalMultipleMetrics.insertSizeMetrics
        }
        
        call baf.Baf {
            input:
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                sampleId = pairInfo.normal,
                tumorFinalBam = pairInfo.tumorFinalBam,
                normalFinalBam = pairInfo.normalFinalBam,
                germlineVcf = unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf[germlineGetIndex.index]
        }
        
        call conpair.Conpair {
            input:
                finalTumorBam = pairInfo.tumorFinalBam,
                finalNormalBam = pairInfo.normalFinalBam,
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                pairName = pairInfo.pairId,
                referenceFa = referenceFa,
                markerBedFile = markerBedFile,
                markerTxtFile = markerTxtFile
        }

        call calling.Calling {
            input:
                mantaJsonLog = mantaJsonLog,
                lancetJsonLog = lancetJsonLog,
                mutectJsonLog = mutectJsonLog,
                mutectJsonLogFilter = mutectJsonLogFilter,
                svabaJsonLog = svabaJsonLog,
                strelkaJsonLog = strelkaJsonLog,
                configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                pairInfo = pairInfo,
                listOfChroms = listOfChroms,
                listOfChromsFull = listOfChromsFull,
                referenceFa = referenceFa,
                callRegions = callRegions,
                bwaReference = bwaReference,
                dbsnpIndels = dbsnpIndels,
                chromBedsWgs = chromBedsWgs,
                readLength = readLength,
                coordReadLength = coordReadLength,
                uniqCoords = uniqCoords,
                bicseq2ConfigFile = bicseq2ConfigFile,
                bicseq2SegConfigFile = bicseq2SegConfigFile,
                tumorMedianInsertSize = tumorGetInsertSize.insertSize,
                normalMedianInsertSize = normalGetInsertSize.insertSize,
                chromFastas = chromFastas,
                bsGenome = bsGenome,
                ponTarGz = ponTarGz,
                gridssAdditionalReference = gridssAdditionalReference
        }

        call msi.Msi {
            input:
                normal = pairInfo.normal,
                pairName = pairInfo.pairId,
                mantisBed = mantisBed,
                intervalListBed = intervalListBed,
                referenceFa = referenceFa,
                tumorFinalBam = pairInfo.tumorFinalBam,
                normalFinalBam = pairInfo.normalFinalBam
        }

        PairRawVcfInfo pairRawVcfInfo = object {
            pairId : pairInfo.pairId,
            filteredMantaSV : Calling.filteredMantaSV,
            strelka2Snv : Calling.strelka2Snv,
            strelka2Indel : Calling.strelka2Indel,
            mutect2 : Calling.mutect2,
            lancet : Calling.lancet,
            svabaSv : Calling.svabaSv,
            svabaIndel : Calling.svabaIndel,
            gridssVcf : Calling.gridssVcf,
            bicseq2Png : Calling.bicseq2Png,
            bicseq2 : Calling.bicseq2,
            tumor : pairInfo.tumor,
            normal : pairInfo.normal,
            tumorFinalBam : pairInfo.tumorFinalBam,
            normalFinalBam : pairInfo.normalFinalBam

        }

        if (library == 'WGS') {
            call mergeVcf.MergeVcf as wgsMergeVcf {
                input:
                    pairRawVcfInfo = pairRawVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponWGSFile,
                    germFile = germFile

            }
        }

        if (library == 'Exome') {
            call mergeVcf.MergeVcf as exomeMergeVcf {
                input:
                    pairRawVcfInfo = pairRawVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponExomeFile,
                    germFile = germFile

            }
        }

        File mergedVcf = select_first([wgsMergeVcf.mergedVcf, exomeMergeVcf.mergedVcf])
        
        call annotate.Annotate {
            input:
                unannotatedVcf = mergedVcf,
                referenceFa = referenceFa,
                production = production,
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                pairName = pairInfo.pairId,
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
        
        call annotate_cnv_sv.AnnotateCnvSv {
            input:
                tumor = pairRawVcfInfo.tumor,
                normal = pairRawVcfInfo.normal,
                pairName = pairRawVcfInfo.pairId,
                listOfChroms =listOfChroms,
                bicseq2 = pairRawVcfInfo.bicseq2,
                cytoBand = cytoBand,
                dgv = dgv,
                thousandG = thousandG,
                cosmicUniqueBed = cosmicUniqueBed,
                cancerCensusBed = cancerCensusBed, 
                ensemblUniqueBed = ensemblUniqueBed,
                
                filteredMantaSV = pairRawVcfInfo.filteredMantaSV,
                svabaSv = pairRawVcfInfo.svabaSv,
                gridssVcf = pairRawVcfInfo.gridssVcf,
                vepGenomeBuild = vepGenomeBuild,
                gap = gap,
                dgvBedpe = dgvBedpe,
                thousandGVcf = thousandGVcf,
                svPon = svPon,
                cosmicBedPe = cosmicBedPe
            
    }
      
   }

    output {
        # Germline
        Array[IndexedVcf] haplotypecallerFinalFiltered = Germline.haplotypecallerFinalFiltered
        Array[File] filteredHaplotypecallerAnnotatedVcf = filteredGermlineAnnotate.haplotypecallerAnnotatedVcf
        Array[File] haplotypecallerAnnotatedVcf = unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf
        Array[File?] alleleCountsTxt = Baf.alleleCountsTxt
        Array[File] kouramiResult = Kourami.result
        # CNV SV output
        Array[File] cnvAnnotatedFinalBed  = AnnotateCnvSv.cnvAnnotatedFinalBed
        Array[File] cnvAnnotatedSupplementalBed  = AnnotateCnvSv.cnvAnnotatedSupplementalBed
        Array[File] svFinalBedPe = AnnotateCnvSv.svFinalBedPe
        Array[File] svHighConfidenceFinalBedPe = AnnotateCnvSv.svHighConfidenceFinalBedPe
        Array[File] svSupplementalBedPe = AnnotateCnvSv.svSupplementalBedPe
        Array[File] svHighConfidenceSupplementalBedPe = AnnotateCnvSv.svHighConfidenceSupplementalBedPe
        # SNV INDELs
        Array[PairVcfInfo] pairVcfInfos = Annotate.pairVcfInfo
        Array[File] mergedVcfs = mergedVcf
        Array[PairRawVcfInfo] pairRawVcfInfos = pairRawVcfInfo
        Array[File] mantisWxsKmerCountsFinal = Msi.mantisWxsKmerCountsFinal
        Array[File] mantisWxsKmerCountsFiltered = Msi.mantisWxsKmerCountsFiltered
        Array[File] mantisExomeTxt = Msi.mantisExomeTxt
        Array[File] mantisStatusFinal = Msi.mantisStatusFinal

        # Conpair
        Array[File] concordanceAll = Conpair.concordanceAll
        Array[File] concordanceHomoz = Conpair.concordanceHomoz
        Array[File] contamination = Conpair.contamination

    }
}

