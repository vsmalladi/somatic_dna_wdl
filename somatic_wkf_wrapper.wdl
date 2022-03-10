version 1.0

import "./wdl_structs.wdl"
import "somatic_wkf.wdl" as somaticWkf
import "widdleware/tasks/label_and_transfer.wdl" as dataTransfer
 

workflow SomaticDNA {
    input {
        Boolean external = false

        BwaReference bwaReference
        IndexedReference referenceFa
        File adaptersFa
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        File bqsrCallRegions
        File chromLengths
        File hsMetricsIntervals
        File randomIntervals
        Array[sampleInfo]+ normalSampleInfos
        Array[sampleInfo]+ sampleInfos
        Array[PairRelationship]+ listOfPairRelationships

        Boolean trim = true
        Boolean production = true
        Boolean bypassQcCheck = false

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

        IndexedVcf omni
        IndexedVcf hapmap
        IndexedVcf onekG

        IndexedVcf whitelist
        IndexedVcf nygcAf
        IndexedVcf pgx
        IndexedTable rwgsPgxBed
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf chdWhitelistVcf

        # signatures
        File cosmicSigs

        Boolean highMem = false

    }

    call somaticWkf.SomaticDNA as somaticWkf {
        input:
            external = external,

            bwaReference = bwaReference,
            referenceFa = referenceFa,
            adaptersFa = adaptersFa,
            MillsAnd1000G = MillsAnd1000G,
            Indels = Indels,
            dbsnp = dbsnp,
            bqsrCallRegions = bqsrCallRegions,
            chromLengths = chromLengths,
            hsMetricsIntervals = hsMetricsIntervals,
            randomIntervals = randomIntervals,
            normalSampleInfos = normalSampleInfos,
            sampleInfos = sampleInfos,
            listOfPairRelationships = listOfPairRelationships,
            trim = trim,
            production = production,
            bypassQcCheck = bypassQcCheck,

            # For Tumor-Normal QC
            markerBedFile = markerBedFile,
            markerTxtFile = markerTxtFile,

            # calling
            listOfChromsFull = listOfChromsFull,
            listOfChroms = listOfChroms,
            callRegions = callRegions,
            dbsnpIndels = dbsnpIndels,
            chromBedsWgs = chromBedsWgs,
            lancetJsonLog = lancetJsonLog,
            mantaJsonLog = mantaJsonLog,
            strelkaJsonLog = strelkaJsonLog,
            svabaJsonLog = svabaJsonLog,
            mutectJsonLog = mutectJsonLog,
            mutectJsonLogFilter = mutectJsonLogFilter,
            configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,

            # BicSeq2
            readLength = readLength,
            coordReadLength = coordReadLength, 
            uniqCoords = uniqCoords,
            bicseq2ConfigFile = bicseq2ConfigFile,
            bicseq2SegConfigFile = bicseq2SegConfigFile,
            chromFastas = chromFastas,
            lambda = lambda,

            # Gridss
            bsGenome = bsGenome,
            ponTarGz = ponTarGz,
            gridssAdditionalReference = gridssAdditionalReference, 

            # merge callers
            intervalListBed = intervalListBed,

            library = library,
            ponWGSFile = ponWGSFile,
            ponExomeFile = ponExomeFile,
            gnomadBiallelic = gnomadBiallelic,

            germFile = germFile,

            # kourami
            kouramiReference = kouramiReference,
            kouramiFastaGem1Index = kouramiFastaGem1Index,

            # mantis
            mantisBed = mantisBed,
            intervalListBed = intervalListBed,

            # annotation:
            vepGenomeBuild = vepGenomeBuild,
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
            masterMind = masterMind,

            # annotate cnv
            cytoBand = cytoBand,
            dgv = dgv, 
            thousandG = thousandG,
            cosmicUniqueBed = cosmicUniqueBed,
            cancerCensusBed = cancerCensusBed,
            ensemblUniqueBed = ensemblUniqueBed,

            # annotate sv
            vepGenomeBuild = vepGenomeBuild,
            # gap,DGV,1000G,PON,COSMIC
            gap = gap,
            dgvBedpe = dgvBedpe,
            thousandGVcf = thousandGVcf,
            svPon = svPon,
            cosmicBedPe = cosmicBedPe,


            # post annotation
            cosmicCensus = cosmicCensus,

            ensemblEntrez = ensemblEntrez,

            # germline
            excludeIntervalList = excludeIntervalList,
            scatterIntervalsHcs = scatterIntervalsHcs,

            omni = omni,
            hapmap = hapmap,
            onekG = onekG,

            whitelist = whitelist,
            nygcAf = nygcAf,
            pgx = pgx,
            rwgsPgxBed = rwgsPgxBed,
            deepIntronicsVcf = deepIntronicsVcf,
            clinvarIntronicsVcf = clinvarIntronicsVcf,
            chdWhitelistVcf = chdWhitelistVcf,

            # signatures
            cosmicSigs = cosmicSigs,

            highMem = false

    }


    File tmpFile = write_json(somaticWkf.finalOutput)
    call dataTransfer.LabelAndTransfer {
        input:
            workflowOutput = read_string(tmpFile),
            bypassQcCheck = bypassQcCheck,
            allQcPass = somaticWkf.allQcPass
    }

    output {
       FinalWorkflowOutput finalOutput = somaticWkf.finalOutput
    }

}


