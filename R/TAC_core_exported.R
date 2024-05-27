#' load genes of somtic TOP-ART criterion as character vector#'
#' @return A character vector.
#' @examples
#' library(TopArtCalc)
#' 
#' somatic_topart_genes <- load_somatic_topart_genes()
#' 
#' @export
load_somatic_topart_genes <- function(){
  message("Gene names originating from hg19 gencode version 19")
  somTopArtGenes <- c(
    'C1orf86','MTOR','MAD2L2','RPA2','SFPQ','PLK3','RAD54L', 'USP1','GADD45A',
    'CDC7','PRMT6','PARP1','SMC6','GEN1','DNMT3A','FANCL','ERCC3','MCM6',
    'PSMD14','NABP1','SUMO1','INO80D','BARD1','SMARCAL1','STK36','FANCD2',
    'TOP2B', 'WDR48','BAP1','POLQ','MCM2','TOPBP1','ATR','RNF168','POLN','HELQ',
    'FAM175A','SMARCA5','MND1','PAPD7','ERCC8','CDK7','POLK','MSH3','XRCC4',
    'RAD50','HUS1B','MDC1','FANCE','POLH','MCM3','MMS22L','REV3L','HDAC2',
    'MCM9','SHPRH','AP5Z1','RPA3','POLM','HUS1','RFC2','NAMPT','CDK5','XRCC2',
    'ESCO2','WRN','GINS4','POLB','SPIDR','PRKDC','NBN','RAD54B','RRM2B','RAD21',
    'NSMCE2','TONSL','RECQL4','SMARCA2','FANCG','SMC5','FANCC','CDC14B',
    'RAD23B','INIP','SWI5','ABL1','IPMK','PTEN','HELLS','SFR1','SMC3','RRM1',
    'FANCF','SSRP1','DDB1','FEN1','KAT5','MUS81','POLD3','MRE11A','ATM','H2AFX',
    'CHEK1','RAD52','RAD51AP1','NABP2','TIMELESS','NAP1L1','UBE2N','UNG',
    'GTF2H3','BRCA2','LIG4','APEX1','REC8','G2E3','FANCM','MNAT1','ZFYVE26',
    'RAD51B','MLH3','YY1','XRCC3','FAN1','RAD51','INO80','TP53BP1','DUT',
    'MORF4L1','FANCI','BLM','EME2','MEIOB','DNASE1L2','SLX4','USP7','ERCC4',
    'SMG1','PALB2','NSMCE1','SLX1B','SLX1A','USP10','GINS2','FANCA','RPA1',
    'ZSWIM7','TOP3A','ATAD5','LIG3','RAD51D','CDK12','PSMC3IP','BRCA1','EME1',
    'RAD51C','BRIP1','ESCO1','RBBP8','MUM1','XAB2','SWSAP1','C19orf40','ERCC1',
    'LIG1','CCDC155','PNKP','RAD21L1','AP5S1','MCM8','SPO11','CDC45','CHEK2',
    'POLR2F','DMC1','XRCC6','MAPK12','FANCB','UBA1','TEX11','UBE2A','BRCC3') 
  return(somTopArtGenes)
}
#' load genes of germline TOP-ART criterion as character vector
#' @export
#' @examples
#' library(TopArtCalc)
#' 
#' germline_topart_genes <- load_germline_topart_genes()
#' @return A character vector.
load_germline_topart_genes <- function(){
  message("Gene names originating from hg19 gencode version 19")
  germTopArtGenes <- c(
    'FANCL','ERCC3','BARD1','FANCD2','BAP1','ATR','FAM175A','RAD50','FANCE',
    'XRCC2','WRN','NBN','RECQL4','FANCG','FANCC','PTEN','FANCF','MRE11A','ATM',
    'BRCA2','MLH3','RAD51','FANCI','BLM','SLX4','ERCC4','PALB2','FANCA',
    'RAD51D','BRCA1','RAD51C','BRIP1','CHEK2','FANCB')
  return(germTopArtGenes)
}


#' calcuation of TOP-ART score
#' @param purity purity of the sample (numeric value between 0 and 1 indicating 
#' the fraction of relevant sample with control/unrelevant tissue)
#' @param ploidy ploidy of the sample (numeric value)
#' @param sex sex of the sample (character: "male", "female", "m", "f")
#' @param somCna GRanges object containing all genomic regions with annotated 
#' total copynumber and cna_type as metadata columns. The total-copynumber 
#' column should be named "tcn" but also some other commonly used names. 
#' It should contain numeric values or characters that can be converted to 
#' numeric values. The cna_type column must contain the information about 
#' loss of heterozygosity (LOH). Therefore the term "LOH" must be explicitely 
#' mentioned in the column. If a genomic region is not present in the object, 
#' it will be taken as heterozygous with neutral TCN of 2. 
#' @param somSmallVars GRanges object containing all somatic small 
#' variants (SNV and INDEL).
#' Required metadata columns are reference base (ref/REF), 
#' alternative base (alt/ALT),
#' annotation of the gene name (gene/GENE) and the allele-frequency (af/AF). 
#' If the object is not provided the tool assumes there are no somatic small 
#' variants.
#' @param germSmallVars GRanges object containing all germline small 
#' variants (SNV and INDEL).
#' Required metadata columns are reference base (ref/REF), alternative 
#' base (alt/ALT),
#' annotation of the gene name (gene/GENE) and the allele-frequency (af/AF)
#' If the object is not provided the tool assumes there are no germline small 
#' variants.
#' @param geneModel GRanges object containing the gene-annoattion of 
#' the used reference genome with metadata column of the gene name (gene)
#' @param seqMethod charcter: either "WGS" or "WES"
#' @param HRD number of LOH events (wait for final text in manuscript)
#' @param LST number of large scale state transition events
#' @param bamDna path to bam-file
#' @param bamRna optional; path to rna file (bam format)
#' @param vcf character; path to variant call file (.vcf.gz format). 
#' Will be used (if provided)
#' for extended SNP phasing if variants on the same gene are too far away from
#' each other for direct haplotype phasing
#' @param haploBlocks GRanges object containing haploblocks. Haploblocks are
#' defined as genomic regions in which SNPs are phased to a specific allele.
#' For example a haploblock could be chr1:1000-10000. This would mean that every
#' genotype annotation in the format "1|0" or "0|1" of a SNP in this region will 
#' be used to phase somatic variants and define their genotype
#' @param SAMPLE_ID optional ID of sample can be added
#' @param yapsaTable table conatining the mutational signatures
#' @param FILE_REFGEN currently not usable.. required if YAPSA is incorporated
#' @param TOTAL_SNVS_YAPSA number of SNVs used to create yapsa signatures
#' @param GR_SOM_SNV_YAPSA currently not usable.. required if YAPSA is incorporated
#' @param GR_FULL_VARS currently not usable.. required if YAPSA is incorporated
#' @param PRINT_OUTPUT should output be printed into the concole?
#' @param OUTPUT_DIR output directory
#' @param MAN_SOM_GENES if somatic genes of relevance deviate from original TOP-ART gene set
#' @param MAN_GERM_GENES if germline genes of relevance deviate from original TOP-ART gene set 
#' @param DEBUG debugging output
#' @param FILTER_SOM_SMALL_VARS should somatic variants be filtered for functional regions? - as defined for TOP-ART
#' @param showReadDetail exports an seperate table containg read info for reads used for phasing (see ZygosityPredictor)
#' @examples
#'  
#'   
#' ## we use the example dataset from the Bioconductor package ZygosiytPredictor
#' library(ZygosityPredictor)
#' 
#' bamfile <- system.file("extdata", "ZP_example.bam", 
#'                        package = "ZygosityPredictor")
#'                        
#' ## load ZygosityPredictor variant datasets
#' data("GR_GERM_SMALL_VARS")
#' data("GR_SCNA")
#' data("GR_SOM_SMALL_VARS")
#' data("GR_GENE_MODEL")
#' 
#' ## create mutational signatures table as it is produced by bioconductor package
#' ## YAPSA
#' yapsaTable <- data.frame(
#'   sig=c("AC1","AC3","AC8"),
#'   exposure=c(4566, 2123, 598),
#'   lower=c(3899, 1902, 430),
#'   upper=c(5698, 2455, 789)
#' )
#' 
#' library(TopArtCalc)
#' ## calculate TOP-ART score
#' topart_run <- calculate_TOP_ART_score(
#'   seqMethod="WGS",
#'   HRD=10,
#'   LST=4,
#'   purity=0.9,
#'   ploidy=2,
#'   sex="female",
#'   ## files
#'   bamDna=bamfile,
#'   ## granges
#'   somCna=GR_SCNA,   # meta cols required: tcn, cna_type
#'   somSmallVars=GR_SOM_SMALL_VARS,   # meta cols required: gene, af, ref, alt
#'   germSmallVars=GR_GERM_SMALL_VARS,  # meta cols required: gene, af, ref, alt, ACMG_class
#'   geneModel=GR_GENE_MODEL,     # meta cols required: gene
#'   yapsaTable=yapsaTable
#'   )
#' @export
calculate_TOP_ART_score <- function(SAMPLE_ID=NULL,
                                    seqMethod,
                                    HRD,
                                    LST,
                                    purity,
                                    ploidy,
                                    sex,
                                    
                                    ## files
                                    bamDna,
                                    bamRna=NULL,
                                    vcf=NULL,
                                    
                                    ## granges
                                    somCna,                   # meta cols required: tcn, cna_type
                                    somSmallVars=NULL,   # meta cols required: gene, af, ref, alt
                                    germSmallVars=NULL,  # meta cols required: gene, af, ref, alt, ACMG_class
                                    geneModel,     # meta cols required: gene
                                    haploBlocks=NULL,
                                    
                                    ## yapsa options
                                    yapsaTable=NULL,
                                    FILE_REFGEN=NULL,
                                    TOTAL_SNVS_YAPSA=NULL,
                                    GR_SOM_SNV_YAPSA=NULL,
                                    GR_FULL_VARS=NULL,
                                    ## options
                                    PRINT_OUTPUT=F,
                                    OUTPUT_DIR=NULL,
                                    MAN_SOM_GENES=NULL,
                                    MAN_GERM_GENES=NULL,
                                    DEBUG=F,
                                    FILTER_SOM_SMALL_VARS=T,
                                    showReadDetail=F
){
  ## for manual runs
  # DEBUG=F
  # FILTER_SOM_SMALL_VARS=F
  # STORE_OUTPUT=F
  # MAN_SOM_GENES=NULL
  # MAN_GERM_GENES=NULL
  
  adapted_input <- check_input_data(SAMPLE_ID,
                                    seqMethod,
                                    HRD,
                                    LST,
                                    purity,
                                    sex,
                                    bamDna,
                                    bamRna,
                                    somCna,              
                                    somSmallVars,   
                                    germSmallVars,
                                    geneModel, 
                                    yapsaTable,
                                    FILE_REFGEN,
                                    TOTAL_SNVS_YAPSA,
                                    GR_SOM_SNV_YAPSA,
                                    PRINT_OUTPUT,
                                    STORE_OUTPUT=FALSE,
                                    OUTPUT_DIR,
                                    MAN_SOM_GENES,
                                    MAN_GERM_GENES,
                                    DEBUG,
                                    FILTER_SOM_SMALL_VARS)
  message("input checked")
  if(adapted_input$proceed==T){
    SOM_TA_GENES=adapted_input$SOM_TA_GENES
    GERM_TA_GENES=adapted_input$GERM_TA_GENES
    ## dont remove thios otherwise it will not work for germine score
    somCna=adapted_input$somCna
    germSmallVars=adapted_input$germSmallVars
    somSmallVars=adapted_input$somSmallVars
    GR_SOM_TA_GENES <- 
      geneModel[get_metacols(geneModel)[,"gene"] %in% SOM_TA_GENES]
    TBL_SOM_SNV_YAPSA <- adapted_input$TBL_SOM_SNV_YAPSA
    ## score calculation
    message("everything prepared for subscore calculation")
    g2_criterion <- 
      estimate_germline_score(
        germSmallVars, 
        GERM_TA_GENES
        )
    message("germline score done")
    g1_criterion <-
      estimate_somatic_score(
        somCna, 
        somSmallVars, 
        purity, 
        sex,
        g2_criterion$gr_germ, 
        DEBUG,
        SOM_TA_GENES,
        GR_SOM_TA_GENES,
        bamDna,
        bamRna,
        ploidy,
        vcf,
        haploBlocks,
        showReadDetail,
        OUTPUT_DIR
        ) 
    message("somatic score done")
    p1_criterion <- 
      estimate_signature_score(
        yapsaTable, 
        seqMethod,
        adapted_input$run_yapsa,
        FILE_REFGEN,
        TBL_SOM_SNV_YAPSA,
        TOTAL_SNVS_YAPSA
        )
    p2_criterion <- 
      estimate_rearrangement_score(
        as.numeric(HRD), 
        as.numeric(LST)
        )
    total_score <- 
      p1_criterion$signature_score+
      p2_criterion$rearrangement_score+
      g1_criterion$somatic_score+
      g2_criterion$germline_score
    ## output creation
    overview <- 
      create_overview_table(
        SAMPLE_ID, 
        seqMethod,
        sex, 
        purity, 
        total_score,
        g2_criterion, 
        g1_criterion,
        p1_criterion, 
        p2_criterion
        )
    output_message <- 
      create_output_message(
        overview, 
        g1_criterion$G1_tbl_eval_per_variant, 
        g1_criterion$G1_tbl_eval_per_gene, 
        g2_criterion$G2_message,
        g1_criterion$G1_tbl_phasing_info 
    )
    ### data export
    if(PRINT_OUTPUT==T){
      cat(output_message)
    }
    return(
      list(
        total_score=total_score, 
        output_message=output_message,
        overview=overview,
        G1_eval_per_variant=g1_criterion$G1_tbl_eval_per_variant, 
        G1_eval_per_gene=g1_criterion$G1_tbl_eval_per_gene, 
        G2_message=g2_criterion$G2_message,
        G1_phasing_info=g1_criterion$G1_tbl_phasing_info,
        uncovered_input=g1_criterion$uncovered_variants,
        ext_snp_phasing=g1_criterion$ext_snp_phasing,
        read_detail=g1_criterion$read_detail
        )
    )
  } else {
    message("cannot calculate TOP-ART score")
    return(NULL)
  }
}