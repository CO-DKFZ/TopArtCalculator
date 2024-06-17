
#' @keywords internal
#' @importFrom dplyr filter rowwise mutate select between
#' @importFrom stringr %>%
estimate_signature_score <- function(yapsaTable=NULL, 
                                     seqMethod,
                                     run_yapsa=F,
                                     FILE_REFGEN=NULL,
                                     TBL_SOM_SNV_YAPSA=NULL,
                                     TOTAL_SNVS_YAPSA=NULL){
  exposure <- sig <- lower <- upper <- conf <- NULL
  if(run_yapsa==T){
    if(nrow(TBL_SOM_SNV_YAPSA)==0){
      ## if no somatic small variants are provided... assuming no variants present
      ## signature score --> 0
      meta_info <- list(AC3_exp=0, AC3_conf=NA, 
                        signature_score=0,
                        total_snvs_yapsa=0,
                        P1_message="no input snvs and no YAPSA result provided -> assuming there are no somatic SNVs -> 0 pts") %>%
        return()
    # } else {
    #   yapsaTable <- run_YAPSA_signature_analysis(FILE_REFGEN, seqMethod, TBL_SOM_SNV_YAPSA)
    }
  }
  #rel_norm_type <- ifelse(seqMethod=='WGS', 'Valid_abs', 'Valid_norm')
  df_signatures <- yapsaTable %>%
    filter(#norm_type==rel_norm_type,
      !exposure==0,
      sig %in% paste0("AC", c(1:30))) %>% 
    rowwise() %>%
    mutate(conf=ifelse(between(0, lower, upper), "low", "high")) %>%
    select(sig, conf, exposure)
  if(is.null(TOTAL_SNVS_YAPSA)){
    TOTAL_SNVS_YAPSA <- df_signatures %>% 
      pull(exposure) %>%
      sum()
  }
  AC3_data <- df_signatures[which(df_signatures$sig=="AC3"),]
  if(nrow(AC3_data)==0){
    signature_score <- 0
    message <- "AC3 not detected"
  } else if(TOTAL_SNVS_YAPSA<25){
    signature_score <- 0
    message <- "number of SNVs not sufficient (below 25)"
  } else if(AC3_data$conf=="low"){
    signature_score <- 1
    message <- "AC3 detected -> zero in confidence interval"
  } else {
    signature_score <- 2
    message <- "AC3 detected -> zero not in confidence interval"
  }
  if(signature_score==0){
    meta_info <- list(AC3_exp=0, AC3_conf=NA, 
                      signature_score=0)
  } else {
    meta_info <- AC3_data %>% 
      select(AC3_exp=exposure, AC3_conf=conf) %>%
      as.list() %>%
      append(list(signature_score=signature_score))
  }
  return(meta_info %>% 
           append(list(P1_message=paste("mutational signature analysis based on",
                                        TOTAL_SNVS_YAPSA,"SNVs ->",
                                        message, "->", signature_score, "pts"),
                       total_snvs_yapsa=TOTAL_SNVS_YAPSA
           )))
}

#' calc rearrangement score
#' @keywords internal
#' @return A List.
#' @importFrom dplyr case_when
estimate_rearrangement_score <- function(HRD_LOH, LST){
  rearrangement_score <- case_when(
    HRD_LOH+LST>20 ~ 2,
    HRD_LOH+LST>10 ~ 1,
    TRUE ~ 0
  )
  message = paste("HRD_LOH:", HRD_LOH, "LST:", LST, "->", rearrangement_score, "pts")
  return(list(rearrangement_score=rearrangement_score, P2_message=message, HRD_LOH=HRD_LOH, LST=LST))
}

#' calc germline score
#' @keywords internal
#' @return A List.
#' @importFrom GenomicRanges elementMetadata ranges
#' @importFrom dplyr as_tibble select
#' @importFrom stringr %>%
estimate_germline_score <- function(germSmallVars, 
                                    GERM_TA_GENES){
  . <- gene <- ACMG_class <- NULL
  if(is.null(germSmallVars)){
    germline_score <- 0
    gr_germ <- NULL
    G2_message=NULL
  } else {
    if(length(GenomicRanges::ranges(germSmallVars))==0){  
      germline_score <- 0
      gr_germ <- NULL
      G2_message=NULL
    } else {
      gr_germ_raw <- germSmallVars %>%
        .[(GenomicRanges::elementMetadata(.)[,"ACMG_class"] %in% c("4","5"))] %>%
        .[(GenomicRanges::elementMetadata(.)[,"gene"] %in% GERM_TA_GENES)]
      if(length(GenomicRanges::ranges(gr_germ_raw))==0){
        germline_score <- 0
        gr_germ <- NULL
        G2_message <- NULL
      } else {
        germline_score <- 1
        G2_message <- gr_germ_raw %>% 
          as_tibble() %>%
          select(relevant_genes=gene,
                 ACMG_class)
        gr_germ <- gr_germ_raw
      }  
    }
  }  
  return(list(germline_score=germline_score,
              gr_germ=gr_germ,
              G2_message=G2_message
  ))
}

#' calc somatic score
#' @keywords internal
#' @return A List.
#' @importFrom dplyr left_join mutate group_by summarize
#' @importFrom ZygosityPredictor predict_zygosity
estimate_somatic_score <- function(somCna, 
                                   somSmallVars, 
                                   purity, 
                                   sex, 
                                   gr_germ, 
                                   DEBUG, 
                                   SOM_TA_GENES,
                                   GR_SOM_TA_GENES,
                                   bamDna,
                                   bamRna,
                                   PLOIDY,
                                   FILE_VCF,
                                   haploBlocks,
                                   showReadDetail,
                                   OUTPUT_DIR){
  # gr_germ <- g2_criterion$gr_germ
  # purity=PURITY
  # ploidy = PLOIDY
  # sex=SEX
  # somCna=somCna
  # somSmallVars=somSmallVars
  # germSmallVars=gr_germ
  # geneModel=GR_SOM_TA_GENES
  # bamDna=bamDna
  # bamRna=bamRna
  # includeHomoDel=TRUE
  # includeIncompleteDel=FALSE
  # showReadDetail=FALSE
  # byTcn=FALSE
  # vcf=raw_vcf
  # assumeSomCnaGaps=FALSE
  # colnameTcn=NULL
  # colnameCnaType=NULL
  # printLog = T
  # distCutOff <- 5000
  # assumeSomCnaGaps=FALSE
  # detailed_genewise_output <- capture.output(
  #print("initialize ZP")
  #pred_zyg <- ZygosityPredictor::predict_zygosity(
  gene <- origin <- is_som <- NULL
  pred_zyg <- predict_zygosity(
    purity=purity, 
    ploidy = PLOIDY,
    sex=sex,
    somCna=somCna, 
    somSmallVars=somSmallVars[which(somSmallVars$gene %in% GR_SOM_TA_GENES$gene)], 
    germSmallVars=gr_germ, 
    geneModel=GR_SOM_TA_GENES,
    bamDna=bamDna,
    bamRna=bamRna,
    includeHomoDel=TRUE,
    includeIncompleteDel=FALSE,
    showReadDetail=showReadDetail,
    byTcn=TRUE,
    vcf=FILE_VCF,
    assumeSomCnaGaps=TRUE,
    colnameTcn=NULL,
    colnameCnaType=NULL,
    printLog = TRUE,
    haploBlocks=haploBlocks,
    verbose=F
    ## we dont provide logdir as the TopARTCalc an ZP are not fully the same in its scoring
    #logDir=OUTPUT_DIR
  )
  #  )
  message("ZP done")
  ## define somatic score
  if(is.null(pred_zyg$eval_per_gene)){
    ## if no gene has any affections
    somatic_score=0
    corr_eval_per_gene <- NULL
  } else {
    
    var_origin <- pred_zyg$eval_per_variant %>%
      group_by(gene) %>%
      summarize(is_som=paste(origin, collapse=", ") %>% str_match("somatic") %>% as.character())
    
    corr_eval_per_gene <- pred_zyg$eval_per_gene %>%
      left_join(var_origin, by="gene") %>%
      mutate(score=case_when(
        status == "all_copies_affected" ~ 2,
        is.na(is_som) ~ 0,
        TRUE ~ 1
      )) %>% select(-is_som)
    somatic_score <- max(corr_eval_per_gene$score)
  }
  return(list(somatic_score=somatic_score, 
              #G1_tbl_eval_per_variant=pred_zyg$eval_per_variant,
              G1_tbl_eval_per_gene=corr_eval_per_gene,
              #G1_tbl_phasing_info=pred_zyg$phasing_info,
              #uncovered_variants=pred_zyg$uncovered_input,
              #ext_snp_phasing=pred_zyg$detailed_phasing_info,
              ZygosityPredictor_result=pred_zyg
              #read_detail=pred_zyg$readpair_info#,
              #detailed_genewise_output=detailed_genewise_output
  )
  )
} 
#' @keywords internal
#' @importFrom dplyr as_tibble
#' @importFrom purrr compact
create_overview_table <- function(SAMPLE_ID, seqMethod,sex, purity, total_score,
                                  g2_criterion, g1_criterion,
                                  p1_criterion, p2_criterion){
  list(
    sample_id=ifelse(!is.null(SAMPLE_ID),
                     SAMPLE_ID, "default"),
    seqMethod=seqMethod,
    sex=sex,
    purity=purity,
    top_art_score=total_score,
    germline_score=g2_criterion$germline_score,
    germline_variants=g2_criterion$G2_message$relevant_genes %>% paste(collapse=", "),
    somatic_score=g1_criterion$somatic_score,
    genes_all_cp_aff=create_gene_overview(g1_criterion$G1_tbl_eval_per_gene, 2),
    genes_wt_cp_left=create_gene_overview(g1_criterion$G1_tbl_eval_per_gene, 1)
  ) %>%
    append(p1_criterion) %>%
    append(p2_criterion) %>%
    compact() %>%
    as_tibble()
}



# G1_tbl_eval_per_variant=g1_criterion$G1_tbl_eval_per_variant
# G1_tbl_eval_per_gene=g1_criterion$G1_tbl_eval_per_gene
# G1_phasing_info=g1_criterion$G1_tbl_phasing_info
# G2_message=g2_criterion$G2_message
#' @keywords internal
#' @importFrom knitr kable
#' @importFrom dplyr as_tibble filter select tibble between rowwise arrange pull mutate_all left_join mutate_at bind_rows desc
#' @importFrom tidyr gather
#' @importFrom stringr str_match str_replace_all %>% str_split str_detect str_replace
#' @importFrom purrr compact map_chr
#' @importFrom utils head
create_output_message <- function(overview, G1_tbl_eval_per_variant, G1_tbl_eval_per_gene, 
                                  G2_message, G1_phasing_info){
  key <- score <- value <- mes <- . <- gene <- gene_info <- phasing_info <- chr <- pos <- af <- tcn <- cna_type <- aff_cp <- pre_info <- NULL
  
  header <- overview[1:5] %>% 
    gather() %>% 
    mutate_all(.funs = insert_tabs) %>%
    mutate(mes= paste(key, value, sep="  ")) %>% 
    pull(mes) %>% 
    paste(collapse="\n") 
  P1_out=paste0("signature score: ", overview$signature_score,
                "\nscoring:\n  ", overview$P1_message) %>% str_replace_all(" ->", "\n  ->")
  P2_out=paste0("rearrangement score: ", overview$rearrangement_score,
                "\nscoring:\n  ", overview$P2_message) %>% str_replace_all(" ->", "\n  ->")
  G2_out1 <- paste0("germline score: ", overview$germline_score) 
  G1_out1 <- paste0("somatic score: ", overview$somatic_score,"\nevaluation per variant:\n  ")
  if(!is.null(G2_message)){
    G2_out2 <- paste(knitr::kable(head(G2_message, n=nrow(G2_message))),collapse = "\n  ") %>% 
      paste0("\nscoring:\n  ", .)
  } else {
    G2_out2 <- NULL
  }
  if(!is.null(G1_tbl_eval_per_gene)){
    if(!is.null(G1_phasing_info)){
      G1_phasing_info <- G1_phasing_info# %>% dplyr::rename(wt_cp=left_wt_cp, no_ovlp=no_overlap, none_rw=none_raw, classes=class_comb) %>%
      #relocate(1:13, DNA_rds, RNA_rds, classes, info)
      #  if(nrow(G1_phasing_info)!=0){
      phasing_list <-lapply(unique(G1_phasing_info$gene), function(GENE){
        tbl_gene <- G1_phasing_info %>% filter(gene==GENE) %>% #select(-gene) %>%
          apply(.,1,paste, collapse="\t", simplify=F) %>% 
          unlist() %>%
          c(paste(names(G1_phasing_info) %>% 
                    .[which(.!="name")],collapse="\t"),.) %>%
          paste(collapse="\n\t") %>%
          ## change things in final message
          str_replace_all("wt-copies left", "wt-cp") %>%
          str_replace_all("muts on diff reads", "at diff")%>% 
          str_replace_all("unclear: phasing -> ", "")%>%
          str_replace_all("unclear: ", "")
        return(c(gene=GENE, phasing_info=tbl_gene))
      }) %>% 
        bind_rows()  
    } else {
      phasing_list <- tibble(gene=NA, phasing_info=NA)
    }
    sorting <- G1_tbl_eval_per_gene %>%
      arrange(desc(score)) %>%
      pull(gene) %>% c(.,G1_tbl_eval_per_variant$gene) %>% unique()
    G1_out3 <- knitr::kable(head(G1_tbl_eval_per_gene %>% mutate(gene=factor(gene, levels=sorting)) %>% 
                                   arrange(gene), n=nrow(G1_tbl_eval_per_gene))) %>% 
      as.list() %>% unlist() %>%
      tibble(gene_info=.) %>%
      mutate(gene=str_split(gene_info, "\\|") %>% map_chr(.,2) %>% str_replace_all(" ","")) %>%
      left_join(phasing_list) %>%
      select(-gene) %>%
      mutate(
        gene_info=ifelse(str_detect(gene_info,"unclear"),
                         paste0("\033[31m", gene_info, "\033[39m"),
                         gene_info),
        mes=ifelse(!is.na(phasing_info),
                   paste(gene_info, phasing_info, sep="\n\t"),
                   paste0(gene_info,""))) %>%
      pull(mes) %>%
      paste0(collapse="\n  ") %>% 
      paste0("\n\nevaluation per gene:\n  ",.)
  } else {
    G1_out3 <- NULL
    sorting <- NULL
  }
  if(!is.null(G1_tbl_eval_per_variant)){
    if(is.null(sorting)){
      sorting <- G1_tbl_eval_per_variant$gene
    }
    G1_out2 <- 
      paste(
        knitr::kable(
          head(
            G1_tbl_eval_per_variant %>% 
              mutate(
                gene=
                  factor(gene, 
                         levels=sorting
                         )
                ) %>%
              arrange(gene) %>%
              select(
                gene, 
                class,
                chr, 
                pos, 
                af=af, 
                tcn, 
                cna_type, 
                aff_cp, 
                pre_info) %>%
              mutate(
                class=
                  str_replace(
                    class,
                    "nonsynonymous", 
                    "nonsyn") %>% 
                  str_replace("frameshift", "fs") %>%
                                                  str_replace("deletion", "del") %>%
                                                  str_replace("insertion", "ins"),
                                                pre_info=str_replace(pre_info, "LOH detected", "LOH") %>% 
                                                  str_replace("somatic-variant", "som") %>%
                                                  str_replace("germline-variant", "germ") %>%
                                                  str_replace("left wt-copies", "wt-cp") %>%
                                                  str_replace("all copies affected", "all aff") %>%
                                                  str_replace("not all affected", "not all aff") %>%
                                                  str_replace("variant lost in tumor", "lost in tumor")
                                         ) %>%
                                         mutate_at(.vars=c("af", "tcn", "aff_cp"), .funs=as.numeric) %>%
                                         mutate_at(.vars=c("af", "tcn", "aff_cp"), .funs=round, digits=2), 
                                       n=nrow(G1_tbl_eval_per_variant))),collapse = "\n  ") %>%
      paste0("\nevaluation per variant:\n  ",.)
  } else {
    G1_out2 <- NULL
  }
  G1_out1 <- paste0("somatic score: ", overview$somatic_score)
  output_message_raw <- list(
    list(header, P1_out, P2_out,G2_out1) %>% compact() %>%
      unlist() %>%
      paste(collapse="\n\n"),   
    G2_out2, "\n\n",
    G1_out1,
    G1_out2,
    G1_out3
  ) %>% compact() %>% unlist() %>% paste(collapse="") %>%
    paste0(.,"\n")
  
  return(output_message_raw)
}
#' 
insert_tabs <- function(COL){
  mx <- max(nchar(COL)) 
  lapply(COL, function(inp){
    paste0(inp, paste0(rep(" ", mx-nchar(inp)), collapse = ""))
  }) %>% c(recursive=T)
}
#' @keywords internal
#' @importFrom dplyr filter pull na_if
create_gene_overview <- function(G1_tbl_eval_per_gene, stf){
  score <- gene <- NULL
  if(!is.null(G1_tbl_eval_per_gene)){
    G1_tbl_eval_per_gene %>% filter(score==stf) %>% pull(gene) %>% paste(collapse=", ") %>%
      na_if("") %>%
      return()
  } else {
    return(NA)
  }
}
#' @keywords internal
#' @importFrom GenomicRanges elementMetadata
get_metacols <- function(obj){
  GenomicRanges::elementMetadata(obj) %>%
    return()
}
#' @keywords internal
#' @importFrom stringi stri_remove_na
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges elementMetadata
#' @importFrom dplyr as_tibble filter select tibble between rowwise all_of
#' @importFrom stringr str_match str_replace_all %>%
check_input_data <- function(SAMPLE_ID=NULL,
                             seqMethod,
                             HRD_LOH,
                             LST,
                             purity,
                             sex,
                             bamDna,
                             bamRna,
                             somCna,               # meta cols required: tcn, cna_type
                             somSmallVars=NULL,    # meta cols required: gene, af, ref, alt
                             germSmallVars=NULL,   # meta cols required: gene, af, ref, alt, ACMG_class
                             geneModel, # meta cols required: gene
                             
                             yapsaTable=NULL,
                             FILE_REFGEN=NULL,
                             TOTAL_SNVS_YAPSA=NULL,
                             GR_SOM_SNV_YAPSA=NULL,
                             
                             PRINT_OUTPUT=F,
                             STORE_OUTPUT=F,
                             OUTPUT_DIR=NULL,
                             MAN_SOM_GENES=NULL,
                             MAN_GERM_GENES=NULL,
                             DEBUG=F,
                             FILTER_SOM_SMALL_VARS){
  . <- score <- gene <- width <- seqnames <- start <- ref <- alt <- NULL
  complete_input <- c("seqMethod", "HRD_LOH", "LST", "purity", "sex", "bamDna",
                      "somCna", "geneModel", "yapsaTable") %>%
    as.character() %>%
    lapply(function(x){exists(as.character(x))}) %>%
    tibble(layer=names(.), exists=.) %>%
    filter(exists==F)
  
  if(nrow(complete_input)!=0){
    message(paste("following inputs are required but missing:", 
                  paste(complete_input$layer, collapse = "; "), "\n  aborting..."))
    return(list(proceed=F))
  } else if(!file.exists(bamDna)){
    ## muss ja net immer da sein... also könnte man auch so amchen dasses nur warnt
    message("Input bamDna seems to not exists\n  aborting...")
    return(list(proceed=F))
    # } else if(!file.exists(bamRna)){  
    #  hier auch noch sagen was passiert falls nicht da 
    #} else if(!i){
  } else {
    ## from here... all basically required inputs are given... now check if they are correct
    if(!seqMethod %in% c("WES", "WGS")){
      message("input seqMethod must be a character: either \'WES\' or \'WGS\'\n  aborting...")
      return(list(proceed=F))
    }
    
    if(is.na(as.numeric(HRD_LOH))){
      message(paste("input HRD_LOH must be numeric or a character that can be converted to numeric;\n  ", HRD_LOH, 
                    "can not be converted to numeric\n  aborting..."))
      return(list(proceed=F))
    }
    if(is.na(as.numeric(LST))){
      message(paste("input LST must be numeric or a character that can be converted to numeric;\n ", LST, 
                    "can not be converted to numeric\n  aborting..."))
      return(list(proceed=F))
    }
    allowed_sex <- c("male", "m", "female", "f") %>% c(.,toupper(.))
    if(!sex %in% allowed_sex){
      message(paste("input sex must be one of", paste(allowed_sex, collapse = "\', \'") %>% paste0("\'", ., "\'"), "\n  aborting..."))
      return(list(proceed=F))
    }
    if(is.na(as.numeric(purity))){
      message(paste("input purity must be numeric or a character that can be converted to numeric;\n  ", purity, 
                    "can not be converted to numeric\n  aborting..."))
      return(list(proceed=F))
    } else if(!between(as.numeric(purity), 0, 1)){
      message(paste("input purity must be a numerci value between 0 and 1;\n  abortin..."))
      return(list(proceed=F))
    }
    
    # ### check germline small variant input
    if(!is.null(germSmallVars)){
      if(!class(germSmallVars)[1]=="GRanges"){
        message(paste("input germSmallVars must be a GRanges object; given input appears to be:",
                      paste(unlist(class(germSmallVars)), collapse = ";")))
        return(list(proceed=F))
      } else if(!(
        #("gene" %in% names(GenomicRanges::elementMetadata(germSmallVars))|
        # "GENE" %in% names(GenomicRanges::elementMetadata(germSmallVars)))&
                  "ACMG_class" %in% names(GenomicRanges::elementMetadata(germSmallVars))#&
                  # ("af"   %in% names(GenomicRanges::elementMetadata(germSmallVars))|
                  #  "AF"   %in% names(GenomicRanges::elementMetadata(germSmallVars)))&
                  # ("ref"  %in% names(GenomicRanges::elementMetadata(germSmallVars))|
                  #  "REF"  %in% names(GenomicRanges::elementMetadata(germSmallVars)))&
                  # ("alt"  %in% names(GenomicRanges::elementMetadata(germSmallVars))|
                  #  "ALT"  %in% names(GenomicRanges::elementMetadata(germSmallVars)))
        )){
        message("input germSmallVars requires the following metadata columns: \'ACMG_class\'  and \n  aborting...")
        return(list(proceed=F))
      } else {
        col_gene <- str_match(names(GenomicRanges::elementMetadata(germSmallVars)), "gene|GENE|Gene") %>%
          stringi::stri_remove_na()
        col_af <- str_match(names(GenomicRanges::elementMetadata(germSmallVars)), "af|AF") %>%
          stringi::stri_remove_na()
        col_ref <- str_match(names(GenomicRanges::elementMetadata(germSmallVars)), "ref|REF|Ref") %>%
          stringi::stri_remove_na()
        col_alt <- str_match(names(GenomicRanges::elementMetadata(germSmallVars)), "alt|ALT|Alt") %>%
          stringi::stri_remove_na()
        num_ACMG_class <- GenomicRanges::elementMetadata(germSmallVars)[,"ACMG_class"] %>%
          str_replace_all("Likely Pathogenic|likely pathogenic|LP", "4") %>%
          str_replace_all("Likely Benign|likely benign|LB", "2") %>%
          str_replace_all("Pathogenic|pathogenic|P", "5") %>%
          str_replace_all("Benign|benign|B", "1") %>%
          str_replace_all("Uncertain Significance|uncertain significance|US|VUS", "3")
        GenomicRanges::elementMetadata(germSmallVars)[,"gene"] <- GenomicRanges::elementMetadata(germSmallVars)[,col_gene]
        GenomicRanges::elementMetadata(germSmallVars)[,"maf"] <- GenomicRanges::elementMetadata(germSmallVars)[,col_af]
        GenomicRanges::elementMetadata(germSmallVars)[,"ref"] <- GenomicRanges::elementMetadata(germSmallVars)[,col_ref]
        GenomicRanges::elementMetadata(germSmallVars)[,"alt"] <- GenomicRanges::elementMetadata(germSmallVars)[,col_alt]
        GenomicRanges::elementMetadata(germSmallVars)[,"ACMG_class"] <- num_ACMG_class
      }
    }
    ################
    ## YAPSA check
    ################
    if(!is.null(yapsaTable)){
      ## if a yapsa table is given
      required_cols <- c("sig", "exposure", "lower", "upper")
      tbl_check_presence <- tibble(col=required_cols) %>%
        rowwise() %>%
        mutate(exists=ifelse(col %in% names(yapsaTable),
                             T, F)) %>%
        filter(exists==F)
      
      if(nrow(tbl_check_presence)!=0){
        message(paste("input yapsaTable must contain the following columns:", paste(required_cols, collapse = "; "),
                      "\n the following are missing:", paste(tbl_check_presence$col, collapse = ", ")))
        
        ## check if enough alternative info is given to perform signature analysis
        run_yapsa=T
      } else {
        run_yapsa=F      
        if(is.null(TOTAL_SNVS_YAPSA)){
          message("input TOTAL_SNV_YAPSA is missing... taking sum of exposure column from yapsaTable")
        }
        TBL_SOM_SNV_YAPSA <- NULL
      }
    } else {
      run_yapsa=T
    }
    if(run_yapsa==T){
      ## no yapsa table given
      if(is.null(FILE_REFGEN)){
        message("input FILE_REFGEN (indexed fasta file) is required to run mutational signature analysis")
        return(list(proceed=F)) 
      } else if(file.exists(FILE_REFGEN)){
        message(paste("input FILE_REFGEN (indexed fasta file) is required to run mutational signature analysis but seems to not exist:",
                      FILE_REFGEN))
        return(list(proceed=F)) 
      }
      if(is.null(GR_SOM_SNV_YAPSA)){
        message("input GR_SOM_SNV_YAPSA not provided")
        if(is.null(somSmallVars)){
          message("\nand no somatic small variants provided (somSmallVars) !! \n... assuming there are no somatic variants in the sample\n !! SIGNATURE SCORE WILL BE ZERO !!")
          TBL_SOM_SNV_YAPSA <- tibble()
        } else {
          message("extracting SNVs from input somSmallVars")
          TBL_SOM_SNV_YAPSA <- as_tibble(somSmallVars) %>%
            filter(width==1) %>%
            select(CHROM=seqnames, POS=start, REF=ref, ALT=alt)
        }
        
      } else {
        if(("ref"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA))|
            "REF"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA)))&
           ("alt"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA))|
            "ALT"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA)))){
          col_ref <- str_match(names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA)), "ref|REF") %>%
            stringi::stri_remove_na()
          col_alt <- str_match(names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA)), "alt|ALT") %>%
            stringi::stri_remove_na()        
          TBL_SOM_SNV_YAPSA <- as_tibble(GR_SOM_SNV_YAPSA) %>%
            filter(width==1) %>%
            select(CHROM=seqnames, POS=start, REF=all_of(col_ref), ALT=all_of(col_alt))
        } else {
          message("input GR_SOM_SNV_YAPSA requires the following meta data columns: REF/ref; ALT/alt")
          return(list(proceed=F))
        }
        
        
      }
    }
    ##### ref genome
    if(!class(geneModel)[1]=="GRanges"){
      message(paste("input geneModel must be a GRanges object; given input appears to be:", 
                    paste(unlist(class(geneModel)), collapse = ";")))
      return(list(proceed=F))
    } else if(!("gene" %in% names(GenomicRanges::elementMetadata(geneModel)))){
      message("input geneModel requires the following metadata columns: \'gene\'\n  aborting...")
      return(list(proceed=F))
    }     
    ## if the input needs to be filtered for functional variants (exonic)
    if(FILTER_SOM_SMALL_VARS==T&!is.null(somSmallVars)){
      somSmallVars <- IRanges::subsetByOverlaps(somSmallVars, 
                                                     geneModel)
    }
    ## warnings  
    if(is.null(SAMPLE_ID)&STORE_OUTPUT==T){
      cat("No sample identification; storing output with default sample ID\n  set SAMPLE_ID=\'your_sample_id\'")
    }
    if(is.null(somSmallVars)){
      message("no somatic small variants provided (somSmallVars=NULL); \n  assuming there are no relevant variants present")
    }
    if(is.null(germSmallVars)){
      message("no germline variants provided (germSmallVars=NULL); \n  assuming there are no relevant variants present")
    }
    if(is.null(MAN_SOM_GENES)){
      SOM_TA_GENES <- load_somatic_topart_genes()
    } else {
      #cat("Using manual list of genes for somatic score")
      SOM_TA_GENES <- MAN_SOM_GENES
    }
    if(is.null(MAN_GERM_GENES)){
      GERM_TA_GENES <- load_germline_topart_genes()
    } else {
      GERM_TA_GENES <- MAN_GERM_GENES
    }
    return(
      list(
        proceed=T,
        SOM_TA_GENES=        SOM_TA_GENES,
        GERM_TA_GENES=       GERM_TA_GENES,
        somCna=              somCna,
        germSmallVars=  germSmallVars,
        somSmallVars=   somSmallVars,
        run_yapsa=           run_yapsa,
        TBL_SOM_SNV_YAPSA=   TBL_SOM_SNV_YAPSA
      )
    )
  }  
}



# run_YAPSA_signature_analysis <- function(FILE_REFGEN, seqMethod, TBL_SOM_SNV_YAPSA){
#   #require(YAPSA)
#   # load data for signatures and cutoffs from the package
#   data(sigs)
#   data(cutoffs)
#   # create the SNV mutational cataloge
#   mutation_catalogue_list <- create_mutation_catalogue_from_df(
#     this_df = TBL_SOM_SNV_YAPSA,
#     this_seqnames.field = "CHROM",
#     this_refGenome = Rsamtools::FaFile(FILE_REFGEN),
#     this_wordLength = 3,
#     this_rownames = rownames(AlexCosmicValid_sig_df),
#     this_verbose = 0)
#   mutation_catalogue_df <- as_tibble(mutation_catalogue_list$matrix)
#   # perform signature anaylysis
#   if(seqMethod == "WES"){
#     sample_col <- "Valid_norm"
#     cor_list <- targetCapture_cor_factors[[targetCapture]]
#     corrected_catalogue_df <- normalizeMotifs_otherRownames(mutation_catalogue_df,
#                                                             cor_list$rel_cor)
#     CosmicValid_LCDlist <- LCD_complex_cutoff(
#       in_mutation_catalogue_df = corrected_catalogue_df,
#       in_signatures_df = AlexCosmicValid_sig_df,
#       in_cutoff_vector = cutoffCosmicValid_rel_df[6, ],
#       in_filename = NULL,
#       in_method = "abs",
#       in_sig_ind = AlexCosmicValid_sigInd_df)
#   } else {
#     sample_col <- "Valid_abs"
#     corrected_catalogue_df <- mutation_catalogue_df
#     CosmicValid_LCDlist <- LCD_complex_cutoff(
#       in_mutation_catalogue_df = corrected_catalogue_df,
#       in_signatures_df = AlexCosmicValid_sig_df,
#       in_cutoff_vector = cutoffCosmicValid_abs_df[6, ],
#       in_filename = NULL,
#       in_method = "abs",
#       in_sig_ind = AlexCosmicValid_sigInd_df)
#   }
#   ## create confidence intervals
#   tbl_mutsig_conf <- variateExp(
#     in_catalogue_df = corrected_catalogue_df,
#     in_sig_df = AlexCosmicValid_sig_df[, rownames(CosmicValid_LCDlist$exposures), drop = FALSE],
#     in_exposures_df = CosmicValid_LCDlist$exposures,
#     in_sigLevel = 0.025, in_delta = 0.4) %>%
#     mutate(sample=sample_col)
#   return(tbl_mutsig_conf)
# }

