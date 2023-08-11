### Functions DEP
#test_diff_to_all_other <- function(se, control){
#  colnames(colData(se)) <- c("type","replicate", "label", "ID" ) # maybehere
#  colData(se)$condition <- ifelse(colData(se)$type == control, control, paste0("not_",control))
#  test_diff(se, type = "control", control = control)
#}

add_rejections_SH <- function(...){
  dep <- add_rejections(...)
  rowData(dep)$name <- rowData(dep)$HGNC_Symbol
  return(dep)
}

#plot_volcano_SH <- function(..., proteins_of_interest = FALSE, additional_title = FALSE){
#  p <- plot_volcano(...)
#  p <- p + 
#    geom_point( data = filter(p$data, name %in% proteins_of_interest ),
#                aes(x =x, y = y), color = "red")
#  if(additional_title != FALSE){
#    p <- p + ggtitle(paste(p$labels$title,additional_title ) )
#  }
#  return(p)
#}

#GGPlotly_Wrapper_Volcano <- function(pp){
#  p <- plotly::plot_ly(data = pp$data, x= ~x, y= ~y, color = ~significant, 
#                                     text = ~name, colors = c("lightgrey", "black")) %>%
#    plotly::layout(title = pp$labels$title,  xaxis = list(title = "log2 Fold change"), 
#           yaxis = list(title = "-log10 p-value") )
#                      
#  return(p)
#}

GGPlotly_Volcano <- function(dep, contrast = "guess", proteins_of_interest = FALSE, 
                             additional_title = FALSE, label_size = 1.8, add_names = TRUE){
  if(contrast %in% colnames(rowData(dep))){
    d <- rowData(dep) %>% as_tibble() 
    title <-  d %>% select(all_of(vars(paste0(contrast, "_diff") ) ) ) %>% colnames() 
    title <- str_remove(title, "_diff")
    d <- d %>% 
      select(Annotated_Sequence, HGNC_Symbol, `log2 Fold change` = vars(paste0(contrast, "_diff") ), p = contains("p.val"), significant) %>%
      mutate( `-log10 p-value` = ( - log10(p) ), peptide = paste(HGNC_Symbol,"\n",  Annotated_Sequence) )
  }else{
  d <- rowData(dep) %>% as_tibble() 
  title <-  d %>% select(all_of( contains("diff") ) ) %>% colnames() 
  title <- str_remove(title, "_diff")
  d <- d %>% select(Annotated_Sequence, HGNC_Symbol, `log2 Fold change` = contains("diff"), p = contains("p.val"), significant) %>%
    mutate( `-log10 p-value` = ( - log10(p) ) )
  }
  
  p <- 
    ggplot( d , aes(`log2 Fold change`,`-log10 p-value`, color = significant, label = HGNC_Symbol, group = Annotated_Sequence ) ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.2) +
    geom_point(size = 1) +
    scale_color_manual(values = c("lightgrey", "black")) + 
    geom_point( data = filter(d, HGNC_Symbol %in% proteins_of_interest ),
                size = 1,
                aes(x =`log2 Fold change`, y = `-log10 p-value`), color = "red") +
    theme_bw()
    if(add_names == TRUE){
      p <- p+
        geom_text( data = filter(d, significant == TRUE ),
                   aes(x = `log2 Fold change`, y = `-log10 p-value`, label = HGNC_Symbol), 
                   size = label_size, nudge_x = 0.04, nudge_y = 0.04) }
  if(additional_title != FALSE){
    p <- p + ggtitle(paste(title, additional_title ) )
  }
  plotly::ggplotly(p) 
}

Run_GSEA <- function(DEP_result, comparison){
  ranked_list <- rowData(DEP_result) %>% 
    as_tibble() %>%
    select(HGNC_Symbol, FC = all_of( comparison) ) %>% 
    mutate(absFC = abs(FC) ) %>%
    group_by(HGNC_Symbol) %>%
    arrange(desc(absFC) ) %>%
    slice(1) %>%
    ungroup %>%
    arrange(FC)
  
  ranked_vector <- unlist(ranked_list[,2])
  names(ranked_vector) <- mapIds(org.Hs.eg.db, 
                                 keys = ranked_list$HGNC_Symbol,
                                 column = "ENTREZID", keytype = "SYMBOL")
  
  ranked_vector <- ranked_vector[!is.na(names(ranked_vector)) ]
  
  pathways <- reactomePathways(names(ranked_vector))
  fgseaRes <- fgsea(pathways, ranked_vector, maxSize=500)
  
  #print(head(fgseaRes))
  if(nrow(fgseaRes  )>1 ){
    DT::datatable(fgseaRes) %>% DT::formatRound(c(2,3,4,5,6), 3)
  }
  
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  ggplot() + theme_void()
  plotGseaTable(pathways[topPathways], ranked_vector, fgseaRes, 
                gseaParam=0.5)
}

Return_DEP_Hits_Plots <- function(data, DEP_result, comparison){
  plot.new()
  Run_GSEA(DEP_result = DEP_result, comparison = all_of( comparison) )
  if(nrow(rowData(DEP_result) %>% 
          as_tibble() %>%
          filter(significant ==TRUE)) > 0 ){
    hits <- tibble()
    hits <- rowData(DEP_result) %>% 
      as_tibble() %>%
      filter(significant ==TRUE) %>%
      as.data.frame()
    hits_all <- hits
    if(nrow(hits) > 100){
      col_padj <- which(grepl("_p.adj", colnames(hits)) )
      hits <- hits[order(hits[,col_padj]),]
      hits <- hits[1:100,]
    }
    hits_mapped <- string_db$map( hits, "HGNC_Symbol", removeUnmappedRows = TRUE )
    hits_mapped_col <- string_db$add_diff_exp_color( hits_mapped, logFcColStr= comparison )
    payload_id <- string_db$post_payload(hits_mapped_col$STRING_id, 
                                         colors=hits_mapped_col$color )
    string_db$plot_network( hits_mapped_col,  payload_id=payload_id)
    
    if(nrow(hits_all) >1){
      data <- data  %>%  
        filter(Annotated_Sequence %in% hits_all$Annotated_Sequence)
      rownames_data <- paste0(data$HGNC_Symbol, "_", data$Annotated_Sequence)
      data <- data  %>%  
        select(-HGNC_Symbol, -Annotated_Sequence,           
               -Master.Protein.Descriptions, -Master.Protein.Accessions) %>%
        as.matrix() 
      colnames(data) <- gsub("log2FC_Xenograft_", "", colnames(data))
      data <- data[, colnames(DEP_result)] 
      rownames(data) <- rownames_data
      data %>%
        pheatmap::pheatmap(col = colorRampPalette(c("navy", "white", "firebrick3"))(50), scale= "row", 
                           fontsize_row = 5, annotation_col = as.data.frame(colData(DEP_result)) %>% select(condition) )
      message("Note: Row-scaling applied for this heatmap")
      
      DT::datatable( hits_all %>% 
                       mutate(
                         Protein = paste0("<a href='https://www.uniprot.org/uniprotkb/", 
                                          Master.Protein.Accessions, "/entry'>", 
                                          HGNC_Symbol, "</a>" )) %>%
                       select(Annotated_Sequence, Protein, contains("diff"), contains("_p.adj"),
                              Master.Protein.Descriptions), 
                     class = 'cell-border stripe', rownames = FALSE, escape = F)  %>% 
        DT::formatRound(c(3,4), 5)
    }else{DT::datatable(hits_all, class = 'cell-border stripe', rownames = FALSE)}
  }else{
    protpless0.05 <- tibble()
    pless0.05 <- (rowData(DEP_result) %>% as_tibble() %>% select(p =contains("p.val")) %>% unlist ) < 0.05
    protpless0.05 <- rowData(DEP_result) %>% 
      as_tibble() %>%
      .[pless0.05,] %>%
      as.data.frame()
    if(nrow(protpless0.05)> 0){
      DT::datatable( protpless0.05 %>% 
                       mutate(
                         Protein = paste0("<a href='https://www.uniprot.org/uniprotkb/", 
                                          Master.Protein.Accessions, "/entry'>", 
                                          HGNC_Symbol, "</a>" )) %>%
                       select(Annotated_Sequence, Protein, contains("diff"), contains("_p.adj"),
                              Master.Protein.Descriptions), 
                     class = 'cell-border stripe', rownames = FALSE, escape = F)  %>% 
        DT::formatRound(c(3,4), 5)
    }
  }
}

# does not work yet
Plot_Enrichment_Single_Pathway <- function(DEP_result, comparison, pw){
  ranked_list <- rowData(DEP_result) %>% 
    as_tibble() %>%
    select(HGNC_Symbol, FC = all_of( comparison) ) %>% 
    mutate(absFC = abs(FC) ) %>%
    group_by(HGNC_Symbol) %>%
    arrange(desc(absFC) ) %>%
    slice(1) %>%
    ungroup %>%
    arrange(FC)
  
  ranked_vector <- unlist(ranked_list[,2])
  names(ranked_vector) <- mapIds(org.Hs.eg.db, 
                                 keys = ranked_list$HGNC_Symbol,
                                 column = "ENTREZID", keytype = "SYMBOL")
  
  ranked_vector <- ranked_vector[!is.na(names(ranked_vector)) ]
  
  pathways <- reactomePathways(names(ranked_vector))
  fgseaRes <- fgsea(pathways, ranked_vector, maxSize=500)
  
  #ggplot() + theme_void()
  plotEnrichment(pathways[pw],
                 ranked_vector) + labs(title=pw)
}


### Function Barplots
plot_barplot_phos <- function(dat, prot, add_title = "", specific_peptide = FALSE){
  if(prot %in% dat$HGNC_Symbol){
    if(specific_peptide != FALSE){
      dat <- dat %>% filter(Annotated_Sequence == specific_peptide)
    }
    dat %>%
      pivot_longer(cols = 3:(ncol(dat)-2), names_to = "sample", values_to = "log2_FC") %>%
      mutate(sample = gsub( "log2FC_", "", sample ) ) %>%
      separate( sample , into = c("xenograft", "treatment", "day", "replicate", "set"), sep = "_", remove = F) %>%
      filter(HGNC_Symbol == prot) %>%
      mutate(AA_1to6 = substr(Annotated_Sequence, 1,6 ) ) %>%
      mutate(treatment = as.factor(treatment)) %>%
      mutate(treatment = factor(treatment, levels = c("ctrl", "E", "EC", "EBC") )) %>%
      ggplot(aes( sample, log2_FC, fill = treatment, group = AA_1to6 )) +
      geom_col(position = "dodge") +
      theme_bw() +
      ggtitle(paste(add_title, prot)) +
      scale_fill_manual(values=PGPalette[c(5,1,2, 4)]) +
      facet_wrap(~Annotated_Sequence) +
      theme(axis.text.x = element_text(angle = 90))
  }else{
    print(paste(prot, "was not detected in dataset", deparse(substitute(dat)) ))
  }
}

plot_barplot_prot <- function(dat, prot, add_title = "", specific_peptide = FALSE){
  if(prot %in% dat$HGNC_Symbol){
    if(specific_peptide != FALSE){
      dat <- dat %>% filter(Annotated_Sequence == specific_peptide)
    }
    dat %>%
      pivot_longer(cols = 2:(ncol(dat)-2), names_to = "sample", values_to = "log2_FC") %>%
      mutate(sample = gsub( "log2FC_", "", sample ) ) %>%
      separate(sample, into = c( "type", "replicate" ), sep = "-", remove = FALSE) %>%
      filter(HGNC_Symbol == prot) %>%
      mutate(treatment = as.factor(treatment)) %>%
      mutate(treatment = factor(treatment, levels = c("ctrl", "E", "EC", "EBC") )) %>%
      ggplot(aes( sample, log2_FC, fill = type)) +
      geom_col(position = "dodge") +
      theme_bw() +
      ggtitle(paste(add_title, prot)) +
      scale_fill_manual(values=PGPalette[c(5,1,2, 4)]) +
      theme(axis.text.x = element_text(angle = 90))
  }else{
    print(paste(prot, "was not detected in dataset", deparse(substitute(dat)) ))
  }
}

### Function Boxplots
plot_boxplot_phos <- function(dat, prot, add_title = "", specific_peptide = FALSE){
  if(prot %in% dat$HGNC_Symbol){
    if(specific_peptide != FALSE){
      dat <- dat %>% filter(Annotated_Sequence == specific_peptide)
    }
    dat %>%
      pivot_longer(cols = 3:(ncol(dat)-2), names_to = "sample", values_to = "log2_FC") %>%
      mutate(sample = gsub( "log2FC_", "", sample ) ) %>%
      #mutate(prep = unlist(prep_l[sample]) ) %>%
      separate( sample , into = c("xenograft", "treatment", "day", "replicate", "set"), sep = "_", remove = F) %>%
      filter(HGNC_Symbol == prot) %>%
      mutate(treatment = as.factor(treatment)) %>%
      mutate(treatment = factor(treatment, levels = c("ctrl", "E", "EC", "EBC") )) %>%
      mutate(AA_1to6 = substr(Annotated_Sequence, 1,6 ) ) %>%
      ggplot(aes( treatment, log2_FC, fill = treatment, group = AA_1to6 )) +
      geom_boxplot(aes(group = treatment, fill = treatment), outlier.size = 0 ) +
      ggbeeswarm::geom_beeswarm() +
      theme_bw() +
      ggtitle(paste(add_title, prot)) +
      scale_fill_manual(values=PGPalette[c(5,1,2, 4)]) +
      #facet_wrap(~Annotated_Sequence+day) +
      facet_wrap(~Annotated_Sequence) +
      theme(axis.text.x = element_text(angle = 90))
  }else{
    print(paste(prot, "was not detected in dataset", deparse(substitute(dat)) ))
  }
}


Test_PhosphoData <- function(pdata, comparison, comparison_base, sig_level = 0.05){
  message(paste("Testing non-log transformed changes with t-test and adjusting p-value using FDR. Significance level =", sig_level))
  stopifnot(c("Annotated_Sequence", "HGNC_Symbol") %in% colnames(pdata) )
  comp_pdata <- pdata %>%
    pivot_longer(where(is.numeric), names_to = "sample", values_to = "log2FC" ) %>%
    mutate(sample = gsub( "log2FC_", "", sample ) ) %>%
    separate( sample , into = c("xenograft", "treatment", "day", "replicate", "set"), 
              sep = "_", remove = F) %>% 
    filter(treatment %in% c( comparison, comparison_base)) 
  tested_pdata <- sapply(unique(comp_pdata$Annotated_Sequence), function(peptide){
    peptide_data <- comp_pdata %>% 
      filter(Annotated_Sequence == peptide)
    x_comp <- peptide_data %>% filter(treatment == comparison) %>% .$ log2FC
    x_base <- peptide_data %>% filter(treatment == comparison_base) %>% .$ log2FC
    tested <- t.test(2^x_comp, 2^x_base)
    p <- tested$p.value
    diff <- mean(x_comp) - mean(x_base) 
    return(c("p" = p, "diff" = diff, "Annotated_Sequence" = peptide) )
  } )
  tested_pdata <- tested_pdata %>% t %>% as_tibble() %>%
    mutate(p = as.numeric(p), diff = as.numeric(diff)) %>%
    mutate( `-log10 p-value` = ( - log10(p) ) ) %>%
    mutate(p_adj = p.adjust(p, method = "fdr"), significant = p_adj < sig_level) 
  left_join( pdata, tested_pdata , by = "Annotated_Sequence" ) 
}

GGPlotly_Volcano_Test <- function(tested_pdata, proteins_of_interest = FALSE, 
                             additional_title = FALSE, label_size = 1.8, add_names = TRUE){
    d <- tested_pdata
  
  p <- 
    ggplot( d , aes(diff,`-log10 p-value`, color = significant, label = HGNC_Symbol, 
                    group = Annotated_Sequence ) ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.2) +
    geom_point(size = 1) +
    scale_color_manual(values = c("lightgrey", "black")) + 
    geom_point( data = filter(d, HGNC_Symbol %in% proteins_of_interest ),
                size = 1,
                aes(x = diff, y = `-log10 p-value`), color = "red") +
    theme_bw()
  if(add_names == TRUE){
    p <- p+
      geom_text( data = filter(d, significant == TRUE ),
                 aes(x = diff, y = `-log10 p-value`, label = HGNC_Symbol), 
                 size = label_size, nudge_x = 0.04, nudge_y = 0.04) }
  if(additional_title != FALSE){
    p <- p + ggtitle(paste(title, additional_title ) )
  }
  plotly::ggplotly(p) 
}

KMeans_Find_Nr_Clusters_elbow <- function(mat_kmean, c_max = 20){
  kclusts <- 
    dplyr::tibble(n_clusts = 1:c_max) %>%
    dplyr::mutate(
      kclust = purrr::map(
        n_clusts,
        ~kmeans(mat_kmean, centers = .x, iter.max = 10, nstart = 5)
      ),
      augmented = purrr::map(kclust, broom::augment, mat_kmean),
      tidied = purrr::map(kclust, broom::tidy),
      glanced = purrr::map(kclust, broom::glance)
    ) %>%
    dplyr::select(-kclust)
  
  point_assignments <- kclusts %>%
    dplyr::select(n_clusts, augmented) %>%
    tidyr::unnest(augmented)
  cluster_info <- kclusts %>%
    dplyr::select(n_clusts, tidied) %>% 
    tidyr::unnest(tidied)
  model_stats <- kclusts %>%
    dplyr::select(n_clusts, glanced) %>% 
    tidyr::unnest(glanced)
  
  # Elbow chart
  ggplot2::ggplot(data = model_stats, aes(n_clusts, tot.withinss)) +
    ggplot2::geom_line() +
    #ggplot2::scale_x_continuous(limits = c(1, 12), breaks = seq(1, 12, 1)) +
    geom_point() +
    theme_bw() +
    ggplot2::ggtitle("Total within sum of squares, by # clusters")
}
