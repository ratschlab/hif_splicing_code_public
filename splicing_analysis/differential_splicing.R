
library("data.table")
library("RColorBrewer")
library("DESeq2")
library("ggplot2")
library("MASS")
library("pheatmap")
library("plyr")

get_PSI_values <- function(res_table, psi_table_file){

  # format PSI table
  psi_table = data.frame(fread(psi_table_file))

  psi_table_val = psi_table[,grep("psi", colnames(psi_table), value=T)]
  hyp_siCTRL_psi = rowMeans(psi_table_val[,grep("Hypoxia_siControl", colnames(psi_table_val))], na.rm=T)
  other_psi = rowMeans(psi_table_val[,grep("Hypoxia_siControl", colnames(psi_table_val), invert=T)], na.rm=T)

  delta_PSI = hyp_siCTRL_psi - other_psi

  hyp_siCTRL_psi = rowMeans(psi_table_val[,grep("Hypoxia_siControl", colnames(psi_table_val))], na.rm=T)
  hyp_siARNT_psi = rowMeans(psi_table_val[,grep("Hypoxia_siARNT", colnames(psi_table_val))], na.rm=T)
  norm_siCTRL_psi = rowMeans(psi_table_val[,grep("Normoxia_siControl", colnames(psi_table_val))], na.rm=T)
  norm_siARNT_psi = rowMeans(psi_table_val[,grep("Normoxia_siARNT", colnames(psi_table_val))], na.rm=T)

  relative_delta_PSI = (hyp_siCTRL_psi-norm_siCTRL_psi) - (hyp_siARNT_psi-norm_siARNT_psi)

  psi_table = data.frame(event_id = psi_table$event_id, delta_PSI, relative_delta_PSI)

  res_table = join(psi_table, res_table, type = "right")

  return(res_table)
}

get_event_coor <- function(res_table, outname){

  spladder_file = paste(indir, "/merge_graphs_", outname, "_C3.confirmed.icgc.txt", sep="")
  spladder_res = data.frame(fread(spladder_file))
  spladder_res = spladder_res[,c(1,3,4,5)]

  total_res = merge(spladder_res, res_table, by="event_id")

  return(total_res)

}

filter_threshhold_counts <- function(countsTable, val){
  #ensure that every row has atleast n counts in every column

  countsTableFiltered = countsTable[rowSums(countsTable>=val) == ncol(countsTable),]

  return(countsTableFiltered)
}

read_in_spladder_counts <- function(countTable_file1){

	# read in the table and get the columns with the counts
	countsTable = data.frame(fread(countTable_file1))

	gene_info = countsTable[,c(1,2)]

	col_interest = c(grep("gene_exp[.]", colnames(countsTable)), grep("event_count[.]", colnames(countsTable)))
	countsTable = countsTable[,col_interest]
	countsTable = ceiling(countsTable)

	colnames(countsTable) = sapply(strsplit(colnames(countsTable), "[.]sorted"),function(x) x[1])

	countsTable = cbind(gene_info, countsTable)
	return(countsTable)
}

filter_and_plot <- function(countsTable, design, plot_colors, RESULTS_FOLDER, plot_main, run_splice_only=TRUE){

	#get plot labels
	plot_labels = colnames(countsTable)

	#normalize
  thresh = 10
  if(run_splice_only){
    thresh = 1
  }
	countsTableFiltered = filter_threshhold_counts(countsTable, thresh)
  print("##################THRESH##################")
  print(thresh)
  print("####################################")

	#num_genes = nrow(norm_counts)

	#plot the boxplots
	exp_idx = which(design$count_type == "exp")
	event_idx = which(!design$count_type == "exp")

	pdf(paste(RESULTS_FOLDER, "exp_log2.pdf", sep=""))
	boxplot(log2(countsTableFiltered[,exp_idx]+1), las=2, main="Raw log2 Values of Expression", col=plot_colors[exp_idx])
	dev.off()

	if(length(event_idx) > 0){
		pdf(paste(RESULTS_FOLDER, "event_log2.pdf", sep=""))
			boxplot(log2(countsTableFiltered[,event_idx]+1), las=2, main="Raw log2 Values of Expression", col=plot_colors[event_idx])
		dev.off()
		plot_pca(countsTableFiltered[,event_idx], paste(RESULTS_FOLDER, "event_pca.pdf", sep=""), plot_colors[event_idx], plot_labels[event_idx], plot_main)
		hierarchical_clust(countsTableFiltered[,event_idx], paste(RESULTS_FOLDER, "event_samp_clust.pdf", sep=""), plot_main)
	}


	#plot the pca of samples
	plot_pca(countsTableFiltered[,exp_idx], paste(RESULTS_FOLDER, "exp_pca.pdf", sep=""), plot_colors[exp_idx], plot_labels[exp_idx], plot_main)

	# plot the clustering
	hierarchical_clust(countsTableFiltered[,exp_idx], paste(RESULTS_FOLDER, "exp_samp_clust.pdf", sep=""), plot_main)


	return(countsTableFiltered)

}

translate_ENS_to_HGNC <- function(in_table, translation_file){
	in_table$id = remove_period(in_table$id)
  annot_table = data.frame(fread(translation_file,header=T), check.names=F)
  colnames(annot_table) = c("id", "ensembl_trans_id", "hgnc_symbol")
  annot_table = unique(annot_table[,c(1,3)])

  in_table = merge(annot_table, in_table, by="id", all.y=T)
  in_table$hgnc_symbol[is.na(in_table$hgnc_symbol)] = ""
  in_table$hgnc_symbol[in_table$hgnc_symbol == ""] = in_table$id[in_table$hgnc_symbol == ""]
  return(in_table)

}


write_results_individial_analysis <- function(res, countsTable, design, gene_event_table, outname, RESULTS_FOLDER, psi_table_file, translation_file){

  # translate gene ids
	res = merge(gene_event_table, res, by="event_id")
	res_trans = translate_ENS_to_HGNC(res, translation_file)
	res_trans = res_trans[ order(res[,c("p_adj")]), ]


	# plot heat map of splicing counts
	event_count_table <- countsTable[,grep("event", colnames(countsTable))]
  event_count_table = data.frame(t(scale(t(event_count_table))))
  sample_ids = colnames(event_count_table)
  event_count_table$event_id = row.names(event_count_table)

  event_count_table = join(res_trans, event_count_table)
  event_count_table = event_count_table[order(event_count_table$p_adj, decreasing=F),]

  # get PSI vslues
  event_count_table = get_PSI_values(event_count_table, psi_table_file)

  # get coordinates
  event_count_table = get_event_coor(event_count_table, outname)
  event_count_table = event_count_table[order(event_count_table$p_adj, decreasing=F), ]


  #write out table
  write.table(event_count_table, paste(RESULTS_FOLDER, "/", outname, ".tsv", sep=""), quote=F, col.name=T, row.name=F, sep="\t")



	annot_info <- design[design$count_type=="event",c("o2_status", "hif_status")]
  row.names(annot_info) = sample_ids


	pdf(paste(RESULTS_FOLDER, "/", outname, "_heatmap_total.pdf", sep=""))
		try(pheatmap(event_count_table[,sample_ids], annotation_col=annot_info, cluster_cols=TRUE, show_rownames=FALSE))
	dev.off()

  event_count_table_sig = event_count_table[event_count_table$p_adj < 0.1,]
  event_count_table_sig = event_count_table_sig[order(abs(event_count_table_sig$estimate), decreasing=T),]
  print(dim(event_count_table_sig))

	pdf(paste(RESULTS_FOLDER, "/", outname, "_heatmap_gene_names.pdf", sep=""))
		try(pheatmap(event_count_table_sig[1:50,sample_ids], annotation_col=annot_info,
              cluster_cols=TRUE, labels_row=event_count_table_sig$hgnc_symbol[1:50],
              method="ward.D2"))
	dev.off()

	pdf(paste(RESULTS_FOLDER, "/", outname, "_heatmap_event_ids.pdf", sep=""))
		try(pheatmap(event_count_table_sig[1:50,sample_ids], annotation_col=annot_info,
              cluster_cols=TRUE, labels_row=event_count_table_sig$event_id[1:50]))
	dev.off()

	pdf(paste(RESULTS_FOLDER, "/", outname, "_heatmap_sig.pdf", sep=""))
		try(pheatmap(event_count_table_sig[,sample_ids],
              annotation_col=annot_info, cluster_cols=TRUE, show_rownames=FALSE))
	dev.off()


	# plot the clustering
	hierarchical_clust(event_count_table_sig[,sample_ids],
                      paste(RESULTS_FOLDER, "/", outname, "_clust_sig.pdf", sep=""),
                      outname)




	pdf(paste(RESULTS_FOLDER, "/", outname, "_fitted.pdf", sep=""))
	for(idx in 1:50){
		plot_counts = data.frame(count=unlist(event_count_table_sig[idx,sample_ids]),annot_info)

		curr_gene_name = event_count_table_sig[which(event_count_table_sig$event_id == gene_event_table$event_id[idx]),"hgnc_symbol"]
		curr_event = event_count_table_sig[which(event_count_table_sig$event_id == gene_event_table$event_id[idx]),"event_id"]

		curr_plot <- ggplot(plot_counts, aes(x=o2_status, y=count, color=hif_status, group=hif_status)) +
			geom_point() + stat_smooth(se=FALSE,method="loess") + ggtitle(paste(curr_gene_name, "\n", curr_event, sep=""))
		print(curr_plot)
	}
	dev.off()


	return(event_count_table)


}

get_nb_corrected_value <- function(res, linear_model_df){


  coef = summary(res)$coefficients
  intercept = coef[1,1]
  GE_est = coef[2,1]
  SF_est = coef[3,1]
  O2_est = coef[4,1]
  hif_est = coef[5,1]
  int_est = coef[6,1]

  GE_stat = linear_model_df$GE
  SF_stat = linear_model_df$sizeFactor
  o2_stat = c(rep(1,3), rep(0,6), rep(1,3))
  hif_stat = c(rep(0,6), rep(1,6))
  intersect = o2_stat*hif_stat

  gene_proportion = intercept + GE_est*GE_stat + SF_est*SF_stat # + O2_est*o2_stat + hif_est*hif_stat
  corrected_value =  linear_model_df$Y - exp(gene_proportion)

  return(corrected_value)

}


do_individual_DE_analysis <- function(design, plot_colors, countsTable, RESULTS_FOLDER, plot_main, gene_event_table, outname, psi_table_file, translation_file, run_splice_only=TRUE){


  	#now that we have the conditions we are interested in comparing, lets normalize
  	#read table
  	rownames(countsTable) <- countsTable$event_id
  	countsTable <- countsTable[,-c(1,2)]

  	# filter low counts and qc-plot
  	countsTable_filtered = filter_and_plot(countsTable, design, plot_colors, RESULTS_FOLDER, plot_main, run_splice_only)

    # make DESeq object so we can get estimated dispersions
    dds_event <- DESeqDataSetFromMatrix(countData = countsTable_filtered[,grep("exp", colnames(countsTable_filtered))],
      colData = design[design$count_type=="event",],
      design = ~ o2_status + hif_status + o2_status:hif_status)

    dds_event = estimateSizeFactors(dds_event)
    dds_event <- DESeq(dds_event)
    event_dispersions = dispersions(dds_event)

    print("finished expression modeling")

    # now we have dispersions, we can run our own model
    total_res = data.frame(estimate=rep(0, nrow(countsTable_filtered)),
                          p_val=rep(-1, nrow(countsTable_filtered)), stringsAsFactors=F)

    corrected_splice_counts = countsTable_filtered[,grep("event", colnames(countsTable_filtered))]

    if(run_splice_only){
      for(idx in 1:nrow(countsTable_filtered)){
        linear_model_df = data.frame(Y = unlist(countsTable_filtered[idx,grep("event", colnames(countsTable_filtered))]),
                                     design[design$count_type=="event",],
                                     GE=unlist(countsTable_filtered[idx,grep("exp", colnames(countsTable_filtered))]),
                                     sizeFactor=sizeFactors(dds_event))

        res = tryCatch(glm.nb(Y ~ GE + sizeFactor + o2_status + hif_status + o2_status:hif_status,
                      data=linear_model_df),error=function(x) "error")
        if(res == "error"){
          estimate = NA
          p_val = NA
          total_res[idx,]= data.frame(estimate, p_val)

          # get corrected value
          corrected_splice_counts[idx,] = NA

        }else{
          estimate = summary(res)$coefficients[6,1]
          p_val = summary(res)$coefficients[6,4]
          total_res[idx,]= data.frame(estimate, p_val)

          # get corrected value
          corr_value = get_nb_corrected_value(res, linear_model_df)
          corrected_splice_counts[idx,] = corr_value
        }

      }
      total_res$p_adj = p.adjust(total_res$p_val, method="BH")

    }else{
      # run splice and expr
      col_data=grep("exp|event", colnames(countsTable_filtered))
      dds_full <- DESeqDataSetFromMatrix(countData = countsTable_filtered[,col_data],
        colData = design,
        design = ~ o2_status + hif_status + count_type + o2_status:hif_status)

      dds_full = estimateSizeFactors(dds_full)
      dds_full <- DESeq(dds_full)

      total_res = results( dds_full, name = c("o2_statusnorm.hif_statusnoHIF") )
      colnames(total_res) = c("baseMean", "estimate", "se", "stat", "p_val", "p_adj")
      total_res = data.frame(total_res)
    }

    print("finished test")
    total_res$event_id = names(dds_event)
    total_res = total_res[,c("event_id", "estimate", "p_val", "p_adj")]
    total_res = total_res[order(total_res$p_adj,decreasing=F),]

    total_res_trans = write_results_individial_analysis(total_res,
                      corrected_splice_counts, design, gene_event_table,
                      outname, RESULTS_FOLDER, psi_table_file, translation_file)

    print("finished writing pt1")

    # filter out anything that is differential on expression in hyp
    genes_intersect_hyp = c()
    expression_interaction_res = results( dds_event, name = c("o2_statusnorm.hif_statusnoHIF") )
    expression_res = results( dds_event, name = c("o2_status_norm_vs_hyp") )
    expression_res = data.frame(expression_res)
    print(expression_res["exon_skip_926",])
    if(nrow(expression_res[expression_res$padj < 0.05 & abs(expression_res$log2FoldChange) > 1,]) > 0){
      expression_res = data.frame(expression_res[expression_res$padj < 0.05 & abs(expression_res$log2FoldChange) > 1,])
      genes_exp_pass = row.names(expression_res)
      genes_intersect = intersect(genes_exp_pass, total_res_trans$event_id[total_res_trans$p_adj < 0.1])

    }
    if(nrow(expression_interaction_res[expression_interaction_res$padj < 0.05 & abs(expression_interaction_res$log2FoldChange) > 1,]) > 0){
      expression_interaction_res = data.frame(expression_interaction_res[expression_interaction_res$padj < 0.05 & abs(expression_interaction_res$log2FoldChange) > 1,])
      genes_exp_pass = row.names(expression_interaction_res)
      genes_interaction_intersect = intersect(genes_exp_pass, total_res_trans$event_id[total_res_trans$p_adj < 0.1])
    }
    if(run_splice_only){
      genes_intersect = union(genes_intersect, genes_interaction_intersect)
    }

    total_res_trans_noRNA = total_res_trans[which(!total_res_trans$event_id %in% genes_intersect & total_res_trans$p_adj < 0.1), ]
    total_res_trans_noRNA = total_res_trans_noRNA[abs(total_res_trans_noRNA$delta_PSI) >= 0.05 | abs(total_res_trans_noRNA$relative_delta_PSI) >= 0.05, ]
    if(nrow(total_res_trans_noRNA) > 0){
      write.table(total_res_trans_noRNA,
                  paste(RESULTS_FOLDER, "/", outname, "_only.tsv", sep=""),
                  sep="\t", quote=F, row.name=F)

    }
    print("finished writing pt2")

    return(total_res_trans)
}

read_and_process_spladder_counts <- function(indir, outname){

  	countTable_file = paste(indir, "/testing_ctrl_vs_arnt/test_results_extended_C3_", outname, ".tsv", sep="")
  	curr_counts_table = read_in_spladder_counts(countTable_file)
  	gene_event_table = curr_counts_table[,c(1,2)]
  	colnames(gene_event_table) = c("event_id", "id")

  	exp_norm_ctrl = paste(rep("exp_norm_ctrl", 3), c("A", "B", "C"), sep="_")
  	exp_norm_noHIF = paste(rep("exp_norm_noHIF", 3), c("A", "B", "C"), sep="_")
  	event_norm_ctrl = paste(rep("event_norm_ctrl", 3), c("A", "B", "C"), sep="_")
  	event_norm_noHIF = paste(rep("event_norm_noHIF", 3), c("A", "B", "C"), sep="_")
  	exp_hyp_ctrl = paste(rep("exp_hyp_ctrl", 3), c("A", "B", "C"), sep="_")
  	exp_hyp_noHIF = paste(rep("exp_hyp_noHIF", 3), c("A", "B", "C"), sep="_")
  	event_hyp_ctrl = paste(rep("event_hyp_ctrl", 3), c("A", "B", "C"), sep="_")
  	event_hyp_noHIF = paste(rep("event_hyp_noHIF", 3), c("A", "B", "C"), sep="_")
  	col_names = c("event_id", "gene_id", exp_norm_ctrl, exp_hyp_ctrl, exp_hyp_noHIF, exp_norm_noHIF,
                  event_norm_ctrl, event_hyp_ctrl, event_hyp_noHIF, event_norm_noHIF)


  	colnames(curr_counts_table) = col_names

    return(list(curr_counts_table, gene_event_table))
}

run <- function(outname, indir, RESULTS_FOLDER, translation_file, run_splice_only=TRUE){

  res = read_and_process_spladder_counts(indir, outname)
  curr_counts_table = res[[1]]
  gene_event_table = res[[2]]

	o2_status = c(rep("norm", 3), rep("hyp", 6), rep("norm", 3), rep("norm", 3), rep("hyp", 6), rep("norm", 3))
	hif_status = c(rep("hif", 6), rep("noHIF", 6), rep("hif", 6), rep("noHIF", 6))
  count_type = c(rep("exp", 12), rep("event", 12))

	design = data.frame(cbind(o2_status, hif_status, count_type))

	pallete = brewer.pal(8,"Dark2")
	plot_colors = rep(pallete, each=3)
	plot_main = "Splicing Test"

  psi_table_file = paste(indir, "merge_graphs_", outname, "_C3.confirmed.txt", sep="")


	res = do_individual_DE_analysis(design, plot_colors, curr_counts_table, RESULTS_FOLDER, plot_main, gene_event_table, outname, psi_table_file, translation_file, run_splice_only)

  return(res)

}

#plotting PCA of samples
plot_pca <-function(counting_table, plot_name, plot_colors, plot_labels, plot_main, takeLog=TRUE){

  ### plot_name is the name of the PDF, include the path is working directory is not set
  ### plot_colors is a vector of colors consisting of how each sample's color (this is usually based on column names
  ### counting_table is the normalized gene counts in a data frame samples are the columns, genes are the rows
  log.samples <- t(counting_table)
  if(takeLog){
    log.samples <- log2(t(counting_table)+1)
  }
  ir.pca <- prcomp(log.samples)
  vars <- apply(ir.pca$x, 2, var)
  props <- vars / sum(vars)
  pc1_var_importance = as.character(format(props[1], digits=3))
  pc2_var_importance = as.character(format(props[2], digits=3))
  pc3_var_importance = as.character(format(props[3], digits=3))

  plot_main = paste(plot_main, 'PC1:', pc1_var_importance, 'PC2:', pc2_var_importance, 'PC3:', pc3_var_importance, sep = ' ')

  pdf(plot_name)
  pairs(ir.pca$x[,1:3], panel = function(x,y) text(x,y, labels=plot_labels,cex=0.5,col=plot_colors), main=plot_main)
  dev.off()

  print(ir.pca$x[,1:2])


  return(ir.pca)

}

#look at the sample clustering
hierarchical_clust <-function(counting_table, plot_name, plot_main){

  ### plot_name is the name of the PDF, include the path is working directory is not set
  ### plot_colors is a vector of colors consisting of how each sample's color (this is usually based on column names
  ### counting_table is the normalized gene counts in a data frame samples are the columns, genes are the rows


  #do ehirarcical clustering
  d <- dist(as.matrix(t(counting_table)))   # find distance matrix
  hc <- hclust(d)                # apply hirarchical clustering

  pdf(plot_name)
  plot(hc, main=plot_main)
  dev.off()


}


get_gene_pos <- function(RESULTS_FOLDER, indir, event_types){

	for(curr_event in event_types){

		spladder_file = paste(indir, "/merge_graphs_", curr_event, "_C3.confirmed.icgc.txt", sep="")
		test_file = paste(RESULTS_FOLDER, curr_event, ".tsv", sep="")

		spladder_res = data.frame(fread(spladder_file))
		spladder_res = spladder_res[,c(1,3,4,5,6)]

		test_res = data.frame(fread(test_file))

		total_res = merge(spladder_res, test_res, by="event_id")
		total_res = total_res[order(total_res$padj, decreasing=F), ]

		total_res_alt_reg = strsplit(total_res$alt_region_coordinates, ":")

		write.table(total_res, paste(RESULTS_FOLDER, curr_event, "_event_coor.tsv", sep=""), sep="\t", quote=F, row.name=F)

	}
}


get_all_splice_test_results <- function(RESULTS_FOLDER, run_splice_only, indir){

  curr_psi_file = PSI_TCGA_FILE_old
  if(run_splice_only){
    curr_psi_file = PSI_TCGA_FILE_newest
  }

  es_res = data.frame(fread(paste(RESULTS_FOLDER, "/exon_skip.tsv", sep="")))
  a3_res = data.frame(fread(paste(RESULTS_FOLDER, "/alt_3prime.tsv", sep="")))
  a5_res = data.frame(fread(paste(RESULTS_FOLDER, "/alt_5prime.tsv", sep="")))

  total_event = rbind(es_res, a3_res, a5_res)

  # get events in TCGA and significant
  events_tcga = data.frame(fread(curr_psi_file))
  events_tcga = events_tcga$event_id
  events_tcga = intersect(events_tcga, total_event$event_id)
  total_event = total_event[which(total_event$event_id %in% events_tcga),]

  # get expression info
  es_res = read_and_process_spladder_counts(indir, "exon_skip")[[1]]
  a3_res = read_and_process_spladder_counts(indir, "alt_3prime")[[1]]
  a5_res = read_and_process_spladder_counts(indir, "alt_5prime")[[1]]
  total_exp = rbind(es_res, a3_res, a5_res)
  total_exp = total_exp[,c(1,2,grep("exp_", colnames(total_exp)))]

  total_res = join(total_exp, total_event, type="inner")

  return(total_res)

}

make_heatmap_plot <- function(RESULTS_FOLDER, run_splice_only, indir){


    total_res = get_all_splice_test_results(RESULTS_FOLDER, run_splice_only, indir)

    exp_ids = grep("exp_norm|exp_hyp", colnames(total_res), value=T)
    splice_ids = grep("event_norm|event_hyp", colnames(total_res), value=T)

    total_res_exp = t(scale(t(total_res[,exp_ids])))
    total_res_spl = t(scale(t(total_res[,splice_ids])))

    total_res = cbind(total_res[,c("event_id", "hgnc_symbol")], total_res_exp, total_res_spl)

    sample_ids = c(exp_ids, splice_ids)
    if(run_splice_only){
      sample_ids = splice_ids
    }
    nx_idx = grep("norm",sample_ids)
    hif_idx = grep("noHIF",sample_ids)

    nx_col = rep("Hypoxic", length(sample_ids))
    nx_col[nx_idx] = "Normoxic"

    hif_col = rep("HIF", length(sample_ids))
    hif_col[hif_idx] = "noHIF"

    design = data.frame(nx_col, hif_col)

  	annot_info_col <- design
    row.names(annot_info_col) = sample_ids

    annot_info_row = rep("NA", nrow(total_res))
    annot_info_row[grep("exon_skip", total_res$event_id)] = "exon_skip"
    annot_info_row[grep("intron_retention", total_res$event_id)] = "intron_retention"
    annot_info_row[grep("alt_3prime", total_res$event_id)] = "alt_3prime"
    annot_info_row[grep("alt_5prime", total_res$event_id)] = "alt_5prime"

    annot_info_row = data.frame(annot_info_row)
    total_res_names = paste(total_res$hgnc_symbol, " (", total_res$event_id, ")", sep="")
    row.names(total_res) = total_res_names
    row.names(annot_info_row) = row.names(total_res)



  	pdf(paste(RESULTS_FOLDER, "/heatmap_sig_all.pdf", sep=""), width=10, height=10)
  		try(pheatmap(total_res[,sample_ids],
                annotation_col=annot_info_col, annotation_row=annot_info_row,
                clustering_method="ward.D2", cluster_cols=TRUE, show_rownames=FALSE))
  	dev.off()

}


# get the parameters from the runner script
args = commandArgs(trailingOnly=TRUE)
indir = as.character(args[1])   ##### where splice results are located
RESULTS_FOLDER = as.character(args[2])   ##### where you want all results to be written to
run_splice_only = as.character(args[3])  ##### TRUE -> run splicing independent of expression
                                        ##### FALSE -> run splicing and expression
translation_file = as.character(args[4])  ##### file for translating gene_ids


if(run_splice_only){
  ############################ RUN ONLY SPLICE ############################
  es_res = run(outname="exon_skip", indir=indir, RESULTS_FOLDER=RESULTS_FOLDER, translation_file, run_splice_only=TRUE)
  ir_res = run(outname="intron_retention", indir=indir, RESULTS_FOLDER=RESULTS_FOLDER, translation_file, run_splice_only=TRUE)

  a3_res = run(outname="alt_3prime", indir=indir, RESULTS_FOLDER=RESULTS_FOLDER, translation_file, run_splice_only=TRUE)
  a5_res = run(outname="alt_5prime", indir=indir, RESULTS_FOLDER=RESULTS_FOLDER, translation_file, run_splice_only=TRUE)

  make_heatmap_plot(RESULTS_FOLDER, run_splice_only=TRUE, indir)

}else{
  ############################ RUN BOTH SPLICE and EXPR ############################

  es_res = run(outname="exon_skip", indir=indir, RESULTS_FOLDER=RESULTS_FOLDER, translation_file, run_splice_only=FALSE)
  ir_res = run(outname="intron_retention", indir=indir, RESULTS_FOLDER=RESULTS_FOLDER, translation_file, run_splice_only=FALSE)

  a3_res = run(outname="alt_3prime", indir=indir, RESULTS_FOLDER=RESULTS_FOLDER, translation_file, run_splice_only=FALSE)
  a5_res = run(outname="alt_5prime", indir=indir, RESULTS_FOLDER=RESULTS_FOLDER, translation_file, run_splice_only=FALSE)

  make_heatmap_plot(RESULTS_FOLDER, run_splice_only=FALSE, indir)

}
