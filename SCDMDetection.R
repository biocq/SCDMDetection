get_DEGs<-function(expression_matrix, samples, sample_origins = NULL, beta = 2, gamma = 0.05)
{
	
	if(is.null(sample_origins))
	{
		sample_origins = rep("tumor", ncol(expression_matrix))
		sample_origins[substr(colnames(expression_matrix), nchar(colnames(expression_matrix)[1])-1, nchar(colnames(expression_matrix)[1])) == "11"] = "normal"	
		if(length(which(sample_origins == "normal")) < 1)
		{	
			print("没有正常样本，不能进行差异表达分析，终止运行。")
			return()
		}
	}
	DEGs = list()
	for(sample in samples)
	{
		diff_genes = get_diff_expressed_genes(expression_matrix, sample, sample_origins, beta, gamma)
		DEGs[[sample]] = diff_genes 
	}
	return(DEGs)
}

get_diff_expressed_genes<-function(expression_matrix, sample, sample_origins, beta = 2, gamma = 0.05)
{
	if(length(which(sample_origins == "normal")) < 1)
	{	
			print("没有正常样本，不能进行差异表达分析，终止运行。")
			return()
	}
	df = data.frame(expression_matrix[, sample])
	healthy_tissues = which(sample_origins == "normal")
	df = cbind(df, expression_matrix[, healthy_tissues])
	colnames(df)[1] = paste(sample)
	coldata = data.frame(rep("healthy", ncol(df)))
	colnames(coldata) = "condition"
	rownames(coldata) = colnames(df)
	coldata$condition = factor(coldata$condition, levels = c("healthy", "tumor"))
	coldata[1, "condition"] = "tumor"
	#countData is a data frame with genes on rows and patients on columns
	#coldata is a data frame with patients on the rows and "tumor"/"healthy" on the first row
	diff_expression_count_matrix = DESeqDataSetFromMatrix(countData = round(df), 
    colData = coldata, 
    design = ~ condition)
	diff_expression_count_matrix = diff_expression_count_matrix[rowSums(counts(diff_expression_count_matrix)) > 1, ]
	diff_expression_count_matrix$condition = factor(diff_expression_count_matrix$condition, levels = c("tumor", "healthy"))
	diff_expression_count_matrix = DESeq(diff_expression_count_matrix, betaPrior = TRUE)
	res = results(diff_expression_count_matrix, alpha = gamma, lfcThreshold = 1, altHypothesis = "greaterAbs")
	res = res[!is.na(res$log2FoldChange), ]
	res = res[!is.na(res$pvalue), ]
	diff_genes = res[res$pvalue<gamma & abs(res$log2FoldChange) > beta, 2]
	names(diff_genes) = rownames(res)[res$pvalue<gamma & abs(res$log2FoldChange) > beta]
	return(abs(diff_genes))
}

get_enriched_pathways<-function(diff_genes, pathwayDB, expression_matrix_genes, pathwayDB_nodes, delta = 0.05)
{
	enrichment_pvalue = c()
	enrichment_pvalue = sort(unlist(lapply(pathwayDB_nodes, function(x) get_enrichment_pvalue(x, diff_genes, expression_matrix_genes))))
	bins = list()
	if(length(which(p.adjust(enrichment_pvalue, method = "BH") < delta)) > 0)
	{
		enrichment_pvalue = p.adjust(enrichment_pvalue, method = "BH")
	}
	pathways_to_include = names(enrichment_pvalue)[enrichment_pvalue < delta]
	if(length(pathways_to_include) == 0){ return(NULL) }
	for(i in 1:length(pathways_to_include))
	{
		pathway_nodes = pathwayDB_nodes[[pathways_to_include[i]]]
		DEGs_in_pathway = intersect(names(diff_genes), pathway_nodes)
		bins[[i]] = DEGs_in_pathway
	}
	names(bins) = pathways_to_include
	return(bins)
}

get_enrichment_pvalue<-function(pathway_nodes, diff_genes, expression_matrix_genes)
{
	target_set = intersect(pathway_nodes, expression_matrix_genes)
	if(length(target_set) == 0)
	{
		return(1) 
	}
	background_set = setdiff(expression_matrix_genes, target_set)
	drawn = names(diff_genes)
	special_drawn = intersect(target_set, drawn)
	if(length(special_drawn) == 0)
	{
		return(1) 
	}
	return(phyper(length(special_drawn), length(target_set), length(background_set), length(drawn), lower.tail = F))	
}

get_pathway_network<-function(pathway_graph, original_network)
{
	#delete duplicated edges
	matches = match_df(data.frame(src = original_network[, 1], dest = original_network[, 2]), as.data.frame(rbind(pathway_graph[, c(1, 2)], pathway_graph[, c(2, 1)])))
	if(nrow(matches) > 0) { 
		original_network = original_network[-as.numeric(rownames(matches)), ] 
	}
	#add the edges to the background network
	#all edges have constant score of 0.1
	original_network = rbind(original_network, matrix(c(pathway_graph[, 1], pathway_graph[, 2], 
					rep(0.1, nrow(pathway_graph))), ncol = 3, nrow = nrow(pathway_graph)))
	#delete self-loops
	original_network = original_network[which(original_network[, 1] != original_network[, 2]), ]
	#use only edgeg where one of the ends (or both) is in the pathway. 
	pathway_nodes = unique(c(as.character(pathway_graph[, 1]), as.character(pathway_graph[, 2])))
	pathway_network = original_network[which(original_network[, 1] %in% pathway_nodes | original_network[, 2] %in% pathway_nodes), ]
	return(pathway_network)
}

run_single_PCSF<-function(pathway_network, driver_gene, original_network, pathway_prizes, alpha)
{
	#add edges that touch the driver and add negative prizes for these nodes
	pathway_nodes = unique(c(pathway_network[, 1], pathway_network[, 2]))
	lines_to_add = which(original_network[, 1] == driver_gene | original_network[, 2] == driver_gene)
	#dont add edges that already touch pathway nodes because we already added them
	pathway_network = rbind(pathway_network, matrix(c(original_network[lines_to_add, 1], 
				original_network[lines_to_add, 2], as.numeric(original_network[lines_to_add, 3])), nrow = length(lines_to_add), ncol = 3))
	if(length(which(duplicated(pathway_network))) > 0) { 
		pathway_network = pathway_network[-which(duplicated(pathway_network)), ] 
	}
	pathway_network_igraph = graph_from_data_frame(as.data.frame(pathway_network[, c(1, 2)]), directed = F)
	#if the putative driver cannot reach any gene in the pathway there is no need to calculate influence score
	if(all(is.infinite(distances(pathway_network_igraph, v = driver_gene, to = pathway_nodes))))
	{
		return(c(0, 0))
	}
	E(pathway_network_igraph)$weight = as.numeric(pathway_network[, 3])
	degs = degree(pathway_network_igraph, v = V(pathway_network_igraph))
	terminals = rep(0, length(V(pathway_network_igraph)))
	names(terminals) = V(pathway_network_igraph)$name
	#DEGs have positive prizes
	terminals[names(pathway_prizes)] = pathway_prizes
	#Steiner nodes have negative penalties as a function of their degree and alpha
	terminals[which(!names(terminals) %in% names(pathway_prizes))] = -degs[which(!names(terminals) %in% names(pathway_prizes))]^alpha
	#assign a high prize to the driver so it will be contained in the solution
	driver_prize = sum(pathway_prizes) * 10
	terminals[driver_gene] = driver_prize
	res = PCSF(pathway_network_igraph, terminals = terminals, w = max(round(sum(pathway_prizes)/10), 20), b = 1, mu = 0)
 	if(length(V(res)) == 0) {
 	  score = c(0, 0)
	}else if(is_connected(res) & (driver_gene %in% V(res)$name)) {
		score = c(sum(V(res)$prize) - sum(E(res)$weight) - driver_prize -degs[driver_gene]^alpha, sum(V(res)$prize) - driver_prize)
	}else { score = c(-1, -1) }
	return(score)
}


PRODIGY<-function(mutated_genes, expression_matrix, network = NULL, sample, diff_genes = NULL, alpha = 0.05, pathway_list = NULL, 
			num_of_cores = 1, sample_origins = NULL, write_results = T, results_folder = "./", 
			beta = 2, gamma = 0.05, delta = 0.05)
{
	
	#if no network is specified, the network derived from STRING is used as in the original publication
	if(is.null(network))
	{
		data(STRING_network)
		network = STRING_network
	}
	network[, "score"] = sapply(as.numeric(network[, "score"]), function(x) min(x, 0.8))
	network[, "score"] = 1-as.numeric(network[, "score"])
	mutated_genes = mutated_genes[mutated_genes %in% unique(c(network[, 1], network[, 2]))]
	if(length(mutated_genes) < 1) {
	 	print("No mutated gene in large PPI network, aborting")
		return() 
	}
	original_network = network
	network = graph_from_data_frame(network, directed = F)
	if(is.null(pathway_list))
 {
		print("no pathway list. using Reactome as default")
		pathway_list = get_pathway_list_from_graphite(source = "reactome", minimal_number_of_nodes = 10, num_of_cores = num_of_cores)
	}
	#get differentially expressed genes list
	if(is.null(diff_genes))
	{
		expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[, 1], network[, 2]))), ]
		#check if sample has expression data
		if(!(sample %in% colnames(expression_matrix)))
		{
			print("sample is missing expression information, aborting")
			return()
		}
		#if sample_origins = NULL we assume the matrices follow the TCGA convention (tumors suffix is "-01" and normals are "-11")
		if(is.null(sample_origins))
		{
			sample_origins = rep("tumor", ncol(expression_matrix))
			sample_origins[substr(colnames(expression_matrix), nchar(colnames(expression_matrix)[1])-1, nchar(colnames(expression_matrix)[1])) == "11"] = "normal"	
		}
		if(length(which(sample_origins == "normal")) < 1)
		{	
			print("没有正常样本，不能进行差异表达分析，终止运行。")
			return()
		}
		print("分析差异表达基因...")	
		diff_genes = get_diff_expressed_genes(expression_matrix, sample, sample_origins, beta, gamma)
	}
 pathway_list_genes = setNames(lapply(pathway_list, function(x) get_pathway_genes_from_network(x)), names(pathway_list))
	#get enriched pathways. Bins contain the DEGs belonging to each pathway
	bins = get_enriched_pathways(diff_genes, pathwayDB = NULL, rownames(expression_matrix), pathwayDB_nodes = pathway_list_genes, delta)
	if(length(bins) == 0){
		print("没有富集的通路，终止运行。")
		return() 
	}
	print(paste("检查", length(mutated_genes), "个突变和", length(bins), "个通路。"))
	Influence_matrix = matrix(ncol = length(mutated_genes))
	colnames(Influence_matrix) = mutated_genes
	Fraction_of_DEGs_matrix = matrix(ncol = length(mutated_genes))
	colnames(Fraction_of_DEGs_matrix) = mutated_genes
	#Main part of the algorithm: calculate influence scores for every mutation and every 
	#deregulated pathway
	for(i in 1:length(bins)){
		pathway_name = names(bins)[i]
		pathway_id = pathway_name
		seed = bins[[i]]
		pathway_network = get_pathway_network(pathway_list[[pathway_name]], original_network)
		if(is.null(nrow(pathway_network))){
		 next
		}
		#prize nodes are DEGs belong to the pathway
		pathway_prizes = diff_genes[bins[[pathway_name]]]
		max_score = sum(diff_genes[seed])
		res = matrix(unlist(mclapply(mutated_genes, function(x) 
			run_single_PCSF(pathway_network, x, original_network, pathway_prizes, alpha), mc.cores = num_of_cores)), ncol = 2, byrow = T)
		if(nrow(res) < length(mutated_genes)){ next }
		influence_scores = res[, 1]/max_score
		fraction_of_DEGs_values = res[, 2]/max_score
		Influence_matrix = rbind(Influence_matrix, influence_scores)
		rownames(Influence_matrix)[nrow(Influence_matrix )] = pathway_name 
		Fraction_of_DEGs_matrix = rbind(Fraction_of_DEGs_matrix , fraction_of_DEGs_values)
		rownames(Fraction_of_DEGs_matrix)[nrow(Fraction_of_DEGs_matrix )] = pathway_name 
		gc()
	}
	if(write_results){
		write.table(Influence_matrix[-1, ], file = paste(results_folder, sample, "_influence_scores.txt", sep = "")
					, quote = F, col.names = T, row.names = T, sep = "\t")
	}
	return(Influence_matrix[-1, ])
}

get_pathway_genes_from_network<-function(pathway)
{
 return(unique(c(pathway[, 1], pathway[, 2])))
}


get_network_from_gene_set<-function(network, gene_set)
{
 return(network[network[, 1] %in% gene_set & network[, 2] %in% gene_set, ])
}


get_pathway_list_from_graphite<-function(source = "reactome", minimal_number_of_nodes = 10, num_of_cores = 1)
{

	list_of_pathways <- graphite::pathways("hsapiens", source)
	num_of_nodes = lapply(list_of_pathways, function(x) length(nodes(x)))
	pathway_names = names(list_of_pathways)[unlist(num_of_nodes) > minimal_number_of_nodes]
	if(num_of_cores == 1)
	{
		pathway_list <- setNames(lapply(pathway_names, function(pathway_name) get_single_graphite_pathway(pathway_name, list_of_pathways)), pathway_names)	
	} else {
		pathway_list <- setNames(mclapply(pathway_names, function(pathway_name) get_single_graphite_pathway(pathway_name, list_of_pathways), mc.cores = num_of_cores), pathway_names)	
	}
	pathway_list <- pathway_list[sapply(pathway_list, length) > 0]
	return(pathway_list)
}

get_single_graphite_pathway<-function(pathway_name, list_of_pathways)
{
	pathway_edges = tryCatch({
		return(graphite::edges(convertIdentifiers(list_of_pathways[[which(names(list_of_pathways) == pathway_name)]], "symbol"))[, c("src", "dest")])
	}, 
	error = function(x){
		return(data.frame())
	}, 
	warning = function(x){
		return(data.frame())
	})
}


PRODIGY_cohort<-function(snv_matrix = NULL, expression_matrix, network = NULL, samples = NULL, DEGs = NULL, alpha = 0.05, pathway_list = NULL, 
			num_of_cores = 1, sample_origins = NULL, write_results = F, results_folder = "./", beta = 2, gamma = 0.05, delta = 0.05, mutation_list = NULL)
{

	if(is.null(samples))
	{
		print("no samples, aborting")
		return()
	}
	all_patients_scores = list()
	if(is.null(sample_origins))
	{
		sample_origins <- rep("tumor", ncol(expression_matrix))
		sample_origins[substr(colnames(expression_matrix), nchar(colnames(expression_matrix)[1])-1, nchar(colnames(expression_matrix)[1])) == "11"] = "normal"	
	}
	if(is.null(network))
	{
		data(STRING_network)
		network <- STRING_network
	}
	#run PRODIGY for all samples
	for(sample in samples)
	{
	 print(sample)
	 if(!is.null(DEGs) & (sample %in% names(DEGs))) 
	 {	 
			diff_genes <- DEGs[[sample]] 
	 } else { diff_genes = NULL }
			if(!is.null(mutation_list))
			{
				sample_mutations <- mutation_list[[sample]]
			} else {
				sample_mutations <- names(snv_matrix[snv_matrix[, sample] == 1, sample])
			}
			all_patients_scores[[sample]] <- PRODIGY(sample_mutations, expression_matrix, 
			  network, sample, diff_genes, alpha = alpha, pathway_list = pathway_list, 
	      num_of_cores = num_of_cores, sample_origins = sample_origins,
			  write_results = write_results, results_folder = results_folder, 
			  beta = beta, gamma = gamma, delta = delta)
	}
	return(all_patients_scores)
}

analyze_PRODIGY_results<-function(all_patients_scores)
{

	if(class(all_patients_scores) == "matrix"){
	  all_patients_scores = list(all_patients_scores)
	}
	for(i in 1:length(all_patients_scores)){	
		all_patients_scores[[i]][which(is.na(all_patients_scores[[i]]))] <- 0
		all_patients_scores[[i]][all_patients_scores[[i]] < 0] <- 0
	}
	Prodigy_rankings <- list()
	for(i in 1:length(all_patients_scores)) {
		if(is.null(all_patients_scores[[i]])){
		  Prodigy_rankings[[i]] <- c();
		  next
		}
		pathways_to_take <- c()
		#single pathway
		if(is.null(nrow(all_patients_scores[[i]])))
		{
			ranking <- sort(all_patients_scores[[i]], decreasing = T)
		#less than 4 mutations, no pathway filtering
		} else if(ncol(all_patients_scores[[i]]) < 4) {
			ranking <- sort(apply(all_patients_scores[[i]], 2, function(x) sum(sort(x, decreasing = T))), decreasing = T)	
		#filter pathways
		} else {
			pathways_to_take <- names(which(apply(all_patients_scores[[i]], 1, function(x) length(which(x > 0))) < ncol(all_patients_scores[[i]])/2))
			if(length(pathways_to_take) < 2)
			{
				ranking <- sort(all_patients_scores[[i]][pathways_to_take, ], decreasing = T)
				#if all pathways are filtered, abort pathway filtering
				if(all(ranking == 0))
				{
					pathways_to_take <- rownames(all_patients_scores[[i]])
					ranking <- sort(apply(all_patients_scores[[i]][pathways_to_take, ], 2, function(x) sum(sort(x, decreasing = T))), decreasing = T)
				}
			} else {
				ranking <- sort(apply(all_patients_scores[[i]][pathways_to_take, ], 2, function(x) sum(sort(x, decreasing = T))), decreasing = T)
				if(all(ranking == 0))
				{
					pathways_to_take <- rownames(all_patients_scores[[i]])
					ranking <- sort(apply(all_patients_scores[[i]][pathways_to_take, ], 2, function(x) sum(sort(x, decreasing = T))), decreasing = T)
				} else {
					pathways_to_take <- pathways_to_take[which(apply(all_patients_scores[[i]][pathways_to_take, ], 1, sum)!= 0)]
				}
			}
		}
		ranking = ranking[ranking > 0]
		if(length(ranking) == 0) {
		  Prodigy_rankings[[i]] = c();
		  next
		}
		# single pathway was used, no bimodel distribution 
		if(length(pathways_to_take) < 2){
		  Prodigy_rankings[[i]] = names(ranking);
		  next
		}
		# check bimodel distribution
		possible_error = tryCatch(
		{
				bimodel_dist <- normalmixEM(ranking, k = 2)
				unimodel_dist <- MASS::fitdistr(ranking, "normal")
		}, 
		error = function(cond) {
				cond
		})
		if(inherits(possible_error, "error")) {
		  Prodigy_rankings[[i]] = names(ranking)
		} else {
			# check which distribution is more likely
			if(bimodel_dist$loglik > unimodel_dist$loglik)
			{
					Prodigy_rankings[[i]] <- names(ranking[bimodel_dist$posterior[, which(bimodel_dist$mu == max(bimodel_dist$mu))] >0.5])
			} else {
					Prodigy_rankings[[i]] <- names(ranking)
			}
		}
	}
	names(Prodigy_rankings) <- names(all_patients_scores)
	return(Prodigy_rankings)
}

calculate_precision<-function(ranked_list, phenotype_genes)
{
	precision_vector <- c()
	for(i in 1:length(ranked_list))
	{
		precision_vector <- c(precision_vector, (length(intersect(ranked_list[1:i], phenotype_genes)))/
						i)
	}
	precision_vector[c(length(ranked_list):100)] <- precision_vector[length(ranked_list)]
	return(precision_vector)
}


calculate_recall <- function(ranked_list, phenotype_genes)
{
	recall_vector = c()
	for(i in 1:length(ranked_list))
	{
		recall_vector <- c(recall_vector, (length(intersect(ranked_list[1:i], phenotype_genes)))/
						length(phenotype_genes))
	}
	recall_vector[c(length(ranked_list):100)] <- recall_vector[length(ranked_list)]
	return(recall_vector)
}

check_performances <- function(ranked_genes_lists, snv_matrix, gold_standard_drivers, plot_figure = T, col = "red")
{

	precision_matrices <- list()
	recall_matrices <- list()
	f1_matrices <- list()	
	for(i in 1 : length(ranked_genes_lists))
	{
		patient_snp <- rownames(snv_matrix)[which(snv_matrix[, paste(names(ranked_genes_lists)[i], sep = "")] == 1)]
		if(length(intersect(gold_standard_drivers, patient_snp)) < 1) 
		{
			 print(paste("no known drivers with SNV mutations for patient", names(ranked_genes_lists)[i]))
			 next
		}
		if(length(intersect(gold_standard_drivers, ranked_genes_lists[[i]])) > 0)
		{
			precision_matrices[[names(ranked_genes_lists[i])]] <- calculate_precision(ranked_genes_lists[[i]][1:min(20, length(ranked_genes_lists[[i]]))], gold_standard_drivers)[1:20]
			recall_matrices[[names(ranked_genes_lists[i])]] <- calculate_recall(ranked_genes_lists[[i]][1:min(20, length(ranked_genes_lists[[i]]))], intersect(gold_standard_drivers, patient_snp))[1:20]
			curr_f1 = (2*precision_matrices[[names(ranked_genes_lists[i])]] * recall_matrices[[names(ranked_genes_lists[i])]])/(precision_matrices[[names(ranked_genes_lists[i])]]+recall_matrices[[names(ranked_genes_lists[i])]])
			curr_f1[is.nan(curr_f1)] <- 0
			f1_matrices[[names(ranked_genes_lists[i])]] <- curr_f1		
		}
	}
	precision_matrix <- do.call(rbind, precision_matrices) 
	recall_matrix <- do.call(rbind, recall_matrices)
	f1_matrix <- do.call(rbind, f1_matrices)		
	colnames(precision_matrix) <- 1:ncol(precision_matrix)
	colnames(recall_matrix) <- 1:ncol(recall_matrix)
	colnames(f1_matrix) <- 1:ncol(f1_matrix)
	df_means = data.frame(x = seq(1, 20, 1))
	if(class(precision_matrix)[1] == "matrix")
	{
		curr_precisions <- apply(precision_matrix, 2, function(x) mean(x))
	} else {
	  curr_precisions <- precision_matrix
	 }
	df_means = cbind(df_means, curr_precisions, rep("PRODIGY", 20))
	names(df_means)[c(2, 3)] <- c("value", "Algorithm")
	precision_plot <- ggplot(data = df_means, aes(x = x, y = value, colour = Algorithm))+ scale_colour_manual(values = col)+ geom_point(size = 1) + geom_line(size = 1) + labs(x = "", y = "Average precision") + ylim(c(0, ceiling(max(df_means[, "value"])*10)/10)) + theme(axis.title.y = element_text(size = 9), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 9))
	
	df_means <- data.frame(x = seq(1, 20, 1))
	if(class(recall_matrix)[1] == "matrix")
	{
		curr_recalls <- apply(recall_matrix, 2, function(x) mean(x))
	} else {
	  curr_recalls <- recall_matrix
	 }
	df_means <- cbind(df_means, curr_recalls, rep("PRODIGY", 20))
	names(df_means)[c(2, 3)] <- c("value", "Algorithm")
	recall_plot <- ggplot(data = df_means, aes(x = x, y = value, colour = Algorithm))+ scale_colour_manual(values = col)+ geom_point(size = 1) + geom_line(size = 1) + labs(x = "", y = "Average recall") + ylim(c(0, ceiling(max(df_means[, "value"])*10)/10)) + theme(axis.title.y = element_text(size = 9), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 9))


	df_means <- data.frame(x = seq(1, 20, 1))
	if(class(f1_matrix)[1] == "matrix"){
		curr_f1 <- apply(f1_matrix, 2, function(x) mean(x))
	} else {
	  curr_f1 <- f1_matrix
	 }
	df_means <- cbind(df_means, curr_f1, rep("PRODIGY", 20))
	names(df_means)[c(2, 3)] <- c("value", "Algorithm")
	f1_plot <- ggplot(data = df_means, aes(x = x, y = value, colour = Algorithm))+ scale_colour_manual(values = col)+ geom_point(size = 1) + geom_line(size = 1) + labs(x = "", y = "Average F1") + ylim(c(0, ceiling(max(df_means[, "value"])*10)/10)) + theme(axis.title.y = element_text(size = 9), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 9))
	if(plot_figure){
		plot_total <- plot_grid(precision_plot, recall_plot, f1_plot, NULL, labels = c("A", "B", "C"), align = "h", nrow = 2, label_size = 12, ncol = 2)
		save_plot("基于金标准癌基因的性能图.png", plot_total)
	}
	return(list("precision" = precision_matrix, "recall" = recall_matrix, "f1" = f1_matrix))
}

install_pkgs <- function(){
if(!require("MASS")) install.packages("MASS")
if(!require("igraph")) install.packages("igraph")
if(!require("ff")) install.packages("ff")
if(!require("plyr")) install.packages("plyr")
if(!require("mixtools")) install.packages("mixtools")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("cowplot")) install.packages("cowplot")
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("DESeq2")) BiocManager::install("DESeq2")
if(!require("graphite")) BiocManager::install("graphite")
if(!require("biomaRt")) BiocManager::install("biomaRt")

if(!require("devtools")) install.packages("devtools")
if(!require("PCSF")) devtools::install_github("IOR-Bioinformatics/PCSF", repos = BiocManager::repositories(), 
    dependencies = TRUE, type = "source", force = TRUE)
}
main <- function(){
 install_pkgs()
 
 Prefix <- "D:/Documents/论文/PRODIGY/"
 setwd(Prefix)
 # Load SNP+expression data derived from TCGA
 load(paste0(Prefix, "COAD_SNV.RData"))
 cat(paste0("加载了", dim(snv_matrix)[1],"基因和",dim(snv_matrix)[2],"个样本","的体细胞突变数据。\n"), file = stdout())
 load(paste0(Prefix, "COAD_Expression.RData"))
 cat(paste0("加载了", dim(expression_matrix)[1],"基因和",dim(expression_matrix)[2],"个样本","的基因表达数据。\n"), file = stdout())
 # Load STRING network data 
 load(paste0(Prefix, "STRING_network.RData"))
 cat(paste0("加载了", dim(STRING_network)[1],"个蛋白质相互作用网络。\n"), file = stdout())
 load(paste0(Prefix, "gold_standard_drivers.RData"))
 cat("加载了249个金标准癌基因。\n", file = stdout())
 network <- STRING_network
 # Take samples for which SNP and expression is available 
 samples <- intersect(colnames(expression_matrix), colnames(snv_matrix))[1:5]
 # Get differentially expressed genes (DEGs) for all samples
 expression_matrix <- expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[, 1], network[, 2]))), ]
 DEGs <- get_DEGs(expression_matrix, samples, sample_origins = NULL, beta = 2, gamma = 0.05)
 # Identify sample origins (tumor or normal)
 sample_origins <- rep("tumor", ncol(expression_matrix))
 sample_origins[substr(colnames(expression_matrix), nchar(colnames(expression_matrix)[1])-1, nchar(colnames(expression_matrix)[1])) == "11"] = "normal"	
 list_of_pathways <- get_pathway_list_from_graphite(source = "reactome", minimal_number_of_nodes = 1, num_of_cores = 1)

 cat("开始运行PRODIGY。\n", file = stdout())
 all_patients_scores <- PRODIGY_cohort(snv_matrix = snv_matrix, expression_matrix = expression_matrix, network = network, samples = samples, DEGs = DEGs, alpha = 0.05, 
 			pathway_list = list_of_pathways, num_of_cores = 1, sample_origins = sample_origins, write_results = F, results_folder = "./", 
 			beta = 2, gamma = 0.05, delta = 0.05)
 cat("分析PRODIGY结果，得到每个样本的癌基因排序。\n", file = stdout())
 results <- analyze_PRODIGY_results(all_patients_scores) 
 result2 <- NULL
 for(i in 1:length(results)){
   for(j in 1:length(results[[i]])){
     result2 <- rbind(result2, data.frame("样本" = names(results)[[i]], "基因" = results[[i]][j], "样本内排序" = j))
   }
 }
 write.table(result2, file = paste0(Prefix, "样本内基因排序.txt"), quote = F, sep = "\t", row.names = FALSE) #Export
 rm(result2)
 perf <- check_performances(ranked_genes_lists = results, snv_matrix = snv_matrix, gold_standard_drivers)

 perf2 <- data.frame("指标" = "精度", "索引" = colnames(perf[[1]]) , "值" = t(perf[[1]]))
 perf2 <- rbind(perf2, data.frame("指标" = "召回率", "索引" = colnames(perf[[2]]) ,"值" = t(perf[[2]])))
 perf2 <- rbind(perf2, data.frame("指标" = "F1", "索引" = colnames(perf[[3]]) , "值" = t(perf[[3]])))
 write.table(perf2, file = paste0(Prefix, "性能.txt"), quote = F, sep = "\t", row.names = FALSE) #Export
 rm(perf2)
 cat("程序运行结束。\n", file = stdout())
}

main()