###########################################
#                                         #
#          Load needed libraries          #
#                                         #
###########################################

library("taRifx")
library("matrixStats")
library("ggplot2")
library("Rtsne")
library("fpc")
library("factoextra")
library("monocle")
library("viridis")
library("gplots")
library("RColorBrewer")
library("destiny")
library("slingshot")
library("rgl")
library("scatterplot3d")
library("made4")
library("pheatmap")
library("matrixStats")
library("statmod")

library("FactoMineR")
library("jackstraw")

library("ReactomePA")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GOSemSim")
library("arulesViz")

library("ggpubr")


# Window size for 3D DM graphs
r3dDefaults$windowRect <- c(0,50, 700, 700)

###########################################
#                                         #
#            Var Gene Selection           #
#                                         #
###########################################

# Code taken and adapted from doi:10.1038.nmeth.2645
# Select the highly variable genes based on the squared coefficient of variation and the mean gene expression and return the RPKM matrix the the HVG

getMostVarGenes <- function(
		data=data,				# RPKM matrix
		fitThr=1.5, 			# Threshold above the fit to select the HGV
		minMeanForFit=1			# Minimum mean gene expression level
		){
	# Remove genes expressed in no cells
	data_no0 <- as.matrix(
		data[rowSums(data)>0,]
	)
	# Compute the mean expression of each genes
	meanGeneExp <- rowMeans(data_no0)
	names(meanGeneExp)<- rownames(data_no0)

	# Compute the squared coefficient of variation
	varGenes <- rowVars(data_no0)
	cv2 <- varGenes / meanGeneExp^2
 
 	# Select the genes which the mean expression is above the expression threshold minMeanForFit
	useForFit <- meanGeneExp >= minMeanForFit

	# Compute the model of the CV2 as a function of the mean expression using GLMGAM
	fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meanGeneExp[useForFit] ), cv2[useForFit] )
	a0 <- unname( fit$coefficients["a0"] )
	a1 <- unname( fit$coefficients["a1tilde"])
 
 	# Get the highly variable gene counts and names
	fit_genes <- names(meanGeneExp[useForFit])
	cv2_fit_genes <- cv2[useForFit]
	fitModel <- fit$fitted.values
	names(fitModel) <- fit_genes
	HVGenes <- fitModel[cv2_fit_genes>fitModel*fitThr]
	print(length(HVGenes))

	# Plot the result
	plot_meanGeneExp <- log10(meanGeneExp)
	plot_cv2 <- log10(cv2)
	plotData <-  data.frame(
		x=plot_meanGeneExp[useForFit],
		y=plot_cv2[useForFit],
		fit=log10(fit$fitted.values),
		HVGenes=log10((fit$fitted.values*fitThr))
	)
	p <- ggplot(plotData, aes(x,y)) +
	geom_point(size=0.1) +
	geom_line(aes(y=fit), color="red") +
	geom_line(aes(y=HVGenes), color="blue") +
	theme_bw() +
	labs(x = "Mean expression (log10)", y="CV2 (log10)")+
	ggtitle(paste(length(HVGenes), " selected genes", sep="")) +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title = element_text(size =16 ,face="bold"),
		legend.position= "none",
		plot.title = element_text(size=18, face="bold", hjust = 0.5),
		aspect.ratio=1
	)+
	scale_color_manual(
		values=c("#595959","#5a9ca9")
	)
	print(p)

	# Return the RPKM matrix containing only the HVG
	HVG <- data_no0[rownames(data_no0) %in% names(HVGenes),]
	return(HVG)
}


###########################################
#                                         #
#              RtSNE Analysis             #
#                                         #
###########################################

# Compute the t-SNE
run_tSNE <- function(
	pca=pca, 					# PCA computed with FactoMineR package
	pc=20, 						# Number of PCs to use in the t-SNE, default is 20
	iter=2000, 					# Number of iteration for th rt-SNE, default is 2000
	perplexity=30				# Set the t-SNE perplexity parameter if needed, default is 30
	){
	set.seed(1)
	# Run t-SNE
	rtsne_out <- Rtsne(
		pca$ind$coord[,pc], 
		pca=FALSE, 
		max_iter=iter, 
		verbose=TRUE
	)
	tSNE <- data.frame(
		rtsne_out$Y
	)
	return(tSNE)
}

# Compute and plot the t-SNE
run_plot_tSNE <- function(
	pca=pca, 
	pc=pc, 
	iter=2000,
	perplexity=30, 
	conditions=conditions, 
	colours=colours, 
	title=FALSE
	){
	tsne <- run_tSNE(
		pca, 
		1:pc, 
		iter,  
		perplexity
	)
	tSNE <- data.frame(
		tsne,
		conditions
	)
	colnames(tSNE)<- c("tSNE_1", "tSNE_2", "cond")
	if (title==TRUE){
		g<-ggplot(tSNE, aes(tSNE_1,tSNE_2)) +
		geom_point(shape = 21, size = 2.5, stroke=0.5, aes(fill=cond, color=cond),alpha=6/10) +
		theme_bw() +
		scale_fill_manual(
			values=colours,
			name=""
			) +
			scale_color_manual(
				values=colours,
				name=""
			) +
			xlab("t-SNE 1") +
			ylab("t-SNE 2") +
		ggtitle(paste("t-SNE plot (", pc," PCs)", sep="")) +
		theme(
			axis.text=element_text(size=16),
			axis.title=element_text(size=16),
			legend.text = element_text(size =16),
			legend.title = element_text(size =16 ,face="bold"),
			plot.title = element_text(size=18, face="bold", hjust = 0.5),
			aspect.ratio=1#,
			#legend.position="none"
		)
	} else {
		g<-ggplot(tSNE, aes(tSNE_1,tSNE_2)) +
		geom_point(shape = 21, size = 2.5, stroke=0.5, aes(fill=cond, color=cond),alpha=6/10) +
		theme_bw() +
		scale_fill_manual(
			values=colours,
			name=""
			) +
			scale_color_manual(
				values=colours,
				name=""
			) +
		xlab("t-SNE 1") +
		ylab("t-SNE 2") +
		theme(
			axis.text=element_text(size=16),
			axis.title=element_text(size=16),
			legend.text = element_text(size =16),
			legend.title = element_text(size =16 ,face="bold"),
			plot.title = element_text(size=18, face="bold", hjust = 0.5),
			aspect.ratio=1#,
			#legend.position="none"
		)
	}
	print(g)
	return(tSNE)
}


# Plot the pre-computed t-SNE
plot_tSNE <- function(tsne=tsne, conditions=conditions, colours=colours){

	tSNE <- data.frame(
		tSNE_1=tsne[,1],
		tSNE_2=tsne[,2],
		c(conditions)
	)

	colnames(tSNE)<- c("tSNE_1", "tSNE_2", "cond")

	print(head(tSNE))


	g<-ggplot(tSNE, aes(tSNE_1, tSNE_2)) +
	geom_point(shape = 21, size = 2, stroke=0.5, aes(fill=cond, color=cond),alpha=6/10) +
	theme_bw() +
	scale_fill_manual(
		values=colours,
		name=""
	) +
	scale_color_manual(
		values=colours,
		name=""
	) +
	# ggtitle(paste("t-SNE plot (", pc," PCs)", sep="")) +
	xlab("t-SNE 1") +
	ylab("t-SNE 2") +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title = element_text(size =16 ,face="bold"),
		plot.title = element_text(size=18, face="bold", hjust = 0.5),
		aspect.ratio=1#,
		#legend.position="none"
	)
	print(g)
	return(tSNE)
}


# t-SNE colored by gene expression level
tsne_gene_per_cell <- function(tsne_result, rpkm){
	gene_per_cell <- rpkm
	gene_per_cell[gene_per_cell>0] <- 1
	gene_per_cell <- colSums(gene_per_cell)

	colnames(tsne_result)<- c("tSNE_1", "tSNE_2")
	warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
	mypalette <- c(warm(20))

	cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494'))
	warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c'))
	mypalette <- c(rev(cold(11)), warm(10))


	p <- ggplot(tsne_result, aes(tSNE_1, tSNE_2)) +
	geom_point(shape = 21, stroke=0.25, aes(colour=as.numeric(gene_per_cell), fill=as.numeric(gene_per_cell)), size = 2) +
	# scale_fill_gradient2(high="darkred", low="yellow")+
	# scale_fill_gradientn(colours = mypalette)+
	# scale_colour_gradientn(colours = mypalette)+
	scale_fill_viridis_c()+
	scale_colour_viridis_c()+
	theme_bw() +
	# ggtitle(gene) +
	xlab("t-SNE 1") +
	ylab("t-SNE 2") +
	theme(
		plot.title = element_text(size=28, face="bold.italic", hjust = 0),
		# axis.text=element_text(size=12),
		# axis.title=element_text(size=16),
		axis.title=element_blank(),
		axis.text=element_blank(),
		legend.text = element_text(size =12),
		legend.title=element_blank(),
		aspect.ratio=1,
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
	)
	# print(p)

}



# t-SNE colored by gene expression level
tsne_gene_exp <- function(tsne_result, gene, rpkm){
	colnames(tsne_result)<- c("tSNE_1", "tSNE_2")
	warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
	mypalette <- c(warm(20))

	cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494'))
	warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c'))
	mypalette <- c(rev(cold(11)), warm(10))

	gene_exp <- rpkm[gene,,drop=FALSE]
	gene_exp <- log(gene_exp+1)
	density_exp <- tsne_result[rownames(tsne_result) %in% names(gene_exp[,gene_exp>rowMeans(gene_exp),drop=FALSE]),]
	# print(density_exp)

	title <- grobTree(textGrob(gene, x=0.05, y=0.93, , hjust = 0, gp=gpar(fontsize=32, fontface="bold.italic")))

	p <- ggplot(tsne_result, aes(tSNE_1, tSNE_2)) +
	geom_point(shape = 21, stroke=0.25, aes(colour=as.numeric(log(rpkm[gene,]+1)), fill=as.numeric(log(rpkm[gene,]+1))), size = 2) +
	geom_density_2d(data=density_exp, aes(x=tSNE_1, y=tSNE_2), bins = 5, colour="black") +
	# scale_fill_gradient2(high="darkred", low="yellow")+
	scale_fill_gradientn(colours = mypalette)+
	scale_colour_gradientn(colours = mypalette)+
	# scale_fill_viridis_c()+
	theme_bw() +
	# ggtitle(gene) +
	annotation_custom(title) +
	xlab("t-SNE 1") +
	ylab("t-SNE 2") +
	theme(
		plot.title = element_text(size=28, face="bold.italic", hjust = 0),
		# axis.text=element_text(size=12),
		# axis.title=element_text(size=16),
		axis.title=element_blank(),
		axis.text=element_blank(),
		legend.text = element_text(size =12),
		legend.title=element_blank(),
		aspect.ratio=1,
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		legend.position="none",
		axis.ticks = element_blank()
	)
	# print(p)

}

violin_gene_exp <- function(gene, rpkm, conditions=conditions, colours=colours, test=TRUE){
	exp <- as.numeric(log(rpkm[gene,]+1))

	gene_exp <- data.frame(
		cell=colnames(rpkm),
		clusters=conditions,
		gene=exp
		)

	if (test==TRUE){
		# Perform pairwise comparisons
		test <- compare_means(gene ~ clusters, data = gene_exp, method = "kruskal.test")
		# print(test)
		# print(compare_means(gene ~ clusters, data = gene_exp, , method = "wilcox.test", ref.group = ".all."))

		p <- ggboxplot(gene_exp, x = "clusters", y = "gene", color = "white")+
		geom_violin(scale = "width", width=0.7, adjust = .5,aes(fill=clusters)) +
		stat_summary(fun.y=mean, geom="point", shape=21, size=3,  stroke = 1, fill="white")+

		# geom_jitter(size=0.3)+  	  
		geom_hline(yintercept = mean(gene_exp$gene), linetype = 2)+

		  # geom_boxplot(fill="white", outlier.shape = NA, width = 0.2)+
		scale_fill_manual(
			values=colours
		) +	  
		  theme_bw()+
		  ggtitle(gene)+
		  expand_limits(y = c(0, max(gene_exp$gene)+1.5)) +
		# stat_compare_means(method = "kruskal.test", label.y = max(gene_exp$gene)+1)+      # Add global p-value
		stat_compare_means(
		  label = "p.signif", 
		  method = "wilcox.test", 
		  ref.group = ".all.", 
		  label.y = max(gene_exp$gene)+0.75, 
		  size=6
		)+ # Pairwise comparison against all
		theme(
			plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
			axis.text=element_text(size=16),
			# axis.title=element_text(size=16),
			axis.title=element_blank(),
			legend.text = element_text(size =16),
			legend.title=element_blank(),
			aspect.ratio=0.5,
			legend.position="none",
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
		)

	}else{
		p <- ggplot(gene_exp, aes(clusters, gene))+
		geom_violin(scale = "width", width=0.7, adjust = .5,aes(fill=clusters)) +
		stat_summary(fun.y=mean, geom="point", shape=21, size=3,  stroke = 1, fill="white")+
		# geom_jitter(size=0.3)+ 
		scale_fill_manual(
			values=colours
		) +	  
		  theme_bw()+
		  ggtitle(gene)+
		theme(
			plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
			axis.text=element_text(size=16),
			# axis.title=element_text(size=16),
			axis.title=element_blank(),
			legend.text = element_text(size =16),
			legend.title=element_blank(),
			aspect.ratio=0.5,
			legend.position="none"
			# panel.grid.major = element_blank(),
			# panel.grid.minor = element_blank()
		)
	}


	print(p)
}


violin_gene_exp_sex <- function(gene, rpkm, conditions=conditions, colours=colours, test=TRUE){
	exp <- as.numeric(log(rpkm[gene,]+1))

	gene_exp <- data.frame(
		cell=colnames(rpkm),
		clusters=conditions,
		gene=exp,
		sex= sapply(strsplit(colnames(rpkm), "_"), `[`, 2)
		)

	if (test==TRUE){
		# Perform pairwise comparisons
		test <- compare_means(gene ~ clusters, data = gene_exp, method = "kruskal.test")
		# print(test)
		# print(compare_means(gene ~ clusters, data = gene_exp, , method = "wilcox.test", ref.group = ".all."))

		p <- ggboxplot(gene_exp, x = "clusters", y = "gene", color = "white")+
		geom_violin(scale = "width", width=0.7, adjust = .5,aes(fill=clusters)) +
		stat_summary(fun.y=mean, geom="point", shape=21, size=3,  stroke = 1, fill="white")+

		# geom_jitter(size=0.3)+  	  
		geom_hline(yintercept = mean(gene_exp$gene), linetype = 2)+

		  # geom_boxplot(fill="white", outlier.shape = NA, width = 0.2)+
		scale_fill_manual(
			values=colours
		) +	  
		  theme_bw()+
		  ggtitle(gene)+
		  expand_limits(y = c(0, max(gene_exp$gene)+1.5)) +
		# stat_compare_means(method = "kruskal.test", label.y = max(gene_exp$gene)+1)+      # Add global p-value
		stat_compare_means(
		  label = "p.signif", 
		  method = "wilcox.test", 
		  ref.group = ".all.", 
		  label.y = max(gene_exp$gene)+0.75, 
		  size=6
		)+ # Pairwise comparison against all
		theme(
			plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
			axis.text=element_text(size=16),
			# axis.title=element_text(size=16),
			axis.title=element_blank(),
			legend.text = element_text(size =16),
			legend.title=element_blank(),
			aspect.ratio=0.5,
			legend.position="none",
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
		)

	}else{

		p <- ggplot(gene_exp, aes(x=clusters, y=gene, fill=sex, dodge=clusters))+
		geom_violin(scale = "width", width=0.7, aes(fill=sex, color=sex), position=position_dodge(0.8), adjust = .5)+
		# geom_jitter()+
		# geom_boxplot(fill="white", outlier.shape = NA, width = 0.15, aes(fill=Sex, color=Sex), position=position_dodge(0.8))+
		geom_point(position=position_jitterdodge(dodge.width=0.8), size=0.3, aes(color=sex))+
		stat_summary(fun.y=median, geom="point", shape=21, position=position_dodge(0.8), size=3,  stroke = 1, aes(color=sex), fill="white")+
		scale_color_manual(
			values=rev(c("#5676f0", "#de62b7"))
		) +
		scale_fill_manual(
			values=rev(c("#85b1f4", "#f4aede"))
		) +
		ggtitle(gene)+
		ylab("Expression level")+
		ylim(c(0,4.5))+
		xlab("Embryonic stages")+
		theme_bw()+
		# coord_flip()+
		# scale_x_discrete(limits=rev(c("E10.5", "E11.5", "E12.5", "E13.5", "E16.5")))+
		# facet_grid(. ~ Type)+
		scale_x_discrete(limits=c("E10.5", "E11.5", "E12.5", "E13.5", "E16.5", "P6"))+
		# facet_grid(Type ~ .)+
		theme(
			axis.text=element_text(size=16),
			axis.title=element_text(size=16),
			legend.text = element_text(size =16),
			legend.title=element_text(size =16, face="bold"),
			plot.title = element_text(size=22, face="bold.italic", hjust = 0.5),
			aspect.ratio=0.5,
			legend.position="bottom"
			# strip.text.x = element_text(size = 16)
			# strip.text.y = element_text(size = 16)
		)

		# print(p)

	}


	print(p)
}


###########################################
#                                         #
#                 Plot PCA                #
#                                         #
###########################################

# Plot pre-computed PCA (with prcomp)
plot_pca <- function(pca=pca, pc=pc, conditions=conditions, colours=colours, title=""){
		# PCs <- 1:pc
		PCs <- paste("PC",1:pc, sep="")
		percent_var_explained <- (pca$sdev^2 / sum(pca$sdev^2))*100
		cond <- factor(conditions)
		col <- factor(conditions)
		levels(col) <- colours
		col <- as.vector(col)

		scores <- as.data.frame(pca$x)
		PCs.combinations <- combn(PCs,2)
		g <- apply(
			PCs.combinations,
			2,
			function(combination)
			{
				p1 <- ggplot(scores, aes_string(x=combination[1], y=combination[2])) +
				geom_point(shape = 21, size = 2, stroke=0.5, aes(fill=cond, colour=cond), alpha=6/10) +
				theme_bw() +
				scale_fill_manual(
					values=colours,
					name=""
				) +
				scale_color_manual(
					values=colours,
					name=""
				) +
				xlab(paste(combination[1], " ", "(",round(percent_var_explained[as.numeric(gsub("PC", "", combination[1]))], digit=3),"%)", sep=""))+
				ylab(paste(combination[2], " ", "(",round(percent_var_explained[as.numeric(gsub("PC", "", combination[2]))], digit=3),"%)", sep=""))+
				ggtitle(title)+
				theme(
					axis.text=element_text(size=16),
					axis.title=element_text(size=16),
					legend.text = element_text(size =16),
					legend.title = element_text(size =16 ,face="bold"),
					plot.title = element_text(size=18, face="bold", hjust = 0.5),
					aspect.ratio=1,
					legend.position="none"
				)
				print(p1)
			}
		)
}

# Plot pre-computed PCA (with prcomp) in 3D
plot_pca_3D <- function(pca=pca, pc=c(1:3), conditions=conditions, colours=colours){
	cond <- factor(conditions)
	col <- factor(conditions)
	levels(col) <- colours
	col <- as.vector(col)

	PCs <- paste("PC",pc, sep="")


	data <- data.frame(
		pca$x[,PCs[1]], 
		pca$x[,PCs[2]], 
		pca$x[,PCs[3]]
	)

	colnames(data) <- PCs

	plot3d(
		data,
		col=col,
		size=6.5
	)

	plot3d(
		data,
		size=8,
		add = TRUE
	)

}



###########################################
#                                         #
#             Plot screeplot              #
#                                         #
###########################################


# Screeplot of a pre-computed PCA (with prcomp)
plot_percent_var <- function(pca=pca, pc=pc){
	percent_var_explained <- (pca$sdev^2 / sum(pca$sdev^2))*100
	percent_var_explained <- data.frame(
		PC=1:length(percent_var_explained),
		percent_Var=percent_var_explained
	)

	sub_percent_var_explained <- percent_var_explained[1:pc,]

	p <- ggplot(sub_percent_var_explained, aes(x=PC, y=percent_Var)) + 
		geom_col()+
		theme_bw() +
		xlab("PCs") +
		ylab("% Variance") +
		ggtitle("Screeplot")+
		theme(
			axis.text=element_text(size=16),
			axis.title=element_text(size=16),
			legend.text = element_text(size =16),
			legend.title = element_text(size =16 ,face="bold"),
			plot.title = element_text(size=18, face="bold", hjust = 0.5),
			aspect.ratio=0.5
		)

	print(p)
}


###########################################
#                                         #
#        Clustering with DBSCAN           #
#                                         #
###########################################

# Cell clustering from t-SNE with DBSCAN
# Obsolete
plot_clusters <- function(tsne=tsne, dist=3.5){
	set.seed(123)
	db <-dbscan(tsne, dist, MinPts = 2)
	g <- fviz_cluster(
		db, 
		tsne, 
		stand = FALSE, 
		geom = "point", 
		ellipse=FALSE, 
		ggtheme=theme_bw(), 
		xlab="tSNE_1", 
		ylab="tSNE_2", 
		pointsize=2.5, 
		shape=19, 
		show.clust.cent=FALSE) +
		ggtitle(paste("DBSCAN clustering (dist=", dist, ")", sep="")) +
		theme(
			axis.text=element_text(size=16),
			axis.title=element_text(size=16),
			legend.text = element_text(size =16),
			legend.title = element_text(size =16 ,face="bold"),
			plot.title = element_text(size=18, face="bold", hjust = 0.5),
			aspect.ratio=1
		)

	print(g)
	return(db$cluster)
}

###########################################
#                                         #
#                   DE                    #
#                                         #
###########################################

# Create the Monocle object with normalized read counts (size-factor)
prepare_for_monocle <- function(count_matrix=count_matrix, stages=stages){
	# Prepare tables for monocle object
	expr_matrix <- as.matrix(count_matrix)
	sample_sheet <- data.frame(cells=names(count_matrix), stages=stages)
	rownames(sample_sheet)<- names(count_matrix)
	gene_annotation <- as.data.frame(rownames(count_matrix))
	rownames(gene_annotation)<- rownames(count_matrix)
	colnames(gene_annotation)<- "genes"
	pd <- new("AnnotatedDataFrame", data = sample_sheet)
	fd <- new("AnnotatedDataFrame", data = gene_annotation)

	# Create a CellDataSet from the relative expression levels
	HSMM <- newCellDataSet(
		as(expr_matrix, "sparseMatrix"),
		phenoData = pd,
		featureData = fd,
		lowerDetectionLimit=0.5,
		expressionFamily=negbinomial.size()
	)
	
	HSMM <- detectGenes(HSMM, min_expr = 5)
	# HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 5, ]
	HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 10, ]

	HSMM <- estimateSizeFactors(HSMM)
	HSMM <- estimateDispersions(HSMM)

	return(HSMM)
}

# Same funcion as previous? (I don't remmember...)
prepare_for_DE <- function(count_matrix=count_matrix, clustering=clustering, stages=stages){
	# Prepare tables for monocle object
	expr_matrix <- as.matrix(count_matrix)
	sample_sheet <- data.frame(cells=names(count_matrix), stages=stages, cellType=clustering)
	rownames(sample_sheet)<- names(count_matrix)
	gene_annotation <- as.data.frame(rownames(count_matrix))
	rownames(gene_annotation)<- rownames(count_matrix)
	colnames(gene_annotation)<- "genes"
	pd <- new("AnnotatedDataFrame", data = sample_sheet)
	fd <- new("AnnotatedDataFrame", data = gene_annotation)

	# Create a CellDataSet from the relative expression levels
	HSMM <- newCellDataSet(
		as(expr_matrix, "sparseMatrix"),
		phenoData = pd,
		featureData = fd,
		lowerDetectionLimit=0.5,
		expressionFamily=negbinomial.size()
	)
	
	HSMM <- detectGenes(HSMM, min_expr = 5)
	# HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 5, ]
	HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 10, ]

	HSMM <- estimateSizeFactors(HSMM)
	HSMM <- estimateDispersions(HSMM)

	return(HSMM)
}

# Compute DE genes
findDEgenes <- function(HSMM=HSMM, qvalue=qvalue){
	diff_test_res <- differentialGeneTest(
		HSMM,
		fullModelFormulaStr="~cellType",
		cores = 3
	)

	sig_genes_0.05 <- subset(diff_test_res, qval < 0.05)
	sig_genes_0.01 <- subset(diff_test_res, qval < 0.01)

	print(paste(nrow(sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
	print(paste(nrow(sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))

	diff_test_res <- subset(diff_test_res, qval< qvalue)

	return(diff_test_res)
}

plot_DE_heatmap <- function(rpkm_matrix=rpkm_matrix, diff_test_res=diff_test_res, qvalue=qvalue, condition=condition) {
	sig_genes <- subset(diff_test_res, qval < qvalue)
	dim(sig_genes)
	de_genes <- as.matrix(log(rpkm_matrix[rownames(rpkm_matrix) %in% rownames(sig_genes),]+1))
	plot_heatmap(de_genes, condition)
}



get_up_reg_clusters <- function(count, clustering, DE_genes){
	cluster_nb <- unique(clustering)
	mean_per_cluster <- vector()
	DE_genes <- DE_genes[order(rownames(DE_genes)),]
	count <- count[order(rownames(count)),]
	count_de_genes <- count[rownames(count) %in% DE_genes$genes,]
	print(dim(count_de_genes))
	for (clusters in cluster_nb) {
		# print(head(count_de_genes[,
		# 		colnames(count_de_genes) %in% names(clustering[clustering==clusters])
		# 	]))
		mean <- rowMeans(
			as.matrix(count_de_genes[,
				colnames(count_de_genes) %in% names(clustering[clustering==clusters])
			])
		)
		names(mean) <- clusters
		mean_per_cluster <- cbind(
			mean_per_cluster,
			mean
		)
	}
	colnames(mean_per_cluster) <- cluster_nb
	up_reg_cluster <- colnames(mean_per_cluster)[apply(mean_per_cluster,1,which.max)]
	de_genes_table <- data.frame(
		DE_genes,
		mean_per_cluster,
		cluster=up_reg_cluster
	)

	return(de_genes_table)
}

get_top_up_reg_clusters <- function(DE_genes, gene_nb){
	cluster_nb <- unique(DE_genes$cluster)
	cluster_nb <- cluster_nb[order(cluster_nb)]
	DE_genes <- DE_genes[order(DE_genes$qval),]
	DE_genes <- DE_genes[DE_genes$qval<1.10e-2,]

	top_DE_genes <- vector()
	for (clusters in cluster_nb) {
		genes <- rownames(DE_genes[DE_genes$cluster==clusters,])[1:gene_nb]
		top_DE_genes <- c(top_DE_genes, genes)
	}

	return(top_DE_genes)
}

###########################################
#                                         #
#              Draw heatmaps              #
#                                         #
###########################################


plot_heatmap <- function(matrix=matrix, clusters=clusters, stages=stages){
	annotation_col <- data.frame(
		Cell_Clusters=clusters,
		Stages=stages
	)

	cell_cluster_colors <- c("#bd66d4", "#77b0f3", "#8dcf38", "#ffd65c", "#fb7072")
	names(cell_cluster_colors) <- unique(clusters)

	annotation_colors <- list(
		Stages=c(P6="#b0f137", E16.5="#fc8d59", E13.5="#fee090", E12.5="#e0f3f8", E11.5="#91bfdb", E10.5="#4575b4"),
		Cell_Clusters=cell_cluster_colors
	)

	pheatmap(
		matrix,  
		show_colnames=FALSE, 
		show_rownames=FALSE, 
		clustering_method="ward.D2",
		annotation_col=annotation_col,
		annotation_colors=annotation_colors,
		color=viridis(100)
	)

}


plot_heatmap_2 <- function(matrix=matrix, clusters=clusters, stages=stages, rowbreaks, colbreaks, cluster_color, stage_color){
	annotation_col <- data.frame(
		Cell_Clusters=clusters,
		Stages=stages
	)

	annotation_colors <- list(
		Stages=stage_color,
		Cell_Clusters=cluster_color
	)

	# Color palette for the heatmap
	cold <- colorRampPalette(c('#41b6c4','#253494','#081d58', '#081d58', '#081d58', '#081d58', '#081d58', '#081d58', '#081d58'))
	warm <- colorRampPalette(c('#fecc5c','#e31a1c','#800026','#800026','#800026','#800026','#800026','#800026'))
	mypalette <- c(rev(cold(20)), warm(20))
	# mypalette <- c(rev(cold(15)), warm(16))
	breaksList = seq(0, 5, by = 0.5)


	pheatmap(
		matrix, 
		# scale="row",
		show_colnames=FALSE, 
		# show_rownames=FALSE, 
		cluster_cols=FALSE,
		# cluster_rows=FALSE,
		clustering_method="ward.D",
		annotation_col=annotation_col,
		annotation_colors=annotation_colors,
		color=viridis(10),
		# color=mypalette,
		# gaps_row=rowbreaks,
		gaps_col=colbreaks,
		border_color = FALSE
		# breaks=breaksList
	)
}


###########################################
#                                         #
#             Diffusion map               #
#                                         #
###########################################


run_diffMap <- function(data=data, condition=condition, sigma="local"){
	destinyObj <- as.ExpressionSet(as.data.frame(t(data)))
	destinyObj$condition <- factor(condition)
	dm <- DiffusionMap(destinyObj, sigma, rotate = TRUE)
	return(dm)
}

plot_eigenVal <- function(dm=dm){
	linepad <- .5
	plot(
		eigenvalues(dm), 
		ylim = 0:1, 
		pch = 20, 
		xlab ='Diffusion component (DC)', 
		ylab ='Eigenvalue'
	)
}

plot_dm_2D <- function(dm=dm, dc=2, condition=condition, colours=colours){
	DCs <- paste("DC",1:dc, sep="")

	dm_eigen <- data.frame(
			dm@eigenvectors
		)

	DCs.combinations <- combn(DCs,2)
	g <- apply(
		DCs.combinations,
		2,
		function(combination)
		{
			p1 <- ggplot(dm_eigen, aes_string(x=combination[1], y=combination[2])) +
			geom_point(shape = 21, size = 2.5, stroke=0.5, aes(fill=condition)) +
			theme_bw() +
			scale_fill_manual(
				values=colours,
				name=""
			) +
			xlab(combination[1])+
			ylab(combination[2])+
			ggtitle("Diffusion Map")+
			theme(
				axis.text=element_text(size=16),
				axis.title=element_text(size=16),
				legend.text = element_text(size =16),
				legend.title = element_text(size =16 ,face="bold"),
				plot.title = element_text(size=18, face="bold", hjust = 0.5),
				aspect.ratio=1
			)
			print(p1)
		}
	)

}


plot_dm_3D <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours){
	cond <- factor(condition)
	col <- factor(condition)
	levels(col) <- colours
	col <- as.vector(col)
	DCs <- paste("DC",dc, sep="")

	data <- data.frame(
		dm@eigenvectors[,DCs[1]], 
		dm@eigenvectors[,DCs[2]], 
		dm@eigenvectors[,DCs[3]]
	)
	colnames(data) <- DCs

	plot3d(
		data,
		col=col,
		size=6.5,
		box = FALSE
	)

	# legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")

	# plot3d(
	# 	data,
	# 	size=7.5,
	# 	add = TRUE
	# )

}




plot_dm_3D_transparent <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours){
	cond <- factor(condition)
	col <- factor(condition)
	levels(col) <- colours
	col <- as.vector(col)
	DCs <- paste("DC",dc, sep="")

	data <- data.frame(
		dm@eigenvectors[,DCs[1]], 
		dm@eigenvectors[,DCs[2]], 
		dm@eigenvectors[,DCs[3]]
	)
	colnames(data) <- DCs

	plot3d(
		data,
		col=col,
		size=6.5,
		alpha=0,
		box = FALSE
	)

	legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")

}




###########################################
#                                         #
#          Lineage prediction             #
#                                         #
###########################################

get_lineage <- function(dm=dm, dim=c(0:20), condition=condition, start=start, end=end, shrink.method="cosine"){

	data <- data.frame(
		dm@eigenvectors[,dim]
	)
	crv <- slingshot(
		data, 
		condition, 
		start.clus = start, 
		end.clus=end,
		maxit=100000,
		shrink.method=shrink.method
		# shrink.method="tricube"
	)

	return(crv)
}


# Thanks to Guillaume Devailly for this function
rankKeepNA <- function(x) {
    return(
        ifelse(
            is.na(x),
            NA,
            rank(
            	x, 
            	na.last = TRUE, 
            	ties.method="random"
            )
        )
    )
}


get_pseudotime <- function(pseudotime, wthres=wthres, ranked=TRUE){
	pseudoT <- list()
	for(lineage in 1:length(pseudotime@curves))local({
		curve <- pseudotime@curves[[lineage]]
		lambda <- curve$lambda
		weight <- curve$w
		ps <- curve$lambda
		ps[weight < wthres] <- NA
		if (ranked==TRUE){
			ps <- rankKeepNA(ps)
		}
		pseudoT[[lineage]] <<- ps
	})
	df <- t(do.call("rbind",pseudoT))
	colnames(df) <- names(pseudotime@curves)
	return(df)
}



plot_gene_per_lineage <- function(
	rpkm_matrix=rpkm_matrix, 
	pseudotime=pseudotime, 
	geneName=geneName, 
	stages=stages, 
	clusters=clusters, 
	stage_colors=stage_colors,
	cluster_colors=cluster_colors
	){

	myplots <- list()
	total_pseudotime <- vector()
	for (lineages in 1:ncol(pseudotime)){
		lineage <- as.vector(pseudotime[,lineages])
		total_pseudotime <- c(total_pseudotime, lineage)
		total_pseudotime <- na.omit(total_pseudotime)
	}
	max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
	max_pseudotime <- max(total_pseudotime)

	pseudotime_data <- data.frame(
		pseudotime=numeric(),
		lineage=numeric(),
		stages=character(),
		clusters=character(),
		gene=numeric()
	)

	for (lineages in 1:ncol(pseudotime)){

		i <- lineages

		lineage <- pseudotime[,lineages]
		sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)

		data <- data.frame(
			pseudotime=lineage,
			lineage=paste("Lineage ",lineages, sep=""),
			stages=stages,
			clusters=clusters,
			gene=t(sub_data)
			)
		colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene")

		pseudotime_data <- rbind(pseudotime_data, data)

	}


	#
	# Old plot I keep for record in case I need it
	#

	# p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
	# scale_shape_manual(values = 21:25) +
	# geom_point(size = 3, stroke=0.7, aes(shape=condition1, fill=condition2), color="white", na.rm = TRUE)+
	# geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.75)+
	# ylab("log(RPKM+1)") +
	# theme_bw() +
	# ggtitle(geneName)+
	# # scale_color_manual(
	# # 	values=colours
	# # ) +
	# scale_fill_manual(
	# 	values=colours
	# ) +
	# guides(
	# 	fill = guide_legend(override.aes = list(color=colours, fill=colours))
 #    ) +
	# theme(
	# 	axis.text=element_text(size=16),
	# 	axis.title=element_text(size=16),
	# 	legend.text = element_text(size =16),
	# 	legend.title=element_blank(),
	# 	plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
	# 	aspect.ratio=0.5#,
	# 	# legend.position="bottom"
	# ) +
	# facet_grid(.~lineage) +
	# theme(
	# 	strip.text.x = element_text(size = 16)
	# )

	# pseudotime_data[pseudotime_data=="Lineage 1"] <- "Sertoli lineage"
	# pseudotime_data[pseudotime_data=="Lineage 2"] <- "Leydig lineage"
	# pseudotime_data[pseudotime_data=="Lineage 3"] <- "Progenitor lineage"

	# p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
	# # geom_point(shape=21, size = 3, stroke=0.7, aes(fill=clusters), color="white", na.rm = TRUE)+
	# # geom_point(shape=108, size = 8,  aes(y=-1, color=stages), na.rm = TRUE)+
	# # geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.75)+
	# geom_smooth(aes(group=lineage, color=lineage, fill=lineage), na.rm = TRUE, method="loess", span=0.5)+
	# ylab("log(RPKM+1)") +
	# theme_bw() +
	# ggtitle(geneName)+
	# scale_color_manual(
	# 	values=cluster_colors
	# ) +
	# scale_fill_manual(
	# 	values=cluster_colors
	# ) +
	# # ylim(0, max(pseudotime_data$gene))+
	# # guides(
	# # 	fill = guide_legend(override.aes = list(color=colours, fill=colours))
 # #    ) +
	# theme(
	# 	axis.text=element_text(size=16),
	# 	axis.title=element_text(size=16),
	# 	legend.text = element_text(size =16),
	# 	legend.title=element_blank(),
	# 	plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
	# 	aspect.ratio=0.5,
	# 	legend.position="bottom"
	# ) +
	# # facet_grid(.~lineage) +
	# theme(
	# 	strip.text.x = element_text(size = 16)
	# )

	p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
	geom_point(shape=21, size = 2.5, stroke=0.5, aes(fill=clusters), color="white", na.rm = TRUE)+
	geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE)+
	geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.5)+
	ylab("log(RPKM+1)") +
	theme_bw() +
	ggtitle(geneName)+
	scale_color_manual(
		values=stage_colors
	) +
	scale_fill_manual(
		values=cluster_colors
	) +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title=element_blank(),
		plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
		aspect.ratio=0.8,
		legend.position="bottom"
	) +
	facet_grid(.~lineage) +
	theme(
		strip.text.x = element_text(size = 16)
	)
	print(p)
}


plot_smoothed_gene_per_lineage <- function(
	rpkm_matrix=rpkm_matrix, 
	pseudotime=pseudotime, 
	lin=lin,
	geneName=geneName, 
	stages=stages, 
	clusters=clusters, 
	stage_colors=stage_colors,
	cluster_colors=cluster_colors,
	lineage_colors=lineage_colors
	){
	max_time <- max(pseudotime[,c(1,2)], na.rm = TRUE)
	pseudotime <- pseudotime[,lin]

	lineage_colors <- lineage_colors[lin]

	myplots <- list()
	total_pseudotime <- vector()

	if (length(lin)==1){

		lineage <- pseudotime
		total_pseudotime <- c(total_pseudotime, lineage)
		total_pseudotime <- na.omit(total_pseudotime)
		max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
		max_pseudotime <- max(total_pseudotime)

		pseudotime_data <- data.frame(
			pseudotime=numeric(),
			lineage=numeric(),
			stages=character(),
			clusters=character(),
			gene=numeric()
		)


		lineage <- pseudotime
		sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)

		data <- data.frame(
			pseudotime=lineage,
			lineage=paste("Lineage ",lin, sep=""),
			stages=stages,
			clusters=clusters,
			gene=t(sub_data)
		)

		colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene")

		# print(lin)

		pseudotime_data <- rbind(pseudotime_data, data)

		p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
		geom_point(shape=21, size = 3,  aes(fill=clusters), color="white", na.rm = TRUE)+
		geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE)+
		geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.5)+
		ylab("log(RPKM+1)") +
		theme_bw() +
		ggtitle(geneName)+
		# scale_color_manual(
		# 	values=cluster_colors
		# ) +
		scale_fill_manual(
			values=cluster_colors
		) +
		scale_color_manual(
			values=stage_colors
		) +
		expand_limits(x = c(0,max_time))+
		expand_limits(y = c(0,2))+
		theme(
			axis.text=element_text(size=16),
			axis.title=element_text(size=16),
			legend.text = element_text(size =16),
			legend.title=element_blank(),
			plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
			aspect.ratio=0.5,
			legend.position="bottom",
			strip.text.x = element_text(size = 16)
		)


	} else {

		for (lineages in 1:ncol(pseudotime)){
			lineage <- as.vector(pseudotime[,lineages])
			total_pseudotime <- c(total_pseudotime, lineage)
			total_pseudotime <- na.omit(total_pseudotime)
		}
	
		max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
		max_pseudotime <- max(total_pseudotime)

		pseudotime_data <- data.frame(
			pseudotime=numeric(),
			lineage=numeric(),
			stages=character(),
			clusters=character(),
			gene=numeric()
		)


		for (lineages in 1:ncol(as.data.frame(pseudotime))){
			lineage <- pseudotime[,lineages]
			sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)

			data <- data.frame(
				pseudotime=lineage,
				lineage=paste("Lineage ",lineages, sep=""),
				stages=stages,
				clusters=clusters,
				gene=t(sub_data)
			)

			colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene")

			pseudotime_data <- rbind(pseudotime_data, data)
		}

		p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
		geom_smooth(aes(group=lineage, color=lineage, fill=lineage), na.rm = TRUE, method="loess", span=0.5)+
		ylab("log(RPKM+1)") +
		theme_bw() +
		ggtitle(geneName)+
		scale_color_manual(
			values=lineage_colors
		) +
		scale_fill_manual(
			values=lineage_colors
		) +
		expand_limits(y = c(0,2))+
		theme(
			axis.text=element_text(size=16),
			axis.title=element_text(size=16),
			legend.text = element_text(size =16),
			legend.title=element_blank(),
			plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
			aspect.ratio=0.5,
			legend.position="bottom",
			strip.text.x = element_text(size = 16)
		)

	}

	print(p)
}


plot_smoothed_gene_per_lineage_sex <- function(
	rpkm_matrix=rpkm_matrix, 
	pseudotime=pseudotime, 
	lin=lin,
	geneName=geneName, 
	stages=stages, 
	clusters=clusters, 
	stage_colors=stage_colors,
	cluster_colors=cluster_colors,
	lineage_colors=lineage_colors
	){
	max_time <- max(pseudotime[,c(1,2)], na.rm = TRUE)
	pseudotime <- pseudotime[,lin]

	lineage_colors <- lineage_colors[lin]

	myplots <- list()
	total_pseudotime <- vector()

	if (length(lin)==1){

		# lineage <- pseudotime
		# total_pseudotime <- c(total_pseudotime, lineage)
		# total_pseudotime <- na.omit(total_pseudotime)
		# max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
		# max_pseudotime <- max(total_pseudotime)

		lineage <- as.vector(pseudotime)
		total_pseudotime <- c(total_pseudotime, lineage)
		total_pseudotime <- na.omit(total_pseudotime)
	
		max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
		max_pseudotime <- max(total_pseudotime)

		pseudotime_data <- data.frame(
			pseudotime=numeric(),
			lineage=numeric(),
			stages=character(),
			clusters=character(),
			gene=numeric(),
			sex= character()
		)


		lineage <- pseudotime
		sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)
		data <- data.frame(
			pseudotime=lineage,
			lineage=paste("Lineage ",lin, sep=""),
			stages=stages,
			clusters=clusters,
			gene=t(sub_data),
			sex=sapply(strsplit(clusters, "_"), `[`, 1)
		)

		colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene", "sex")
		pseudotime_data <- rbind(pseudotime_data, data)

		print(head(pseudotime_data))

		pseudotime_data<-pseudotime_data[!(pseudotime_data[,"lineage"]=="Lineage 1" & pseudotime_data[,"sex"]=="XX"),]
		pseudotime_data<-pseudotime_data[!(pseudotime_data[,"lineage"]=="Lineage 2" & pseudotime_data[,"sex"]=="XY"),]

		rect1 <- data.frame(xmin=240, xmax=290, ymin=-Inf, ymax=Inf)
		rect2 <- data.frame(xmin=380, xmax=410, ymin=-Inf, ymax=Inf)

		p3 <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
		# geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
		geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
		geom_point(shape=21, size = 2,  aes(fill=clusters), color="white", stroke = 0.3, na.rm = TRUE, position="dodge", alpha=0.8)+
		geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE, position="dodge")+
		geom_smooth(aes(group=sex, color=sex), size=2, se=FALSE, na.rm = TRUE, method="loess", span=0.5)+
		ylab("log(RPKM+1)") +
		theme_bw()+
		ggtitle(geneName)+
		# scale_color_manual(
		# 	values=cluster_colors
		# ) +
		scale_fill_manual(
			values=cluster_colors
		) +
		scale_color_manual(
			values=stage_colors
		) +
		expand_limits(x = c(0,max_time))+
		expand_limits(y = c(0,2))+
		theme(
			axis.text=element_text(size=12),
			axis.title=element_text(size=12),
			legend.text = element_text(size =12),
			legend.title=element_blank(),
			plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
			aspect.ratio=0.5,
			legend.position="none",
			strip.text.x = element_text(size = 12)
		)


	} else {

		for (lineages in 1:ncol(pseudotime)){
			lineage <- as.vector(pseudotime[,lineages])
			total_pseudotime <- c(total_pseudotime, lineage)
			total_pseudotime <- na.omit(total_pseudotime)
		}
	
		max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
		max_pseudotime <- max(total_pseudotime)

		pseudotime_data <- data.frame(
			pseudotime=numeric(),
			lineage=numeric(),
			stages=character(),
			clusters=character(),
			gene=numeric(),
			sex= character()
		)


		for (lineages in 1:ncol(as.data.frame(pseudotime))){
			lineage <- pseudotime[,lineages]
			sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)

			data <- data.frame(
				pseudotime=lineage,
				lineage=paste("Lineage ",lineages, sep=""),
				stages=stages,
				clusters=clusters,
				gene=t(sub_data),
				sex=sapply(strsplit(clusters, "_"), `[`, 1)
			)

			colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene", "sex")

			pseudotime_data <- rbind(pseudotime_data, data)
		}


		pseudotime_data<-pseudotime_data[!(pseudotime_data[,"lineage"]=="Lineage 1" & pseudotime_data[,"sex"]=="XX"),]
		pseudotime_data<-pseudotime_data[!(pseudotime_data[,"lineage"]=="Lineage 2" & pseudotime_data[,"sex"]=="XY"),]

		rect1 <- data.frame(xmin=240, xmax=290, ymin=-Inf, ymax=Inf)
		rect2 <- data.frame(xmin=380, xmax=410, ymin=-Inf, ymax=Inf)

		p3 <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
		geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
		geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
		geom_point(shape=21, size = 2,  aes(fill=clusters), color="white", stroke = 0.3, na.rm = TRUE, position="dodge", alpha=0.8)+
		geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE, position="dodge")+
		geom_smooth(aes(group=lineage, color=lineage), size=2, se=FALSE, na.rm = TRUE, method="loess", span=0.5)+
		ylab("log(RPKM+1)") +
		theme_bw()+
		ggtitle(geneName)+
		# scale_color_manual(
		# 	values=cluster_colors
		# ) +
		scale_fill_manual(
			values=cluster_colors
		) +
		scale_color_manual(
			values=stage_colors
		) +
		expand_limits(x = c(0,max_time))+
		expand_limits(y = c(0,2))+
		theme(
			axis.text=element_text(size=12),
			axis.title=element_text(size=12),
			legend.text = element_text(size =12),
			legend.title=element_blank(),
			plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
			aspect.ratio=0.5,
			legend.position="none",
			strip.text.x = element_text(size = 12)
		)

	}
	# print(p1)
	# print(p3)
	return(p3)
}




plot_aligned_lineage_sex <- function(
	rpkm_matrix=rpkm_matrix, 
	pseudotime=pseudotime, 
	lin=lin,
	stage_colors=stage_colors
	){

	sex <- sapply(strsplit(rownames(pseudotime), "_"), `[`, 2)

	lin_female <- pseudotime[grep("XX",rownames(pseudotime)), lin]
	lin_male <- pseudotime[grep("XY",rownames(pseudotime)), lin]
	lin_pseudotime <- c(lin_female, lin_male)

	cells <- c(names(lin_female), names(lin_male))
	sex <- sapply(strsplit(cells, "_"), `[`, 2)
	stages <- sapply(strsplit(cells, "_"), `[`, 1)


	pseudotime_data <- data.frame(
		pseudotime=lin_pseudotime,
		stages=stages,
		sex= sex
	)

	rownames(pseudotime_data) <- cells

	print(head(pseudotime_data))

	# rect1 <- data.frame(xmin=240, xmax=290, ymin=-Inf, ymax=Inf)
	# rect2 <- data.frame(xmin=380, xmax=410, ymin=-Inf, ymax=Inf)

	g <- ggplot(pseudotime_data, aes(x=pseudotime, y=sex))+
	# geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
	geom_point(size = 4,  aes(color=stages), stroke = 0.3, na.rm = TRUE, position="dodge", alpha=0.8)+
	# geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE, position="dodge")+
	# geom_smooth(aes(group=sex, color=sex), size=2, se=FALSE, na.rm = TRUE, method="loess", span=0.5)+
	ylab("pseudotime") +
	theme_bw()+

	scale_color_manual(
		values=stage_colors
	) +
	theme(
		axis.text=element_text(size=12),
		axis.title=element_text(size=12),
		legend.text = element_text(size =12),
		legend.title=element_blank(),
		plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
		aspect.ratio=0.05,
		legend.position="none",
		strip.text.x = element_text(size = 12)
	)

	print(g)
}




compare_lineage <- function(pseudotime=pseudotime, condition=condition){
	myLineages <- list()
	for (lineages in 1:length(pseudotime))local({
		i <- lineages
		lineage <- pseudotime[[lineages]]$pseudotime
		myLineages[[i]] <<- lineage
	})
	df <- t(do.call("rbind",myLineages))
	df[df>=0] <- 1
	df[is.na(df)] <- 0
	sdf <- as.character(rowSums(df))
	return(sdf)

}


###########################################
#                                         #
#            Loess smoothing              #
#                                         #
###########################################

smooth_gene_exp <- function(data=data, pseudotime=pseudotime, span=0.75){
	smooth_data <- data
	for (gene in 1:nrow(data)){
		gene_exp <- t(data[gene,])
		smooth <- loess(formula=gene_exp~pseudotime, span=span)
		smooth_data[gene,] <- predict(smooth, newdata=pseudotime)
	}
	return(smooth_data)
}


get_var_genes_pseudotime <- function(rpkm, count, pseudotime, lineageNb=0, clusters){


	lineage_pseudotime_all <- pseudotime[,lineageNb]
	lineage <- pseudotime[!is.na(pseudotime[,lineageNb]),lineageNb,drop=FALSE]

	rpkm_exp_lineage <- rpkm[,rownames(lineage)]
	count_exp_lineage <- count[,rownames(lineage)]
	clusters_in_lineage <- clusters[names(clusters) %in% rownames(lineage)]

	print("Prepare for DE...")
	genes_pseudoT <- prepare_for_DE (
		count, 
		lineage_pseudotime_all, 
		lineage_pseudotime_all
	)
	print("Done.")

	genes_pseudoT$Pseudotime <- lineage_pseudotime_all

	print("Compute DE genes...")

	#Select cells from the lineage of interest
	genes_pseudoT <- genes_pseudoT[,rownames(lineage)]
	# Re-detect how many cells express each genes in the subset of cells
	genes_pseudoT <- detectGenes(genes_pseudoT, min_expr = 5)
	# Remove genes expressed in less than 10 cells
	genes_pseudoT <- genes_pseudoT[fData(genes_pseudoT)$num_cells_expressed >= 10, ]

	DE_genes_pseudoT <- differentialGeneTest(
		genes_pseudoT, 
		fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
		cores = 3
	)
	print("Done.")

	return(DE_genes_pseudoT)
}


get_var_genes_pseudotime_special <- function(rpkm, count, pseudotime, lineageNb=0, clusters){

	if (lineageNb > 0) {
		lineage_pseudotime_all <- pseudotime[,lineageNb]
		lineage <- pseudotime[!is.na(pseudotime[,lineageNb]),lineageNb,drop=FALSE]
	} else {
		lineage_pseudotime_all <- pseudotime
		lineage <- pseudotime[!is.na(pseudotime),drop=FALSE]
	}


	rpkm_exp_lineage <- rpkm[,names(lineage)]
	count_exp_lineage <- count[,names(lineage)]
	clusters_in_lineage <- clusters[names(clusters) %in% names(lineage)]

	print("Prepare for DE...")
	genes_pseudoT <- prepare_for_DE (
		count, 
		lineage_pseudotime_all, 
		lineage_pseudotime_all
	)
	print("Done.")

	genes_pseudoT$Pseudotime <- lineage_pseudotime_all

	print("Compute DE genes...")

	#Select cells from the lineage of interest
	genes_pseudoT <- genes_pseudoT[,names(lineage)]
	# Re-detect how many cells express each genes in the subset of cells
	genes_pseudoT <- detectGenes(genes_pseudoT, min_expr = 5)
	# Remove genes expressed in less than 10 cells
	genes_pseudoT <- genes_pseudoT[fData(genes_pseudoT)$num_cells_expressed >= 10, ]

	DE_genes_pseudoT <- differentialGeneTest(
		genes_pseudoT, 
		fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
		cores = 3
	)
	print("Done.")

	return(DE_genes_pseudoT)
}


get_gene_clustering <- function(DE_genes_pseudoT, rpkm, pseudotime, lineageNb, cell_clusters, clusterNb=10, qvalue=0.1){
	lineage <- pseudotime[!is.na(pseudotime[,lineageNb]),lineageNb,drop=FALSE]
	clusters_in_lineage <- cell_clusters[names(cell_clusters) %in% rownames(lineage)]

	print(paste("Select genes with q-value<", qvalue, "...", sep=""))
	sig_gene_pseudoT <- row.names(subset(DE_genes_pseudoT, qval < qvalue))
	sig_var_gene_pseudoT <- rpkm[sig_gene_pseudoT,rownames(lineage)]

	print(paste(length(sig_gene_pseudoT), " selected genes.", sep=""))

	print("Smooth gene expression...")
	smooth <- smooth_gene_exp(
		log(sig_var_gene_pseudoT+1), 
		lineage[,1], 
		span=0.5
	)
	print("Done.")

	smooth <- smooth[ , order(lineage[,1])]

	print("Cluster the genes...")


	# cell_cluster_palette <- c(
	# 	C1="#560047", 
	# 	C2="#a53bad", 
	# 	C3="#eb6bac", 
	# 	C4="#ffa8a0"
	# )

	cell_cluster_palette <- c(
		C2="#aeff01", 
		C4="#009900", 
		C1="#457cff", 
		C3="#00c5ec",
		C6="#8dfaff",
		C5="#a38cff"
	)


	gene_clustering <- pheatmap(
		smooth, 
		scale="row", 
		cutree_rows=clusterNb,
		clustering_method="ward.D",
		silent=TRUE
	)

	print("Plot heatmap...")

	clusters <- cutree(gene_clustering$tree_row, k = clusterNb)
	clustering <- data.frame(clusters)
	clustering[,1] <- as.character(clustering[,1])
	colnames(clustering) <- "Gene_Clusters"

	gene_cluster_palette <- c(
	  '#a6cee3',
	  '#1f78b4',
	  '#b2df8a',
	  '#33a02c',
	  '#fb9a99',
	  '#e31a1c',
	  '#fdbf6f',
	  '#ff7f00',
	  '#cab2d6',
	  '#6a3d9a',
	  '#ffff99',
	  '#b15928',
	  '#49beaa',
	  '#611c35',
	  '#2708a0'
	  )
	gene_cluster_colors <- gene_cluster_palette[1:max(clusters)]
	names(gene_cluster_colors) <- 1:max(clusters)


	annotation_row <- clustering

	annotation_col <- data.frame(
		Cell_Clusters=clusters_in_lineage,
		Stages=sapply(strsplit(names(clusters_in_lineage), "_"), `[`, 1)
	)

	annotation_colors <- list(
		Stages=c(
			E10.5="#2754b5", 
			E11.5="#8a00b0", 
			E12.5="#d20e0f", 
			E13.5="#f77f05", 
			E16.5="#f9db21",
			P6="#b0f137"
		),
		Cell_Clusters=cell_cluster_palette,
		Gene_Clusters=gene_cluster_colors
	)

	# Color palette for the heatmap
	cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58', '#081d58'))
	warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
	mypalette <- c(rev(cold(20)), warm(21))
	# mypalette <- c(rev(cold(15)), warm(16))
	breaksList = seq(-10, 10, by = 0.5)

	pheatmap(
		smooth, 
		scale="row", 
		cluster_cols=FALSE, 
		show_colnames=FALSE, 
		show_rownames=FALSE, 
		cutree_rows=clusterNb, 
		clustering_method="ward.D",
		annotation_col=annotation_col,
		annotation_row=annotation_row,
		annotation_colors=annotation_colors,
		annotation_names_row=FALSE,
		color=mypalette
	)

	return(clustering)
}
