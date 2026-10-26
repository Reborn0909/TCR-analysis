



argv <- commandArgs(trailingOnly=TRUE)
if (length(argv) %% 2L) stop('Arguments must be --key value pairs')
opt <- setNames(as.list(argv[seq(2,length(argv),2)]),sub('^--','',argv[seq(1,length(argv),2)]))
pick <- function(k, default=NULL) if (!is.null(opt[[k]])) opt[[k]] else default
input <- pick('input'); if (is.null(input)) stop('--input is required')
out <- pick('output','results'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
seed <- as.integer(pick('seed','20261004')); B <- as.integer(pick('iterations','2000'))
if (B<1000L) stop('At least 1000 downsampling iterations are required')
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(edgeR); library(ggplot2)})
set.seed(seed)
save_csv <- function(d,n) write.csv(d,file.path(out,paste0('main_',n,'.csv')),row.names=FALSE)
save_plot <- function(p,n,w=8,h=5) {
 ggsave(file.path(out,paste0('main_',n,'.pdf')),p,width=w,height=h)
 ggsave(file.path(out,paste0('main_',n,'.png')),p,width=w,height=h,dpi=220,device=grDevices::png)
}
cat('Start:',format(Sys.time(),tz='Asia/Shanghai'),'seed=',seed,' iterations=',B,'\n')
obj <- readRDS(input); md <- obj[[]]
stopifnot(ncol(obj)==nrow(md),!anyDuplicated(rownames(md)))
md$cell <- rownames(md)
sample_order <- c('H-CD8Trp-1','E-CD8Trp-1','E-CD8DP-1','H-CD8Trp-2','E-CD8Trp-2','E-CD8DP-2')
groups <- c('H-Trp','E-Trp','E-DP')
stopifnot(setequal(unique(md$sample),sample_order))
map <- unique(md[c('sample','group')]); map <- map[match(sample_order,map$sample),]
map$donor <- paste0('Donor',sub('.*-','',map$sample))
map$donor_evidence <- 'Historical user-executed code: DeepSeek 0orvarra385ufmj58b message key 413; sample suffix 1/2'
map$SS_Pre_Post_label <- c('SS','Pre-sel','Post-sel')[match(map$group,groups)]
map$label_status <- 'Historical user code key 413 plus terminology keys 415/416; biological workflow checked against local PPT'
map$experimental_strategy <- c('Unexpanded input','Sorted before expansion','Expanded before sorting')[match(map$group,groups)]
axis_labels <- c('H-Trp'='SS\n(H-Trp)','E-Trp'='Pre-sel\n(E-Trp)','E-DP'='Post-sel\n(E-DP)')
if (!is.null(pick('mapping'))) {
 supplied <- read.csv(pick('mapping'),stringsAsFactors=FALSE)
 stopifnot(all(c('sample','donor','group') %in% names(supplied)), !anyDuplicated(supplied$sample), setequal(supplied$sample,map$sample))
 map <- supplied[match(sample_order,supplied$sample),]
}
save_csv(map,'sample_mapping')
md$donor <- map$donor[match(md$sample,map$sample)]
md$group <- factor(md$group,levels=groups)
qc <- do.call(rbind,lapply(split(md,md$sample),function(d)data.frame(sample=d$sample[1],group=as.character(d$group[1]),donor=d$donor[1],post_QC_cells=nrow(d),TCR_assigned_cells=sum(!is.na(d$CTstrict)),missing_TCR_as_single=sum(is.na(d$CTstrict)&d$cloneSize_real==1),median_genes=median(d$nFeature_RNA),median_UMI=median(d$nCount_RNA),median_mito_pct=median(d$percent.mt))))
save_csv(qc,'qc_and_TCR_assignment_audit'); print(qc)

eligible <- !is.na(md$CTstrict) & nzchar(md$CTstrict) & md$CTstrict!='NA'
eligibility <- 'Post-QC cells with assigned TCR; productive status NOT independently verified'
productive_status <- rep(NA_integer_,nrow(map)); contig_stats <- list()
if (!is.null(pick('contig-root'))) {
 root <- pick('contig-root')
 all_files <- list.files(root,pattern='filtered_contig_annotations\\.csv(\\.gz)?$',recursive=TRUE,full.names=TRUE)
 if (!length(all_files)) stop('No filtered_contig_annotations.csv files found under --contig-root')
 productive_cells <- character()
 for (i in seq_len(nrow(map))) {
  s <- map$sample[i]; f <- all_files[grepl(s,all_files,fixed=TRUE)]
  if (length(f)!=1L) stop('Need one unambiguous contig file for ',s, '; found ',length(f))
  z <- read.csv(f,stringsAsFactors=FALSE)
  stopifnot(all(c('barcode','productive','chain') %in% names(z)))
  tf <- function(v) toupper(as.character(v)) %in% c('TRUE','T','1')
  k <- tf(z$productive) & z$chain %in% c('TRA','TRB')
  if ('is_cell' %in% names(z)) k <- k & tf(z$is_cell)
  if ('high_confidence' %in% names(z)) k <- k & tf(z$high_confidence)
  cells <- unique(paste0(s,'_',z$barcode[k]))
  productive_status[i] <- sum(md$sample==s & md$cell %in% cells)
  productive_cells <- c(productive_cells,cells)
  contig_stats[[i]] <- data.frame(sample=s,file=f,n_contigs=nrow(z),n_productive_contigs=sum(k),productive_barcodes=length(cells),productive_postQC_cells=productive_status[i],has_reads='reads' %in% names(z),has_UMIs='umis' %in% names(z))
 }
 eligible <- eligible & md$cell %in% productive_cells
 eligibility <- 'Post-QC assigned TCR cells with original productive TRA/TRB contig (is_cell/high_confidence when available)'
 save_csv(do.call(rbind,contig_stats),'original_contig_audit')
}
save_csv(data.frame(sample=map$sample,productive_TCR_bearing_cells=productive_status,status=eligibility),'productive_status')

shannon <- function(x) {p<-as.numeric(table(x)); p<-p/sum(p); -sum(p*log(p))}
tcr_summary <- list(); draws <- list(); trends <- list()
for (definition in c('CTgene','CTstrict')) {
 e <- eligible & !is.na(md[[definition]]) & nzchar(md[[definition]])
 ns <- table(factor(md$sample[e],levels=sample_order)); n_down <- min(ns)
 if (n_down<1) stop('An eligible TCR sample is empty')
 cat('TCR:',definition,'n_downsample=',n_down,'\n')
 for (s in sample_order) {
  x <- md[[definition]][e & md$sample==s]; m <- match(s,map$sample)
  sim <- vapply(seq_len(B),function(i){if(length(x)==n_down)shannon(x) else shannon(x[sample.int(length(x),n_down,replace=FALSE)])},numeric(1))
  key <- paste(definition,s,sep='_')
  tcr_summary[[key]] <- data.frame(sample=s,group=map$group[m],donor=map$donor[m],clonotype_definition=definition,eligible_TCR_cells=length(x),unique_clonotypes=length(unique(x)),raw_Shannon=shannon(x),n_downsample=n_down,iterations=B,seed=seed,downsample_median=median(sim),sampling_lower_2.5pct=unname(quantile(sim,.025)),sampling_upper_97.5pct=unname(quantile(sim,.975)),eligibility=eligibility)
  draws[[key]] <- data.frame(sample=s,clonotype_definition=definition,iteration=seq_len(B),Shannon=sim)
 }
}
ts <- do.call(rbind,tcr_summary); td<-do.call(rbind,draws)
save_csv(ts,'TCR_equal_cell_downsample_summary');save_csv(td,'TCR_equal_cell_downsample_iterations')
for (definition in c('CTgene','CTstrict')) for (donor in unique(map$donor)) {
 d <- ts[ts$clonotype_definition==definition & ts$donor==donor,]; d<-d[match(groups,d$group),]
 for (j in 1:2) trends[[paste(definition,donor,j)]]<-data.frame(donor=donor,clonotype_definition=definition,comparison=paste(d$group[j+1],'vs',d$group[j]),raw_delta=d$raw_Shannon[j+1]-d$raw_Shannon[j],downsample_delta=d$downsample_median[j+1]-d$downsample_median[j],direction_preserved=sign(d$raw_Shannon[j+1]-d$raw_Shannon[j])==sign(d$downsample_median[j+1]-d$downsample_median[j]))
}
save_csv(do.call(rbind,trends),'TCR_donor_trends')
ts$group<-factor(ts$group,levels=groups)
p<-ggplot(ts,aes(group,downsample_median,group=donor,color=donor))+geom_line()+geom_point(size=2)+geom_errorbar(aes(ymin=sampling_lower_2.5pct,ymax=sampling_upper_97.5pct),width=.08)+facet_wrap(~clonotype_definition)+scale_x_discrete(labels=axis_labels)+theme_bw(base_size=11)+labs(x='Experimental condition',y='Shannon diversity, natural log',title='Equal-cell TCR resampling by donor',subtitle=paste0(B,' draws; bars = sampling percentiles, not donor confidence intervals'),caption=eligibility)
save_plot(p,'TCR_equal_cell_downsample')

legacy_path <- pick('legacy-metadata',file.path(out,'main_legacy_minimal_metadata.rds'))
if (file.exists(legacy_path)) {
 leg <- readRDS(legacy_path); score_cols <- grep('(Score1$|GBD_Pathogenic1$|CD8_Treg1$|GBD_Effector1$)',names(leg),value=TRUE)
 summaries<-list()
 for (s in sample_order) for (sc in score_cols) {
  ix<-leg$sample==s; z<-leg[[sc]][ix]; m<-match(s,map$sample)
  summaries[[paste(s,sc)]]<-data.frame(sample=s,group=map$group[m],donor=map$donor[m],signature=sc,n_legacy_cells=sum(ix),n_scores=sum(is.finite(z)),mean_score=mean(z,na.rm=TRUE),median_score=median(z,na.rm=TRUE),source='F:/downstream/seurat_minimal_for_sharing.rds; historical QC subset, not recalculated on new QC')
 }
 ls<-do.call(rbind,summaries);save_csv(ls,'legacy_signature_donor_summary')
 ls$group<-factor(ls$group,levels=groups)
 save_plot(ggplot(ls,aes(group,mean_score,group=donor,color=donor))+geom_line()+geom_point()+facet_wrap(~signature,scales='free_y',ncol=4)+theme_bw(base_size=9)+labs(x='Original condition',y='Sample mean of historical cell score',title='Historical signatures summarized by donor (old QC object)',subtitle='Two donor trajectories; no cell-level P values; historical gene-set definitions not recovered'), 'legacy_signature_donor_summary',12,10)
}

counts <- LayerData(obj,assay='RNA',layer='counts')
stopifnot(identical(colnames(counts),rownames(md)),all(counts@x>=0),all(counts@x==round(counts@x)))
cell_index <- match(md$sample,sample_order)
indicator <- sparseMatrix(i=seq_len(nrow(md)),j=cell_index,x=1,dims=c(nrow(md),length(sample_order)))
pb <- as.matrix(counts %*% indicator);colnames(pb)<-sample_order
stopifnot(all(abs(colSums(pb)-tapply(md$nCount_RNA,factor(md$sample,levels=sample_order),sum))<.1))
save_csv(data.frame(gene=rownames(pb),pb,check.names=FALSE),'pseudobulk_raw_UMI_counts')
cat('Pseudobulk libraries:\n');print(colSums(pb))

memory_genes<-c('TCF7','CCR7','IL7R','LEF1');exhaust_genes<-c('TOX','HAVCR2','ENTPD1','LAG3')
requested<-c(memory_genes,exhaust_genes)
stopifnot(all(requested %in% rownames(obj)))


recipe <- list(
 CD8_Trp_50=c('CCL5','GZMK','CD74','GZMA','NKG7','CST7','DUSP2','HLA-DRB1','HLA-DPA1','HLA-DPB1','COTL1','CD8A','HLA-DRA','IL32','CMC1','LYAR','KLRK1','CD8B','HLA-DQA1','HLA-DRB5','HLA-C','HCST','GZMM','HLA-DQB1','ZFP36L2','HLA-A','B2M','HLA-B','CXCR3','HLA-DMA','CLIC1','EOMES','CCL4','SRGN','SH2D1A','CD99','PLAAT4','TIGIT','CTSW','KLRG1','CD69','PSMB9','CD2','FYN','FABP5','CXCR4','GUK1','PFN1','CKLF','CCL4L2'),
 CD8_Treg_associated=c('FOXP3','IL2RA','CTLA4','TIGIT','IKZF2','LAG3'),
 Stemness=c('TCF7','LEF1','SELL','CCR7','CD27','IL7R','CD28'),
 Exhaustion7=c('PDCD1','HAVCR2','LAG3','TIGIT','TOX','ENTPD1','CTLA4'),
 Pathogenic_cytotoxic=c('GZMB','PRF1','NKG7','GNLY','FGFBP2','CCL5'),
 Costimulation=c('CD28','TNFRSF4','TNFRSF9','TNFRSF18','ICOS','CD27'),
 Migration=c('CXCR3','CXCR4','CCR5','CCR7','SELL','CXCR5'),
 TissueResident=c('ITGAE','ITGA1','CXCR6','CD69','CD103','CD49A'),
 Glycolysis=c('HK2','GPI','PFKL','ALDOA','GAPDH','LDHA','PKM','ENO1'),
 OXPHOS=c('NDUFA1','SDHB','UQCRC1','COX5A','ATP5A1','ATP5B','CYCS'))
alias <- c(CD103='ITGAE',CD49A='ITGA1',ATP5A1='ATP5F1A',ATP5B='ATP5F1B')
resolve <- function(g) {if(g %in% rownames(obj)) g else if(g %in% names(alias))unname(alias[g]) else g}
gene_audit<-do.call(rbind,lapply(names(recipe),function(k){r<-vapply(recipe[[k]],resolve,character(1));data.frame(signature=k,original_gene=recipe[[k]],resolved_gene=r,present_in_newQC=r %in% rownames(obj),duplicate_after_resolution=duplicated(r),source='DeepSeek 0orvarra385ufmj58b message key 384; historical recipe recomputed') }))
save_csv(gene_audit,'recomputed_signature_gene_audit')
features <- lapply(recipe,function(g) unique(intersect(vapply(g,resolve,character(1)),rownames(obj))))
stopifnot(all(lengths(features)>=3))
all_features<-c(list(memory4=memory_genes,exhaustion4=exhaust_genes),features)
obj<-AddModuleScore(obj,features=all_features,assay='RNA',name='ReviewerModule',nbin=24,ctrl=100,seed=seed,search=FALSE)
logdata<-LayerData(obj,assay='RNA',layer='data')[unique(unlist(all_features)),,drop=FALSE]
recomputed<-list()
for(i in seq_along(all_features))for(s in sample_order){ix<-md$sample==s;m<-match(s,map$sample);sc<-obj[[paste0('ReviewerModule',i)]][,1][ix];recomputed[[paste(i,s)]]<-data.frame(sample=s,group=map$group[m],donor=map$donor[m],signature=names(all_features)[i],n_cells=sum(ix),n_genes=length(all_features[[i]]),mean_control_adjusted_score=mean(sc),median_control_adjusted_score=median(sc),mean_log1p_expression=mean(Matrix::colMeans(logdata[all_features[[i]],ix,drop=FALSE])),source=if(i<=2)'Reviewer-specified 4-gene list' else 'Historical recipe key384 recomputed on newQC; alias resolution audited')}
rs<-do.call(rbind,recomputed);save_csv(rs,'recomputed_signature_donor_summary')
rs$group<-factor(rs$group,levels=groups)
save_plot(ggplot(rs,aes(group,mean_control_adjusted_score,group=donor,color=donor))+geom_line()+geom_point(size=1.7)+facet_wrap(~signature,scales='free_y',ncol=4)+scale_x_discrete(labels=axis_labels)+theme_bw(base_size=10)+labs(x='Experimental condition',y='Sample mean control-adjusted score',title='New QC signatures: each point is one donor-condition',subtitle='Historical recipe key 384, Seurat AddModuleScore; n = 2 donors; no cell-level P values'), 'recomputed_signature_donor_summary',13,10)
score_sample<-list();gene_sample<-list()
for (s in sample_order) {
 ix<-md$sample==s;m<-match(s,map$sample)
 score_sample[[s]]<-data.frame(sample=s,group=map$group[m],donor=map$donor[m],n_cells=sum(ix),memory4_mean_log1p=mean(Matrix::colMeans(logdata[memory_genes,ix,drop=FALSE])),exhaustion4_mean_log1p=mean(Matrix::colMeans(logdata[exhaust_genes,ix,drop=FALSE])),memory4_control_adjusted_mean=mean(obj$ReviewerModule1[ix]),exhaustion4_control_adjusted_mean=mean(obj$ReviewerModule2[ix]))
 for (g in requested) gene_sample[[paste(s,g)]]<-data.frame(sample=s,group=map$group[m],donor=map$donor[m],gene=g,category=if(g %in% memory_genes)'memory/progenitor' else 'exhaustion-associated',n_cells=sum(ix),mean_log1p_normalized_RNA=mean(logdata[g,ix]),percent_expressing=100*mean(logdata[g,ix]>0))
}
ss<-do.call(rbind,score_sample);gs<-do.call(rbind,gene_sample)
save_csv(ss,'requested_memory_exhaustion_donor_scores'); save_csv(gs,'requested_eight_genes_donor_expression')
save_csv(data.frame(module=c(rep('memory4',4),rep('exhaustion4',4)),gene=requested,score_definition='Mean RNA log1p normalized expression; additionally Seurat AddModuleScore nbin=24 ctrl=100 seed=20261004'),'requested_score_gene_sets')
writeLines(capture.output(sessionInfo()),file.path(out,'main_R_sessionInfo.txt'))
rm(obj,counts,indicator,logdata);gc()
gs$group<-factor(gs$group,levels=groups)
save_plot(ggplot(gs,aes(group,mean_log1p_normalized_RNA,group=donor,color=donor))+geom_line()+geom_point(size=2)+facet_wrap(~gene,scales='free_y',ncol=4)+scale_x_discrete(labels=axis_labels)+theme_bw(base_size=11)+labs(x='Experimental condition',y='Sample mean log1p normalized RNA',title='Requested genes: each point is one donor-condition',subtitle='n = 2 donors; descriptive expression, no cell-level P values'), 'requested_eight_genes_donor_expression',12,6)
sl<-do.call(rbind,lapply(setdiff(names(ss),c('sample','group','donor','n_cells')),function(k)data.frame(ss[c('sample','group','donor')],signature=k,score=ss[[k]])))
sl$group<-factor(sl$group,levels=groups)
save_plot(ggplot(sl,aes(group,score,group=donor,color=donor))+geom_line()+geom_point(size=2)+facet_wrap(~signature,scales='free_y',ncol=2)+scale_x_discrete(labels=axis_labels)+theme_bw(base_size=11)+labs(x='Experimental condition',y='Donor-condition score',title='Memory and exhaustion-associated scores',subtitle='Two donors; mean log expression and control-adjusted module scores; no significance stars'), 'requested_memory_exhaustion_donor_scores',10,7)
delta<-list()
for(k in unique(sl$signature))for(d in unique(sl$donor)) {
 z<-sl[sl$signature==k & sl$donor==d,];z<-z[match(groups,z$group),]
 for(j in 1:2)delta[[paste(k,d,j)]]<-data.frame(signature=k,donor=d,comparison=paste(groups[j+1],'vs',groups[j]),paired_change=z$score[j+1]-z$score[j])
}
save_csv(do.call(rbind,delta),'requested_score_paired_changes')

info<-map;info$group<-factor(info$group,levels=groups); info$donor<-factor(info$donor)
design<-model.matrix(~donor+group,data=info);rownames(design)<-info$sample
stopifnot(qr(design)$rank==ncol(design), nrow(design)-qr(design)$rank==2L)
save_csv(data.frame(sample=rownames(design),design,check.names=FALSE),'pseudobulk_design_matrix')
y<-DGEList(pb); keep<-filterByExpr(y,design=design); y<-y[keep,,keep.lib.sizes=FALSE]; y<-calcNormFactors(y,method='TMM')
y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
save_csv(data.frame(sample=colnames(y),group=info$group,donor=info$donor,postQC_cells=qc$post_QC_cells[match(colnames(y),qc$sample)],total_raw_UMIs=colSums(pb),filtered_library_size=y$samples$lib.size,TMM_norm_factor=y$samples$norm.factors,genes_tested=sum(keep),n_independent_donors=length(unique(info$donor)),residual_df=2),'pseudobulk_library_summary')
contrasts<-list('E_Trp_vs_H_Trp'=c(0,0,1,0),'E_DP_vs_E_Trp'=c(0,0,-1,1),'E_DP_vs_H_Trp'=c(0,0,0,1))
de_summary<-list(); donor_lfc<-list();lcpm<-cpm(y,log=TRUE,prior.count=2)
for (nm in names(contrasts)) {
 q<-glmQLFTest(fit,contrast=contrasts[[nm]])
 z<-topTags(q,n=Inf,sort.by='PValue')$table;z$gene<-rownames(z)
 pairs<-switch(nm,E_Trp_vs_H_Trp=c('E-Trp','H-Trp'),E_DP_vs_E_Trp=c('E-DP','E-Trp'),E_DP_vs_H_Trp=c('E-DP','H-Trp'))
 for(d in levels(info$donor)) {
  a<-which(info$donor==d & info$group==pairs[1]); b<-which(info$donor==d & info$group==pairs[2])
  z[[paste0(d,'_logCPM_delta')]]<-lcpm[z$gene,a]-lcpm[z$gene,b]
 }
 z$both_donors_same_direction<-sign(z$Donor1_logCPM_delta)==sign(z$Donor2_logCPM_delta)
 z$inference_status<-'Exploratory paired edgeR QL; 2 independent donors; limited population inference'
 z<-z[,c('gene',setdiff(names(z),'gene'))];save_csv(z,paste0('pseudobulk_DE_',nm))
 de_summary[[nm]]<-data.frame(comparison=nm,genes_tested=nrow(z),FDR_lt_0.05=sum(z$FDR<.05),FDR_lt_0.05_absLFC_gt1=sum(z$FDR<.05 & abs(z$logFC)>1),FDR_lt_0.05_same_direction=sum(z$FDR<.05 & z$both_donors_same_direction),n_donors=2,residual_df=2)
 z$display<-ifelse(z$FDR<.05 & abs(z$logFC)>1,'FDR < 0.05 and |logFC| > 1','Other');z$neglogFDR <- -log10(pmax(z$FDR,1e-300))
 p<-ggplot(z,aes(logFC,neglogFDR,color=display))+geom_point(size=.7,alpha=.65)+scale_color_manual(values=c('FDR < 0.05 and |logFC| > 1'='#b34e3d','Other'='grey65'))+theme_bw(base_size=11)+labs(x=paste(pairs[1],'vs',pairs[2],'log2 fold change'),y='-log10 FDR',title=paste('Donor-paired pseudobulk:',nm),subtitle='Raw UMI sums, TMM, edgeR QL ~ donor + group; exploratory with two donors')+theme(legend.position='bottom')
 save_plot(p,paste0('pseudobulk_volcano_',nm),8,6)
}
ds<-do.call(rbind,de_summary);save_csv(ds,'pseudobulk_DE_summary');print(ds)
gene_pb<-do.call(rbind,lapply(requested,function(g)data.frame(sample=sample_order,group=info$group,donor=info$donor,gene=g,raw_UMI=pb[g,],logCPM=if(g %in% rownames(lcpm))as.numeric(lcpm[g,]) else NA_real_)))
save_csv(gene_pb,'requested_eight_genes_pseudobulk_expression')
cat('Completed:',format(Sys.time(),tz='Asia/Shanghai'),'\n');print(sessionInfo())
