

args <- commandArgs(trailingOnly=TRUE);out<-if(length(args))args[1] else 'results'
suppressPackageStartupMessages({library(edgeR);library(ggplot2)})
save_csv<-function(d,n)write.csv(d,file.path(out,paste0('main_',n,'.csv')),row.names=FALSE)
read_result<-function(n)read.csv(file.path(out,paste0('main_',n,'.csv')),check.names=FALSE,stringsAsFactors=FALSE)
save_plot<-function(p,n,w=8,h=5){ggsave(file.path(out,paste0('main_',n,'.pdf')),p,width=w,height=h);ggsave(file.path(out,paste0('main_',n,'.png')),p,width=w,height=h,dpi=220,device=grDevices::png)}
groups<-c('H-Trp','E-Trp','E-DP');axis_labels<-c('H-Trp'='SS\n(H-Trp)','E-Trp'='Pre-sel\n(E-Trp)','E-DP'='Post-sel\n(E-DP)')
cat('Finalization start:',format(Sys.time(),tz='Asia/Shanghai'),'\n')
info<-read_result('sample_mapping'); sample_order<-info$sample
qc<-read_result('qc_and_TCR_assignment_audit')
pbtab<-read_result('pseudobulk_raw_UMI_counts');pb<-as.matrix(pbtab[,-1]);rownames(pb)<-pbtab$gene
stopifnot(identical(colnames(pb),sample_order),nrow(pb)>0L,all(pb>=0),all(pb==round(pb)))
info$group<-factor(info$group,levels=groups);info$donor<-factor(info$donor)
design<-model.matrix(~donor+group,data=info);rownames(design)<-info$sample
stopifnot(qr(design)$rank==ncol(design),nrow(design)-qr(design)$rank==2L)
save_csv(data.frame(sample=rownames(design),design,check.names=FALSE),'pseudobulk_design_matrix')
y<-DGEList(pb);keep<-filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y,method='TMM')
y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
save_csv(data.frame(sample=colnames(y),group=info$group,donor=info$donor,postQC_cells=qc$post_QC_cells[match(colnames(y),qc$sample)],total_raw_UMIs=colSums(pb),filtered_library_size=y$samples$lib.size,TMM_norm_factor=y$samples$norm.factors,genes_tested=sum(keep),n_independent_donors=length(unique(info$donor)),residual_df=2),'pseudobulk_library_summary')
contrasts<-list('E_Trp_vs_H_Trp'=c(0,0,1,0),'E_DP_vs_E_Trp'=c(0,0,-1,1),'E_DP_vs_H_Trp'=c(0,0,0,1))
de_summary<-list();lcpm<-cpm(y,log=TRUE,prior.count=2);de_list<-list()
for(nm in names(contrasts)){
 q<-glmQLFTest(fit,contrast=contrasts[[nm]]);z<-topTags(q,n=Inf,sort.by='PValue')$table;z$gene<-rownames(z)
 pairs<-switch(nm,E_Trp_vs_H_Trp=c('E-Trp','H-Trp'),E_DP_vs_E_Trp=c('E-DP','E-Trp'),E_DP_vs_H_Trp=c('E-DP','H-Trp'))
 for(d in levels(info$donor)){a<-which(info$donor==d & info$group==pairs[1]);b<-which(info$donor==d & info$group==pairs[2]);z[[paste0(d,'_logCPM_delta')]]<-lcpm[z$gene,a]-lcpm[z$gene,b]}
 z$both_donors_same_direction<-sign(z$Donor1_logCPM_delta)==sign(z$Donor2_logCPM_delta)
 z$inference_status<-'Exploratory paired edgeR QL; 2 independent donors; limited population inference'
 z<-z[,c('gene',setdiff(names(z),'gene'))];save_csv(z,paste0('pseudobulk_DE_',nm));de_list[[nm]]<-z
 de_summary[[nm]]<-data.frame(comparison=nm,genes_tested=nrow(z),FDR_lt_0.05=sum(z$FDR<.05),FDR_lt_0.05_absLFC_gt1=sum(z$FDR<.05 & abs(z$logFC)>1),FDR_lt_0.05_same_direction=sum(z$FDR<.05 & z$both_donors_same_direction),n_donors=2,residual_df=2)
 z$display<-ifelse(z$FDR<.05 & abs(z$logFC)>1,'FDR < 0.05 and |logFC| > 1','Other');z$neglogFDR<--log10(pmax(z$FDR,1e-300))
 p<-ggplot(z,aes(logFC,neglogFDR,color=display))+geom_point(size=.7,alpha=.65)+scale_color_manual(values=c('FDR < 0.05 and |logFC| > 1'='#b34e3d','Other'='grey65'))+theme_bw(base_size=11)+labs(x=paste(pairs[1],'vs',pairs[2],'log2 fold change'),y='-log10 FDR',title=paste('Donor-paired pseudobulk:',nm),subtitle='Raw UMI sums, TMM, edgeR QL ~ donor + group; exploratory with two donors')+theme(legend.position='bottom')
 save_plot(p,paste0('pseudobulk_volcano_',nm),8,6)
}
ds<-do.call(rbind,de_summary);save_csv(ds,'pseudobulk_DE_summary');print(ds)
requested<-c('TCF7','CCR7','IL7R','LEF1','TOX','HAVCR2','ENTPD1','LAG3')
gene_pb<-do.call(rbind,lapply(requested,function(g)data.frame(sample=sample_order,group=info$group,donor=info$donor,gene=g,raw_UMI=pb[g,],logCPM=if(g %in% rownames(lcpm))as.numeric(lcpm[g,]) else NA_real_)))
save_csv(gene_pb,'requested_eight_genes_pseudobulk_expression')

ts<-read_result('TCR_equal_cell_downsample_summary');td<-read_result('TCR_equal_cell_downsample_iterations')
stopifnot(length(unique(ts$n_downsample))==1L,unique(ts$n_downsample)==min(ts$eligible_TCR_cells),length(unique(ts$iterations))==1L,all(ts$iterations>=1000),length(unique(ts$seed))==1L)
for(i in seq_len(nrow(ts))){z<-td$Shannon[td$sample==ts$sample[i] & td$clonotype_definition==ts$clonotype_definition[i]];stopifnot(length(z)==ts$iterations[i],abs(median(z)-ts$downsample_median[i])<1e-12,abs(quantile(z,.025)-ts$sampling_lower_2.5pct[i])<1e-12,abs(quantile(z,.975)-ts$sampling_upper_97.5pct[i])<1e-12)}
ts$group<-factor(ts$group,levels=groups)
save_plot(ggplot(ts,aes(group,downsample_median,group=donor,color=donor))+geom_line()+geom_point(size=2)+geom_errorbar(aes(ymin=sampling_lower_2.5pct,ymax=sampling_upper_97.5pct),width=.08)+facet_wrap(~clonotype_definition)+scale_x_discrete(labels=axis_labels)+theme_bw(base_size=11)+labs(x='Experimental condition',y='Shannon diversity, natural log',title='Equal-cell TCR resampling by donor',subtitle=paste0(unique(ts$iterations),' draws; bars are sampling percentiles, not donor confidence intervals'),caption=unique(ts$eligibility)), 'TCR_equal_cell_downsample',9,5)
rs<-read_result('recomputed_signature_donor_summary');rs$group<-factor(rs$group,levels=groups)
save_plot(ggplot(rs,aes(group,mean_control_adjusted_score,group=donor,color=donor))+geom_line()+geom_point(size=1.7)+facet_wrap(~signature,scales='free_y',ncol=4)+scale_x_discrete(labels=axis_labels)+theme_bw(base_size=10)+labs(x='Experimental condition',y='Sample mean control-adjusted score',title='New QC signatures: each point is one donor-condition',subtitle='Historical recipe key 384, Seurat AddModuleScore; n = 2 donors; no cell-level P values'), 'recomputed_signature_donor_summary',13,10)
gs<-read_result('requested_eight_genes_donor_expression');gs$group<-factor(gs$group,levels=groups)
save_plot(ggplot(gs,aes(group,mean_log1p_normalized_RNA,group=donor,color=donor))+geom_line()+geom_point(size=2)+facet_wrap(~gene,scales='free_y',ncol=4)+scale_x_discrete(labels=axis_labels)+theme_bw(base_size=11)+labs(x='Experimental condition',y='Sample mean log1p normalized RNA',title='Requested genes: each point is one donor-condition',subtitle='n = 2 donors; descriptive expression, no cell-level P values'), 'requested_eight_genes_donor_expression',12,6)
ss<-read_result('requested_memory_exhaustion_donor_scores')
sl<-do.call(rbind,lapply(setdiff(names(ss),c('sample','group','donor','n_cells')),function(k)data.frame(ss[c('sample','group','donor')],signature=k,score=ss[[k]])));sl$group<-factor(sl$group,levels=groups)
save_plot(ggplot(sl,aes(group,score,group=donor,color=donor))+geom_line()+geom_point(size=2)+facet_wrap(~signature,scales='free_y',ncol=2)+scale_x_discrete(labels=axis_labels)+theme_bw(base_size=11)+labs(x='Experimental condition',y='Donor-condition score',title='Memory and exhaustion-associated scores',subtitle='Two donors; mean log expression and control-adjusted module scores; no significance stars'), 'requested_memory_exhaustion_donor_scores',10,7)
delta<-list()
for(k in unique(sl$signature))for(d in unique(sl$donor)){z<-sl[sl$signature==k & sl$donor==d,];z<-z[match(groups,z$group),];for(j in 1:2)delta[[paste(k,d,j)]]<-data.frame(signature=k,donor=d,comparison=paste(groups[j+1],'vs',groups[j]),paired_change=z$score[j+1]-z$score[j])}
save_csv(do.call(rbind,delta),'requested_score_paired_changes')
ls<-read_result('legacy_signature_donor_summary');ls$group<-factor(ls$group,levels=groups)
ls$source<-'Historical seurat_minimal_for_sharing.rds; old QC subset, not recalculated on new QC'
save_csv(ls,'legacy_signature_donor_summary')
save_plot(ggplot(ls,aes(group,mean_score,group=donor,color=donor))+geom_line()+geom_point()+facet_wrap(~signature,scales='free_y',ncol=4)+scale_x_discrete(labels=axis_labels)+theme_bw(base_size=9)+labs(x='Experimental condition',y='Sample mean of historical cell score',title='Historical signatures summarized by donor (old QC object)',subtitle='Two donor trajectories; no cell-level P values; use newly recomputed signatures for revised figures'), 'legacy_signature_donor_summary',13,10)

focus<-c('IL2RA','FOXP3','CTLA4','TIGIT','PDCD1','IKZF2','TCF7','CD28','CD27','HLA-DRA','HLA-DRB1','CD8A','GZMB','PRF1','CCL5',requested)
focus<-unique(intersect(focus,rownames(lcpm)))
heat<-function(genes,name,title){m<-lcpm[genes,,drop=FALSE];v<-t(scale(t(m)));v[!is.finite(v)]<-0;long<-data.frame(gene=rep(rownames(v),ncol(v)),sample=rep(colnames(v),each=nrow(v)),z=as.vector(v));long$gene<-factor(long$gene,levels=rev(genes));long$sample<-factor(long$sample,levels=sample_order);save_csv(long,paste0(name,'_source'));save_plot(ggplot(long,aes(sample,gene,fill=z))+geom_tile(color='white')+scale_fill_gradient2(low='#386b92',mid='white',high='#b5453e',name='Row z-score')+theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())+labs(x='Each column = one donor-condition',y=NULL,title=title,subtitle='TMM-normalized logCPM; each gene standardized across six pseudobulk libraries'),name,9,max(5,length(genes)*.24))}
heat(focus,'predefined_key_genes_pseudobulk_heatmap','Predefined genes: donor-level expression display')
chosen<-list();gene_union<-character()
for(nm in names(de_list)[1:2]){z<-de_list[[nm]];z<-z[z$FDR<.05 & abs(z$logFC)>1 & z$both_donors_same_direction,];z<-head(z,20);if(nrow(z)){chosen[[nm]]<-data.frame(comparison=nm,z[c('gene','logFC','FDR','Donor1_logCPM_delta','Donor2_logCPM_delta')]);gene_union<-union(gene_union,z$gene)}}
if(length(gene_union)){save_csv(do.call(rbind,chosen),'G_pseudobulk_selected_candidates');heat(gene_union,'G_pseudobulk_selected_candidates_heatmap','G candidate replacement: selected from paired pseudobulk DE')}

stopifnot(nrow(gs)==nrow(info)*8,nrow(rs)==nrow(info)*12,nrow(ss)==nrow(info),all(ss$n_cells==qc$post_QC_cells[match(ss$sample,qc$sample)]),sum(qc$TCR_assigned_cells)+sum(qc$missing_TCR_as_single)==sum(qc$post_QC_cells))
stopifnot(all(read_result('TCR_donor_trends')$direction_preserved))
cat('CHECKS PASS: pseudobulk integer matrix; paired design rank 4 / df2; TCR median/percentiles match all',nrow(td),'draws; all 8 trend directions retained;',sum(qc$post_QC_cells),'QC cells;',sum(qc$TCR_assigned_cells),'assigned TCR cells;',sum(qc$missing_TCR_as_single),'missing excluded;',nrow(gs),'gene summaries;',nrow(rs),'score summaries.\n')
writeLines(capture.output(sessionInfo()),file.path(out,'main_finalization_R_sessionInfo.txt'))
cat('FINALIZATION_COMPLETED:',format(Sys.time(),tz='Asia/Shanghai'),'\n');print(sessionInfo())
