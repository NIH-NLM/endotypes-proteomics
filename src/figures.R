# Paths -- cohorts/ is committed and persists, data/run_artifacts/ is regenerable.
# See src/paths.R, which is the only place those locations are written down.
.f    <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.here <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else "src"
source(file.path(.here, "paths.R"))

suppressMessages({library(WGCNA); library(ComplexHeatmap); library(circlize)})
options(stringsAsFactors=FALSE); set.seed(42)
w<-readRDS(art("wgcna_A.rds")); e<-readRDS(art("eigengenes_A.rds"))
sig<-readRDS(art("sig_modules_A.rds")); en<-readRDS(art("endotypes_A.rds"))
pr<-readRDS(art("projection_healthy.rds"))
X<-w$X; m<-w$meta; mods<-w$mods; ME<-e$ME; endo<-en$endo; ifn<-pr$ifn
ann <- m[,c("Disease_activity","SLEDAI_2K","C3_level","Sm_status","Ro_60_status","Age_group")]
keep <- union(ifn, sig)
sel <- unlist(lapply(keep, function(k){ g<-colnames(X)[mods==k]
  g[order(-abs(cor(X[,g,drop=FALSE], ME[[paste0("ME",k)]])))][1:min(14,length(g))] }))
ordm <- keep[order(sapply(keep,function(k) sum(mods==k)))]
Z <- scale(as.matrix(X[,sel]))
png(art("fig1_blocks.png"), width=max(1600,58*length(sel)), height=1700, res=150)
draw(Heatmap(Z, name="z-score", col=colorRamp2(c(-2,0,2),c("#2166AC","white","#B2182B")),
  row_split=endo, column_split=factor(mods[match(sel,colnames(X))], levels=ordm),
  cluster_rows=TRUE, cluster_row_slices=FALSE, cluster_columns=TRUE,
  clustering_method_rows="ward.D2", clustering_method_columns="ward.D2",
  show_row_dend=TRUE, row_dend_width=unit(25,"mm"), show_row_names=FALSE,
  column_names_side="top", column_names_gp=gpar(fontsize=7), column_dend_side="top",
  column_title_gp=gpar(fontsize=8), column_title_rot=90,
  left_annotation=rowAnnotation(df=ann, annotation_name_gp=gpar(fontsize=7)),
  row_title="endotype",
  column_title=sprintf("cohort A (n=%d) -- clinically associated modules; interferon module = %s",nrow(X),ifn)))
invisible(dev.off()); cat("fig1_blocks.png --",length(sel),"proteins\n")

# the interferon module alone, patients ordered by its eigenprotein
g <- colnames(X)[mods==ifn]; ei <- ME[[paste0("ME",ifn)]]
o <- order(-ei)
png(art("fig2_interferon_module.png"), width=1500, height=1700, res=150)
draw(Heatmap(scale(as.matrix(X[o,g])), name="z-score",
  col=colorRamp2(c(-2,0,2),c("#2166AC","white","#B2182B")),
  cluster_rows=FALSE, cluster_columns=TRUE, clustering_method_columns="ward.D2",
  show_row_names=FALSE, column_names_side="top", column_names_gp=gpar(fontsize=9),
  left_annotation=rowAnnotation(df=m[o,c("Sm_status","Ro_60_status","dsDNA_status","C3_level","SLEDAI_2K")],
                                annotation_name_gp=gpar(fontsize=8)),
  row_title="patients, ordered by interferon eigenprotein (high at top)",
  column_title=sprintf("the interferon module '%s' -- %d proteins",ifn,length(g))))
invisible(dev.off()); cat("fig2_interferon_module.png\n")

# healthy volunteers projected
P<-pr$projection; tab<-pr$table
top <- head(tab$module[order(tab$hv_mean)], 6); top <- union(top, ifn)
df <- do.call(rbind, lapply(top, function(k){ s<-P[[k]]$sle
  rbind(data.frame(m=k,g="SLE",z=(s-mean(s))/sd(s)),
        data.frame(m=k,g="healthy",z=(P[[k]]$new-mean(s))/sd(s))) }))
df$m <- factor(df$m, levels=top)
png(art("fig3_healthy_projection.png"), width=1800, height=1100, res=150)
par(mar=c(9,4.5,3,1))
boxplot(z~g+m, data=df, las=2, col=c("#2166AC","#B2182B"), xlab="",
        ylab="eigenprotein, SLE SD units", cex.axis=0.75,
        main="healthy volunteers projected onto SLE-defined modules")
abline(h=0, lty=2); legend("topright", c("healthy","SLE"), fill=c("#2166AC","#B2182B"), bty="n")
invisible(dev.off()); cat("fig3_healthy_projection.png\n")
