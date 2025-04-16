library(dplyr)
library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
setwd('/public/home/b20213040320/AAA_highadaption/plastic/lung')

lung<-readRDS('/public/home/b20213040320/AAA_highadaption/lung_adult.rds')
DefaultAssay(lung)<-'RNA'


Idents(lung)<-lung$time

wk1_p <- FindMarkers(lung, ident.1 = "1wk", ident.2 = "plain", min.pct = 0.25,logfc.threshold = 0.25)
wk2_p <- FindMarkers(lung, ident.1 = "2wk", ident.2 = "plain", min.pct = 0.25,logfc.threshold = 0.25)
wk3_p <- FindMarkers(lung, ident.1 = "3wk", ident.2 = "plain", min.pct = 0.25,logfc.threshold = 0.25)
mon8_p <- FindMarkers(lung, ident.1 = "8 mon", ident.2 = "plain", min.pct = 0.25,logfc.threshold = 0.25)

wk1_ts <- FindMarkers(lung, ident.1 = "Tibetan", ident.2 = "1wk", min.pct = 0.25,logfc.threshold = 0.25)
wk2_ts <- FindMarkers(lung, ident.1 = "Tibetan", ident.2 = "2wk", min.pct = 0.25,logfc.threshold = 0.25)
wk3_ts <- FindMarkers(lung, ident.1 = "Tibetan", ident.2 = "3wk", min.pct = 0.25,logfc.threshold = 0.25)
mon8_ts <- FindMarkers(lung, ident.1 = "Tibetan", ident.2 = "8 mon", min.pct = 0.25,logfc.threshold = 0.25)

wk1_p <-subset(wk1_p,wk1_p$p_val_adj<0.05)
wk2_p <-subset(wk2_p,wk2_p$p_val_adj<0.05)
wk3_p <-subset(wk3_p,wk3_p$p_val_adj<0.05)
mon8_p <-subset(mon8_p,mon8_p$p_val_adj<0.05)
wk1_ts <-subset(wk1_ts,wk1_ts$p_val_adj<0.05)
wk2_ts <-subset(wk2_ts,wk2_ts$p_val_adj<0.05)
wk3_ts <-subset(wk3_ts,wk3_ts$p_val_adj<0.05)
mon8_ts <-subset(mon8_ts,mon8_ts$p_val_adj<0.05)

wk1_p_up <-subset(wk1_p,wk1_p$avg_log2FC>0)
wk1_p_down <-subset(wk1_p,wk1_p$avg_log2FC<0)
wk2_p_up <-subset(wk2_p,wk2_p$avg_log2FC>0)
wk2_p_down <-subset(wk2_p,wk2_p$avg_log2FC<0)
wk3_p_up <-subset(wk3_p,wk3_p$avg_log2FC>0)
wk3_p_down <-subset(wk3_p,wk3_p$avg_log2FC<0)
mon8_p_up <-subset(mon8_p,mon8_p$avg_log2FC>0)
mon8_p_down <-subset(mon8_p,mon8_p$avg_log2FC<0)
wk1_ts_up <-subset(wk1_ts,wk1_ts$avg_log2FC>0)
wk1_ts_down <-subset(wk1_ts,wk1_ts$avg_log2FC<0)
wk2_ts_up <-subset(wk2_ts,wk2_ts$avg_log2FC>0)
wk2_ts_down <-subset(wk2_ts,wk2_ts$avg_log2FC<0)
wk3_ts_up <-subset(wk3_ts,wk3_ts$avg_log2FC>0)
wk3_ts_down <-subset(wk3_ts,wk3_ts$avg_log2FC<0)
mon8_ts_up <-subset(mon8_ts,mon8_ts$avg_log2FC>0)
mon8_ts_down <-subset(mon8_ts,mon8_ts$avg_log2FC<0)

##reversing
wk1_rev1<-as.data.frame(intersect(rownames(wk1_p_up),rownames(wk1_ts_down))) 
wk1_rev2<-as.data.frame(intersect(rownames(wk1_p_down),rownames(wk1_ts_up)))
colnames(wk1_rev1)<-'gene'
colnames(wk1_rev2)<-'gene'
wk1_rev<-rbind(wk1_rev1,wk1_rev2)
wk2_rev1<-as.data.frame(intersect(rownames(wk2_p_up),rownames(wk2_ts_down))) 
wk2_rev2<-as.data.frame(intersect(rownames(wk2_p_down),rownames(wk2_ts_up)))
colnames(wk2_rev1)<-'gene'
colnames(wk2_rev2)<-'gene'
wk2_rev<-rbind(wk2_rev1,wk2_rev2)
wk3_rev1<-as.data.frame(intersect(rownames(wk3_p_up),rownames(wk3_ts_down))) 
wk3_rev2<-as.data.frame(intersect(rownames(wk3_p_down),rownames(wk3_ts_up)))
colnames(wk3_rev1)<-'gene'
colnames(wk3_rev2)<-'gene'
wk3_rev<-rbind(wk3_rev1,wk2_rev2)
mon8_rev1<-as.data.frame(intersect(rownames(mon8_p_up),rownames(mon8_ts_down))) 
mon8_rev2<-as.data.frame(intersect(rownames(mon8_p_down),rownames(mon8_ts_up)))
colnames(mon8_rev1)<-'gene'
colnames(mon8_rev2)<-'gene'
mon8_rev<-rbind(mon8_rev1,mon8_rev2)


##reinforcing
wk1_rein1<-as.data.frame(intersect(rownames(wk1_p_up),rownames(wk1_ts_up))) 
wk1_rein2<-as.data.frame(intersect(rownames(wk1_p_down),rownames(wk1_ts_down)))
colnames(wk1_rein1)<-'gene'
colnames(wk1_rein2)<-'gene'
wk1_rein<-rbind(wk1_rein1,wk1_rein2)
wk2_rein1<-as.data.frame(intersect(rownames(wk2_p_up),rownames(wk2_ts_up))) 
wk2_rein2<-as.data.frame(intersect(rownames(wk2_p_down),rownames(wk2_ts_down)))
colnames(wk2_rein1)<-'gene'
colnames(wk2_rein2)<-'gene'
wk2_rein<-rbind(wk2_rein1,wk2_rein2)
wk3_rein1<-as.data.frame(intersect(rownames(wk3_p_up),rownames(wk3_ts_up))) 
wk3_rein2<-as.data.frame(intersect(rownames(wk3_p_down),rownames(wk3_ts_down)))
colnames(wk3_rein1)<-'gene'
colnames(wk3_rein2)<-'gene'
wk3_rein<-rbind(wk3_rein1,wk2_rein2)
mon8_rein1<-as.data.frame(intersect(rownames(mon8_p_up),rownames(mon8_ts_up))) 
mon8_rein2<-as.data.frame(intersect(rownames(mon8_p_down),rownames(mon8_ts_down)))
colnames(mon8_rein1)<-'gene'
colnames(mon8_rein2)<-'gene'
mon8_rein<-rbind(mon8_rein1,mon8_rein2)

parametric_bootstrap <- function(Lo, SEo, Lp, SEp, La, SEa, n_boot=10000) {
  results <- replicate(n_boot, {
    sim_o <- rnorm(1, Lo, SEo)
    sim_p <- rnorm(1, Lp, SEp)
    sim_a <- rnorm(1, La, SEa)
    PC <- sim_p - sim_o
    GC <- sim_a - sim_p
    if ((PC > 0 & GC > 0) | (PC < 0 & GC < 0)) {
      "reinforcing"
    } else if ((PC > 0 & GC < 0) | (PC < 0 & GC > 0)) {
      "reversing"
    } else {
      "neither"
    }
  })
  
  tibble(
    RI_count = sum(results == "reinforcing"),
    RV_count = sum(results == "reversing"),
    classification = case_when(
      RV_count >= 0.95*n_boot ~ "reversing",
      RI_count >= 0.95*n_boot ~ "reinforcing",
      TRUE ~ "neither"
    )
  )
}

wk1_p$genename<-rownames(wk1_p)
wk1_ts$genename<-rownames(wk1_ts)
genename<-rbind(wk1_p,wk1_ts)
genename<-genename$gene

metadata <-lung@meta.data
df <-unique(genename) 
data<-t(LayerData(lung, layer = "counts", features = df))
metadata$cell<-rownames(metadata)
genes <- colnames(data)
gene_groups <- split(genes, ceiling(seq_along(genes)/500))  
results <- lapply(gene_groups, function(g) {
  as.data.frame(data[, g]) %>%
    tibble::rownames_to_column("cell") %>%
    pivot_longer(-cell, names_to = "gene", values_to = "expression") %>%
    left_join(metadata, by = "cell") %>%
    group_by(gene, time) %>%
    summarise(
      mean_expr = mean(expression),
      se_expr = sd(expression)/sqrt(n()), 
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = time,
      values_from = c(mean_expr, se_expr),
      names_sep = "_"
    )
})
stage_stats <- bind_rows(results)

wk1 <- stage_stats %>%
  mutate(
    boot_result = pmap(
      list(
        mean_expr_plain, se_expr_plain,
        mean_expr_1wk, se_expr_1wk,
        mean_expr_Tibetan, se_expr_Tibetan
      ),
      ~ parametric_bootstrap(..1, ..2, ..3, ..4, ..5, ..6)
    )
  ) %>%
  unnest(boot_result)

wk2 <- stage_stats %>%
  mutate(
    boot_result = pmap(
      list(
        mean_expr_plain, se_expr_plain,
        mean_expr_2wk, se_expr_2wk,
        mean_expr_Tibetan, se_expr_Tibetan
      ),
      ~ parametric_bootstrap(..1, ..2, ..3, ..4, ..5, ..6)
    )
  ) %>%
  unnest(boot_result)

wk3 <- stage_stats %>%
  mutate(
    boot_result = pmap(
      list(
        mean_expr_plain, se_expr_plain,
        mean_expr_3wk, se_expr_3wk,
        mean_expr_Tibetan, se_expr_Tibetan
      ),
      ~ parametric_bootstrap(..1, ..2, ..3, ..4, ..5, ..6)
    )
  ) %>%
  unnest(boot_result)

colnames(stage_stats)<-gsub('8 mon','mon8',colnames(stage_stats))
mon8 <- stage_stats %>%
  mutate(
    boot_result = pmap(
      list(
        mean_expr_plain, se_expr_plain,
        mean_expr_mon8, se_expr_mon8,
        mean_expr_Tibetan, se_expr_Tibetan
      ),
      ~ parametric_bootstrap(..1, ..2, ..3, ..4, ..5, ..6)
    )
  ) %>%
  unnest(boot_result)

save(wk1, wk2,wk3,mon8, wk1_rev, wk2_rev,wk3_rev,mon8_rev,file = "lung_all_compare.RData")

#####celltype
wk1_gene<-rbind(wk1_p,ts_wk1)
wk2_gene<-rbind(wk2_p,ts_wk2)
wk3_gene<-rbind(wk3_p,ts_wk3)
mon8_gene<-rbind(mon8_p,ts_mon8)
load("D:/360MoveData/Users/Admin/Desktop/PNAS revise/gene/new/lung_celltype_deg.RData")
lung<-readRDS('/public/home/b20213040320/AAA_highadaption/lung_adult.rds')

DefaultAssay(lung)<-'RNA'
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
cell<-list()
metadata<-list()
gene<-list()
data<-list()
gene_groups<-list()
results<-list()
stage_stats<-list()
wk1_result<-list()
wk2_result<-list()
wk3_result<-list()
mon8_result<-list()
wk1_deg<-list()
wk2_deg<-list()
wk3_deg<-list()
mon8_deg<-list()
for (i in 1:length(df)){
  tryCatch({
    cell[[i]]<-subset(lung,idents = df[[i]])
    metadata[[i]] <-cell[[i]]@meta.data
    wk1_deg[[i]]<-subset(wk1_gene,wk1_gene$class%in%df[[i]])
    wk2_deg[[i]]<-subset(wk2_gene,wk2_gene$class%in%df[[i]])
    wk3_deg[[i]]<-subset(wk3_gene,wk3_gene$class%in%df[[i]])
    mon8_deg[[i]]<-subset(mon8_gene,mon8_gene$class%in%df[[i]])
    gene[[i]]<-unique(wk1_deg[[i]]$gene) 
    data[[i]]<-t(LayerData(cell[[i]], layer = "counts", features = wk1_deg[[i]]$gene))
    metadata[[i]]$cell<-rownames(metadata[[i]])
    gene_groups[[i]] <- split(gene[[i]], ceiling(seq_along(gene[[i]])/500))
    results[[i]] <- lapply(gene_groups[[i]], function(g) {
      as.data.frame(data[[i]][, g]) %>%
        tibble::rownames_to_column("cell") %>%
        pivot_longer(-cell, names_to = "gene", values_to = "expression") %>%
        left_join(metadata[[i]], by = "cell") %>%
        group_by(gene, time) %>%
        summarise(
          mean_expr = mean(expression),
          se_expr = sd(expression)/sqrt(n()), 
          .groups = "drop"
        ) %>%
        pivot_wider(
          names_from = time,
          values_from = c(mean_expr, se_expr),
          names_sep = "_"
        )
    })
    stage_stats[[i]] <- bind_rows(results[[i]])
    colnames(stage_stats[[i]])<-gsub('8 mon','mon8',colnames(stage_stats[[i]]))
    wk1_result[[i]] <- stage_stats[[i]] %>%
      mutate(
        boot_result = pmap(
          list(
            mean_expr_plain, se_expr_plain,
            mean_expr_1wk, se_expr_1wk,
            mean_expr_Tibetan, se_expr_Tibetan
          ),
          ~ parametric_bootstrap(..1, ..2, ..3, ..4, ..5, ..6)
        )
      ) %>%
      unnest(boot_result)
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
save(wk1_result,file = 'lung_wk1_result.rds')

DefaultAssay(lung)<-'RNA'
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
cell<-list()
metadata<-list()
gene<-list()
data<-list()
gene_groups<-list()
results<-list()
stage_stats<-list()
wk1_result<-list()
wk2_result<-list()
wk3_result<-list()
mon8_result<-list()
wk1_deg<-list()
wk2_deg<-list()
wk3_deg<-list()
mon8_deg<-list()
for (i in 1:length(df)){
  tryCatch({
    cell[[i]]<-subset(lung,idents = df[[i]])
    metadata[[i]] <-cell[[i]]@meta.data
    wk1_deg[[i]]<-subset(wk1_gene,wk1_gene$class%in%df[[i]])
    wk2_deg[[i]]<-subset(wk2_gene,wk2_gene$class%in%df[[i]])
    wk3_deg[[i]]<-subset(wk3_gene,wk3_gene$class%in%df[[i]])
    mon8_deg[[i]]<-subset(mon8_gene,mon8_gene$class%in%df[[i]])
    gene[[i]]<-unique(wk2_deg[[i]]$gene) 
    data[[i]]<-t(LayerData(cell[[i]], layer = "counts", features = wk2_deg[[i]]$gene))
    metadata[[i]]$cell<-rownames(metadata[[i]])
    gene_groups[[i]] <- split(gene[[i]], ceiling(seq_along(gene[[i]])/500))
    results[[i]] <- lapply(gene_groups[[i]], function(g) {
      as.data.frame(data[[i]][, g]) %>%
        tibble::rownames_to_column("cell") %>%
        pivot_longer(-cell, names_to = "gene", values_to = "expression") %>%
        left_join(metadata[[i]], by = "cell") %>%
        group_by(gene, time) %>%
        summarise(
          mean_expr = mean(expression),
          se_expr = sd(expression)/sqrt(n()), 
          .groups = "drop"
        ) %>%
        pivot_wider(
          names_from = time,
          values_from = c(mean_expr, se_expr),
          names_sep = "_"
        )
    })
    stage_stats[[i]] <- bind_rows(results[[i]])
    colnames(stage_stats[[i]])<-gsub('8 mon','mon8',colnames(stage_stats[[i]]))
    wk2_result[[i]] <- stage_stats[[i]] %>%
      mutate(
        boot_result = pmap(
          list(
            mean_expr_plain, se_expr_plain,
            mean_expr_2wk, se_expr_2wk,
            mean_expr_Tibetan, se_expr_Tibetan
          ),
          ~ parametric_bootstrap(..1, ..2, ..3, ..4, ..5, ..6)
        )
      ) %>%
      unnest(boot_result)
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
save(wk2_result,file = 'lung_wk2_result.rds')


DefaultAssay(lung)<-'RNA'
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
cell<-list()
metadata<-list()
gene<-list()
data<-list()
gene_groups<-list()
results<-list()
stage_stats<-list()
wk1_result<-list()
wk2_result<-list()
wk3_result<-list()
mon8_result<-list()
wk1_deg<-list()
wk2_deg<-list()
wk3_deg<-list()
mon8_deg<-list()
for (i in 1:length(df)){
  tryCatch({
    cell[[i]]<-subset(lung,idents = df[[i]])
    metadata[[i]] <-cell[[i]]@meta.data
    wk1_deg[[i]]<-subset(wk1_gene,wk1_gene$class%in%df[[i]])
    wk2_deg[[i]]<-subset(wk2_gene,wk2_gene$class%in%df[[i]])
    wk3_deg[[i]]<-subset(wk3_gene,wk3_gene$class%in%df[[i]])
    mon8_deg[[i]]<-subset(mon8_gene,mon8_gene$class%in%df[[i]])
    gene[[i]]<-unique(wk3_deg[[i]]$gene) 
    data[[i]]<-t(LayerData(cell[[i]], layer = "counts", features = wk3_deg[[i]]$gene))
    metadata[[i]]$cell<-rownames(metadata[[i]])
    gene_groups[[i]] <- split(gene[[i]], ceiling(seq_along(gene[[i]])/500))
    results[[i]] <- lapply(gene_groups[[i]], function(g) {
      as.data.frame(data[[i]][, g]) %>%
        tibble::rownames_to_column("cell") %>%
        pivot_longer(-cell, names_to = "gene", values_to = "expression") %>%
        left_join(metadata[[i]], by = "cell") %>%
        group_by(gene, time) %>%
        summarise(
          mean_expr = mean(expression),
          se_expr = sd(expression)/sqrt(n()), 
          .groups = "drop"
        ) %>%
        pivot_wider(
          names_from = time,
          values_from = c(mean_expr, se_expr),
          names_sep = "_"
        )
    })
    stage_stats[[i]] <- bind_rows(results[[i]])
    colnames(stage_stats[[i]])<-gsub('8 mon','mon8',colnames(stage_stats[[i]]))
    wk3_result[[i]] <- stage_stats[[i]] %>%
      mutate(
        boot_result = pmap(
          list(
            mean_expr_plain, se_expr_plain,
            mean_expr_3wk, se_expr_3wk,
            mean_expr_Tibetan, se_expr_Tibetan
          ),
          ~ parametric_bootstrap(..1, ..2, ..3, ..4, ..5, ..6)
        )
      ) %>%
      unnest(boot_result)
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
save(wk3_result,file = 'lung_wk3_result.rds')

DefaultAssay(lung)<-'RNA'
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
cell<-list()
metadata<-list()
gene<-list()
data<-list()
gene_groups<-list()
results<-list()
stage_stats<-list()
wk1_result<-list()
wk2_result<-list()
wk3_result<-list()
mon8_result<-list()
wk1_deg<-list()
wk2_deg<-list()
wk3_deg<-list()
mon8_deg<-list()
for (i in 1:length(df)){
  tryCatch({
    cell[[i]]<-subset(lung,idents = df[[i]])
    metadata[[i]] <-cell[[i]]@meta.data
    wk1_deg[[i]]<-subset(wk1_gene,wk1_gene$class%in%df[[i]])
    wk2_deg[[i]]<-subset(wk2_gene,wk2_gene$class%in%df[[i]])
    wk3_deg[[i]]<-subset(wk3_gene,wk3_gene$class%in%df[[i]])
    mon8_deg[[i]]<-subset(mon8_gene,mon8_gene$class%in%df[[i]])
    gene[[i]]<-unique(mon8_deg[[i]]$gene) 
    data[[i]]<-t(LayerData(cell[[i]], layer = "counts", features = mon8_deg[[i]]$gene))
    metadata[[i]]$cell<-rownames(metadata[[i]])
    gene_groups[[i]] <- split(gene[[i]], ceiling(seq_along(gene[[i]])/500))
    results[[i]] <- lapply(gene_groups[[i]], function(g) {
      as.data.frame(data[[i]][, g]) %>%
        tibble::rownames_to_column("cell") %>%
        pivot_longer(-cell, names_to = "gene", values_to = "expression") %>%
        left_join(metadata[[i]], by = "cell") %>%
        group_by(gene, time) %>%
        summarise(
          mean_expr = mean(expression),
          se_expr = sd(expression)/sqrt(n()), 
          .groups = "drop"
        ) %>%
        pivot_wider(
          names_from = time,
          values_from = c(mean_expr, se_expr),
          names_sep = "_"
        )
    })
    stage_stats[[i]] <- bind_rows(results[[i]])
    colnames(stage_stats[[i]])<-gsub('8 mon','mon8',colnames(stage_stats[[i]]))
    mon8_result[[i]] <- stage_stats[[i]] %>%
      mutate(
        boot_result = pmap(
          list(
            mean_expr_plain, se_expr_plain,
            mean_expr_mon8, se_expr_mon8,
            mean_expr_Tibetan, se_expr_Tibetan
          ),
          ~ parametric_bootstrap(..1, ..2, ..3, ..4, ..5, ..6)
        )
      ) %>%
      unnest(boot_result)
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
save(mon8_result,file = 'lung_mon8_result.rds')



####figure
#####reversing
sub<-list()
int<-list()
per<-list()
rev<-list()
all<-data.frame()
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
for (i in 1:length(df)){
  tryCatch({
    sub[[i]]<-subset(wk1_result[[i]],wk1_result[[i]]$classification%in%'reversing')
    int[[i]]<-length(intersect(sub[[i]]$gene,reversing_wk1[[i]]$gene))%>%as.data.frame()
    colnames(int[[i]])<-'parametric_bootstrap'
    rev[[i]]<-length(reversing_wk1[[i]]$gene)%>%as.data.frame()
    colnames(rev[[i]])<-'all'
    per[[i]]<-cbind(int[[i]],rev[[i]])
    per[[i]]$False_Positive<-per[[i]]$all-per[[i]]$parametric_bootstrap
    per[[i]]$cell<-df[[i]]
    all<-rbind(all,per[[i]])
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
all$time<-'wk1'
all1<-all

sub<-list()
int<-list()
per<-list()
rev<-list()
all<-data.frame()
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
for (i in 1:length(df)){
  tryCatch({
    sub[[i]]<-subset(wk2_result[[i]],wk2_result[[i]]$classification%in%'reversing')
    int[[i]]<-length(intersect(sub[[i]]$gene,reversing_wk2[[i]]$gene))%>%as.data.frame()
    colnames(int[[i]])<-'parametric_bootstrap'
    rev[[i]]<-length(reversing_wk2[[i]]$gene)%>%as.data.frame()
    colnames(rev[[i]])<-'all'
    per[[i]]<-cbind(int[[i]],rev[[i]])
    per[[i]]$False_Positive<-per[[i]]$all-per[[i]]$parametric_bootstrap
    per[[i]]$cell<-df[[i]]
    all<-rbind(all,per[[i]])
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
all$time<-'wk2'
all2<-all

sub<-list()
int<-list()
per<-list()
rev<-list()
all<-data.frame()
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
for (i in 1:length(df)){
  tryCatch({
    sub[[i]]<-subset(wk3_result[[i]],wk3_result[[i]]$classification%in%'reversing')
    int[[i]]<-length(intersect(sub[[i]]$gene,reversing_wk3[[i]]$gene))%>%as.data.frame()
    colnames(int[[i]])<-'parametric_bootstrap'
    rev[[i]]<-length(reversing_wk3[[i]]$gene)%>%as.data.frame()
    colnames(rev[[i]])<-'all'
    per[[i]]<-cbind(int[[i]],rev[[i]])
    per[[i]]$False_Positive<-per[[i]]$all-per[[i]]$parametric_bootstrap
    per[[i]]$cell<-df[[i]]
    all<-rbind(all,per[[i]])
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
all$time<-'wk3'
all3<-all

sub<-list()
int<-list()
per<-list()
rev<-list()
all<-data.frame()
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
for (i in 1:length(df)){
  tryCatch({
    sub[[i]]<-subset(mon8_result[[i]],mon8_result[[i]]$classification%in%'reversing')
    int[[i]]<-length(intersect(sub[[i]]$gene,reversing_mon8[[i]]$gene))%>%as.data.frame()
    colnames(int[[i]])<-'parametric_bootstrap'
    rev[[i]]<-length(reversing_mon8[[i]]$gene)%>%as.data.frame()
    colnames(rev[[i]])<-'all'
    per[[i]]<-cbind(int[[i]],rev[[i]])
    per[[i]]$False_Positive<-per[[i]]$all-per[[i]]$parametric_bootstrap
    per[[i]]$cell<-df[[i]]
    all<-rbind(all,per[[i]])
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
all$time<-'mon8'
all4<-all
reversing<-rbind(all1,all2,all3,all4)


#####reinforcing
sub<-list()
int<-list()
per<-list()
rev<-list()
all<-data.frame()
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
for (i in 1:length(df)){
  tryCatch({
    sub[[i]]<-subset(wk1_result[[i]],wk1_result[[i]]$classification%in%'reinforcing')
    int[[i]]<-length(intersect(sub[[i]]$gene,reinforcing_wk1[[i]]$gene))%>%as.data.frame()
    colnames(int[[i]])<-'parametric_bootstrap'
    rev[[i]]<-length(reinforcing_wk1[[i]]$gene)%>%as.data.frame()
    colnames(rev[[i]])<-'all'
    per[[i]]<-cbind(int[[i]],rev[[i]])
    per[[i]]$False_Positive<-per[[i]]$all-per[[i]]$parametric_bootstrap
    per[[i]]$cell<-df[[i]]
    all<-rbind(all,per[[i]])
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
all$time<-'wk1'
all1<-all

sub<-list()
int<-list()
per<-list()
rev<-list()
all<-data.frame()
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
for (i in 1:length(df)){
  tryCatch({
    sub[[i]]<-subset(wk2_result[[i]],wk2_result[[i]]$classification%in%'reinforcing')
    int[[i]]<-length(intersect(sub[[i]]$gene,reinforcing_wk2[[i]]$gene))%>%as.data.frame()
    colnames(int[[i]])<-'parametric_bootstrap'
    rev[[i]]<-length(reinforcing_wk2[[i]]$gene)%>%as.data.frame()
    colnames(rev[[i]])<-'all'
    per[[i]]<-cbind(int[[i]],rev[[i]])
    per[[i]]$False_Positive<-per[[i]]$all-per[[i]]$parametric_bootstrap
    per[[i]]$cell<-df[[i]]
    all<-rbind(all,per[[i]])
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
all$time<-'wk2'
all2<-all

sub<-list()
int<-list()
per<-list()
rev<-list()
all<-data.frame()
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
for (i in 1:length(df)){
  tryCatch({
    sub[[i]]<-subset(wk3_result[[i]],wk3_result[[i]]$classification%in%'reinforcing')
    int[[i]]<-length(intersect(sub[[i]]$gene,reinforcing_wk3[[i]]$gene))%>%as.data.frame()
    colnames(int[[i]])<-'parametric_bootstrap'
    rev[[i]]<-length(reinforcing_wk3[[i]]$gene)%>%as.data.frame()
    colnames(rev[[i]])<-'all'
    per[[i]]<-cbind(int[[i]],rev[[i]])
    per[[i]]$False_Positive<-per[[i]]$all-per[[i]]$parametric_bootstrap
    per[[i]]$cell<-df[[i]]
    all<-rbind(all,per[[i]])
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
all$time<-'wk3'
all3<-all

sub<-list()
int<-list()
per<-list()
rev<-list()
all<-data.frame()
Idents(lung)<-lung$cell
df<-as.data.frame(table(lung$cell))
df<-as.character(df$Var1)
for (i in 1:length(df)){
  tryCatch({
    sub[[i]]<-subset(mon8_result[[i]],mon8_result[[i]]$classification%in%'reinforcing')
    int[[i]]<-length(intersect(sub[[i]]$gene,reinforcing_mon8[[i]]$gene))%>%as.data.frame()
    colnames(int[[i]])<-'parametric_bootstrap'
    rev[[i]]<-length(reinforcing_mon8[[i]]$gene)%>%as.data.frame()
    colnames(rev[[i]])<-'all'
    per[[i]]<-cbind(int[[i]],rev[[i]])
    per[[i]]$False_Positive<-per[[i]]$all-per[[i]]$parametric_bootstrap
    per[[i]]$cell<-df[[i]]
    all<-rbind(all,per[[i]])
  }, error = function(e) {
    warning(paste("Warning: Iteration", i, "caused an error and was skipped."))
  })
}
all$time<-'mon8'
all4<-all
reinforcing<-rbind(all1,all2,all3,all4)

reinforcing$group<-'reinforcing'
reversing$group<-'reversing'
all<-rbind(reversing,reinforcing)
all$percentage<-all$parametric_bootstrap/all$all
all[is.na(all)] <- 0

all$percentage<-NULL
all$False_Positive<-NULL


library(tidyr)
library(dplyr)

wk1<-subset(all,time=='mon8')
df_long <- wk1 %>%
  pivot_longer(
    cols = c(parametric_bootstrap, all),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(
    cell = factor(cell, levels = unique(cell))  
  )

df_long <- df_long %>%
  mutate(
    cell_num = as.numeric(cell),
    x_pos = ifelse(type == "parametric_bootstrap", cell_num - 0.2, cell_num + 0.2)
  )

library(ggplot2)

ggplot(df_long, aes(x = x_pos, y = value, fill = group)) +
  geom_col(width = 0.3, position = position_stack()) +
  scale_x_continuous(
    breaks = unique(df_long$cell_num),
    labels = levels(df_long$cell)
  ) +
  labs(x = "Cell", y = "Value", fill = "Group") +
  coord_flip() +
  theme_classic()+scale_fill_manual(values = c('#9bd9ff','#ff999a'))

