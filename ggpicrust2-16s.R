pkgs <- c("phyloseq", "ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "BiocManager", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

devtools::install_github("cafferychen777/ggpicrust2")

library(ggpicrust2)
library(tidyverse)
library(GGally)
library(ggprism)
library(patchwork)
library(ggh4x)
library(openxlsx)

data("metacyc_abundance")
data("metadata")

setwd("G:/熊猫宏基因组课题/16S/picrust2fmr79sZ2A8YxxFjHjKt8Uw_result/picrust2fmr79sZ2A8YxxFjHjKt8Uw_result/")
"G:\熊猫宏基因组课题\16S\picrust2fmr79sZ2A8YxxFjHjKt8Uw_result\picrust2fmr79sZ2A8YxxFjHjKt8Uw_result\MetaCyc.pathway.abundance.xls"

metacyc_abundance1 <- read.delim("MetaCyc.pathway.abundance.xls", sep = "\t",header = TRUE)
metadata1 <- read.delim("G:/熊猫宏基因组课题/16S/metadata2.txt", sep = "\t",header = TRUE, fileEncoding = "UTF-16")


# 假设 metadata 中有一列名为 sample_id
# 获取要筛选的列名
selected_columns <- metadata1$sample_name

# 筛选 metacyc_abundance 中的列
filtered_abundance <- metacyc_abundance1[, colnames(metacyc_abundance1) %in% selected_columns]

# 查看筛选后的数据框
head(filtered_abundance)
filtered_abundance$pathway <- metacyc_abundance2$pathway

#差异丰度分析
metacyc_daa_results_df <- pathway_daa(abundance = filtered_abundance %>% 
                                        column_to_rownames("pathway"), metadata = metadata1, 
                                      group = "tooth", daa_method = "LinDA",
                                      select = NULL, p.adjust = "BH", reference = NULL)
#多模型结果比较
methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = filtered_abundance %>% 
                column_to_rownames("pathway"), metadata = metadata1, group = "tooth", daa_method = method)
})

method_names <- c("ALDEx2_Welch's t test","ALDEx2_Wilcoxon rank test","DESeq2", "edgeR")
# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = method_names)
#kegg通路注释
data("kegg_abundance")
data("metadata")
"G:\熊猫宏基因组课题\16S\picrust2fmr79sZ2A8YxxFjHjKt8Uw_result\picrust2fmr79sZ2A8YxxFjHjKt8Uw_result\Sample.KO.abundance.xls"
kegg_abundance2 <- read.table("Sample.KO.abundance.xls", sep = "\t", header = TRUE, fill = TRUE, quote = "")
#metacyc_abundance1 <- read.delim("MetaCyc.pathway.abundance.xls", sep = "\t",header = TRUE)
kegg_abundance1 <- kegg_abundance1[, colnames(kegg_abundance1) %in% selected_columns]
metadata = metadata, group = "tooth", daa_method = "LinDA")
ko_abundance1<-kegg_abundance1
data("ko_abundance")
data("metadata")
ko_abundance1$NAME <-kegg_abundance2$function.
ko_abundance1 <- ko_abundance1[, c("NAME", setdiff(names(ko_abundance1), "NAME"))]
kegg_abundance <- ko2kegg_abundance(data = ko_abundance)
kegg_abundance1 <- ko2kegg_abundance(data = ko_abundance1)

daa_results_df <- pathway_daa(kegg_abundance1, metadata = metadata1,
                              group = "tooth", daa_method = "LinDA")
daa_annotated_results_df <- pathway_annotation(pathway = "KO",
                                               daa_results_df = daa_results_df, ko_to_kegg = TRUE)

p <- pathway_errorbar(abundance = kegg_abundance1,
                      daa_results_df = daa_annotated_results_df,
                      Group = metadata1$tooth,
                      ko_to_kegg = TRUE,
                      p_values_threshold = 0.05,
                      order = "pathway_class",
                      select = NULL,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "pathway_name")


data("metacyc_abundance")
data("metadata")
metacyc_daa_results_df <- pathway_daa(abundance = filtered_abundance %>% 
                                        column_to_rownames("pathway"), metadata = metadata1,
                                      group = "tooth", daa_method = "LinDA")

metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc",
                                                       daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)

p1 <- pathway_errorbar(abundance = filtered_abundance %>% 
                        column_to_rownames("pathway"),
                      daa_results_df = metacyc_daa_annotated_results_df,
                      Group = metadata1$tooth,
                      ko_to_kegg = FALSE,
                      p_values_threshold = 0.05,
                      order = "group",
                      select = NULL,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "description")

# Perform differential abundance analysis
metacyc_daa_results_df <- pathway_daa(
  abundance = filtered_abundance %>% column_to_rownames("pathway"),
  metadata = metadata1,
  group = "tooth",
  daa_method = "LinDA"
)

# Annotate the results

annotated_metacyc_daa_results_df <- pathway_annotation(
  pathway = "MetaCyc",
  daa_results_df = metacyc_daa_results_df,
  ko_to_kegg = FALSE
)

feature_with_p_0.05 <- metacyc_daa_results_df %>% 
  filter(p_adjust < 0.05)

library(ggplot2)
library(pheatmap)

# 假设 filtered_abundance 和其他数据框已经准备好了

# 创建热图
p <- pathway_heatmap(
  abundance = filtered_abundance %>% 
    right_join(
      annotated_metacyc_daa_results_df %>% select(all_of(c("feature","description"))),
      by = c("pathway" = "feature")
    ) %>% 
    filter(pathway %in% feature_with_p_0.05$feature) %>% 
    select(-"pathway") %>% 
    column_to_rownames("description"),
  metadata = metadata1, 
  group = "tooth"
)

# 美化热图，去掉 x 轴
p + 
  theme(
    text = element_text(size = 8),  # 设置字体大小
    axis.text.x = element_blank(),  # 去掉 x 轴标签
    axis.ticks.x = element_blank(),  # 去掉 x 轴刻度线
    panel.grid.major = element_line(colour = "grey", size = 0.5),  # 添加主网格线
    panel.grid.minor = element_blank()  # 去掉次网格线
  ) +
  geom_tile(color = "white")  # 每个单元格的边缘线

filtered_abundance <- filtered_abundance[, c("pathway", setdiff(names(filtered_abundance), "pathway"))]
pathway_pca(abundance = metacyc_abundance1  %>% column_to_rownames("pathway"),
            metadata = metadata1, group = "tooth")
