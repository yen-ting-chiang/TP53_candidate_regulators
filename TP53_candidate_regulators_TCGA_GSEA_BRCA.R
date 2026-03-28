################################################################################
# Target Genes TCGA GSEA Analysis — limma-voom 版本 (多 Dataset 批次分析)
#
# 分析流程:
#   1. 準備 MSigDB gene sets (過濾 C2:CGP 僅包含 p53 相關)
#   2. 進入 Dataset 迴圈 (e.g., brca_tcga_gdc, coad_tcga_gdc, ...)
#      - 讀取 read counts (data_mrna_seq_read_counts.txt)
#      - 判定 TP53 mutation status (data_mutations.txt)
#      - TP53 Mutant vs Wild Type 差異表現與 GSEA 分析
#      - 針對 144 個 candidate genes 執行 limma 差異表現分析 與 GSEA (僅 TP53 WT)
#      - 統整各個 Gene Set 的 GSEA 結果
#
# 輸出至各自的 [dataset]_results_limma_[gene_list]/ 資料夾
################################################################################

# ---- 0. 環境設定 ----
set.seed(42)

library(limma)
library(edgeR)
library(fgsea)
library(msigdbr)
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)

# 設定正確的 working directory
setwd("C:/Users/danny/Documents/R_project/TP53_candidate_regulators")

# ---- 0.1 設定要執行的 Datasets ----
# 您可以在這裡自由註解或新增想要跑的 datasets 目錄名稱
datasets_to_run <- c(
  "brca_tcga_gdc",
  "coad_tcga_gdc",
  "difg_tcga_gdc",
  "gbm_tcga_gdc",
  "hcc_tcga_gdc",
  "luad_tcga_gdc",
  "paad_tcga_gdc"
)
# 如果只要跑少數幾個，可以改成如下：
# datasets_to_run <- c("brca_tcga_gdc")

# 讀取 Candidate Genes
target_genes_file <- "p53_candidate_regulators_144_genes.csv"
if (!file.exists(target_genes_file)) {
  stop(sprintf("找不到目標基因名單檔案: %s", target_genes_file))
}
candidate_genes_df <- read.csv(target_genes_file, stringsAsFactors = FALSE)
candidate_genes <- candidate_genes_df$Gene_name

# 取得檔名主體作為資料夾命名後綴
target_genes_file_base <- sub("\\.[^.]+$", "", target_genes_file)

# ---- 1. 準備 MSigDB 與自定義 Gene Sets (僅需執行一次) ----
cat("\n=== 準備 GSEA gene sets (全域環境) ===\n")

# Hallmark 僅取 HALLMARK_P53_PATHWAY
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
hallmark_list <- hallmark_list["HALLMARK_P53_PATHWAY"]
cat("Hallmark gene sets 取用數:", length(hallmark_list), "\n")

# C2:CGP 僅取特定五個
c2cgp_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
c2cgp_list_all <- split(c2cgp_df$gene_symbol, c2cgp_df$gs_name)
target_c2cgp_names <- c("FISCHER_DIRECT_P53_TARGETS_META_ANALYSIS",
                        "ONGUSAHA_TP53_TARGETS",
                        "PEREZ_TP53_TARGETS",
                        "WU_APOPTOSIS_BY_CDKN1A_VIA_TP53",
                        "SCIAN_CELL_CYCLE_TARGETS_OF_TP53_AND_TP73_DN")
c2cgp_list <- c2cgp_list_all[target_c2cgp_names]
cat("C2:CGP gene sets 取用數:", length(c2cgp_list), "\n")

# 自定義 gene set "p53_direct_target_genes_Fischer_2017"
fischer_file <- "p53_direct_target_genes_Fischer_2017.csv"
if (!file.exists(fischer_file)) {
  stop(sprintf("找不到自定義 gene set 檔案: %s", fischer_file))
}
fischer_df <- read.csv(fischer_file, stringsAsFactors = FALSE)
# 去除空白的 gene symbols
fischer_genes <- fischer_df$Gene_Symbol[fischer_df$Gene_Symbol != ""]
custom_list <- list(p53_direct_target_genes_Fischer_2017 = fischer_genes)
cat("Custom gene sets 取用數:", length(custom_list), "\n")

# 將所有 gene sets 合併成一個 unified_gene_sets
unified_gene_sets <- c(hallmark_list, c2cgp_list, custom_list)
cat("統整後總 Gene Sets 數量:", length(unified_gene_sets), "\n")

# ---- 2. 定義 limma-voom + GSEA 分析函式 (全域環境) ----

run_limma_continuous <- function(counts, target_expr,
                                 target_entrez_id,
                                 gene_id_mapping,
                                 design_mat = NULL,
                                 coef_name = NULL) {
  # 移除 target gene 自身
  counts_no_target <- counts[rownames(counts) != target_entrez_id, ]
  
  # 建立 DGEList 並計算 normalization factors
  dge <- DGEList(counts = counts_no_target)
  dge <- calcNormFactors(dge, method = "TMM")
  
  # 過濾低表現基因 (使用 edgeR::filterByExpr)
  if (is.null(design_mat)) {
    design_mat <- model.matrix(~ target_expr)
    coef_target <- 2  # target coefficient (第 2 欄)
  } else {
    if (is.null(coef_name)) {
      stop("使用自定義 design matrix 時必須指定 coef_name")
    }
    coef_target <- coef_name
  }
  
  keep <- filterByExpr(dge, design = design_mat)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat(sprintf("    低表現基因過濾 (filterByExpr): 保留 %d / %d 基因\n",
              sum(keep), length(keep)))
  
  # voom 轉換
  v <- voom(dge, design = design_mat, plot = FALSE)
  
  # Linear modeling
  fit <- lmFit(v, design_mat)
  fit <- eBayes(fit)
  
  # 取得結果
  res <- topTable(fit, coef = coef_target, number = Inf, sort.by = "none")
  
  # 加入 Entrez Gene Id
  res$Entrez_Gene_Id <- rownames(res)
  
  # 加入 Gene Symbol
  res$Gene_Symbol <- gene_id_mapping$Gene_Symbol[
    match(res$Entrez_Gene_Id, gene_id_mapping$Entrez_Gene_Id)
  ]
  
  return(res)
}

run_gsea <- function(limma_res, gene_set_list, collection_name,
                     analysis_label, out_dir) {
  if (length(gene_set_list) == 0) {
     cat(sprintf("    [%s] 由於 gene sets 為空，跳過此 GSEA (collection: %s)\n", analysis_label, collection_name))
     return(NULL)
  }

  # 建立 ranked gene list (使用 moderated t-statistic)
  ranked <- limma_res[!is.na(limma_res$Gene_Symbol) & !is.na(limma_res$t), ]
  
  # 若有重複 gene symbol, 取 |t| 最大者
  ranked <- ranked[order(abs(ranked$t), decreasing = TRUE), ]
  ranked <- ranked[!duplicated(ranked$Gene_Symbol), ]
  
  # 建立 named vector
  ranks <- ranked$t
  names(ranks) <- ranked$Gene_Symbol
  ranks <- sort(ranks, decreasing = TRUE)
  
  cat(sprintf("    Ranked genes: %d\n", length(ranks)))
  
  # 執行 fgsea
  suppressWarnings({
    gsea_res <- fgsea(
      pathways = gene_set_list,
      stats = ranks,
      minSize = 15,
      maxSize = 1500
    )
  })
  
  # 整理結果
  gsea_res <- gsea_res[order(gsea_res$pval), ]
  
  # 將 leadingEdge 從 list 轉為字串
  gsea_out <- as.data.frame(gsea_res)
  gsea_out$leadingEdge <- sapply(gsea_out$leadingEdge, function(x) {
    paste(x, collapse = "; ")
  })
  
  # 輸出 GSEA 結果
  result_file <- file.path(out_dir,
                           sprintf("GSEA_%s_%s_results.csv", analysis_label, collection_name))
  write.csv(gsea_out, result_file, row.names = FALSE)
  cat(sprintf("    結果已儲存: %s\n", result_file))
  
  # 輸出 leading edge genes (僅 FDR < 0.25 的 pathways)
  sig_pathways <- gsea_res[gsea_res$padj < 0.25 & !is.na(gsea_res$padj), ]
  if (nrow(sig_pathways) > 0) {
    le_list <- lapply(seq_len(nrow(sig_pathways)), function(i) {
      data.frame(
        Pathway = sig_pathways$pathway[i],
        NES = sig_pathways$NES[i],
        padj = sig_pathways$padj[i],
        LeadingEdge_Gene = sig_pathways$leadingEdge[[i]],
        stringsAsFactors = FALSE
      )
    })
    le_df <- do.call(rbind, le_list)
    le_file <- file.path(out_dir,
                         sprintf("GSEA_%s_%s_leading_edge.csv", analysis_label, collection_name))
    write.csv(le_df, le_file, row.names = FALSE)
    cat(sprintf("    Leading edge genes 已儲存: %s (FDR<0.25 pathways: %d)\n",
                le_file, nrow(sig_pathways)))
  } else {
    cat(sprintf("    無 FDR<0.25 的顯著 pathways (collection: %s)\n",
                collection_name))
  }
  
  # 統計摘要
  cat(sprintf("    [%s - %s] Total: %d | FDR<0.05: %d | FDR<0.25: %d\n",
              analysis_label, collection_name,
              nrow(gsea_res),
              sum(gsea_res$padj < 0.05, na.rm = TRUE),
              sum(gsea_res$padj < 0.25, na.rm = TRUE)))
  
  return(gsea_res)
}

# ==============================================================================
# ---- 3. 循序處理多個 Dataset 迴圈 ----
# ==============================================================================

for (dataset_dir in datasets_to_run) {
  cat("\n")
  cat("****************************************************************\n")
  cat(sprintf("*** 開始分析 Dataset: %s ***\n", dataset_dir))
  cat("****************************************************************\n")
  
  # 取得 dataset prefix (e.g., "brca" from "brca_tcga_gdc")
  dataset_prefix <- sub("_tcga_gdc$", "", dataset_dir)
  
  # 根據使用者需求建立動態總輸出資料夾
  output_dir <- sprintf("%s_results_limma_%s", dataset_prefix, target_genes_file_base)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat(sprintf(">> 結果將輸出至目錄: %s/\n", output_dir))
  
  # ---- 3.1 讀取 Read Counts ----
  cat("\n=== 讀取 read counts 資料 ===\n")
  counts_file <- file.path(dataset_dir, "data_mrna_seq_read_counts.txt")
  if (!file.exists(counts_file)) {
    cat(sprintf("!! 警告: 找不到檔案 %s，跳過此 dataset。\n", counts_file))
    next
  }
  
  counts_raw <- fread(counts_file, header = TRUE, sep = "\t", data.table = FALSE)
  
  entrez_col_idx <- which(colnames(counts_raw) == "Entrez_Gene_Id")
  if (length(entrez_col_idx) == 0) {
    cat(sprintf("!! 警告: 找不到 Entrez_Gene_Id 欄位，跳過 %s\n", dataset_dir))
    next
  }
  
  entrez_col <- counts_raw[[entrez_col_idx]]
  annotation_cols <- seq_len(entrez_col_idx)
  count_cols <- counts_raw[, -annotation_cols]
  
  valid_rows <- !is.na(entrez_col)
  entrez_col <- entrez_col[valid_rows]
  count_cols <- count_cols[valid_rows, ]
  
  counts_mat <- rowsum(as.matrix(count_cols), group = entrez_col)
  counts_mat <- round(counts_mat)
  storage.mode(counts_mat) <- "integer"
  cat("Read counts 維度:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")
  
  # ---- 3.2 Entrez ID 轉換為 Gene Symbol ----
  entrez_ids <- rownames(counts_mat)
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = entrez_ids,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")
  id_mapping <- data.frame(
    Entrez_Gene_Id = entrez_ids,
    Gene_Symbol = as.character(gene_symbols),
    stringsAsFactors = FALSE
  )
  
  # ---- 3.3 判定 TP53 Mutation Status ----
  cat("\n=== 判定 TP53 mutation status ===\n")
  mut_data_file <- file.path(dataset_dir, "data_mutations.txt")
  if (!file.exists(mut_data_file)) {
    cat(sprintf("!! 警告: 找不到 %s，未能分析此 dataset TP53 狀態。\n", mut_data_file))
    next
  }
  
  mut_data <- fread(mut_data_file, header = TRUE, sep = "\t", data.table = FALSE, skip = "Hugo_Symbol")
  tp53_mut <- mut_data[mut_data$Hugo_Symbol == "TP53", ]
  
  protein_altering_classes <- c(
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
    "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del",
    "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site"
  )
  
  tp53_protein_altering <- tp53_mut[tp53_mut$Variant_Classification %in% protein_altering_classes, ]
  tp53_mutant_barcodes_raw <- unique(tp53_protein_altering$Tumor_Sample_Barcode)
  tp53_mutant_samples <- unique(substr(tp53_mutant_barcodes_raw, 1, 16))
  
  cat("TP53 protein-altering mutations 數量:", nrow(tp53_protein_altering), "\n")
  cat("TP53 mutant samples 數量:", length(tp53_mutant_samples), "\n")
  
  sample_ids <- colnames(counts_mat)
  sample_ids_trunc <- substr(sample_ids, 1, 16)
  tp53_status <- ifelse(sample_ids_trunc %in% tp53_mutant_samples, "Mutant", "WT")
  names(tp53_status) <- sample_ids
  
  cat("TP53 Mutant:", sum(tp53_status == "Mutant"), "/ WT:", sum(tp53_status == "WT"), "\n")
  if (sum(tp53_status == "Mutant") == 0) {
    warning(sprintf("[%s] 無任何 TP53 mutant sample 匹配! 請檢查 barcode", dataset_dir))
  }
  
  # ---- 3.4 分析: TP53 Mutant vs Wild Type 差異表現與 GSEA 分析 ----
  cat("\n========================================\n")
  cat("=== 分析 5: TP53 Mutant vs Wild Type ===\n")
  cat("========================================\n")
  
  out_dir_mt_vs_wt <- file.path(output_dir, sprintf("%s_TP53_mt_vs_wt", dataset_prefix))
  if (!dir.exists(out_dir_mt_vs_wt)) dir.create(out_dir_mt_vs_wt, recursive = TRUE)
  
  tp53_factor_group <- factor(tp53_status, levels = c("WT", "Mutant"))
  print(table(tp53_factor_group))
  
  if (sum(tp53_status == "Mutant") >= 3 && sum(tp53_status == "WT") >= 3) {
    design_mt_vs_wt <- model.matrix(~ tp53_factor_group)
    colnames(design_mt_vs_wt) <- c("Intercept", "Mutant_vs_WT")
    
    dge_mt_wt <- DGEList(counts = counts_mat)
    dge_mt_wt <- calcNormFactors(dge_mt_wt, method = "TMM")
    
    keep_mt_wt <- filterByExpr(dge_mt_wt, design = design_mt_vs_wt)
    dge_mt_wt <- dge_mt_wt[keep_mt_wt, , keep.lib.sizes = FALSE]
    
    v_mt_wt <- voom(dge_mt_wt, design = design_mt_vs_wt, plot = FALSE)
    fit_mt_wt <- lmFit(v_mt_wt, design_mt_vs_wt)
    fit_mt_wt <- eBayes(fit_mt_wt)
    
    res_mt_wt <- topTable(fit_mt_wt, coef = "Mutant_vs_WT", number = Inf, sort.by = "none")
    res_mt_wt$Entrez_Gene_Id <- rownames(res_mt_wt)
    res_mt_wt$Gene_Symbol <- id_mapping$Gene_Symbol[match(res_mt_wt$Entrez_Gene_Id, id_mapping$Entrez_Gene_Id)]
    
    limma_file_mt_wt <- file.path(out_dir_mt_vs_wt, "limma_TP53_mutant_vs_wt_results.csv")
    write.csv(res_mt_wt, limma_file_mt_wt, row.names = FALSE)
    
    cat("\n  --- GSEA: TP53 Mutant vs WT - Unified Gene Sets ---\n")
    gsea_mt_wt_unified <- run_gsea(res_mt_wt, unified_gene_sets, "Unified_GeneSets",
                                   "TP53_mutant_vs_wt", out_dir_mt_vs_wt)
  } else {
    cat(sprintf("WARNING: [%s] TP53 Mutant 或 WT 分組樣本數量不足 (<3)，跳過此分析\n", dataset_dir))
  }
  
  # ---- 3.5 迴圈分析開始: 針對每個候選基因 (僅 TP53 Wild Type 樣本) ----
  cat("\n================================================================\n")
  cat(sprintf("=== 針對 %d 個候選基因進行分析 (僅 TP53 Wild Type) ===\n", length(candidate_genes)))
  cat("================================================================\n")
  
  wt_samples <- sample_ids[tp53_status == "WT"]
  if (length(wt_samples) < 20) {
    cat(sprintf("!! ERROR: [%s] TP53 Wild Type samples 數量不足 (<20)，跳過所有的連續變數 GSEA 分析\n", dataset_dir))
    next
  }
  
  counts_wt <- counts_mat[, wt_samples]
  
  dge_wt <- DGEList(counts = counts_wt)
  dge_wt <- calcNormFactors(dge_wt, method = "TMM")
  logcpm_wt <- cpm(dge_wt, log = TRUE, prior.count = 1)
  
  for (target_gene in candidate_genes) {
    cat(sprintf("\n-- 處理目標基因: %s (WT) --\n", target_gene))
    
    gene_out_dir <- file.path(output_dir, target_gene)
    if (!dir.exists(gene_out_dir)) dir.create(gene_out_dir, recursive = TRUE)
    
    target_entrez <- mapIds(org.Hs.eg.db,
                           keys = target_gene,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")
                           
    if (is.na(target_entrez) || !target_entrez %in% rownames(counts_wt)) {
      cat(sprintf("    !! 警告: %s (Entrez ID %s) 不在 read counts 資料中! 跳過該基因。\n", target_gene, target_entrez))
      next
    }
    
    target_logcpm <- logcpm_wt[target_entrez, ]
    names(target_logcpm) <- colnames(counts_wt)
    target_scaled <- scale(target_logcpm)[, 1]
    
    # limma WT samples
    cat(sprintf("    > limma-voom: %s (WT samples)\n", target_gene))
    limma_res_wt <- run_limma_continuous(
      counts = counts_wt,
      target_expr = target_scaled,
      target_entrez_id = target_entrez,
      gene_id_mapping = id_mapping
    )
    
    limma_file_wt <- file.path(gene_out_dir, sprintf("limma_%s_wt_results.csv", target_gene))
    write.csv(limma_res_wt, limma_file_wt, row.names = FALSE)
    
    # GSEA WT samples (Unified Gene Sets)
    cat(sprintf("    > GSEA: %s (WT) - Unified Gene Sets\n", target_gene))
    gsea_wt_unified <- run_gsea(limma_res_wt, unified_gene_sets, "Unified_GeneSets",
                                sprintf("%s_wt", target_gene), gene_out_dir)
  }
  
  # ---- 3.6 結果統整: 針對各個 Gene Set 彙整所有候選基因的 GSEA 結果 ----
  cat("\n================================================================\n")
  cat("=== 開始統整各個 Gene Set 在所有候選基因的 GSEA 結果 ===\n")
  cat("================================================================\n")
  
  summary_dir <- file.path(output_dir, "GSEA_Summary_by_GeneSet")
  if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)
  
  all_pathways <- names(unified_gene_sets)
  
  for (pw in all_pathways) {
    pw_summary_list <- list()
    
    for (target_gene in candidate_genes) {
      res_file <- file.path(output_dir, target_gene, sprintf("GSEA_%s_wt_Unified_GeneSets_results.csv", target_gene))
      
      if (file.exists(res_file)) {
        res_df <- read.csv(res_file, stringsAsFactors = FALSE)
        pw_row <- res_df[res_df$pathway == pw, ]
        
        if (nrow(pw_row) > 0) {
          extracted_info <- data.frame(
            Target_Gene = target_gene,
            pval = pw_row$pval,
            padj = pw_row$padj,
            log2err = pw_row$log2err,
            ES = pw_row$ES,
            NES = pw_row$NES,
            size = pw_row$size,
            leadingEdge = pw_row$leadingEdge,
            stringsAsFactors = FALSE
          )
          pw_summary_list[[length(pw_summary_list) + 1]] <- extracted_info
        }
      }
    }
    
    if (length(pw_summary_list) > 0) {
      pw_summary_df <- do.call(rbind, pw_summary_list)
      safe_pw_name <- gsub("[^A-Za-z0-9_.-]", "_", pw)
      out_file <- file.path(summary_dir, sprintf("GSEA_summary_%s.csv", safe_pw_name))
      write.csv(pw_summary_df, out_file, row.names = FALSE)
      cat(sprintf("  -> 已輸出 %s 的統整結果至: %s\n", pw, out_file))
    }
  }

  
  # ---- 3.7 分析摘要 (針對個別 Dataset) ----
  cat("\n================================================================\n")
  cat(sprintf("=== [%s] 分析完成摘要 ===\n", dataset_dir))
  cat("================================================================\n")
  cat(sprintf("總 samples: %d\n", ncol(counts_mat)))
  cat(sprintf("TP53 Mutant: %d | TP53 Wild Type: %d\n",
              sum(tp53_status == "Mutant"), sum(tp53_status == "WT")))
  cat(sprintf("統整後 GSEA Gene Sets 數量: %d sets\n", length(unified_gene_sets)))
  cat(sprintf("所有結果已輸出至: %s/\n", output_dir))
  cat("================================================================\n")
}

cat("\n\n################################################################\n")
cat("### 所有選定之 Dataset 分析完成! ###\n")
cat("################################################################\n")

# ==============================================================================
# ---- 4. 尋找共同顯著基因: NES > 1.5 & padj < 0.05 交集分析 (獨立執行區塊) ----
# ==============================================================================

cat("\n\n################################################################\n")
cat("### 開始針對各 Dataset 進行交集分析 ###\n")
cat("################################################################\n")

for (dataset_dir in datasets_to_run) {
  dataset_prefix <- sub("_tcga_gdc$", "", dataset_dir)
  output_dir <- sprintf("%s_results_limma_%s", dataset_prefix, target_genes_file_base)
  summary_dir <- file.path(output_dir, "GSEA_Summary_by_GeneSet")
  
  if (!dir.exists(summary_dir)) {
    cat(sprintf("\n!!警告: 找不到 %s 的統整資料夾 (%s)，跳過此 dataset 的交集分析。\n", dataset_dir, summary_dir))
    next
  }
  
  cat(sprintf("\n=== [%s] 顯著候選基因交集分析 ===\n", dataset_dir))
  
  target_summary_files <- c(
    "p53_direct_target_genes_Fischer_2017" = "GSEA_summary_p53_direct_target_genes_Fischer_2017.csv",
    "FISCHER_DIRECT" = "GSEA_summary_FISCHER_DIRECT_P53_TARGETS_META_ANALYSIS.csv",
    "HALLMARK_P53" = "GSEA_summary_HALLMARK_P53_PATHWAY.csv"
  )
  
  sig_genes_list <- list()
  
  for (set_name in names(target_summary_files)) {
    file_path <- file.path(summary_dir, target_summary_files[[set_name]])
    if (file.exists(file_path)) {
      df <- read.csv(file_path, stringsAsFactors = FALSE)
      
      # 確保轉為數值後篩選: NES > 1.5 且 padj < 0.05
      valid_rows <- !is.na(df$NES) & !is.na(df$padj) & 
                    (as.numeric(df$NES) > 1.5) & 
                    (as.numeric(df$padj) < 0.05)
      sig_df <- df[valid_rows, ]
      
      sig_genes_list[[set_name]] <- sig_df$Target_Gene
      cat(sprintf("  [%s] 符合條件之基因數: %d\n", set_name, length(sig_genes_list[[set_name]])))
    } else {
      cat(sprintf("  !!警告: 找不到檔案 %s，請確認該 Gene Set 有成功統整。\n", target_summary_files[[set_name]]))
      sig_genes_list[[set_name]] <- character(0)
    }
  }
  
  # 取交集
  if (length(sig_genes_list) == 3) {
    intersect_genes <- Reduce(intersect, sig_genes_list)
    cat(sprintf("  -> 三個 Gene Sets 皆符合條件的交集基因數: %d\n", length(intersect_genes)))
    
    intersect_df <- data.frame(Target_Gene = intersect_genes, stringsAsFactors = FALSE)
    intersect_out_file <- file.path(summary_dir, "Intersection_Significant_Target_Genes.csv")
    write.csv(intersect_df, intersect_out_file, row.names = FALSE)
    cat(sprintf("  -> 交集結果已輸出至: %s\n", intersect_out_file))
  } else {
    cat("  未能成功載入三個指定的 Gene Sets，略過交集分析。\n")
  }
}
cat("\n=== 所有 Dataset 交集分析完成 ===\n")

# ==============================================================================
# ---- 5. 跨 Dataset 統整: Intersection_Significant_Target_Genes 交集結果彙整 獨立執行區塊 ----
# ==============================================================================
cat("\n\n################################################################\n")
cat("### 開始進行跨 Dataset 的顯著基因交集結果彙整 ###\n")
cat("################################################################\n")

cross_dataset_list <- list()

for (dataset_dir in datasets_to_run) {
  dataset_prefix <- sub("_tcga_gdc$", "", dataset_dir)
  output_dir <- sprintf("%s_results_limma_%s", dataset_prefix, target_genes_file_base)
  summary_dir <- file.path(output_dir, "GSEA_Summary_by_GeneSet")
  intersect_file <- file.path(summary_dir, "Intersection_Significant_Target_Genes.csv")
  
  if (file.exists(intersect_file)) {
    intersect_df <- read.csv(intersect_file, stringsAsFactors = FALSE)
    if (nrow(intersect_df) > 0) {
      for (gene in intersect_df$Target_Gene) {
        if (is.null(cross_dataset_list[[gene]])) {
          cross_dataset_list[[gene]] <- c(toupper(dataset_prefix))
        } else {
          cross_dataset_list[[gene]] <- c(cross_dataset_list[[gene]], toupper(dataset_prefix))
        }
      }
    }
  } else {
    cat(sprintf("  !!警告: 找不到 %s 的交集結果檔案，跳過此 dataset。\n", dataset_dir))
  }
}

if (length(cross_dataset_list) > 0) {
  cross_summary_df <- data.frame(
    Target_Gene = names(cross_dataset_list),
    Total_Count = sapply(cross_dataset_list, length),
    Datasets = sapply(cross_dataset_list, function(x) paste(x, collapse = ", ")),
    stringsAsFactors = FALSE
  )
  
  # 加入 Category 欄位 (對應回 target_genes_file)
  if ("Category" %in% colnames(candidate_genes_df)) {
    category_map <- setNames(candidate_genes_df$Category, candidate_genes_df$Gene_name)
    cross_summary_df$Category <- category_map[cross_summary_df$Target_Gene]
    
    # 重新排列欄位順序
    cross_summary_df <- cross_summary_df[, c("Target_Gene", "Category", "Total_Count", "Datasets")]
  }
  
  # 依據出現次數由大到小排序
  cross_summary_df <- cross_summary_df[order(cross_summary_df$Total_Count, decreasing = TRUE), ]
  
  cross_out_file <- sprintf("Cross_Dataset_Intersection_Summary_%s.csv", target_genes_file_base)
  write.csv(cross_summary_df, cross_out_file, row.names = FALSE)
  cat(sprintf("\n-> 跨 Dataset 統整結果已輸出至: %s\n", cross_out_file))
} else {
  cat("\n-> 沒有在任何 Dataset 中找到符合條件的顯著基因。\n")
}
cat("\n=== 跨 Dataset 統整完成 ===\n")

# ==============================================================================
# ---- 6. 跨 Dataset 統整: TP53 Mutant vs WT GSEA 結果彙整 獨立執行區塊 ----
# ==============================================================================
cat("\n\n################################################################\n")
cat("### 開始進行跨 Dataset 的 TP53 Mutant vs WT GSEA 結果彙整 ###\n")
cat("################################################################\n")

gsea_mt_wt_list <- list()

for (dataset_dir in datasets_to_run) {
  dataset_prefix <- sub("_tcga_gdc$", "", dataset_dir)
  output_dir <- sprintf("%s_results_limma_%s", dataset_prefix, target_genes_file_base)
  out_dir_mt_vs_wt <- file.path(output_dir, sprintf("%s_TP53_mt_vs_wt", dataset_prefix))
  
  gsea_res_file <- file.path(out_dir_mt_vs_wt, "GSEA_TP53_mutant_vs_wt_Unified_GeneSets_results.csv")
  
  if (file.exists(gsea_res_file)) {
    df <- read.csv(gsea_res_file, stringsAsFactors = FALSE)
    if (nrow(df) > 0) {
      req_cols <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")
      avail_cols <- intersect(req_cols, colnames(df))
      
      if (length(avail_cols) > 0) {
        df_sub <- df[, avail_cols, drop = FALSE]
        dataset_name_upper <- toupper(dataset_prefix)
        df_sub$Dataset <- dataset_name_upper
        
        # 重新排列欄位，將 Dataset 放在最前面
        df_sub <- df_sub[, c("Dataset", avail_cols)]
        
        gsea_mt_wt_list[[dataset_name_upper]] <- df_sub
      }
    }
  } else {
    cat(sprintf("  !!警告: 找不到 %s 的 TP53_mt_vs_wt GSEA 結果檔案，跳過此 dataset。\n", dataset_dir))
  }
}

if (length(gsea_mt_wt_list) > 0) {
  cross_gsea_mt_wt_df <- do.call(rbind, gsea_mt_wt_list)
  
  # 依據 pathway 排序，方便同一個 pathway 在各 dataset 中比較
  cross_gsea_mt_wt_df <- cross_gsea_mt_wt_df[order(cross_gsea_mt_wt_df$pathway, cross_gsea_mt_wt_df$Dataset), ]
  
  cross_gsea_out_file <- sprintf("Cross_Dataset_TP53_mt_vs_wt_GSEA_Summary_%s.csv", target_genes_file_base)
  write.csv(cross_gsea_mt_wt_df, cross_gsea_out_file, row.names = FALSE)
  cat(sprintf("\n-> 跨 Dataset TP53 Mutant vs WT GSEA 統整結果已輸出至: %s\n", cross_gsea_out_file))
} else {
  cat("\n-> 沒有在任何 Dataset 中找到符合的 TP53 Mutant vs WT GSEA 結果。\n")
}
cat("\n=== 跨 Dataset TP53 Mutant vs WT GSEA 統整完成 ===\n")
