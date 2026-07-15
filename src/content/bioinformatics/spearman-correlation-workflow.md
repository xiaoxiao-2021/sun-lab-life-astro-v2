---
title: "Spearman correlation workflow"
description: "基因表达量与免疫评分之间的 Spearman 相关性分析流程。"
date: 2026-07-13
updated: 2026-07-13
category: "Correlation"
tags: ["Spearman", "ssGSEA", "R", "Bioinformatics"]
status: "Building"
draft: false
---


## 1. 这篇笔记解决什么问题？

这篇笔记用于整理基因表达量与免疫评分之间的 Spearman 相关性分析流程。

适用场景：

- 基因表达量 vs 免疫评分
- 差异基因 vs ssGSEA score
- 分组内相关性分析
- 四象限相关性筛选

## 2. 输入数据

需要准备三类数据：

| 数据 | 说明 |
|---|---|
| expression matrix | log2 normalized counts |
| score table | ssGSEA score |
| sample metadata | sample / group 信息 |

## 3. 核心思路

先保证样本名完全匹配，然后对每个基因和目标 score 计算 Spearman 相关性。

## 4. 基础代码

```r
cor.test(
  x = gene_expression,
  y = score,
  method = "spearman"
)
```

## 5. 输出结果

推荐输出字段：

- gene
- score_name
- r
- pvalue
- padj
- direction
- quadrant

## 6. 常见 debug

### 样本名不匹配

可以用：

```r
setdiff(target_df$sample, colnames(feature_mat))
setdiff(colnames(feature_mat), target_df$sample)
```

### p 值显著但 r 不高

说明存在统计显著性，但相关强度可能一般，不能只看 p 值。

## 7. 下一步

- 加入 r / p 阈值开关
- 加入四象限输出
- 批量导出 CSV
- 批量绘制散点图