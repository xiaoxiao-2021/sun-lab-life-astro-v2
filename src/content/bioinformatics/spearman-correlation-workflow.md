---

title: "Spearman correlation workflow"
description: "中文通常称为 Spearman 秩相关分析，是一种用于评估两个变量之间单调相关关系的非参数统计方法"
date: 2026-07-16
updated: 2026-07-16
category: "Correlation"
tags: ["Spearman", "Correlation", "R", "Bioinformatics"]
status: "Finish"
draft: false
------------

## 一、方法概述

### 1.1 这篇笔记解决什么问题？

这篇笔记用于整理一个可重复使用的 **Spearman correlation workflow**，主要包括：

* 输入数据整理
* 样本名称检查
* 样本顺序匹配
* 单次相关性分析
* 多个 feature 的批量相关性分析
* 多个 target 的批量相关性分析
* 分组内相关性分析
* 多重检验校正
* 结果筛选与导出
* 单个相关关系的散点图绘制
* 常见报错与 Debug

Spearman correlation 并不局限于某一种数据类型。

只要两个变量：

1. 来自同一批样本或观测对象
2. 能够按照样本一一对应
3. 数值可以进行排序

就可以考虑使用 Spearman correlation。

常见应用场景包括：

* 基因表达量 vs 免疫评分
* 差异基因表达量 vs ssGSEA score
* 基因表达量 vs 临床指标
* 基因表达量 vs 免疫细胞浸润评分
* 基因 vs 基因
* 蛋白表达量 vs 表型数据
* 微生物丰度 vs 代谢指标
* 分组内相关性分析
* 多个 feature 与多个 target 的批量相关性分析

因此，这个工作流本质上解决的问题是：

```text
一个或多个 feature
        vs
一个或多个 target variable
```

基因表达量与 ssGSEA score 只是其中一种具体的生物信息学应用。

---

### 1.2 Spearman correlation 的基本原理

Spearman correlation，中文通常称为 **Spearman 秩相关分析**。

它不是直接根据原始数值计算相关性，而是先将原始数值转换为**秩次（rank）**，再判断两个变量的秩次是否具有一致的变化趋势。

例如，某个变量的原始数值为：

```text
10, 20, 50, 100
```

转换为秩次后为：

```text
1, 2, 3, 4
```

如果另一个变量的秩次也随着它持续升高，二者就可能具有较强的正相关关系。

因此，Spearman correlation 更关注：

> 两个变量的相对顺序是否一致。

而不是：

> 两个变量的原始数值具体相差多少。

Spearman 相关系数通常记为 `ρ`，在实际结果表中也经常简写为 `r`。

其取值范围为：

```text
-1 ≤ r ≤ 1
```

相关系数的基本含义如下：

* `r` 接近 `1`：较强的正向单调关系
* `r` 接近 `-1`：较强的负向单调关系
* `r` 接近 `0`：没有明显的单调关系，也就是两个变量暂时“不太熟”

正相关表示：

> 一个变量增大时，另一个变量总体上也倾向于增大。

负相关表示：

> 一个变量增大时，另一个变量总体上倾向于减小。

Spearman correlation 关注的是两个变量是否具有一致的上升或下降趋势，并不要求这种关系必须表现为严格的直线。

---

### 1.3 Spearman 与 Pearson 的区别

Spearman correlation 和 Pearson correlation 都可以用于相关性分析，但二者评估的关系不同。

| 方法       | 主要评估内容 | 计算依据 | 对异常值的敏感程度 |
| -------- | ------ | ---- | --------- |
| Pearson  | 线性关系   | 原始数值 | 相对较高      |
| Spearman | 单调关系   | 数值秩次 | 相对较低      |

#### Pearson correlation

Pearson correlation 主要用于评估两个变量之间是否存在线性关系。

通常适用于：

* 两个变量之间的关系近似线性
* 数据中没有明显的极端值
* 数据分布和方差结构相对稳定
* 原始数值大小本身具有重要意义

#### Spearman correlation

Spearman correlation 主要评估两个变量之间是否存在单调关系。

由于它根据秩次进行计算，因此具有以下特点：

* 不要求原始数据服从正态分布
* 不要求两个变量之间必须是严格的线性关系
* 对极端值和异常值相对不敏感
* 可以用于连续变量
* 可以用于有序变量

基因表达数据经常具有以下特点：

* 数据分布偏斜
* 不同基因的表达范围差异较大
* 可能存在极端高表达或极端低表达样本
* 基因表达量与目标变量之间不一定满足严格线性关系

因此，对于基因表达量相关性分析，**Spearman correlation 通常是一个相对稳健的选择**。

但这并不意味着所有基因表达数据都必须使用 Spearman correlation。

如果表达量经过合理标准化和转换，两个变量之间近似线性，并且不存在明显异常值，也可以考虑使用 Pearson correlation。

具体选择哪一种相关性分析方法，需要结合：

```text
数据分布
+
散点图形态
+
异常值情况
+
研究关注的是线性关系还是单调关系
```

---

## 二、数据准备

### 2.1 工作流的输入形式

为了固定输入格式、减少代码分支并简化后续 Debug，本工作流统一使用两个 `data.frame` 作为输入：

```text
feature_df
target_df
```

其中：

```text
feature_df：保存待分析的 feature 数值

target_df：保存 sample、target 和 group 等信息
```

需要注意的是：

> 整个工作流输入的是两个数据框，但 `cor.test()` 真正计算的是从两个数据框中提取出的两个数值向量。

整个过程可以理解为：

```text
两个 data.frame
       ↓
根据 feature、target、sample 和 group 筛选数据
       ↓
提取两个一一对应的 numeric vector
       ↓
运行 Spearman correlation
```

---

### 2.2 feature_df 的数据结构

`feature_df` 用于保存待进行相关性分析的特征数据。

推荐格式为：

* 第一列为 feature 名称
* 其余列为 sample
* 每一行表示一个 feature
* 每个单元格表示该 feature 在对应 sample 中的数值

例如：

| feature | S01 | S02 | S03 | S04 | S05 | S06 |
| ------- | --: | --: | --: | --: | --: | --: |
| Gene_A  | 4.2 | 5.1 | 5.8 | 3.1 | 3.5 | 4.0 |
| Gene_B  | 2.3 | 2.8 | 3.2 | 5.1 | 4.8 | 4.3 |
| Gene_C  | 7.4 | 6.9 | 8.0 | 6.1 | 7.2 | 7.8 |

其中，feature 可以是：

* gene
* protein
* metabolite
* microbial species
* clinical measurement
* questionnaire item
* 其他可以进行排序的变量

在基因表达分析中，`feature_df` 可以保存：

```text
log2 normalized expression
VST-transformed expression
rlog-transformed expression
经过合理转换的 TPM 或 CPM
```

一般不建议直接使用未经标准化和转换的 raw counts 进行相关性分析。

---

### 2.3 target_df 的数据结构

`target_df` 用于保存样本名称、目标变量和分组信息。

推荐格式为：

* 每一行表示一个 sample
* `sample` 列保存样本名称
* 数值列保存待分析的 target
* `group` 列保存样本分组

例如：

| sample | Score_A | Score_B | group |
| ------ | ------: | ------: | ----- |
| S01    |    0.42 |    1.15 | IA    |
| S02    |    0.81 |    0.94 | IA    |
| S03    |    1.10 |    0.75 | IA    |
| S04    |    0.35 |    1.34 | ENEG  |
| S05    |    0.53 |    1.08 | ENEG  |
| S06    |    0.91 |    0.82 | ENEG  |

目标变量可以是：

* ssGSEA score
* immune infiltration score
* pathway score
* clinical index
* laboratory measurement
* phenotype score
* 其他连续变量或有序变量

因此，两个数据框可以概括为：

```text
feature_df
行 = feature
列 = sample
值 = feature value
```

```text
target_df
行 = sample
列 = sample、target、group
```

---

### 2.4 数值矩阵转换为 feature_df

R 中的特征数据也可能保存为数值型矩阵。

例如：

```r
feature_mat
```

显示为：

```text
       S01 S02 S03
Gene_A 4.2 5.1 3.8
Gene_B 2.3 2.8 4.1
```

此时：

```r
colnames(feature_mat)
```

返回：

```text
"S01" "S02" "S03"
```

而：

```r
rownames(feature_mat)
```

返回：

```text
"Gene_A" "Gene_B"
```

也就是说，`Gene_A` 和 `Gene_B` 保存在矩阵的行名中，并不是正式的数据列。

可以使用：

```r
feature_df <- tibble::rownames_to_column(
  as.data.frame(feature_mat),
  var = "feature"
)
```

转换后得到：

| feature | S01 | S02 | S03 |
| ------- | --: | --: | --: |
| Gene_A  | 4.2 | 5.1 | 3.8 |
| Gene_B  | 2.3 | 2.8 | 4.1 |

如果不使用 `tibble`，也可以使用基础 R：

```r
feature_df <- as.data.frame(feature_mat)

feature_df <- data.frame(
  feature = rownames(feature_df),
  feature_df,
  row.names = NULL,
  check.names = FALSE
)
```

数值矩阵本身也可以作为相关性分析的数据来源。

矩阵与数据框的主要区别在于：

* 数据的保存形式不同
* 行、列的提取方式不同
* 数据框可以同时保存不同类型的列

最终进入 `cor.test()` 的仍然是两个数值向量。

为了固定工作流格式，本流程统一要求在正式分析前将特征数据整理为 `feature_df`。

---

### 2.5 target_df 中不同列的作用

`target_df` 中的列可以分为两类。

#### 信息列和筛选列

例如：

```text
sample
group
batch
sex
stage
```

这些列主要用于：

* 样本识别
* 样本匹配
* 样本筛选
* 分组分析
* 结果注释

它们不会自动进入相关性计算。

#### 目标变量列

例如：

```text
Score_A
Score_B
Age
BMI
Clinical_index
```

只有明确指定为目标变量的数值列，才会被提取并参与相关性分析。

例如：

```r
target_col <- "Score_A"
```

程序真正提取的是：

```r
y <- target_df[[target_col]]
```

而不是将整个 `target_df` 传入 `cor.test()`。

因此，`group` 列本身不会影响相关系数的计算。

它只用于决定哪些样本进入某一次分析。

可以概括为：

```text
sample 决定两个数据框如何匹配

group 决定哪些样本参与分析

feature 决定从 feature_df 中提取哪一行

target 决定从 target_df 中提取哪一个数值列
```

---

### 2.6 读取数据

读取 `feature_df`：

```r
feature_df <- read.csv(
  "feature_df.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
```

读取 `target_df`：

```r
target_df <- read.csv(
  "target_df.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
```

这里设置：

```r
check.names = FALSE
```

可以避免 R 自动修改样本名称。

例如，样本名称：

```text
Sample-01
```

在某些情况下可能被自动修改为：

```text
Sample.01
```

这可能导致两个数据框中的样本名称无法匹配。

---

### 2.7 检查输入数据结构

检查两个数据框的基本结构：

```r
str(feature_df)
str(target_df)
```

检查前几行：

```r
head(feature_df)
head(target_df)
```

检查列名：

```r
colnames(feature_df)
colnames(target_df)
```

检查数据维度：

```r
dim(feature_df)
dim(target_df)
```

---

## 三、自动化工作流程

### 3.1 参数设置

首先统一设置工作流参数。

```r
# ============================================================
# 1. 参数设置
# ============================================================

feature_col <- "feature"
sample_col <- "sample"
group_col <- "group"

target_cols <- c(
  "Score_A",
  "Score_B"
)

use_group <- TRUE

selected_groups <- c(
  "IA",
  "ENEG"
)

min_n <- 3

r_cutoff <- 0.5
p_cutoff <- 0.05
padj_cutoff <- 0.05
```

参数说明：

| 参数                | 说明             |
| ----------------- | -------------- |
| `feature_col`     | feature 名称所在列  |
| `sample_col`      | 样本名称所在列        |
| `group_col`       | 分组信息所在列        |
| `target_cols`     | 需要分析的目标变量      |
| `use_group`       | 是否进行分组内相关性分析   |
| `selected_groups` | 需要分析的组         |
| `min_n`           | 允许进行分析的最小有效样本量 |
| `r_cutoff`        | 相关系数筛选阈值       |
| `p_cutoff`        | 原始 p 值筛选阈值     |
| `padj_cutoff`     | 校正后 p 值筛选阈值    |

---

### 3.2 检查必要列是否存在

检查 feature 列：

```r
if (!feature_col %in% colnames(feature_df)) {
  stop("feature_df 中不存在指定的 feature_col。")
}
```

检查 sample 列：

```r
if (!sample_col %in% colnames(target_df)) {
  stop("target_df 中不存在指定的 sample_col。")
}
```

检查 target 列：

```r
missing_target_cols <- setdiff(
  target_cols,
  colnames(target_df)
)

if (length(missing_target_cols) > 0) {
  stop(
    paste0(
      "以下 target 列不存在：",
      paste(missing_target_cols, collapse = ", ")
    )
  )
}
```

如果进行分组分析，还需要检查 group 列：

```r
if (
  use_group &&
  !group_col %in% colnames(target_df)
) {
  stop("target_df 中不存在指定的 group_col。")
}
```

---

### 3.3 检查 feature 和 sample 是否重复

检查 feature 名称：

```r
duplicated_features <- feature_df[[feature_col]][
  duplicated(feature_df[[feature_col]])
]

duplicated_features
```

如果没有重复，结果通常为：

```r
character(0)
```

检查样本名称：

```r
duplicated_samples <- target_df[[sample_col]][
  duplicated(target_df[[sample_col]])
]

duplicated_samples
```

feature 或 sample 出现重复时，需要先明确重复数据的来源，再决定：

* 删除重复记录
* 合并重复记录
* 计算平均值
* 保留其中一条记录

不建议在不了解数据含义的情况下直接去重。

---

### 3.4 获取 feature_df 中的样本列

除了 feature 列以外，其余列都作为样本列：

```r
sample_cols <- setdiff(
  colnames(feature_df),
  feature_col
)
```

查看样本列：

```r
head(sample_cols)
```

检查样本数量：

```r
length(sample_cols)
```

---

### 3.5 检查样本是否匹配

相关性分析最重要的前提之一，是两个变量中的每个数值必须来自同一个样本。

例如：

```text
Gene_A 在 S01 中的表达量
```

必须与：

```text
Score_A 在 S01 中的评分
```

进行配对。

不能仅仅因为两个向量长度相同，就直接进行相关性分析。

检查 `target_df` 中有、但 `feature_df` 中没有的样本：

```r
setdiff(
  target_df[[sample_col]],
  sample_cols
)
```

检查 `feature_df` 中有、但 `target_df` 中没有的样本：

```r
setdiff(
  sample_cols,
  target_df[[sample_col]]
)
```

理想情况下，两段代码都返回：

```r
character(0)
```

---

### 3.6 提取共同样本

```r
common_samples <- intersect(
  sample_cols,
  target_df[[sample_col]]
)
```

检查共同样本数量：

```r
length(common_samples)
```

如果共同样本数量过少，应当停止分析并检查样本名称。

```r
if (length(common_samples) < min_n) {
  stop("两个数据框中的共同样本数量不足。")
}
```

---

### 3.7 统一样本顺序

即使两个数据框包含相同的样本，如果样本顺序不同，也可能得到完全错误的相关性结果。

首先按照共同样本整理 `feature_df`：

```r
feature_df <- feature_df[
  ,
  c(
    feature_col,
    common_samples
  ),
  drop = FALSE
]
```

然后按照相同顺序整理 `target_df`：

```r
target_df <- target_df[
  match(
    common_samples,
    target_df[[sample_col]]
  ),
  ,
  drop = FALSE
]
```

再次确认：

```r
identical(
  colnames(feature_df)[-1],
  target_df[[sample_col]]
)
```

理想情况下返回：

```r
TRUE
```

如果返回 `FALSE`，不能继续进行后续分析。

---

### 3.8 检查分析列是否为数值型

检查 `feature_df` 中的样本列：

```r
feature_numeric_check <- sapply(
  feature_df[common_samples],
  is.numeric
)

feature_numeric_check
```

如果存在非数值型样本列：

```r
names(
  feature_numeric_check[
    !feature_numeric_check
  ]
)
```

检查目标变量：

```r
target_numeric_check <- sapply(
  target_df[target_cols],
  is.numeric
)

target_numeric_check
```

如果存在非数值型 target：

```r
names(
  target_numeric_check[
    !target_numeric_check
  ]
)
```

不建议直接对所有列使用 `as.numeric()`，因为非数值字符可能被转换为 `NA`，从而掩盖原始数据问题。

应当先检查是否存在：

```text
-
unknown
not detected
空字符串
单位字符
```

再决定如何处理。

---

### 3.9 单次 Spearman correlation

假设需要分析：

```text
Gene_A vs Score_A
```

设置 feature 和 target：

```r
feature_name <- "Gene_A"
target_col <- "Score_A"
```

从 `feature_df` 中提取 feature vector：

```r
x <- as.numeric(
  unlist(
    feature_df[
      feature_df[[feature_col]] == feature_name,
      common_samples,
      drop = FALSE
    ],
    use.names = FALSE
  )
)
```

从 `target_df` 中提取 target vector：

```r
y <- target_df[[target_col]]
```

检查两个向量长度：

```r
length(x)
length(y)
```

检查缺失值：

```r
sum(is.na(x))
sum(is.na(y))
```

只保留两个变量都有数值的样本：

```r
complete_idx <- complete.cases(
  x,
  y
)
```

进行 Spearman correlation：

```r
test_result <- cor.test(
  x = x[complete_idx],
  y = y[complete_idx],
  method = "spearman",
  exact = FALSE
)
```

提取相关系数：

```r
r <- unname(
  test_result$estimate
)
```

提取 p 值：

```r
pvalue <- test_result$p.value
```

记录实际参与分析的样本数量：

```r
n <- sum(complete_idx)
```

整理结果：

```r
single_result <- data.frame(
  feature = feature_name,
  target = target_col,
  group = "All",
  r = r,
  pvalue = pvalue,
  n = n
)
```

---

### 3.10 为什么设置 exact = FALSE？

当数据中存在相同数值时，会产生相同秩次，也就是 tied ranks。

在这种情况下，R 可能提示：

```text
Cannot compute exact p-value with ties
```

因此，通常可以设置：

```r
exact = FALSE
```

使用近似方法计算 p 值。

```r
cor.test(
  x,
  y,
  method = "spearman",
  exact = FALSE
)
```

---

### 3.11 自动化批量分析函数

下面将样本匹配、整体分析、分组分析、多 target 分析和结果整理封装为一个函数。

```r
run_spearman_workflow <- function(
    feature_df,
    target_df,
    feature_col = "feature",
    sample_col = "sample",
    target_cols,
    group_col = NULL,
    selected_groups = NULL,
    min_n = 3
) {

  # ----------------------------------------------------------
  # 1. 检查输入对象
  # ----------------------------------------------------------

  if (!is.data.frame(feature_df)) {
    stop("feature_df 必须是 data.frame。")
  }

  if (!is.data.frame(target_df)) {
    stop("target_df 必须是 data.frame。")
  }

  if (!feature_col %in% colnames(feature_df)) {
    stop("feature_df 中不存在指定的 feature_col。")
  }

  if (!sample_col %in% colnames(target_df)) {
    stop("target_df 中不存在指定的 sample_col。")
  }

  missing_target_cols <- setdiff(
    target_cols,
    colnames(target_df)
  )

  if (length(missing_target_cols) > 0) {
    stop(
      paste0(
        "以下 target 列不存在：",
        paste(missing_target_cols, collapse = ", ")
      )
    )
  }

  if (
    !is.null(group_col) &&
    !group_col %in% colnames(target_df)
  ) {
    stop("target_df 中不存在指定的 group_col。")
  }

  # ----------------------------------------------------------
  # 2. 检查重复名称
  # ----------------------------------------------------------

  if (anyDuplicated(feature_df[[feature_col]]) > 0) {
    stop("feature_df 中存在重复的 feature 名称。")
  }

  if (anyDuplicated(target_df[[sample_col]]) > 0) {
    stop("target_df 中存在重复的 sample 名称。")
  }

  # ----------------------------------------------------------
  # 3. 提取共同样本
  # ----------------------------------------------------------

  sample_cols <- setdiff(
    colnames(feature_df),
    feature_col
  )

  common_samples <- intersect(
    sample_cols,
    target_df[[sample_col]]
  )

  if (length(common_samples) < min_n) {
    stop("两个数据框中的共同样本数量不足。")
  }

  # ----------------------------------------------------------
  # 4. 统一样本顺序
  # ----------------------------------------------------------

  feature_df <- feature_df[
    ,
    c(
      feature_col,
      common_samples
    ),
    drop = FALSE
  ]

  target_df <- target_df[
    match(
      common_samples,
      target_df[[sample_col]]
    ),
    ,
    drop = FALSE
  ]

  if (
    !identical(
      colnames(feature_df)[-1],
      target_df[[sample_col]]
    )
  ) {
    stop("样本顺序匹配失败。")
  }

  # ----------------------------------------------------------
  # 5. 检查数据类型
  # ----------------------------------------------------------

  feature_numeric_check <- sapply(
    feature_df[common_samples],
    is.numeric
  )

  if (!all(feature_numeric_check)) {
    stop("feature_df 中存在非数值型样本列。")
  }

  target_numeric_check <- sapply(
    target_df[target_cols],
    is.numeric
  )

  if (!all(target_numeric_check)) {
    stop("target_df 中存在非数值型 target 列。")
  }

  # ----------------------------------------------------------
  # 6. 确定分析组别
  # ----------------------------------------------------------

  if (is.null(group_col)) {

    group_list <- "All"

  } else {

    if (is.null(selected_groups)) {

      group_list <- unique(
        target_df[[group_col]][
          !is.na(target_df[[group_col]])
        ]
      )

    } else {

      group_list <- selected_groups

    }
  }

  # ----------------------------------------------------------
  # 7. 批量计算
  # ----------------------------------------------------------

  result_list <- list()
  result_index <- 1

  for (current_group in group_list) {

    if (is.null(group_col)) {

      group_target_df <- target_df
      group_samples <- target_df[[sample_col]]
      group_label <- "All"

    } else {

      group_keep <- (
        !is.na(target_df[[group_col]]) &
        target_df[[group_col]] == current_group
      )

      group_target_df <- target_df[
        group_keep,
        ,
        drop = FALSE
      ]

      group_samples <- group_target_df[[sample_col]]
      group_label <- as.character(current_group)
    }

    for (target_col in target_cols) {

      y <- group_target_df[[target_col]]

      for (i in seq_len(nrow(feature_df))) {

        feature_name <- feature_df[[feature_col]][i]

        x <- as.numeric(
          unlist(
            feature_df[
              i,
              group_samples,
              drop = FALSE
            ],
            use.names = FALSE
          )
        )

        complete_idx <- complete.cases(
          x,
          y
        )

        n_complete <- sum(complete_idx)

        x_complete <- x[complete_idx]
        y_complete <- y[complete_idx]

        if (
          n_complete < min_n ||
          length(unique(x_complete)) < 2 ||
          length(unique(y_complete)) < 2
        ) {

          current_result <- data.frame(
            feature = feature_name,
            target = target_col,
            group = group_label,
            r = NA_real_,
            pvalue = NA_real_,
            n = n_complete
          )

        } else {

          test_result <- cor.test(
            x = x_complete,
            y = y_complete,
            method = "spearman",
            exact = FALSE
          )

          current_result <- data.frame(
            feature = feature_name,
            target = target_col,
            group = group_label,
            r = unname(
              test_result$estimate
            ),
            pvalue = test_result$p.value,
            n = n_complete
          )
        }

        result_list[[result_index]] <- current_result
        result_index <- result_index + 1
      }
    }
  }

  cor_result <- do.call(
    rbind,
    result_list
  )

  # ----------------------------------------------------------
  # 8. 多重检验校正
  # ----------------------------------------------------------

  correction_group <- interaction(
    cor_result$target,
    cor_result$group,
    drop = TRUE
  )

  cor_result$padj <- ave(
    cor_result$pvalue,
    correction_group,
    FUN = function(x) {
      p.adjust(
        x,
        method = "BH"
      )
    }
  )

  # ----------------------------------------------------------
  # 9. 添加方向
  # ----------------------------------------------------------

  cor_result$direction <- ifelse(
    is.na(cor_result$r),
    NA_character_,
    ifelse(
      cor_result$r > 0,
      "Positive",
      ifelse(
        cor_result$r < 0,
        "Negative",
        "No direction"
      )
    )
  )

  # ----------------------------------------------------------
  # 10. 添加相关强度
  # ----------------------------------------------------------

  cor_result$strength <- cut(
    abs(cor_result$r),
    breaks = c(
      -Inf,
      0.2,
      0.4,
      0.6,
      0.8,
      Inf
    ),
    labels = c(
      "Very weak",
      "Weak",
      "Moderate",
      "Strong",
      "Very strong"
    ),
    right = FALSE
  )

  return(cor_result)
}
```

---

### 3.12 整体相关性分析

如果不进行分组分析：

```r
cor_result_all <- run_spearman_workflow(
  feature_df = feature_df,
  target_df = target_df,
  feature_col = feature_col,
  sample_col = sample_col,
  target_cols = target_cols,
  group_col = NULL,
  selected_groups = NULL,
  min_n = min_n
)
```

此时结果中的：

```text
group = "All"
```

表示使用全部共同样本进行分析。

---

### 3.13 分组内相关性分析

如果需要分析 IA 和 ENEG 组内部的相关性：

```r
cor_result_group <- run_spearman_workflow(
  feature_df = feature_df,
  target_df = target_df,
  feature_col = feature_col,
  sample_col = sample_col,
  target_cols = target_cols,
  group_col = group_col,
  selected_groups = selected_groups,
  min_n = min_n
)
```

需要注意的是：

> `group` 不会直接进入 `cor.test()`。

程序的实际逻辑是：

```text
根据 group 筛选样本
        ↓
提取该组内的 feature vector
        ↓
提取该组内的 target vector
        ↓
计算该组内部的 Spearman correlation
```

例如 IA 组内的分析，真正进入 `cor.test()` 的只有：

```text
x = IA 组样本中的 Gene_A 表达量

y = IA 组样本中的 Score_A
```

`group` 只负责筛选样本，不会作为第三个变量参与相关性计算。

---

### 3.14 使用分析开关

也可以通过 `use_group` 控制是否进行分组分析：

```r
if (use_group) {

  cor_result <- run_spearman_workflow(
    feature_df = feature_df,
    target_df = target_df,
    feature_col = feature_col,
    sample_col = sample_col,
    target_cols = target_cols,
    group_col = group_col,
    selected_groups = selected_groups,
    min_n = min_n
  )

} else {

  cor_result <- run_spearman_workflow(
    feature_df = feature_df,
    target_df = target_df,
    feature_col = feature_col,
    sample_col = sample_col,
    target_cols = target_cols,
    group_col = NULL,
    selected_groups = NULL,
    min_n = min_n
  )
}
```

这样，正式运行时只需要修改：

```r
use_group <- TRUE
```

或者：

```r
use_group <- FALSE
```

---

## 四、结果整理与输出

### 4.1 输出结果结构

自动化函数输出长格式结果表。

| feature | target  | group |     r | pvalue |  n |  padj | direction | strength  |
| ------- | ------- | ----- | ----: | -----: | -: | ----: | --------- | --------- |
| Gene_A  | Score_A | IA    |  0.72 |  0.012 | 18 | 0.084 | Positive  | Strong    |
| Gene_B  | Score_A | IA    | -0.61 |  0.031 | 18 | 0.120 | Negative  | Strong    |
| Gene_A  | Score_A | ENEG  |  0.15 |  0.630 | 16 | 0.890 | Positive  | Very weak |

每一行表示：

```text
一个 feature
×
一个 target
×
一个 group
```

对应的一次 Spearman correlation。

如果不进行分组分析，则：

```text
group = "All"
```

---

### 4.2 多重检验校正

当同时对大量 feature 进行分析时，会产生大量 p 值。

如果仅使用：

```text
pvalue < 0.05
```

随着检验数量增加，出现假阳性结果的概率也会增加。

本工作流使用 Benjamini-Hochberg 方法计算：

```text
padj
```

并且在每个：

```text
target × group
```

内部独立进行多重检验校正。

也就是说：

```text
Score_A × IA
Score_A × ENEG
Score_B × IA
Score_B × ENEG
```

分别进行 p 值校正。

---

### 4.3 结果筛选

使用原始 p 值筛选：

```r
filtered_result_p <- subset(
  cor_result,
  !is.na(r) &
    abs(r) >= r_cutoff &
    pvalue < p_cutoff
)
```

使用校正后的 p 值筛选：

```r
filtered_result_padj <- subset(
  cor_result,
  !is.na(r) &
    abs(r) >= r_cutoff &
    padj < padj_cutoff
)
```

例如：

```r
r_cutoff <- 0.5
p_cutoff <- 0.05
padj_cutoff <- 0.05
```

分别对应：

```text
|r| ≥ 0.5
pvalue < 0.05
padj < 0.05
```

不能只根据 p 值判断结果是否重要。

结果解释通常需要同时考虑：

```text
相关方向
+
相关强度
+
统计显著性
+
有效样本量
+
散点图分布
```

---

### 4.4 导出结果

创建结果文件夹：

```r
dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)
```

导出完整结果：

```r
write.csv(
  cor_result,
  file = "results/spearman_all_results.csv",
  row.names = FALSE
)
```

导出使用原始 p 值筛选的结果：

```r
write.csv(
  filtered_result_p,
  file = "results/spearman_filtered_by_pvalue.csv",
  row.names = FALSE
)
```

导出使用校正后 p 值筛选的结果：

```r
write.csv(
  filtered_result_padj,
  file = "results/spearman_filtered_by_padj.csv",
  row.names = FALSE
)
```

---

### 4.5 按 group 分别导出结果

如果进行了分组分析，可以将不同组分别导出。

```r
split_result <- split(
  cor_result,
  cor_result$group
)
```

批量导出：

```r
invisible(
  lapply(
    names(split_result),
    function(group_name) {

      write.csv(
        split_result[[group_name]],
        file = file.path(
          "results",
          paste0(
            "spearman_",
            group_name,
            ".csv"
          )
        ),
        row.names = FALSE
      )
    }
  )
)
```

---

### 4.6 单个相关关系的散点图

相关性结果不能只看 `r` 和 p 值。

还需要通过散点图判断：

* 是否存在异常值
* 是否由少数样本驱动
* 是否具有总体单调趋势
* 是否存在明显的分组结构

例如绘制：

```text
Gene_A vs Score_A
```

首先准备绘图数据：

```r
library(ggplot2)

feature_name <- "Gene_A"
target_col <- "Score_A"

common_samples <- intersect(
  setdiff(
    colnames(feature_df),
    feature_col
  ),
  target_df[[sample_col]]
)

plot_target_df <- target_df[
  match(
    common_samples,
    target_df[[sample_col]]
  ),
  ,
  drop = FALSE
]

x <- as.numeric(
  unlist(
    feature_df[
      feature_df[[feature_col]] == feature_name,
      common_samples,
      drop = FALSE
    ],
    use.names = FALSE
  )
)

plot_df <- data.frame(
  sample = common_samples,
  feature_value = x,
  target_value = plot_target_df[[target_col]],
  group = plot_target_df[[group_col]]
)
```

绘制散点图：

```r
ggplot(
  plot_df,
  aes(
    x = feature_value,
    y = target_value
  )
) +
  geom_point(
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    title = paste0(
      feature_name,
      " vs ",
      target_col
    ),
    x = feature_name,
    y = target_col
  ) +
  theme_classic()
```

需要注意：

> Spearman correlation 衡量的是单调关系，而 `geom_smooth(method = "lm")` 绘制的是线性拟合线。

因此，这条拟合线只是用于辅助观察总体趋势，并不是 Spearman correlation 的计算过程。

如果不希望显示线性拟合线，可以只绘制散点：

```r
ggplot(
  plot_df,
  aes(
    x = feature_value,
    y = target_value
  )
) +
  geom_point(
    size = 2
  ) +
  labs(
    x = feature_name,
    y = target_col
  ) +
  theme_classic()
```

---

### 4.7 分组散点图

如果需要查看不同 group 的数据分布，可以按照 group 着色：

```r
ggplot(
  plot_df,
  aes(
    x = feature_value,
    y = target_value,
    color = group
  )
) +
  geom_point(
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  theme_classic()
```

也可以按照 group 分面：

```r
ggplot(
  plot_df,
  aes(
    x = feature_value,
    y = target_value
  )
) +
  geom_point(
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  facet_wrap(
    ~ group
  ) +
  theme_classic()
```

分面图可以更直观地观察不同组内部的相关趋势。

---

## 五、结果解释

### 5.1 相关系数的方向

当：

```text
r > 0
```

表示正相关。

当：

```text
r < 0
```

表示负相关。

当：

```text
r ≈ 0
```

表示没有明显的单调关系。

---

### 5.2 相关强度

可以使用 `abs(r)` 描述相关强度。

| `abs(r)`  | 描述          |
| --------- | ----------- |
| 0.00–0.19 | Very weak   |
| 0.20–0.39 | Weak        |
| 0.40–0.59 | Moderate    |
| 0.60–0.79 | Strong      |
| 0.80–1.00 | Very strong |

需要注意的是，这种分级没有完全统一的标准。

不同研究领域、样本量和研究目的可能使用不同阈值，因此不应机械地使用固定分级。

---

### 5.3 p 值显著但 r 不高

这种情况表示统计检验发现相关关系可能不同于零，但相关强度较弱。

这种现象在样本量较大时更容易出现。

因此不能只看 p 值，还需要关注：

* `r` 的大小
* 实际样本量
* 散点图分布
* 是否存在异常值
* 结果是否具有实际或生物学意义

---

### 5.4 r 较高但 p 值不显著

可能原因包括：

* 样本量较小
* 分组后有效样本数量不足
* 个别样本对结果影响较大
* 数据中存在大量相同秩次
* 相关关系不够稳定

此时需要结合散点图和实际样本量进一步判断。

---

### 5.5 不同 group 的相关系数不同

例如：

```text
IA：r = 0.70

ENEG：r = 0.20
```

这只能说明两个组中观察到的相关系数大小不同。

不能仅凭：

```text
0.70 > 0.20
```

就认为两个组的相关性存在统计学差异。

如果研究目的是正式比较两个组的相关系数，还需要进一步使用相关系数差异检验，或者建立包含 group 交互项的统计模型。

---

### 5.6 相关性不等于因果关系

Spearman correlation 只能说明两个变量之间存在单调关联。

即使相关系数很高，也不能直接证明：

```text
变量 A 导致变量 B
```

观察到的相关关系可能来源于：

* A 影响 B
* B 影响 A
* 第三个变量同时影响 A 和 B
* 分组结构造成的表观相关
* 批次效应
* 偶然的统计关联

因此，相关性分析通常用于：

* 发现候选关系
* 提出研究假设
* 筛选潜在 feature
* 辅助后续实验设计

不能单独用于证明因果机制。

---

## 六、常见 Debug

### 6.1 样本名称不匹配

检查：

```r
setdiff(
  target_df[[sample_col]],
  sample_cols
)

setdiff(
  sample_cols,
  target_df[[sample_col]]
)
```

可能原因包括：

* 样本名称前后存在空格
* 大小写不一致
* 使用了不同的连接符
* 一个数据框中缺少部分样本
* 样本名称被 R 自动修改

去掉前后空格：

```r
target_df[[sample_col]] <- trimws(
  target_df[[sample_col]]
)

colnames(feature_df) <- trimws(
  colnames(feature_df)
)
```

---

### 6.2 样本顺序不一致

即使样本集合相同，也必须重新排列。

```r
target_df <- target_df[
  match(
    sample_cols,
    target_df[[sample_col]]
  ),
  ,
  drop = FALSE
]
```

再次检查：

```r
identical(
  sample_cols,
  target_df[[sample_col]]
)
```

应返回：

```r
TRUE
```

---

### 6.3 feature 名称保存在行名中

检查：

```r
rownames(feature_df)
```

如果 feature 名称保存在行名中，需要转换为普通列：

```r
feature_df <- tibble::rownames_to_column(
  feature_df,
  var = feature_col
)
```

否则：

```r
feature_df[[feature_col]]
```

可能返回：

```r
NULL
```

---

### 6.4 数据不是数值型

检查 feature 样本列：

```r
sapply(
  feature_df[sample_cols],
  is.numeric
)
```

检查 target：

```r
sapply(
  target_df[target_cols],
  is.numeric
)
```

如果某一列不是数值型，应先检查其中是否存在字符、单位或异常标记，而不是直接强制转换。

---

### 6.5 存在缺失值

检查：

```r
sum(is.na(x))
sum(is.na(y))
```

只保留完整观测值：

```r
complete_idx <- complete.cases(
  x,
  y
)
```

相关性计算使用：

```r
x[complete_idx]
y[complete_idx]
```

结果中的 `n` 应记录实际参与分析的样本数，而不是原始总样本数。

---

### 6.6 feature 在所有样本中完全相同

例如：

```text
5, 5, 5, 5, 5
```

该变量没有变化，无法计算相关性。

检查：

```r
length(
  unique(x)
)
```

如果结果小于 `2`，说明该变量没有变异。

---

### 6.7 target 在所有样本中完全相同

例如：

```text
1, 1, 1, 1, 1
```

同样无法计算相关性。

检查：

```r
length(
  unique(y)
)
```

---

### 6.8 分组后样本量过少

分组分析会减少实际参与计算的样本数。

例如，整体数据有 50 个样本，但 IA 组只有 6 个样本。

即使得到较高的 `r`，结果也可能不稳定。

因此，分组分析结果中必须保留：

```text
n
```

并结合实际样本量进行解释。

---

### 6.9 出现 tied ranks 提示

提示：

```text
Cannot compute exact p-value with ties
```

通常是因为数据中存在重复数值。

设置：

```r
exact = FALSE
```

即可使用近似方法计算 p 值。

---

## 七、工作流总结

整个 Spearman correlation workflow 可以概括为：

```text
准备 feature_df 和 target_df
            ↓
检查两个数据框的结构
            ↓
识别 feature、sample、target 和 group
            ↓
检查 feature 和 sample 是否重复
            ↓
检查样本名称是否一致
            ↓
提取共同样本
            ↓
统一样本顺序
            ↓
检查数据类型和缺失值
            ↓
根据 group 筛选样本
            ↓
提取 feature numeric vector
            ↓
提取 target numeric vector
            ↓
计算 Spearman correlation
            ↓
整理 r、pvalue 和 n
            ↓
进行多重检验校正
            ↓
添加 direction 和 strength
            ↓
根据 r、pvalue 或 padj 筛选
            ↓
导出结果
            ↓
通过散点图检查数据分布
```

核心逻辑为：

```text
工作流输入：
两个 data.frame

实际统计计算：
两个一一匹配的 numeric vector

sample 的作用：
完成两个数据框之间的样本匹配

group 的作用：
筛选进入某一次分析的样本

feature 的作用：
指定从 feature_df 中提取哪一行

target 的作用：
指定从 target_df 中提取哪一个数值列
```

矩阵和数据框不会改变 Spearman correlation 的统计原理。

它们的差异主要在于数据保存形式和行、列提取方式不同。

最终真正参与分析的始终是两个按照样本一一对应的数值向量。
