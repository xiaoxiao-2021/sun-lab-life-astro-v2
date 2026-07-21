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

这篇笔记用于整理一个可重复使用的 **Spearman correlation workflow**，重点说明：

* Spearman correlation 适合解决什么问题
* 输入数据需要整理成什么格式
* 两个数据框之间如何匹配样本
* `group` 列在分组分析中发挥什么作用
* 工作流内部如何完成相关性计算
* 最终会输出哪些结果
* 如何对结果进行基本解释和 Debug

这篇文章主要用于理解和使用工作流。

完整的自动化 R 代码、示例数据和可复用模板将单独整理在 GitHub repository 中。

---

### 1.2 什么是 Spearman correlation？

Spearman correlation，中文通常称为 **Spearman 秩相关分析**，是一种用于评估两个变量之间单调相关关系的非参数统计方法。

它并不直接根据原始数值之间的距离计算相关性，而是先将数据转换为秩次，再判断两个变量的秩次是否具有一致的变化趋势。

例如，一个变量的原始数值为：

```text
10, 20, 50, 100
```

对应秩次为：

```text
1, 2, 3, 4
```

如果另一个变量也随着第一个变量持续升高，那么二者可能具有正向单调关系。

因此，Spearman correlation 更关注：

> 两个变量的相对顺序是否一致。

而不是：

> 两个变量的原始数值具体相差多少。

---

### 1.3 相关系数如何理解？

Spearman 相关系数通常记为 `ρ`，在结果表中也经常写作 `rho` 或 `r`。

其取值范围为：

```text
-1 ≤ rho ≤ 1
```

基本含义如下：

* `rho` 接近 `1`：较强的正向单调关系
* `rho` 接近 `-1`：较强的负向单调关系
* `rho` 接近 `0`：没有明显的单调关系，它俩 “不熟”

正相关表示：

> 一个变量升高时，另一个变量总体上也倾向于升高。

负相关表示：

> 一个变量升高时，另一个变量总体上倾向于降低。

需要注意的是，Spearman correlation 不要求两个变量之间必须呈现严格的直线关系。

只要二者整体上保持一致的上升或下降趋势，就可能得到较高的 Spearman 相关系数。

---

### 1.4 Spearman 与 Pearson 的区别

Spearman correlation 和 Pearson correlation 都可以用于相关性分析，但二者关注的关系不同。

| 方法       | 主要评估内容 | 计算依据 | 对异常值的敏感程度 |
| -------- | ------ | ---- | --------- |
| Pearson  | 线性关系   | 原始数值 | 相对较高      |
| Spearman | 单调关系   | 数值秩次 | 相对较低      |

Pearson correlation 更适合分析近似线性的关系。

Spearman correlation 则具有以下特点：

* 不要求原始数据服从正态分布
* 不要求两个变量之间严格线性
* 可以分析单调但非线性的关系
* 对极端值和异常值相对稳健
* 可以用于连续变量和有序变量

基因表达数据经常具有分布偏斜、表达范围差异较大以及极端表达值等特点，因此 Spearman correlation 经常用于基因表达量与评分、表型或临床指标之间的相关性分析。

但这并不意味着所有基因表达数据都必须使用 Spearman。

实际选择方法时，还需要结合：

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

### 1.5 适用场景

Spearman correlation 并不局限于生物信息学数据。

只要两个变量来自同一批样本或观测对象，并且能够一一对应，就可以考虑进行相关性分析。

常见场景包括：

* 基因表达量 vs 免疫评分
* 基因表达量 vs ssGSEA score
* 基因表达量 vs 临床指标
* 基因表达量 vs 免疫细胞浸润评分
* 基因 vs 基因
* 蛋白表达量 vs 表型数据
* 微生物丰度 vs 代谢指标
* 分组内部的相关性分析
* 多个 feature 与多个 target 的批量相关性分析

因此，这个工作流本质上处理的是：

```text
一个或多个 feature
        vs
一个或多个 target
```

基因表达量与 ssGSEA score 只是其中一种具体应用。

---

## 二、输入数据准备

### 2.1 工作流需要哪些输入？

为了固定数据格式并减少后续 Debug，本工作流统一使用两个 `data.frame` 作为输入：

```text
feature_df
target_df
```

其中：

```text
feature_df
用于保存待分析的 feature 数值
```

```text
target_df
用于保存 sample、target 和 group 等信息
```

需要注意的是：

> 工作流接收的是两个数据框，但 `cor.test()` 真正计算的是从两个数据框中提取出的两个数值向量。

整体逻辑为：

```text
feature_df + target_df
          ↓
按照 sample 匹配数据
          ↓
根据需要使用 group 筛选样本
          ↓
提取 feature vector 和 target vector
          ↓
计算 Spearman correlation
```

---

### 2.2 feature_df 的数据格式

`feature_df` 用于保存待分析的特征数据。

推荐结构为：

* 第一列为 feature 名称
* 其余列为样本
* 每一行表示一个 feature
* 每个单元格表示该 feature 在对应样本中的数值

例如：

| feature | S01 | S02 | S03 | S04 |
| ------- | --: | --: | --: | --: |
| Gene_A  | 4.2 | 5.1 | 5.8 | 3.1 |
| Gene_B  | 2.3 | 2.8 | 3.2 | 5.1 |
| Gene_C  | 7.4 | 6.9 | 8.0 | 6.1 |

这里的 feature 可以是：

* gene
* protein
* metabolite
* microbial species
* clinical measurement
* 其他可以进行排序的变量

在基因表达分析中，通常可以使用：

```text
log2 normalized expression
VST-transformed expression
rlog-transformed expression
经过合理转换的 TPM 或 CPM
```

一般不建议直接使用未经标准化和转换的 raw counts 进行相关性分析。

---

### 2.3 数值矩阵如何整理？

原始表达数据可能保存为数值矩阵：

```text
       S01 S02 S03
Gene_A 4.2 5.1 3.8
Gene_B 2.3 2.8 4.1
```

此时 `Gene_A` 和 `Gene_B` 可能保存在矩阵的行名中，而不是正式的数据列。

可以在进入工作流前转换为数据框：

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

数值矩阵本身也可以提取行或列进行相关性分析。

矩阵和数据框的区别主要在于数据保存形式和提取方式不同，最终进入 `cor.test()` 的仍然是两个数值向量。

为了简化工作流，本流程统一规定输入对象为数据框。

---

### 2.4 target_df 的数据格式

`target_df` 用于保存样本名称、目标变量和分组信息。

推荐结构为：

* 每一行表示一个样本
* `sample` 列保存样本名称
* 数值列保存需要分析的 target
* `group` 列保存样本分组

例如：

| sample | Score_A | Score_B | group |
| ------ | ------: | ------: | ----- |
| S01    |    0.42 |    1.15 | IA    |
| S02    |    0.81 |    0.94 | IA    |
| S03    |    1.10 |    0.75 | ENEG  |
| S04    |    0.35 |    1.34 | ENEG  |

target 可以是：

* ssGSEA score
* immune infiltration score
* pathway score
* clinical index
* laboratory measurement
* phenotype score
* 其他连续变量或有序变量

---

### 2.5 不同列分别有什么作用？

两个数据框中不同列的作用可以概括为：

| 内容        | 作用                         |
| --------- | -------------------------- |
| `feature` | 指定需要从 `feature_df` 中提取哪个特征 |
| `sample`  | 完成两个数据框之间的样本匹配             |
| `target`  | 指定需要与 feature 计算相关性的变量     |
| `group`   | 筛选参与某一次组内分析的样本             |

其中，`group` 不会直接进入 Spearman correlation 的计算。

例如进行 IA 组内分析时，程序首先筛选：

```text
group == "IA"
```

然后仅提取 IA 组样本对应的：

```text
Gene_A expression vector
Score_A vector
```

最后计算：

```r
cor.test(
  x,
  y,
  method = "spearman",
  exact = FALSE
)
```

因此：

> `group` 决定分析哪些样本，但不作为第三个变量参与相关系数计算。

---

### 2.6 数据分析前的基本要求

在运行相关性分析前，需要确保：

1. `feature_df` 中的样本名称与 `target_df$sample` 可以匹配
2. 两个数据框中的样本顺序一致
3. feature 数值和 target 数值为数值型
4. 每个样本只能对应一条 target 记录
5. feature 名称和 sample 名称不存在未处理的重复
6. 参与分析的变量不是完全相同的常数
7. 缺失值已经被识别和处理

其中最重要的是：

> 两个向量中的相同位置必须代表同一个样本。

不能仅仅因为两个向量长度相同，就直接进行相关性分析。

---

## 三、核心分析逻辑

### 3.1 样本匹配

首先需要检查两个数据框中的样本是否一致。

```r
setdiff(
  target_df$sample,
  colnames(feature_df)[-1]
)

setdiff(
  colnames(feature_df)[-1],
  target_df$sample
)
```

第一段代码用于寻找：

```text
target_df 中存在，但 feature_df 中不存在的样本
```

第二段代码用于寻找：

```text
feature_df 中存在，但 target_df 中不存在的样本
```

理想情况下，两段代码均返回：

```r
character(0)
```

---

### 3.2 统一样本顺序

即使两个数据框中包含相同样本，如果样本顺序不同，也可能得到错误的相关系数。

因此，需要按照样本名称重新排列：

```r
common_samples <- intersect(
  colnames(feature_df)[-1],
  target_df$sample
)

target_df <- target_df[
  match(common_samples, target_df$sample),
  ,
  drop = FALSE
]
```

随后确认：

```r
identical(
  common_samples,
  target_df$sample
)
```

应返回：

```r
TRUE
```

---

### 3.3 提取两个数值向量

假设需要分析：

```text
Gene_A vs Score_A
```

程序需要从 `feature_df` 中提取 `Gene_A` 在所有匹配样本中的数值：

```text
x = Gene_A expression vector
```

再从 `target_df` 中提取：

```text
y = Score_A vector
```

最终，实际参与计算的是：

```text
numeric vector x
numeric vector y
```

而不是两个完整的数据框。

---

### 3.4 单次相关性计算

Spearman correlation 的核心代码为：

```r
cor.test(
  x = x,
  y = y,
  method = "spearman",
  exact = FALSE
)
```

其中：

* `x`：feature vector
* `y`：target vector
* `method = "spearman"`：指定使用 Spearman correlation
* `exact = FALSE`：当数据中存在相同秩次时，使用近似方法计算 p 值

一次计算通常会得到：

```text
rho
pvalue
n
```

其中 `n` 表示实际参与本次计算的完整样本数量。

---

### 3.5 批量分析的基本思路

如果 `feature_df` 中有大量 feature，工作流会依次完成：

```text
提取 Feature_1
→ 与 target 计算 Spearman correlation

提取 Feature_2
→ 与 target 计算 Spearman correlation

提取 Feature_3
→ 与 target 计算 Spearman correlation
```

如果存在多个 target，则进一步循环：

```text
所有 feature vs Score_A
所有 feature vs Score_B
所有 feature vs Score_C
```

最终将每次相关性分析的结果合并为一个长格式结果表。

---

### 3.6 分组分析的基本思路

如果指定了 `group_col`，工作流会先按照 group 筛选样本，再在每个组内部独立计算相关性。

例如：

```text
IA 组：
所有 feature vs Score_A

ENEG 组：
所有 feature vs Score_A
```

每个组的分析均只使用该组内部的样本。

需要注意：

* 分组后样本量会减少
* 小样本条件下相关系数可能不稳定
* 不同组中相关系数大小不同，不等于两组之间存在显著差异

例如：

```text
IA：rho = 0.70
ENEG：rho = 0.20
```

只能说明两个组中观察到的相关系数不同。

不能仅根据 `0.70 > 0.20` 就认为两个组的相关关系存在统计学差异。

---

## 四、参数与输出结果

### 4.1 核心参数

一个通用 Spearman workflow 通常需要设置：

| 参数                | 说明            |
| ----------------- | ------------- |
| `feature_col`     | feature 名称所在列 |
| `sample_col`      | 样本名称所在列       |
| `target_cols`     | 需要分析的目标变量     |
| `group_col`       | 分组信息所在列       |
| `selected_groups` | 需要分析的组        |
| `min_n`           | 最小有效样本数量      |
| `r_cutoff`        | 相关系数筛选阈值      |
| `p_cutoff`        | 原始 p 值筛选阈值    |
| `padj_cutoff`     | 校正后 p 值筛选阈值   |

如果不进行分组分析，可以不指定 `group_col`。

---

### 4.2 推荐输出字段

推荐结果表至少包含：

| 字段          | 说明            |
| ----------- | ------------- |
| `feature`   | 特征名称          |
| `target`    | 目标变量名称        |
| `group`     | 分析使用的组别       |
| `rho`       | Spearman 相关系数 |
| `pvalue`    | 原始 p 值        |
| `padj`      | 多重检验校正后的 p 值  |
| `n`         | 实际参与分析的样本数量   |
| `direction` | 正相关或负相关       |

例如：

| feature | target  | group |   rho | pvalue |  padj |  n | direction |
| ------- | ------- | ----- | ----: | -----: | ----: | -: | --------- |
| Gene_A  | Score_A | IA    |  0.72 |  0.012 | 0.084 | 18 | Positive  |
| Gene_B  | Score_A | IA    | -0.61 |  0.031 | 0.120 | 18 | Negative  |
| Gene_A  | Score_A | ENEG  |  0.15 |  0.630 | 0.890 | 16 | Positive  |

每一行代表：

```text
一个 feature
×
一个 target
×
一个 group
```

对应的一次 Spearman correlation。

如果没有进行分组分析，可以将：

```text
group = "All"
```

---

### 4.3 pvalue 与 padj

当只分析一对变量时，可以直接查看 `pvalue`。

当同时对大量 feature 进行相关性分析时，会产生大量 p 值，也会增加假阳性结果的概率。

因此，批量分析通常需要进行多重检验校正，例如 Benjamini-Hochberg 方法：

```r
p.adjust(
  pvalue,
  method = "BH"
)
```

校正后的结果通常记录为：

```text
padj
```

在多个 target 或多个 group 的分析中，通常可以在每个：

```text
target × group
```

内部进行多重检验校正。

---

### 4.4 结果筛选原则

常见筛选条件包括：

```text
|rho| ≥ 设定阈值
pvalue < 设定阈值
```

或者：

```text
|rho| ≥ 设定阈值
padj < 设定阈值
```

例如：

```text
|rho| ≥ 0.5
pvalue < 0.05
```

但不能只看 p 值。

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

## 五、结果解释

### 5.1 相关方向

当：

```text
rho > 0
```

表示正相关。

当：

```text
rho < 0
```

表示负相关。

当：

```text
rho ≈ 0
```

表示没有明显的单调关系。

---

### 5.2 相关强度

相关强度可以根据 `abs(rho)` 进行描述。

例如：

| `abs(rho)` | 描述          |
| ---------- | ----------- |
| 0.00–0.19  | Very weak   |
| 0.20–0.39  | Weak        |
| 0.40–0.59  | Moderate    |
| 0.60–0.79  | Strong      |
| 0.80–1.00  | Very strong |

这种分级没有完全统一的标准。

不同研究领域、样本量和研究目的可以采用不同阈值，因此不应机械套用。

---

### 5.3 p 值显著但 rho 不高

这表示统计检验发现相关关系可能不同于零，但相关强度较弱。

这种情况在样本量较大时更容易出现。

因此还需要查看：

* `rho` 的大小
* 有效样本量
* 散点图分布
* 是否存在异常值
* 是否具有实际或生物学意义

---

### 5.4 rho 较高但 p 值不显著

可能原因包括：

* 样本量较小
* 分组后有效样本数量不足
* 少数样本对结果影响较大
* 数据中存在大量相同秩次
* 相关关系不够稳定

此时需要结合散点图和实际样本量进一步判断。

---

### 5.5 相关性不等于因果关系

Spearman correlation 只能说明两个变量之间存在单调关联。

即使相关系数很高，也不能直接证明：

```text
变量 A 导致变量 B
```

观察到的相关关系可能来自：

* A 影响 B
* B 影响 A
* 第三个变量同时影响 A 和 B
* 分组结构造成的表观相关
* 批次效应
* 偶然的统计关联

因此，相关性分析更适合用于：

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
  target_df$sample,
  colnames(feature_df)[-1]
)

setdiff(
  colnames(feature_df)[-1],
  target_df$sample
)
```

常见原因包括：

* 样本名称前后存在空格
* 大小写不一致
* 使用了不同连接符
* 一个数据框中缺少部分样本
* 样本名称被 R 自动修改

---

### 6.2 样本顺序不一致

即使样本集合相同，也必须按照样本名称重新排序。

分析前应确认：

```r
identical(
  feature_samples,
  target_df$sample
)
```

结果为：

```r
TRUE
```

后再进行计算。

---

### 6.3 feature 或 target 不是数值型

检查：

```r
str(feature_df)
str(target_df)
```

可能导致列被读取为字符型的内容包括：

```text
-
unknown
not detected
空字符串
单位字符
```

应先检查原始数据，再决定如何转换或处理。

---

### 6.4 存在缺失值

Spearman correlation 应当使用两个变量均有完整数值的样本。

```r
complete_idx <- complete.cases(
  x,
  y
)
```

结果中的 `n` 应记录实际参与分析的样本数，而不是原始总样本数。

---

### 6.5 变量在所有样本中完全相同

例如：

```text
5, 5, 5, 5, 5
```

该变量没有变化，因此无法计算相关性。

可以检查：

```r
length(
  unique(x)
)
```

如果结果小于 `2`，说明该变量没有有效变异。

---

### 6.6 分组后样本量过少

分组分析会减少实际参与分析的样本数量。

即使得到较高的 `rho`，小样本结果也可能不稳定。

因此，结果表中应当保留：

```text
n
```

并结合样本量解释结果。

---

### 6.7 出现 tied ranks 提示

如果出现：

```text
Cannot compute exact p-value with ties
```

通常是因为数据中存在重复数值，从而产生相同秩次。

可以设置：

```r
exact = FALSE
```

使用近似方法计算 p 值。

---

## 七、工作流总结

整个流程可以概括为：

```text
准备 feature_df 和 target_df
            ↓
检查数据格式
            ↓
按照 sample 匹配样本
            ↓
统一样本顺序
            ↓
根据 group 筛选样本
            ↓
提取 feature vector
            ↓
提取 target vector
            ↓
计算 Spearman correlation
            ↓
整理 rho、pvalue、padj 和 n
            ↓
筛选并解释结果
```

核心理解为：

```text
工作流输入：
两个 data.frame

实际统计计算：
两个一一匹配的 numeric vector

sample：
用于匹配样本

group：
用于筛选样本，不参与相关系数计算

feature：
指定从 feature_df 中提取哪一个特征

target：
指定从 target_df 中提取哪一个目标变量
```

---

## 八、完整代码

完整的 R workflow、参数说明与后续示例文件将持续整理在 GitHub repository 中。

- [查看 workflow 项目目录](https://github.com/xiaoxiao-2021/sun-bioinformatics-workflows/tree/main/spearman-correlation)
- [查看核心 R 脚本](https://github.com/xiaoxiao-2021/sun-bioinformatics-workflows/blob/main/spearman-correlation/R/spearman_correlation_workflow.R)
- [R function 运行脚本](https://github.com/xiaoxiao-2021/sun-bioinformatics-workflows/blob/main/spearman-correlation/R/run_spearman_workflow.R)