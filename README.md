# PCA with Partial SVD Implementation

[![R-CMD-check](https://github.com/dereklei12/pca_R/workflows/R-CMD-check/badge.svg)](https://github.com/dereklei12/pca_R/actions)

## 项目简介

改进 mixOmics 中的 PCA 实现,使用 partial SVD (rARPACK::svds) 提升大数据集的计算效率。

## 跨平台测试

本项目使用 GitHub Actions 在以下平台自动测试:

- ✅ Ubuntu (Linux)
- ✅ Windows
- ✅ macOS

## 安装依赖

```r
install.packages(c("testthat", "ggplot2", "dplyr", "tidyr", "knitr"))
install.packages("rARPACK")
BiocManager::install("mixOmics")
```

## 运行测试

```r
library(testthat)
test_dir("test/testthat")
```

## 项目结构

```
.
├── .github/
│   └── workflows/
│       └── R-CMD-check.yml    # GitHub Actions 配置
├── R/
│   ├── pca_pSVD.R             # 核心实现
│   ├── pca_asym.R
│   └── pca_em_woodbury.R
├── test/
│   └── testthat/
│       ├── test-pca_pSVD.R    # 单元测试
│       └── test-pca_woodbury.R
├── Data/
│   └── HaffinaCovidPBMC_30000cells_dense.rds
├── pca_pSVD.Rmd               # 分析报告
└── README.md
```

## 测试覆盖

- ✅ 基本功能测试
- ✅ 与 prcomp 对比
- ✅ 与 mixOmics::pca 对比
- ✅ 符号不变性处理
- ✅ 跨平台数值稳定性
- ✅ 边缘情况处理

## 作者

Kim-Anh Lê Cao & Saritha Kodikara

## 许可

MIT License
