#!/usr/bin/env Rscript

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

suppressPackageStartupMessages({
  library(pak)
})

pak::pkg_install(c(
  "DOSE",
  "ReactomePA",
  "downloader",
  "enrichplot",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "org.Rn.eg.db",
  "org.Sc.sgd.db",
  "org.Mmu.eg.db",
  "msigdbr",
  "pathview",
  "ggridges",
  "ggplot2"
))
