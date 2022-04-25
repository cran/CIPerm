## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup1-------------------------------------------------------------------
library(CIPerm)
x <- c(19, 22, 25, 26)
y <- c(23, 33, 40)

## ----table2-------------------------------------------------------------------
demo <- dset(x, y, returnData = TRUE)
knitr::kable(demo, digits = 2)

## ----pvals.p6-----------------------------------------------------------------
# Difference in means
pval(demo, tail = "Left", value = "m")

# Sum of treatment group
pval(demo, tail = "Left", value = "s")

# Difference in medians
pval(demo, tail = "Left", value = "d")

# Wilcoxon rank sum statistic
pval(demo, tail = "Left", value = "w")

## ----cint.p11-----------------------------------------------------------------
cint(demo, conf.level = 1-2/35, tail = "Left")

## ----setup2-------------------------------------------------------------------
wl1 <- c(1.72, 1.64, 1.74, 1.70, 1.82, 1.82, 1.90, 1.82, 2.08)
wl2 <- c(1.78, 1.86, 1.96, 2.00, 2.00, 1.96)
wl <- dset(wl1, wl2)

al1 <- c(1.24, 1.38, 1.36, 1.40, 1.38, 1.48, 1.38, 1.54, 1.56)
al2 <- c(1.14, 1.20, 1.30, 1.26, 1.28, 1.18)
al <- dset(al1, al2)

## ----pvals.p14----------------------------------------------------------------
pval(wl, tail = "Two", value = "s")
pval(al, tail = "Two", value = "s")

## ----cint.p15-----------------------------------------------------------------
cint(wl, conf.level = .95, tail = "Two")
cint(al, conf.level = .95, tail = "Two")

## ----montecarlo---------------------------------------------------------------
wl <- dset(wl1, wl2, nmc = 999)
al <- dset(al1, al2, nmc = 999)
pval(wl, tail = "Two", value = "s")
pval(al, tail = "Two", value = "s")
cint(wl, conf.level = .95, tail = "Two")
cint(al, conf.level = .95, tail = "Two")

