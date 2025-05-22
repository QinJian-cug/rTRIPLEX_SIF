# install.packages(c("devtools","roxygen2","testthat","knitr"))
# install.packages("rstudioapi")
# rstudioapi::isAvailable("1.0.0")
# devtools::install_github("hadley/devtools")
# 写入r包的内置数据集
library(devtools)
library(roxygen2)
# data_test <- read.csv('./data/Inputvariable_month.csv')
setwd('E:/code/rTRIPLEX_SIF')
Inputvariable <- read.csv('./data/Inputvariable_DBS.csv')
Inputpara <- read.csv('./data/Inputpara_DBS.csv')
RC_PSII <- read.csv('./data/RC_PSII.csv')
onemonth_exam <- read.csv('./data/onemonth_exam.csv')
# usethis::use_data(data_test)
usethis::use_data(Inputpara, overwrite = TRUE)
usethis::use_data(Inputvariable, overwrite = TRUE)
usethis::use_data(RC_PSII, overwrite = TRUE)
usethis::use_data(onemonth_exam, overwrite = TRUE)

usethis::use_testthat() # 测试初始化
devtools::test() # 开始测试
devtools::document() # 更新帮助文档（每次更新R后都需要运行）
usethis::use_vignette("my-vignette")
devtools::check() # 执行检查

