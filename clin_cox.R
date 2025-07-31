cal_index_clin <- function(data = NULL, # ML.Dev.Prog.Sig, 函数计算结果
                        sets = NULL,
                        vars=NULL
){
  if (T) {
    library(survival)
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(ggsci)
    library(data.table)
    Sys.setenv(LANGUAGE = "en") # 显示英文报错信息
    options(stringsAsFactors = FALSE) # 禁止chr转成factor
  }
  message("--- C-index Calculate ---")
result <- data.frame()
index_result <- lapply(1:length(sets), function(x){
  cc <- data.frame(Cindex = sapply(vars,function(y){
    as.numeric(summary(coxph(as.formula(paste0('Surv(OS.time,OS)~',y)),data = data[[x]]))$concordance[1])
    }),
    se=sapply(vars,function(y){
      as.numeric(summary(coxph(as.formula(paste0('Surv(OS.time,OS)~',y)),data = data[[x]]))$concordance[2])
    }),
    Datasets=sets[x]
    ) %>% 
    rownames_to_column(var='Variables')
  return(cc)
})
index_result <- do.call(rbind, index_result)
}