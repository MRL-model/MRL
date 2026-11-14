library(timeROC)
library(survival)
clin_os$risk_score=t(exp_matrix)%*%cof_g[,1]
clin_os$OS.time=as.numeric(
  as.character(clin_oss$OS.time)
)

clin_os$OS=as.numeric(
  as.character(clin_os$OS)
)

clin_os$risk_score <- as.numeric(
  as.character(clin_os$risk_score)
)
clin_os <- clin_os[
  complete.cases(
    clin_os[, c( "OS.time", "OS", "risk_score")]
  ) &
    clin_os$OS.time >= 0 &
    clin_os$OS %in% c(0, 1),
  ,
  drop = FALSE
]
time_points=c(12, 36, 60)#month
tROC <- timeROC(
  T = clin_os$OS.time,
  delta = clin_os$OS,
  marker = clin_os$risk_score,
  cause = 1,
  weighting = "marginal",
  times = time_points,
  ROC = TRUE,
  iid = TRUE
)
