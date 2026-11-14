clin_os$OS=as.numeric(clin_os$OS)
library(timeROC)
tROC <- timeROC(
  T = clin_os$OS.time,
  delta = clin_os$OS,
  marker = clin_os$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)
