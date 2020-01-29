### test for survminer
#### how to use survminer ####
library("survminer")
# Fit survival curves
require("survival")
lung
fit <- survfit(Surv(time, status) ~ 1, data = lung)
# Drawing curves
ggsurvplot(fit, color = "#2E9FDF")
# Fit survival curves
require("survival")
fit<- survfit(Surv(time, status) ~ sex, data = lung)
# Drawing survival curves
x <- ggsurvplot(fit,pval = TRUE)
