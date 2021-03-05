# Test linear mixed model examples

#################################################################################################
# Gaussian examples
#################################################################################################

# https://www.jaredknowles.com/journal/2013/11/25/getting-started-with-mixed-effect-models-in-r
library(lme4) # load library
# load lmm.data data.frame
lmm.data <- read.table("http://bayes.acs.unt.edu:8083/BayesContent/class/Jon/R_SC/Module9/lmm.data.txt",
                       header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)

lmm_lme4 <- lmer(extro ~ open + agree + social + (1|school), data=lmm.data)
# random intercept standard deviation: 9.79
lmm_glmm = glmm(extro ~ open + agree + social + (1|school), data=lmm.data, 
                family = "gaussian", trace = 1, optim_options = optimControl(var_start = 100))
# Results:
## var_start needed to start large to avoid divergence issues
## Took 100 iterations
## fixef and logLik very similar between lme4 and glmm
## Random intercept variance and scale values very different between the two methods


# https://cran.r-project.org/web/packages/multilevelTools/vignettes/lmer-vignette.html
# data(aces_daily, package = "JWileymisc")
# can download data aces_daily from "https://github.com/JWiley/JWileymisc/blob/master/data/aces_daily.rda"
load("C:/Users/hheiling/Documents/GitHub/aces_daily.rda")
aces_lme4 <- lmer(NegAff ~ STRESS + (1 + STRESS | UserID),
                  data = aces_daily)
# Random intercept standard deviation: 0.206
aces_glmm = glmm(NegAff ~ STRESS + (1 + STRESS | UserID),
                 data = aces_daily, family = "gaussian",
                 optim_options = optimControl(var_start = 1)) # default var_start = 1
## Results: fixed effects and log-likelihood match very well; 
## log-like for glmm() is larger (closer to 0, less negative) than lme4 result
## random effect variance of glmm() result is larger than lme4


# ChickWieght dataset automatically provided by R
chick_lme4 <- lmer(weight ~ Time * Diet + (1 + Time | Chick), data=ChickWeight)
# random intercept standard deviation: 10.8 (variance 116)
chick_glmm = glmm(weight ~ Time * Diet + (1 + Time | Chick), data=ChickWeight, 
                  family = "gaussian", trace = 1,
                  optim_options = optimControl(var_start = 121))

# Exam data from library mlmRev, downloaded to personal folder
# Can download directly from "https://github.com/cran/mlmRev/blob/master/data/Exam.rda"
load("C:/Users/hheiling/Documents/GitHub/Exam.rda")
exam_lme4 = lmer(normexam ~ standLRT + schavg + (standLRT | school), data=Exam)
# random intercept standard deviation: 0.278
exam_glmm = glmm(normexam ~ standLRT + schavg + (standLRT | school), data=Exam, family = "gaussian")
# Results: fixed effects very similar; logLik very similar, with glmm logLik slightly less negative
#   random effect variances very similar; residual error in glmm smaller than in lme4 result


# Running glmm() on sleepstudy data
library(lme4)
# sleepstudy data fom lme4 package
data("sleepstudy")

# sleepstudy_sub = sleepstudy[which(sleepstudy$Subject %in% c(308,309,310,330,331,332)),]
fit_lme4 = lmer(formula = Reaction ~ Days + (Days | Subject), data = sleepstudy)
# random intercept standard deviation: 24.74 (variance = 612.1)
fitA = glmm(formula = Reaction ~ Days + (Days | Subject), data = sleepstudy,
            family = "gaussian", covar = "unstructured", 
            optim_options = optimControl(nMC_start = 1000, nMC_max = 10^4,var_start = 600), 
            trace = 1)



##################################################################################################