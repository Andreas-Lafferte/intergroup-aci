#*******************************************************************************************************************
#
# 0. Identification ---------------------------------------------------
# Title: Data anlaysis for Research Paper
# Author: Andreas Laffert            
# Overview: Analysis of ELSOC 2018-2023      
# Date: 04-01-2025           
#
#******************************************************************************************************************

# 1. Packages  -----------------------------------------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse,
               sjmisc, 
               sjPlot,
               here,
               lavaan,
               psych,
               corrplot,
               ggdist,
               sjlabelled,
               semTools,
               gtools,
               summarytools,
               lme4,
               performance)

options(scipen=999)
options(survey.lonely.psu = "certainty")
rm(list = ls())

# 2. Data -----------------------------------------------------------------

load(file = here("input/data/proc/df_study1_wide.RData"))
load(file = here("input/data/proc/df_study1_long.RData"))

glimpse(df_study1_wide)

# 3. Analysis -------------------------------------------------------------


# 3.1 Descriptive ----

df_study1_wide %>% 
  select(starts_with(c("freq", "positiv_contc", "aci")), ess, sex) %>% 
  summarytools::dfSummary()

m0 <- lmer(aci ~ 1 + (1 | idencuesta), data = db_proc)

performance::icc(m0, by_group = T)
# ICC = 0.11 
# Un 11% de la varianza de aci se asocia a diferencias entre 
# individuos, mientras que un 89% de su varianza corresponde al cambio en el tiempo en los mismos individuos
# Es posble usar un RI-CLPM

# Estimation RI-CLPM

between <- ' 
# Componentes between
   RI_x =~ 1*freq_cont_lc1 + 1*freq_cont_lc2 + 1*freq_cont_lc3 + 1*freq_cont_lc4
   RI_y =~ 1*aci1 + 1*aci2 + 1*aci3 + 1*aci4
           '

within  <- '
# Componentes within
   cx1 =~ 1*freq_cont_lc1
   cx2 =~ 1*freq_cont_lc2
   cx3 =~ 1*freq_cont_lc3
   cx4 =~ 1*freq_cont_lc4
   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3
   cy4 =~ 1*aci4

# Constreñir las varianzas del error de medicion a 0
   freq_cont_lc1 ~~ 0*freq_cont_lc1
   freq_cont_lc2 ~~ 0*freq_cont_lc2
   freq_cont_lc3 ~~ 0*freq_cont_lc3
   freq_cont_lc4 ~~ 0*freq_cont_lc4
   aci1 ~~ 0*aci1
   aci2 ~~ 0*aci2
   aci3 ~~ 0*aci3
   aci4 ~~ 0*aci4
          '

efectos  <- '
# Estimar los efectos lagged
    cx2 ~ cx1 + cy1
    cx3 ~ cx2 + cy2
    cx4 ~ cx3 + cy3
    cy2 ~ cx1 + cy1
    cy3 ~ cx2 + cy2
    cy4 ~ cx3 + cy3
           '

covarianzas  <- '
# Covarianza entre los componentes within t=1
    cx1 ~~ cy1
   
# Covarianzas entre los residuos componente within
    cx2 ~~ cy2
    cx3 ~~ cy3
    cx4 ~~ cy4

# Varianzas residuales del componente within
    cx1 ~~ cx1
    cy1 ~~ cy1 
    cx2 ~~ cx2
    cy2 ~~ cy2 
    cx3 ~~ cx3 
    cy3 ~~ cy3 
    cx4 ~~ cx4 
    cy4 ~~ cy4 
    
# Varianza y covarianza entre RI
    RI_x ~~ RI_x
    RI_y ~~ RI_y
    RI_x ~~ RI_y
    
# Correlacion entre los RI y componentes within t=1 
    RI_x ~~ 0*cx1
    RI_x ~~ 0*cy1
    RI_y ~~ 0*cx1
    RI_y ~~ 0*cy1 
           '


fit1 <- lavaan(model = c(between, within, efectos, covarianzas),
                      data = df_study1_wide[,c(1:19)], 
                      estimator = "MLR", missing = "FIML",
                      meanstructure = T, int.ov.free = T)

summary(fit1, fit.measures = T, ci = T, standardized = T)

bwcomp <- ' 
           # Crear los componentes between
            RI_x =~ 1*freq_cont_lc1 + 1*freq_cont_lc2 + 1*freq_cont_lc3 + 1*freq_cont_lc4
            RI_y =~ 1*aci1 + 1*aci2 + 1*aci3 + 1*aci4
            
            # Crear los componentes within
            cx1 =~ 1*freq_cont_lc1
            cx2 =~ 1*freq_cont_lc2
            cx3 =~ 1*freq_cont_lc3
            cx4 =~ 1*freq_cont_lc4

            cy1 =~ 1*aci1
            cy2 =~ 1*aci2
            cy3 =~ 1*aci3
            cy4 =~ 1*aci4

            # constreñir las varianzas del error de medicion a cero
            freq_cont_lc1 ~~ 0*freq_cont_lc1
            freq_cont_lc2 ~~ 0*freq_cont_lc2
            freq_cont_lc3 ~~ 0*freq_cont_lc3
            freq_cont_lc4 ~~ 0*freq_cont_lc4

            aci1 ~~ 0*aci1
            aci2 ~~ 0*aci2
            aci3 ~~ 0*aci3
            aci4 ~~ 0*aci4
          '


varcov <- ' 
           # Estimar la covarianza entre los componentes within t=1
           cx1 ~~ cy1
           
           # Estimar las covarianzas entre los residuos del componente within
           cx2 ~~ cy2
           cx3 ~~ cy3
           cx4 ~~ cy4

           # Estimar las varianzas residuales del componente within
           cx1 ~~ cx1 # Varianzas
           cy1 ~~ cy1 
           cx2 ~~ cx2 # Varianzas residuales
           cy2 ~~ cy2 
           cx3 ~~ cx3 
           cy3 ~~ cy3 
           cx4 ~~ cx4 
           cy4 ~~ cy4 

           # Estimar la varianza y covarianza entre los RI. 
           RI_x ~~ RI_x
           RI_y ~~ RI_y
           RI_x ~~ RI_y
           
           # Fijar la correlacion entre los RI y componentes within t=1 a cero 
           RI_x ~~ 0*cx1
           RI_x ~~ 0*cy1
           RI_y ~~ 0*cx1
           RI_y ~~ 0*cy1 
          '


a1 <- ' # Sin constreñir
           cx2 ~ cx1
           cx3 ~ cx2
           cx4 ~ cx3
           cy2 ~ cy1
           cy3 ~ cy2
           cy4 ~ cy3
          '

a2 <- ' # Constreñido
           cx2 ~ a*cx1
           cx3 ~ a*cx2
           cx4 ~ a*cx3
           cy2 ~ d*cy1
           cy3 ~ d*cy2
           cy4 ~ d*cy3
          '

b1 <- ' # Sin constreñir
           cx2 ~ cx1
           cx3 ~ cx2 
           cx4 ~ cx3
           cy2 ~ cx1 + cy1
           cy3 ~ cx2 + cy2
           cy4 ~ cx3 + cy3
          '

b2 <- ' # Constreñido
           cx2 ~ a*cx1
           cx3 ~ a*cx2
           cx4 ~ a*cx3
           cy2 ~ c*cx1 + d*cy1
           cy3 ~ c*cx2 + d*cy2
           cy4 ~ c*cx3 + d*cy3
          '



df1 <- df_study1_wide %>% 
  select(idencuesta, starts_with(c("freq_cont", "aci")), sex, ess)

models <- c("a1","a2","b1","b2")

fit <- list()
for (i in models){
  fit[[i]] <- lavaan(model = c(bwcomp,get(i),varcov),
                     data = df1, 
                     estimator = "MLR", 
                     missing = "FIML",
                     meanstructure = T, 
                     int.ov.free = T)
}

gofdt <- list()
for (i in names(fit)){
  x <- fitMeasures(fit[[i]])[c("chisq.scaled", "df.scaled",
                               "pvalue.scaled", "cfi.scaled",
                               "tli.scaled", "rmsea.scaled",
                               "srmr_mplus", "aic",
                               "bic", "bic2",
                               "logl", "npar",
                               "scaling.factor.h0")]
  gofdt[[i]] <- setNames(as.numeric(x),
                         c("X2","df",
                           "pvalue","CFI",
                           "TLI","RMSEA",
                           "SRMR","AIC",
                           "BIC","aBIC",
                           "LL","par",
                           "LLcorrectf"))}
gofdt <- data.table(m=names(gofdt),dplyr::bind_rows(gofdt))

kableExtra::kable(gofdt,digits = 2,
                  format = "markdown")


summary(fit$b2, fit.measures = T, standardized = T)

