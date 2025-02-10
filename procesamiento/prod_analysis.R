# 0. Identification ---------------------------------------------------

# Title: Data analysis for research paper on Contact and Classism
# Institution: Centro de Estudios de Conflicto y Cohesión Social (COES)
# Responsible: Researcher

# Executive Summary: This script contains the code to create the analysis code for Contact and Classism
# Date: February 6, 2025

# 1. Packages  -----------------------------------------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse,
               sjmisc, 
               sjPlot,
               here,
               lavaan,
               psych,
               sjlabelled,
               semTools,
               gtools,
               summarytools,
               lme4,
               performance,
               data.table,
               effectsize)

options(scipen=999)
options(survey.lonely.psu = "certainty")
rm(list = ls())

# 2. Data -----------------------------------------------------------------

load(file = here("input/data/proc/df_study1_wide.RData"))
load(file = here("input/data/proc/df_study1_long.RData"))

glimpse(df_study1_wide)

df_study1_wide_or <- df_study1_wide

# 3. Analysis -------------------------------------------------------------

# function models comparation
gof.comp  <- function(data, pairs, measures = c("CFI", "TLI", "RMSEA", "SRMR", "AIC", 
                                               "BIC", "aBIC", "par", "LL")){
  comp <- list()
  for (i in 1:length(pairs)){
    gof <- data
    nest <- pairs[[i]][1]
    full <- pairs[[i]][2]
    delta <- NULL
    for (k in measures){
      delta[paste0(k,"_D")] <- gof[m==nest, get(k)] - gof[m==full, get(k)] }
    par_LLcorf_nest <- gof[m==nest,par]*gof[m==nest,LLcorrectf]
    par_LLcorf_full <- gof[m==full,par]*gof[m==full,LLcorrectf]
    delta["CD"] <- (par_LLcorf_nest-par_LLcorf_full)/delta["par_D"]
    delta["TRd"] <- (-2*delta["LL_D"])/delta["CD"]
    delta["TRd_df"] <- gof[m==full, "par"] - gof[m==nest, "par"]
    delta["TRd_pvalue"] <- pchisq(as.numeric(delta["TRd"]),
                                  as.numeric(delta["TRd_df"]), lower.tail = F)
    comp[[paste0(nest," vs. ",full,sep="")]] <- delta }
  comp <- data.table(comp=names(comp),dplyr::bind_rows(comp))
  return(comp)
}


df_study1_wide <- df_study1_wide_or[,c(1:5,10:19)]

# 3.1 Descriptive -------------------------------------------------------------

df_study1_wide %>% 
  select(-idencuesta) %>% 
  summarytools::dfSummary()

M <- df_study1_wide %>% 
  remove_all_labels() %>% 
  as.data.frame() %>% 
  select(-idencuesta)

M$sex <- as.numeric(M$sex)

sjPlot::tab_corr(M, 
                 na.deletion = "pairwise", 
                 corr.method = "pearson", 
                 triangle = "lower")

m0 <- lmer(aci ~ 1 + (1 | idencuesta), data = db_proc)

performance::icc(m0, by_group = T)
# ICC = 0.11 

# 11% of the variance of aci is associated with differences between individuals, while 89% of its variance corresponds to change over time within the same individuals. 
# individuals, while 89% of its variance corresponds to the change over time in the same individuals.
# Is it possible to use an RI-CLPM ? ....


# 3.1 CLPM without controls --------------------------------------------------------------

within  <- '
# Componentes within
   cx1 =~ 1*freq_cont_lc1
   cx2 =~ 1*freq_cont_lc2
   cx3 =~ 1*freq_cont_lc3
   cx4 =~ 1*freq_cont_lc4
   
   cz1 =~ 1*positiv_contc_lc1
   cz2 =~ 1*positiv_contc_lc2
   cz3 =~ 1*positiv_contc_lc3
   cz4 =~ 1*positiv_contc_lc4
   
   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3
   cy4 =~ 1*aci4

# Constreñir las varianzas del error de medicion a 0
   freq_cont_lc1 ~~ 0*freq_cont_lc1
   freq_cont_lc2 ~~ 0*freq_cont_lc2
   freq_cont_lc3 ~~ 0*freq_cont_lc3
   freq_cont_lc4 ~~ 0*freq_cont_lc4
   
   positiv_contc_lc1 ~~ 0*positiv_contc_lc1
   positiv_contc_lc2 ~~ 0*positiv_contc_lc2
   positiv_contc_lc3 ~~ 0*positiv_contc_lc3
   positiv_contc_lc4 ~~ 0*positiv_contc_lc4
   
   aci1 ~~ 0*aci1
   aci2 ~~ 0*aci2
   aci3 ~~ 0*aci3
   aci4 ~~ 0*aci4
          '


covarianzas  <- '
# Covarianza entre los componentes within t=1 con corr entre predictores
    cx1 ~~ cy1
    cz1 ~~ cy1
    cx1 ~~ cz1

# Covarianzas entre los residuos componente within con corr entre predictores
    cx2 ~~ cy2
    cx3 ~~ cy3
    cx4 ~~ cy4
    
    cz2 ~~ cy2
    cz3 ~~ cy3
    cz4 ~~ cy4
    
    cx2 ~~ cz2
    cx3 ~~ cz3
    cx4 ~~ cz4

# Varianzas residuales del componente within 
    cx1 ~~ cx1
    cz1 ~~ cz1
    cy1 ~~ cy1 
    cx2 ~~ cx2
    cz2 ~~ cz2
    cy2 ~~ cy2 
    cx3 ~~ cx3
    cz3 ~~ cz3
    cy3 ~~ cy3 
    cx4 ~~ cx4 
    cz4 ~~ cz4
    cy4 ~~ cy4 
           '

## Model A: Autoregressive ----

a1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
    
    cz2 ~ cz1 
    cz3 ~ cz2 
    cz4 ~ cz3
  
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '

a2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
    
    cz2 ~ b*cz1 
    cz3 ~ b*cz2 
    cz4 ~ b*cz3
  
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model B: Forward ----

b1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
    
    cz2 ~ cz1 
    cz3 ~ cz2 
    cz4 ~ cz3
  
    cy2 ~ cx1 + cz2 + cy1
    cy3 ~ cx2 + cz3 + cy2
    cy4 ~ cx3 + cz4 + cy3
           '
b2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
    
    cz2 ~ b*cz1 
    cz3 ~ b*cz2 
    cz4 ~ b*cz3
  
    cy2 ~ d*cx1 + e*cz2 + c*cy1
    cy3 ~ d*cx2 + e*cz3 + c*cy2
    cy4 ~ d*cx3 + e*cz4 + c*cy3
           '

## Model C: Backward ----

c1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1
    cx3 ~ cx2 + cz2 + cy2
    cx4 ~ cx3 + cz3 + cy3
    
    cz2 ~ cz1 + cx1 + cy1
    cz3 ~ cz2 + cx2 + cy2
    cz4 ~ cz3 + cx3 + cy3
  
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '
c2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + f*cz2 + g*cy2
    cx2 ~ a*cx1 + f*cz1 + g*cy1
    cx4 ~ a*cx3 + f*cz3 + g*cy3
    
    cz2 ~ b*cz1 + h*cx1 + i*cy1
    cz3 ~ b*cz2 + h*cx2 + i*cy2
    cz4 ~ b*cz3 + h*cx3 + i*cy3
  
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model D: Bidirectional ---- 

d1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1
    cx3 ~ cx2 + cz2 + cy2
    cx4 ~ cx3 + cz3 + cy3
    
    cz2 ~ cz1 + cx1 + cy1
    cz3 ~ cz2 + cx2 + cy2
    cz4 ~ cz3 + cx3 + cy3
  
    cy2 ~ cx1 + cz2 + cy1
    cy3 ~ cx2 + cz3 + cy2
    cy4 ~ cx3 + cz4 + cy3
           '


d2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + f*cz2 + g*cy2
    cx2 ~ a*cx1 + f*cz1 + g*cy1
    cx4 ~ a*cx3 + f*cz3 + g*cy3
    
    cz2 ~ b*cz1 + h*cx1 + i*cy1
    cz3 ~ b*cz2 + h*cx2 + i*cy2
    cz4 ~ b*cz3 + h*cx3 + i*cy3
  
    cy2 ~ d*cx1 + e*cz2 + c*cy1
    cy3 ~ d*cx2 + e*cz3 + c*cy2
    cy4 ~ d*cx3 + e*cz4 + c*cy3
           '

## Estimation CLPM ----

# models
models <- c("a1","a2","b1","b2","c1","c2","d1","d2")

m <- list()
for (i in models){
  m[[i]] <- lavaan(model = c(within,get(i),covarianzas),
                     data = df_study1_wide, 
                     estimator = "MLR", 
                     missing = "FIML",
                     meanstructure = T, 
                     int.ov.free = T)
}

# fit measures
gofdt <- list()
for (i in names(m)){
  x <- fitMeasures(m[[i]])[c("chisq.scaled", "df.scaled",
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

gofdt <- gofdt %>% mutate(
  interpret_CFI = effectsize::interpret_cfi(CFI),
  interpret_RMSEA = effectsize::interpret_rmsea(RMSEA),
) %>% relocate(interpret_CFI, .after = CFI) %>%
  relocate(interpret_RMSEA, .after = RMSEA)


# comparacion modelos

# testear si efectos son iguales en el tiempo

comp1 <- gof.comp(data = gofdt, 
                  pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# testear direccion relationes
comp2 <- gof.comp(data = gofdt, 
                  pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                               c("b2","d2"), c("c2","d2"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# conclusion, modelo d mejor modelo

summary(m[["d1"]], fit.measures = T, ci = T, standardized = T)

summary(m[["d2"]], fit.measures = T, ci = T, standardized = T)


# 3.2 CLPM with controls -------------------------------------------------------------

within  <- '
# Componentes within
   cx1 =~ 1*freq_cont_lc1
   cx2 =~ 1*freq_cont_lc2
   cx3 =~ 1*freq_cont_lc3
   cx4 =~ 1*freq_cont_lc4
   
   cz1 =~ 1*positiv_contc_lc1
   cz2 =~ 1*positiv_contc_lc2
   cz3 =~ 1*positiv_contc_lc3
   cz4 =~ 1*positiv_contc_lc4
   
   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3
   cy4 =~ 1*aci4

# Constreñir las varianzas del error de medicion a 0
   freq_cont_lc1 ~~ 0*freq_cont_lc1
   freq_cont_lc2 ~~ 0*freq_cont_lc2
   freq_cont_lc3 ~~ 0*freq_cont_lc3
   freq_cont_lc4 ~~ 0*freq_cont_lc4
   
   positiv_contc_lc1 ~~ 0*positiv_contc_lc1
   positiv_contc_lc2 ~~ 0*positiv_contc_lc2
   positiv_contc_lc3 ~~ 0*positiv_contc_lc3
   positiv_contc_lc4 ~~ 0*positiv_contc_lc4
   
   aci1 ~~ 0*aci1
   aci2 ~~ 0*aci2
   aci3 ~~ 0*aci3
   aci4 ~~ 0*aci4
          '

covarianzas  <- '
# Covarianza entre los componentes within t=1 con corr entre predictores
    cx1 ~~ cy1
    cz1 ~~ cy1
    cx1 ~~ cz1

# Covarianzas entre los residuos componente within con corr entre predictores
    cx2 ~~ cy2
    cx3 ~~ cy3
    cx4 ~~ cy4
    
    cz2 ~~ cy2
    cz3 ~~ cy3
    cz4 ~~ cy4
    
    cx2 ~~ cz2
    cx3 ~~ cz3
    cx4 ~~ cz4

# Varianzas residuales del componente within 
    cx1 ~~ cx1
    cz1 ~~ cz1
    cy1 ~~ cy1 
    cx2 ~~ cx2
    cz2 ~~ cz2
    cy2 ~~ cy2 
    cx3 ~~ cx3
    cz3 ~~ cz3
    cy3 ~~ cy3 
    cx4 ~~ cx4 
    cz4 ~~ cz4
    cy4 ~~ cy4 
           '

## Model A: Autoregressive ----

a1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
    
    cz2 ~ cz1 
    cz3 ~ cz2 
    cz4 ~ cz3
  
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '

a2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
    
    cz2 ~ b*cz1 
    cz3 ~ b*cz2 
    cz4 ~ b*cz3
  
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model B: Forward ----

b1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
    
    cz2 ~ cz1 
    cz3 ~ cz2 
    cz4 ~ cz3
  
    cy2 ~ cx1 + cz2 + cy1 + sex + ess
    cy3 ~ cx2 + cz3 + cy2 + sex + ess
    cy4 ~ cx3 + cz4 + cy3 + sex + ess
           '
b2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
    
    cz2 ~ b*cz1 
    cz3 ~ b*cz2 
    cz4 ~ b*cz3
  
    cy2 ~ d*cx1 + e*cz2 + c*cy1 + sexdep*sex + essdep*ess
    cy3 ~ d*cx2 + e*cz3 + c*cy2 + sexdep*sex + essdep*ess
    cy4 ~ d*cx3 + e*cz4 + c*cy3 + sexdep*sex + essdep*ess
           '

## Model C: Backward ----

c1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1 + sex + ess
    cx3 ~ cx2 + cz2 + cy2 + sex + ess
    cx4 ~ cx3 + cz3 + cy3 + sex + ess
    
    cz2 ~ cz1 + cx1 + cy1 + sex + ess
    cz3 ~ cz2 + cx2 + cy2 + sex + ess
    cz4 ~ cz3 + cx3 + cy3 + sex + ess
  
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '
c2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + f*cz2 + g*cy2 + sexindepx*sex + sexindepx*ess
    cx2 ~ a*cx1 + f*cz1 + g*cy1 + sexindepx*sex + sexindepx*ess
    cx4 ~ a*cx3 + f*cz3 + g*cy3 + sexindepx*sex + sexindepx*ess
    
    cz2 ~ b*cz1 + h*cx1 + i*cy1 + sexindepz*sex + sexindepz*ess
    cz3 ~ b*cz2 + h*cx2 + i*cy2 + sexindepz*sex + sexindepz*ess
    cz4 ~ b*cz3 + h*cx3 + i*cy3 + sexindepz*sex + sexindepz*ess
  
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model D: Bidirectional ---- 

d1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1 + sex + ess
    cx3 ~ cx2 + cz2 + cy2 + sex + ess
    cx4 ~ cx3 + cz3 + cy3 + sex + ess
    
    cz2 ~ cz1 + cx1 + cy1 + sex + ess
    cz3 ~ cz2 + cx2 + cy2 + sex + ess
    cz4 ~ cz3 + cx3 + cy3 + sex + ess
  
    cy2 ~ cx1 + cz2 + cy1 + sex + ess
    cy3 ~ cx2 + cz3 + cy2 + sex + ess
    cy4 ~ cx3 + cz4 + cy3 + sex + ess
           '


d2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + f*cz2 + g*cy2 + sexindepx*sex + sexindepx*ess
    cx2 ~ a*cx1 + f*cz1 + g*cy1 + sexindepx*sex + sexindepx*ess
    cx4 ~ a*cx3 + f*cz3 + g*cy3 + sexindepx*sex + sexindepx*ess
    
    cz2 ~ b*cz1 + h*cx1 + i*cy1 + sexindepz*sex + sexindepz*ess
    cz3 ~ b*cz2 + h*cx2 + i*cy2 + sexindepz*sex + sexindepz*ess
    cz4 ~ b*cz3 + h*cx3 + i*cy3 + sexindepz*sex + sexindepz*ess
  
    cy2 ~ d*cx1 + e*cz2 + c*cy1 + sexdep*sex + essdep*ess
    cy3 ~ d*cx2 + e*cz3 + c*cy2 + sexdep*sex + essdep*ess
    cy4 ~ d*cx3 + e*cz4 + c*cy3 + sexdep*sex + essdep*ess
           '

## Estimation CLPM ----

# models
models <- c("a1","a2","b1","b2","c1","c2","d1","d2")

m_con <- list()
for (i in models){
  m_con[[i]] <- lavaan(model = c(within,get(i),covarianzas),
                       data = df_study1_wide, 
                       estimator = "MLR", 
                       missing = "FIML",
                       meanstructure = T, 
                       int.ov.free = T)
}

# fit measures
gofdt_con <- list()
for (i in names(m_con)){
  x <- fitMeasures(m_con[[i]])[c("chisq.scaled", "df.scaled",
                                 "pvalue.scaled", "cfi.scaled",
                                 "tli.scaled", "rmsea.scaled",
                                 "srmr_mplus", "aic",
                                 "bic", "bic2",
                                 "logl", "npar",
                                 "scaling.factor.h0")]
  gofdt_con[[i]] <- setNames(as.numeric(x),
                         c("X2","df",
                           "pvalue","CFI",
                           "TLI","RMSEA",
                           "SRMR","AIC",
                           "BIC","aBIC",
                           "LL","par",
                           "LLcorrectf"))}

gofdt_con <- data.table(m=names(gofdt_con),dplyr::bind_rows(gofdt_con))


# comparacion modelos

# testear si efectos son iguales en el tiempo

comp3 <- gof.comp(data = gofdt_con, 
                  pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# testear direccion relationes
comp4 <- gof.comp(data = gofdt_con, 
                  pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                               c("b2","d2"), c("c2","d2"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# Con o sin controles?

lavTestLRT(m[["d2"]], m_con[["d2"]]) # modelo con controles es sign mejor que sin controles


# 3.3 RI- CLPM without controls -------------------------------------------------------------

between <- ' 
# Componentes between
   RI_x =~ 1*freq_cont_lc1 + 1*freq_cont_lc2 + 1*freq_cont_lc3 + 1*freq_cont_lc4
   RI_z =~ 1*positiv_contc_lc1 + 1*positiv_contc_lc2 + 1*positiv_contc_lc3 + 1*positiv_contc_lc4
   RI_y =~ 1*aci1 + 1*aci2 + 1*aci3 + 1*aci4
           '

within  <- '
# Componentes within
   cx1 =~ 1*freq_cont_lc1
   cx2 =~ 1*freq_cont_lc2
   cx3 =~ 1*freq_cont_lc3
   cx4 =~ 1*freq_cont_lc4
   
   cz1 =~ 1*positiv_contc_lc1
   cz2 =~ 1*positiv_contc_lc2
   cz3 =~ 1*positiv_contc_lc3
   cz4 =~ 1*positiv_contc_lc4
   
   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3
   cy4 =~ 1*aci4

# Constreñir las varianzas del error de medicion a 0
   freq_cont_lc1 ~~ 0*freq_cont_lc1
   freq_cont_lc2 ~~ 0*freq_cont_lc2
   freq_cont_lc3 ~~ 0*freq_cont_lc3
   freq_cont_lc4 ~~ 0*freq_cont_lc4
   
   positiv_contc_lc1 ~~ 0*positiv_contc_lc1
   positiv_contc_lc2 ~~ 0*positiv_contc_lc2
   positiv_contc_lc3 ~~ 0*positiv_contc_lc3
   positiv_contc_lc4 ~~ 0*positiv_contc_lc4
   
   aci1 ~~ 0*aci1
   aci2 ~~ 0*aci2
   aci3 ~~ 0*aci3
   aci4 ~~ 0*aci4
          '

covarianzas  <- '
# Covarianza entre los componentes within t=1 con corr entre predictores
    cx1 ~~ cy1
    cz1 ~~ cy1
    cx1 ~~ cz1

# Covarianzas entre los residuos componente within con corr entre predictores
    cx2 ~~ cy2
    cx3 ~~ cy3
    cx4 ~~ cy4
    
    cz2 ~~ cy2
    cz3 ~~ cy3
    cz4 ~~ cy4
    
    cx2 ~~ cz2
    cx3 ~~ cz3
    cx4 ~~ cz4

# Varianzas residuales del componente within 
    cx1 ~~ cx1
    cz1 ~~ cz1
    cy1 ~~ cy1 
    cx2 ~~ cx2
    cz2 ~~ cz2
    cy2 ~~ cy2 
    cx3 ~~ cx3
    cz3 ~~ cz3
    cy3 ~~ cy3 
    cx4 ~~ cx4 
    cz4 ~~ cz4
    cy4 ~~ cy4 
    
# Varianza y covarianza entre RI
    RI_x ~~ RI_x
    RI_z ~~ RI_z
    RI_y ~~ RI_y
    RI_x ~~ RI_y
    RI_x ~~ RI_z
    RI_z ~~ RI_y
    
# Correlacion entre los RI y componentes within t=1 
    RI_x ~~ 0*cx1
    RI_x ~~ 0*cz1
    RI_x ~~ 0*cy1
    
    RI_z ~~ 0*cx1
    RI_z ~~ 0*cz1
    RI_z ~~ 0*cy1
    
    RI_y ~~ 0*cx1
    RI_y ~~ 0*cz1
    RI_y ~~ 0*cy1 
           '

## Model A: Autoregressive ----

a1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
    
    cz2 ~ cz1 
    cz3 ~ cz2 
    cz4 ~ cz3
  
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '

a2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
    
    cz2 ~ b*cz1 
    cz3 ~ b*cz2 
    cz4 ~ b*cz3
  
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model B: Forward ----

b1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
    
    cz2 ~ cz1 
    cz3 ~ cz2 
    cz4 ~ cz3
  
    cy2 ~ cx1 + cz2 + cy1
    cy3 ~ cx2 + cz3 + cy2
    cy4 ~ cx3 + cz4 + cy3
           '
b2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
    
    cz2 ~ b*cz1 
    cz3 ~ b*cz2 
    cz4 ~ b*cz3
  
    cy2 ~ d*cx1 + e*cz2 + c*cy1
    cy3 ~ d*cx2 + e*cz3 + c*cy2
    cy4 ~ d*cx3 + e*cz4 + c*cy3
           '

## Model C: Backward ----

c1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1
    cx3 ~ cx2 + cz2 + cy2
    cx4 ~ cx3 + cz3 + cy3
    
    cz2 ~ cz1 + cx1 + cy1
    cz3 ~ cz2 + cx2 + cy2
    cz4 ~ cz3 + cx3 + cy3
  
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '
c2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + f*cz2 + g*cy2
    cx2 ~ a*cx1 + f*cz1 + g*cy1
    cx4 ~ a*cx3 + f*cz3 + g*cy3
    
    cz2 ~ b*cz1 + h*cx1 + i*cy1
    cz3 ~ b*cz2 + h*cx2 + i*cy2
    cz4 ~ b*cz3 + h*cx3 + i*cy3
  
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model D: Bidirectional ---- 

d1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1
    cx3 ~ cx2 + cz2 + cy2
    cx4 ~ cx3 + cz3 + cy3
    
    cz2 ~ cz1 + cx1 + cy1
    cz3 ~ cz2 + cx2 + cy2
    cz4 ~ cz3 + cx3 + cy3
  
    cy2 ~ cx1 + cz2 + cy1
    cy3 ~ cx2 + cz3 + cy2
    cy4 ~ cx3 + cz4 + cy3
           '


d2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + f*cz2 + g*cy2
    cx2 ~ a*cx1 + f*cz1 + g*cy1
    cx4 ~ a*cx3 + f*cz3 + g*cy3
    
    cz2 ~ b*cz1 + h*cx1 + i*cy1
    cz3 ~ b*cz2 + h*cx2 + i*cy2
    cz4 ~ b*cz3 + h*cx3 + i*cy3
  
    cy2 ~ d*cx1 + e*cz2 + c*cy1
    cy3 ~ d*cx2 + e*cz3 + c*cy2
    cy4 ~ d*cx3 + e*cz4 + c*cy3
           '

## Estimation RI-CLPM ----

# models
models <- c("a1","a2","b1","b2","c1","c2","d1","d2")

m_ri <- list()
for (i in models){
  m_ri[[i]] <- lavaan(model = c(between,within,get(i),covarianzas),
                      data = df_study1_wide, 
                      estimator = "MLR", 
                      missing = "FIML",
                      meanstructure = T, 
                      int.ov.free = T)
}

# fit measures
gofdt_ri <- list()
for (i in names(m_ri)){
  x <- fitMeasures(m_ri[[i]])[c("chisq.scaled", "df.scaled",
                                "pvalue.scaled", "cfi.scaled",
                                "tli.scaled", "rmsea.scaled",
                                "srmr_mplus", "aic",
                                "bic", "bic2",
                                "logl", "npar",
                                "scaling.factor.h0")]
  gofdt_ri[[i]] <- setNames(as.numeric(x),
                            c("X2","df",
                              "pvalue","CFI",
                              "TLI","RMSEA",
                              "SRMR","AIC",
                              "BIC","aBIC",
                              "LL","par",
                              "LLcorrectf"))}

gofdt_ri <- data.table(m=names(gofdt_ri),dplyr::bind_rows(gofdt_ri))

gofdt_ri <- gofdt_ri %>% mutate(
  interpret_CFI = effectsize::interpret_cfi(CFI),
  interpret_RMSEA = effectsize::interpret_rmsea(RMSEA),
) %>% relocate(interpret_CFI, .after = CFI) %>%
  relocate(interpret_RMSEA, .after = RMSEA)


# comparacion modelos

# testear si efectos son iguales en el tiempo

comp5 <- gof.comp(data = gofdt_ri, 
                  pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# testear direccion relationes
comp6 <- gof.comp(data = gofdt_ri, 
                  pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                               c("b2","d2"), c("c2","d2"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# conclusion, modelo d mejor modelo

summary(m_ri[["d1"]], fit.measures = T, ci = T, standardized = T)

summary(m_ri[["d2"]], fit.measures = T, ci = T, standardized = T)


# 3.3 RI- CLPM with controls -------------------------------------------------------------

between <- ' 
# Componentes between
   RI_x =~ 1*freq_cont_lc1 + 1*freq_cont_lc2 + 1*freq_cont_lc3 + 1*freq_cont_lc4
   RI_z =~ 1*positiv_contc_lc1 + 1*positiv_contc_lc2 + 1*positiv_contc_lc3 + 1*positiv_contc_lc4
   RI_y =~ 1*aci1 + 1*aci2 + 1*aci3 + 1*aci4
           '

within  <- '
# Componentes within
   cx1 =~ 1*freq_cont_lc1
   cx2 =~ 1*freq_cont_lc2
   cx3 =~ 1*freq_cont_lc3
   cx4 =~ 1*freq_cont_lc4
   
   cz1 =~ 1*positiv_contc_lc1
   cz2 =~ 1*positiv_contc_lc2
   cz3 =~ 1*positiv_contc_lc3
   cz4 =~ 1*positiv_contc_lc4
   
   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3
   cy4 =~ 1*aci4

# Constreñir las varianzas del error de medicion a 0
   freq_cont_lc1 ~~ 0*freq_cont_lc1
   freq_cont_lc2 ~~ 0*freq_cont_lc2
   freq_cont_lc3 ~~ 0*freq_cont_lc3
   freq_cont_lc4 ~~ 0*freq_cont_lc4
   
   positiv_contc_lc1 ~~ 0*positiv_contc_lc1
   positiv_contc_lc2 ~~ 0*positiv_contc_lc2
   positiv_contc_lc3 ~~ 0*positiv_contc_lc3
   positiv_contc_lc4 ~~ 0*positiv_contc_lc4
   
   aci1 ~~ 0*aci1
   aci2 ~~ 0*aci2
   aci3 ~~ 0*aci3
   aci4 ~~ 0*aci4
          '

covarianzas  <- '
# Covarianza entre los componentes within t=1 con corr entre predictores
    cx1 ~~ cy1
    cz1 ~~ cy1
    cx1 ~~ cz1

# Covarianzas entre los residuos componente within con corr entre predictores
    cx2 ~~ cy2
    cx3 ~~ cy3
    cx4 ~~ cy4
    
    cz2 ~~ cy2
    cz3 ~~ cy3
    cz4 ~~ cy4
    
    cx2 ~~ cz2
    cx3 ~~ cz3
    cx4 ~~ cz4

# Varianzas residuales del componente within 
    cx1 ~~ cx1
    cz1 ~~ cz1
    cy1 ~~ cy1 
    cx2 ~~ cx2
    cz2 ~~ cz2
    cy2 ~~ cy2 
    cx3 ~~ cx3
    cz3 ~~ cz3
    cy3 ~~ cy3 
    cx4 ~~ cx4 
    cz4 ~~ cz4
    cy4 ~~ cy4 
    
# Varianza y covarianza entre RI
    RI_x ~~ RI_x
    RI_z ~~ RI_z
    RI_y ~~ RI_y
    RI_x ~~ RI_y
    RI_x ~~ RI_z
    RI_z ~~ RI_y
    
# Correlacion entre los RI y componentes within t=1 
    RI_x ~~ 0*cx1
    RI_x ~~ 0*cz1
    RI_x ~~ 0*cy1
    
    RI_z ~~ 0*cx1
    RI_z ~~ 0*cz1
    RI_z ~~ 0*cy1
    
    RI_y ~~ 0*cx1
    RI_y ~~ 0*cz1
    RI_y ~~ 0*cy1 
           '

## Model A: Autoregressive ----
a1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
    
    cz2 ~ cz1 
    cz3 ~ cz2 
    cz4 ~ cz3
  
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '

a2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
    
    cz2 ~ b*cz1 
    cz3 ~ b*cz2 
    cz4 ~ b*cz3
  
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model B: Forward ----

b1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
    
    cz2 ~ cz1 
    cz3 ~ cz2 
    cz4 ~ cz3
  
    cy2 ~ cx1 + cz2 + cy1 + sex + ess
    cy3 ~ cx2 + cz3 + cy2 + sex + ess
    cy4 ~ cx3 + cz4 + cy3 + sex + ess
           '
b2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
    
    cz2 ~ b*cz1 
    cz3 ~ b*cz2 
    cz4 ~ b*cz3
  
    cy2 ~ d*cx1 + e*cz2 + c*cy1 + sexdep*sex + essdep*ess
    cy3 ~ d*cx2 + e*cz3 + c*cy2 + sexdep*sex + essdep*ess
    cy4 ~ d*cx3 + e*cz4 + c*cy3 + sexdep*sex + essdep*ess
           '

## Model C: Backward ----

c1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1 + sex + ess
    cx3 ~ cx2 + cz2 + cy2 + sex + ess
    cx4 ~ cx3 + cz3 + cy3 + sex + ess
    
    cz2 ~ cz1 + cx1 + cy1 + sex + ess
    cz3 ~ cz2 + cx2 + cy2 + sex + ess
    cz4 ~ cz3 + cx3 + cy3 + sex + ess
  
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '
c2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + f*cz2 + g*cy2 + sexindepx*sex + sexindepx*ess
    cx2 ~ a*cx1 + f*cz1 + g*cy1 + sexindepx*sex + sexindepx*ess
    cx4 ~ a*cx3 + f*cz3 + g*cy3 + sexindepx*sex + sexindepx*ess
    
    cz2 ~ b*cz1 + h*cx1 + i*cy1 + sexindepz*sex + sexindepz*ess
    cz3 ~ b*cz2 + h*cx2 + i*cy2 + sexindepz*sex + sexindepz*ess
    cz4 ~ b*cz3 + h*cx3 + i*cy3 + sexindepz*sex + sexindepz*ess
  
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model D: Bidirectional ---- 

d1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1 + sex + ess
    cx3 ~ cx2 + cz2 + cy2 + sex + ess
    cx4 ~ cx3 + cz3 + cy3 + sex + ess
    
    cz2 ~ cz1 + cx1 + cy1 + sex + ess
    cz3 ~ cz2 + cx2 + cy2 + sex + ess
    cz4 ~ cz3 + cx3 + cy3 + sex + ess
  
    cy2 ~ cx1 + cz2 + cy1 + sex + ess
    cy3 ~ cx2 + cz3 + cy2 + sex + ess
    cy4 ~ cx3 + cz4 + cy3 + sex + ess
           '


d2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + f*cz2 + g*cy2 + sexindepx*sex + sexindepx*ess
    cx2 ~ a*cx1 + f*cz1 + g*cy1 + sexindepx*sex + sexindepx*ess
    cx4 ~ a*cx3 + f*cz3 + g*cy3 + sexindepx*sex + sexindepx*ess
    
    cz2 ~ b*cz1 + h*cx1 + i*cy1 + sexindepz*sex + sexindepz*ess
    cz3 ~ b*cz2 + h*cx2 + i*cy2 + sexindepz*sex + sexindepz*ess
    cz4 ~ b*cz3 + h*cx3 + i*cy3 + sexindepz*sex + sexindepz*ess
  
    cy2 ~ d*cx1 + e*cz2 + c*cy1 + sexdep*sex + essdep*ess
    cy3 ~ d*cx2 + e*cz3 + c*cy2 + sexdep*sex + essdep*ess
    cy4 ~ d*cx3 + e*cz4 + c*cy3 + sexdep*sex + essdep*ess
           '
## Estimation RI-CLPM ----

# models
models <- c("a1","a2","b1","b2","c1","c2","d1","d2")

m_ri_con <- list()
for (i in models){
  m_ri_con[[i]] <- lavaan(model = c(between,within,get(i),covarianzas),
                         data = df_study1_wide, 
                         estimator = "MLR", 
                         missing = "FIML",
                         meanstructure = T, 
                         int.ov.free = T)
}

# fit measures
gofdt_ri_con <- list()
for (i in names(m_ri_con)){
  x <- fitMeasures(m_ri_con[[i]])[c("chisq.scaled", "df.scaled",
                                   "pvalue.scaled", "cfi.scaled",
                                   "tli.scaled", "rmsea.scaled",
                                   "srmr_mplus", "aic",
                                   "bic", "bic2",
                                   "logl", "npar",
                                   "scaling.factor.h0")]
  gofdt_ri_con[[i]] <- setNames(as.numeric(x),
                               c("X2","df",
                                 "pvalue","CFI",
                                 "TLI","RMSEA",
                                 "SRMR","AIC",
                                 "BIC","aBIC",
                                 "LL","par",
                                 "LLcorrectf"))}

gofdt_ri_con <- data.table(m=names(gofdt_ri_con),dplyr::bind_rows(gofdt_ri_con))

gofdt_ri_con <- gofdt_ri_con %>% mutate(
  interpret_CFI = effectsize::interpret_cfi(CFI),
  interpret_RMSEA = effectsize::interpret_rmsea(RMSEA),
) %>% relocate(interpret_CFI, .after = CFI) %>%
  relocate(interpret_RMSEA, .after = RMSEA)


# comparacion modelos

# testear si efectos son iguales en el tiempo

comp7 <- gof.comp(data = gofdt_ri_con, 
                  pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# testear direccion relationes
comp8 <- gof.comp(data = gofdt_ri_con, 
                  pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                               c("b2","d2"), c("c2","d2"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# conclusion, modelo d mejor modelo

summary(m_ri_con[["d1"]], fit.measures = T, ci = T, standardized = T)

summary(m_ri_con[["d2"]], fit.measures = T, ci = T, standardized = T)

# Con o sin controles?

lavTestLRT(m_ri[["d2"]], m_ri_con[["d2"]]) # modelo con controles es sign mejor que sin controles


# RI-CLPM d15 y d18 sin controles -------------------------------------------

df_study1_wide <- df_study1_wide_or[,c(1,6:9,14:19)]

between <- ' 
# Componentes between
   RI_x =~ 1*freq_neg_cont_lc1 + 1*freq_neg_cont_lc2 + 1*freq_neg_cont_lc3 + 1*freq_neg_cont_lc4
   RI_y =~ 1*aci1 + 1*aci2 + 1*aci3 + 1*aci4
           '

within  <- '
# Componentes within
   cx1 =~ 1*freq_neg_cont_lc1
   cx2 =~ 1*freq_neg_cont_lc2
   cx3 =~ 1*freq_neg_cont_lc3
   cx4 =~ 1*freq_neg_cont_lc4
   
   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3
   cy4 =~ 1*aci4

# Constreñir las varianzas del error de medicion a 0
   freq_neg_cont_lc1 ~~ 0*freq_neg_cont_lc1
   freq_neg_cont_lc2 ~~ 0*freq_neg_cont_lc2
   freq_neg_cont_lc3 ~~ 0*freq_neg_cont_lc3
   freq_neg_cont_lc4 ~~ 0*freq_neg_cont_lc4
   
   aci1 ~~ 0*aci1
   aci2 ~~ 0*aci2
   aci3 ~~ 0*aci3
   aci4 ~~ 0*aci4
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

## Model A: Autoregressive ----
a1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
  
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '

a2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
  
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model B: Forward ----

b1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
  
    cy2 ~ cx1 + cy1
    cy3 ~ cx2 + cy2
    cy4 ~ cx3 + cy3
           '
b2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
  
    cy2 ~ d*cx1 + c*cy1 
    cy3 ~ d*cx2 + c*cy2 
    cy4 ~ d*cx3 + c*cy3 
           '

## Model C: Backward ----

c1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cy1 
    cx3 ~ cx2 + cy2 
    cx4 ~ cx3 + cy3 
    
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '
c2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + g*cy2 
    cx2 ~ a*cx1 + g*cy1 
    cx4 ~ a*cx3 + g*cy3 
    
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model D: Bidirectional ---- 

d1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cy1 
    cx3 ~ cx2 + cy2 
    cx4 ~ cx3 + cy3 
  
    cy2 ~ cx1 + cy1 
    cy3 ~ cx2 + cy2 
    cy4 ~ cx3 + cy3 
           '


d2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + g*cy2 
    cx2 ~ a*cx1 + g*cy1 
    cx4 ~ a*cx3 + g*cy3 
    
    cy2 ~ d*cx1 + c*cy1 
    cy3 ~ d*cx2 + c*cy2 
    cy4 ~ d*cx3 + c*cy3 
           '
## Estimation RI-CLPM ----

# models
models <- c("a1","a2","b1","b2","c1","c2","d1","d2")

m2_ri <- list()
for (i in models){
  m2_ri[[i]] <- lavaan(model = c(between,within,get(i),covarianzas),
                       data = df_study1_wide, 
                       estimator = "MLR", 
                       missing = "FIML",
                       meanstructure = T, 
                       int.ov.free = T)
}


gofdt2_ri <- list()
for (i in names(m2_ri)){
  x <- fitMeasures(m2_ri[[i]])[c("chisq.scaled", "df.scaled",
                                 "pvalue.scaled", "cfi.scaled",
                                 "tli.scaled", "rmsea.scaled",
                                 "srmr_mplus", "aic",
                                 "bic", "bic2",
                                 "logl", "npar",
                                 "scaling.factor.h0")]
  gofdt2_ri[[i]] <- setNames(as.numeric(x),
                             c("X2","df",
                               "pvalue","CFI",
                               "TLI","RMSEA",
                               "SRMR","AIC",
                               "BIC","aBIC",
                               "LL","par",
                               "LLcorrectf"))}

gofdt2_ri <- data.table(m=names(gofdt2_ri),dplyr::bind_rows(gofdt2_ri))

gofdt2_ri <- gofdt2_ri %>% mutate(
  interpret_CFI = effectsize::interpret_cfi(CFI),
  interpret_RMSEA = effectsize::interpret_rmsea(RMSEA),
) %>% relocate(interpret_CFI, .after = CFI) %>%
  relocate(interpret_RMSEA, .after = RMSEA)

# comparacion modelos

# testear si efectos son iguales en el tiempo

comp9 <- gof.comp(data = gofdt2_ri, 
                  pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# testear direccion relationes
comp10 <- gof.comp(data = gofdt2_ri, 
                   pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                                c("b2","d2"), c("c2","d2"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))


summary(m2_ri[["d1"]], fit.measures = T, ci = T, standardized = T)

summary(m2_ri[["d2"]], fit.measures = T, ci = T, standardized = T)


# RI-CLPM d15 y d18 con controles -------------------------------------------

between <- ' 
# Componentes between
   RI_x =~ 1*freq_neg_cont_lc1 + 1*freq_neg_cont_lc2 + 1*freq_neg_cont_lc3 + 1*freq_neg_cont_lc4
   RI_y =~ 1*aci1 + 1*aci2 + 1*aci3 + 1*aci4
           '

within  <- '
# Componentes within
   cx1 =~ 1*freq_neg_cont_lc1
   cx2 =~ 1*freq_neg_cont_lc2
   cx3 =~ 1*freq_neg_cont_lc3
   cx4 =~ 1*freq_neg_cont_lc4
   
   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3
   cy4 =~ 1*aci4

# Constreñir las varianzas del error de medicion a 0
   freq_neg_cont_lc1 ~~ 0*freq_neg_cont_lc1
   freq_neg_cont_lc2 ~~ 0*freq_neg_cont_lc2
   freq_neg_cont_lc3 ~~ 0*freq_neg_cont_lc3
   freq_neg_cont_lc4 ~~ 0*freq_neg_cont_lc4
   
   aci1 ~~ 0*aci1
   aci2 ~~ 0*aci2
   aci3 ~~ 0*aci3
   aci4 ~~ 0*aci4
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

## Model A: Autoregressive ----
a1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
  
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '

a2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
  
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model B: Forward ----

b1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 
    cx4 ~ cx3 
  
    cy2 ~ cx1 + cy1 + sex + ess
    cy3 ~ cx2 + cy2 + sex + ess
    cy4 ~ cx3 + cy3 + sex + ess
           '
b2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
  
    cy2 ~ d*cx1 + c*cy1 + sexdep*sex + essdep*ess
    cy3 ~ d*cx2 + c*cy2 + sexdep*sex + essdep*ess
    cy4 ~ d*cx3 + c*cy3 + sexdep*sex + essdep*ess
           '

## Model C: Backward ----

c1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cy1 + sex + ess
    cx3 ~ cx2 + cy2 + sex + ess
    cx4 ~ cx3 + cy3 + sex + ess
    
    cy2 ~ cy1
    cy3 ~ cy2
    cy4 ~ cy3
           '
c2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + g*cy2 + sexindepx*sex + sexindepx*ess
    cx2 ~ a*cx1 + g*cy1 + sexindepx*sex + sexindepx*ess
    cx4 ~ a*cx3 + g*cy3 + sexindepx*sex + sexindepx*ess
    
    cy2 ~ c*cy1
    cy3 ~ c*cy2
    cy4 ~ c*cy3
           '

## Model D: Bidirectional ---- 

d1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cy1 + sex + ess
    cx3 ~ cx2 + cy2 + sex + ess
    cx4 ~ cx3 + cy3 + sex + ess
    
    cy2 ~ cx1 + cy1 + sex + ess
    cy3 ~ cx2 + cy2 + sex + ess
    cy4 ~ cx3 + cy3 + sex + ess
           '


d2  <- '
# Estimar los efectos constreñidos
    cx3 ~ a*cx2 + g*cy2 + sexindepx*sex + sexindepx*ess
    cx2 ~ a*cx1 + g*cy1 + sexindepx*sex + sexindepx*ess
    cx4 ~ a*cx3 + g*cy3 + sexindepx*sex + sexindepx*ess
    
    cy2 ~ d*cx1 + c*cy1 + sexdep*sex + essdep*ess
    cy3 ~ d*cx2 + c*cy2 + sexdep*sex + essdep*ess
    cy4 ~ d*cx3 + c*cy3 + sexdep*sex + essdep*ess
           '
## Estimation RI-CLPM ----

# models
models <- c("a1","a2","b1","b2","c1","c2","d1","d2")

m2_ri_con <- list()
for (i in models){
  m2_ri_con[[i]] <- lavaan(model = c(between,within,get(i),covarianzas),
                           data = df_study1_wide, 
                           estimator = "MLR", 
                           missing = "FIML",
                           meanstructure = T, 
                           int.ov.free = T)
}


gofdt2_ri_con <- list()
for (i in names(m2_ri_con)){
  x <- fitMeasures(m2_ri_con[[i]])[c("chisq.scaled", "df.scaled",
                                     "pvalue.scaled", "cfi.scaled",
                                     "tli.scaled", "rmsea.scaled",
                                     "srmr_mplus", "aic",
                                     "bic", "bic2",
                                     "logl", "npar",
                                     "scaling.factor.h0")]
  gofdt2_ri_con[[i]] <- setNames(as.numeric(x),
                                 c("X2","df",
                                   "pvalue","CFI",
                                   "TLI","RMSEA",
                                   "SRMR","AIC",
                                   "BIC","aBIC",
                                   "LL","par",
                                   "LLcorrectf"))}

gofdt2_ri_con <- data.table(m=names(gofdt2_ri_con),dplyr::bind_rows(gofdt2_ri_con))

gofdt2_ri_con <- gofdt2_ri_con %>% mutate(
  interpret_CFI = effectsize::interpret_cfi(CFI),
  interpret_RMSEA = effectsize::interpret_rmsea(RMSEA),
) %>% relocate(interpret_CFI, .after = CFI) %>%
  relocate(interpret_RMSEA, .after = RMSEA)

# comparacion modelos

# testear si efectos son iguales en el tiempo

comp11 <- gof.comp(data = gofdt2_ri_con, 
                  pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# testear direccion relationes
comp12 <- gof.comp(data = gofdt2_ri_con, 
                   pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                                c("b2","d2"), c("c2","d2"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))


summary(m2_ri_con[["d1"]], fit.measures = T, ci = T, standardized = T)

summary(m2_ri_con[["d2"]], fit.measures = T, ci = T, standardized = T)

# 4. Save and export ------------------------------------------------------

save(m, m_con, m_ri, m_ri_con, gofdt, gofdt_con, gofdt_ri, gofdt_ri_con,
     m2_ri, m2_ri_con, gofdt2_ri, gofdt2_ri_con, 
     comp1, comp2, comp3, comp4, comp5, comp6, comp7, comp8,
     comp9, comp10, comp11, comp12,
     file = here("output/models/ri_clpm_models.RData"),
     compress = TRUE)

param1a <- data.table(parameterEstimates(m_ri_con[["d1"]]))

param1b <- data.table(parameterEstimates(m_ri_con[["d2"]]))

param2a <- data.table(parameterEstimates(m2_ri_con[["d1"]]))

param2b <- data.table(parameterEstimates(m2_ri_con[["d2"]]))

a <- param1a %>% 
  filter(op == "~" & lhs %in% c("cy2", "cy3", "cy4")) %>% 
  mutate(pvalue=gtools::stars.pval(pvalue),
         ci = paste0("[", round(ci.lower, 3), "-", round(ci.upper, 3), "]")) %>% 
  select(-c(z, ci.lower, ci.upper))


b <- param1b %>% 
  filter(op == "~" & lhs %in% c("cy2")) %>% 
  mutate(pvalue=gtools::stars.pval(pvalue),
         ci = paste0("[", round(ci.lower, 3), "-", round(ci.upper, 3), "]")) %>% 
  select(-c(label, z, ci.lower, ci.upper))

c <- param2a %>% 
  filter(op == "~" & lhs %in% c("cy2", "cy3", "cy4")) %>% 
  mutate(pvalue=gtools::stars.pval(pvalue),
         ci = paste0("[", round(ci.lower, 3), "-", round(ci.upper, 3), "]")) %>% 
  select(-c(z, ci.lower, ci.upper))


d <- param2b %>% 
  filter(op == "~" & lhs %in% c("cy2")) %>% 
  mutate(pvalue=gtools::stars.pval(pvalue),
         ci = paste0("[", round(ci.lower, 3), "-", round(ci.upper, 3), "]")) %>% 
  select(-c(label, z, ci.lower, ci.upper))

library(writexl)

dataset_names <- list(a, b, c, d)

writexl::write_xlsx(dataset_names, path = here("output/comparacion_modelos.xlsx"))

