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

# 3.1 Descriptive -------------------------------------------------------------

df_study1_wide <- df_study1_wide_or[,c(1:9,14:19)]

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


# CLPM without controls -------------------------------------------

within  <- '
# Componentes within
   cx1 =~ 1*freq_cont_lc1
   cx2 =~ 1*freq_cont_lc2
   cx3 =~ 1*freq_cont_lc3
   cx4 =~ 1*freq_cont_lc4

   cz1 =~ 1*freq_neg_cont_lc1
   cz2 =~ 1*freq_neg_cont_lc2
   cz3 =~ 1*freq_neg_cont_lc3
   cz4 =~ 1*freq_neg_cont_lc4

   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3
   cy4 =~ 1*aci4

# Constreñir las varianzas del error de medicion a 0
   freq_cont_lc1 ~~ 0*freq_cont_lc1
   freq_cont_lc2 ~~ 0*freq_cont_lc2
   freq_cont_lc3 ~~ 0*freq_cont_lc3
   freq_cont_lc4 ~~ 0*freq_cont_lc4
   
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
  
    cy2 ~ cx1 + cz1 + cy1
    cy3 ~ cx2 + cz2 + cy2
    cy4 ~ cx3 + cz3 + cy3
           '
b2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
    
    cz2 ~ b*cz1 
    cz3 ~ b*cz2 
    cz4 ~ b*cz3
  
    cy2 ~ d*cx1 + e*cz1 + c*cy1
    cy3 ~ d*cx2 + e*cz2 + c*cy2
    cy4 ~ d*cx3 + e*cz3 + c*cy3
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
    cx2 ~ a*cx1 + f*cz1 + g*cy1
    cx3 ~ a*cx2 + f*cz2 + g*cy2
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
  
    cy2 ~ cx1 + cz1 + cy1
    cy3 ~ cx2 + cz2 + cy2
    cy4 ~ cx3 + cz3 + cy3
           '

d2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 + f*cz1 + g*cy1
    cx3 ~ a*cx2 + f*cz2 + g*cy2
    cx4 ~ a*cx3 + f*cz3 + g*cy3
    
    cz2 ~ b*cz1 + h*cx1 + i*cy1
    cz3 ~ b*cz2 + h*cx2 + i*cy2
    cz4 ~ b*cz3 + h*cx3 + i*cy3
  
    cy2 ~ d*cx1 + e*cz1 + c*cy1
    cy3 ~ d*cx2 + e*cz2 + c*cy2
    cy4 ~ d*cx3 + e*cz3 + c*cy3
           '

## Estimation CLPM ----

# models
models <- c("a1","a2","b1","b2","c1","c2","d1","d2")

m1 <- list()
for (i in models){
  m1[[i]] <- lavaan(model = c(within,get(i),covarianzas),
                    data = df_study1_wide, 
                    estimator = "MLR", 
                    missing = "FIML",
                    meanstructure = T, 
                    int.ov.free = T)
}


gofdt1 <- list()
for (i in names(m1)){
  x <- fitMeasures(m1[[i]])[c("chisq.scaled", "df.scaled",
                              "pvalue.scaled", "cfi.scaled",
                              "tli.scaled", "rmsea.scaled",
                              "srmr_mplus", "aic",
                              "bic", "bic2",
                              "logl", "npar",
                              "scaling.factor.h0")]
  gofdt1[[i]] <- setNames(as.numeric(x),
                          c("X2","df",
                            "pvalue","CFI",
                            "TLI","RMSEA",
                            "SRMR","AIC",
                            "BIC","aBIC",
                            "LL","par",
                            "LLcorrectf"))}

gofdt1 <- data.table(m=names(gofdt1),dplyr::bind_rows(gofdt1))

gofdt1 <- gofdt1 %>% mutate(
  interpret_CFI = effectsize::interpret_cfi(CFI),
  interpret_RMSEA = effectsize::interpret_rmsea(RMSEA),
) %>% relocate(interpret_CFI, .after = CFI) %>%
  relocate(interpret_RMSEA, .after = RMSEA)

# comparacion modelos

# testear si efectos son iguales en el tiempo

comp1 <- gof.comp(data = gofdt1, 
                   pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# testear direccion relationes
comp2 <- gof.comp(data = gofdt1, 
                   pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                                c("b2","d2"), c("c2","d2"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))


summary(m1[["d1"]], fit.measures = T, ci = T, standardized = T)

summary(m1[["d2"]], fit.measures = T, ci = T, standardized = T)


# CLPM with controls -------------------------------------------

within  <- '
# Componentes within
   cx1 =~ 1*freq_cont_lc1
   cx2 =~ 1*freq_cont_lc2
   cx3 =~ 1*freq_cont_lc3
   cx4 =~ 1*freq_cont_lc4

   cz1 =~ 1*freq_neg_cont_lc1
   cz2 =~ 1*freq_neg_cont_lc2
   cz3 =~ 1*freq_neg_cont_lc3
   cz4 =~ 1*freq_neg_cont_lc4

   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3
   cy4 =~ 1*aci4

# Constreñir las varianzas del error de medicion a 0
   freq_cont_lc1 ~~ 0*freq_cont_lc1
   freq_cont_lc2 ~~ 0*freq_cont_lc2
   freq_cont_lc3 ~~ 0*freq_cont_lc3
   freq_cont_lc4 ~~ 0*freq_cont_lc4
   
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
  
    cy2 ~ cx1 + cz1 + cy1 + sex + ess
    cy3 ~ cx2 + cz2 + cy2 + sex + ess
    cy4 ~ cx3 + cz3 + cy3 + sex + ess
           '
b2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 
    cx4 ~ a*cx3 
    
    cz2 ~ b*cz1 
    cz3 ~ b*cz2 
    cz4 ~ b*cz3
  
    cy2 ~ d*cx1 + e*cz1 + c*cy1 + sexdep*sex + essdep*ess
    cy3 ~ d*cx2 + e*cz2 + c*cy2 + sexdep*sex + essdep*ess
    cy4 ~ d*cx3 + e*cz3 + c*cy3 + sexdep*sex + essdep*ess
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
    cx2 ~ a*cx1 + f*cz1 + g*cy1 + sexindepx*sex + sexindepx*ess
    cx3 ~ a*cx2 + f*cz2 + g*cy2 + sexindepx*sex + sexindepx*ess
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

    cy2 ~ cx1 + cz1 + cy1 + sex + ess
    cy3 ~ cx2 + cz2 + cy2 + sex + ess
    cy4 ~ cx3 + cz3 + cy3 + sex + ess
           '

d2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 + f*cz1 + g*cy1 + sexindepx*sex + sexindepx*ess
    cx3 ~ a*cx2 + f*cz2 + g*cy2 + sexindepx*sex + sexindepx*ess
    cx4 ~ a*cx3 + f*cz3 + g*cy3 + sexindepx*sex + sexindepx*ess

    cz2 ~ b*cz1 + h*cx1 + i*cy1 + sexindepz*sex + sexindepz*ess
    cz3 ~ b*cz2 + h*cx2 + i*cy2 + sexindepz*sex + sexindepz*ess
    cz4 ~ b*cz3 + h*cx3 + i*cy3 + sexindepz*sex + sexindepz*ess

    cy2 ~ d*cx1 + e*cz1 + c*cy1 + sexdep*sex + essdep*ess
    cy3 ~ d*cx2 + e*cz2 + c*cy2 + sexdep*sex + essdep*ess
    cy4 ~ d*cx3 + e*cz3 + c*cy3 + sexdep*sex + essdep*ess
           '

## Estimation CLPM ----

# models
models <- c("a1","a2","b1","b2","c1","c2","d1","d2")

m1_con <- list()
for (i in models){
  m1_con[[i]] <- lavaan(model = c(within,get(i),covarianzas),
                       data = df_study1_wide, 
                       estimator = "MLR", 
                       missing = "FIML",
                       meanstructure = T, 
                       int.ov.free = T)
}

# fit measures
gofdt1_con <- list()
for (i in names(m1_con)){
  x <- fitMeasures(m1_con[[i]])[c("chisq.scaled", "df.scaled",
                                 "pvalue.scaled", "cfi.scaled",
                                 "tli.scaled", "rmsea.scaled",
                                 "srmr_mplus", "aic",
                                 "bic", "bic2",
                                 "logl", "npar",
                                 "scaling.factor.h0")]
  gofdt1_con[[i]] <- setNames(as.numeric(x),
                             c("X2","df",
                               "pvalue","CFI",
                               "TLI","RMSEA",
                               "SRMR","AIC",
                               "BIC","aBIC",
                               "LL","par",
                               "LLcorrectf"))}

gofdt1_con <- data.table(m=names(gofdt1_con),dplyr::bind_rows(gofdt1_con))

gofdt1_con <- gofdt1_con %>% mutate(
  interpret_CFI = effectsize::interpret_cfi(CFI),
  interpret_RMSEA = effectsize::interpret_rmsea(RMSEA),
) %>% relocate(interpret_CFI, .after = CFI) %>%
  relocate(interpret_RMSEA, .after = RMSEA)

# comparacion modelos

# testear si efectos son iguales en el tiempo

comp3 <- gof.comp(data = gofdt1_con, 
                  pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# testear direccion relationes
comp4 <- gof.comp(data = gofdt1_con, 
                  pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                               c("b2","d2"), c("c2","d2"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))


summary(m1_con[["d1"]], fit.measures = T, ci = T, standardized = T)

summary(m1_con[["d2"]], fit.measures = T, ci = T, standardized = T)

# constreñido mejor?

lavTestLRT(m1_con[["d1"]], m1_con[["d2"]]) 

# Con o sin controles?

lavTestLRT(m1[["d2"]], m1_con[["d2"]]) # el modelo con conytoles mejora significativamente el ajuste

performance::compare_performance(m1$d2, m1_con$d2)

data.table(parameterEstimates(m1_con[["d2"]])) %>% 
  filter(op == "~" & lhs %in% c("cy2")) %>% 
  mutate(pvalue=gtools::stars.pval(pvalue),
         ci = paste0("[", round(ci.lower, 3), "-", round(ci.upper, 3), "]")) %>% 
  select(-c(label, z, ci.lower, ci.upper))

data.table(parameterEstimates(m1_con[["d1"]])) %>% 
  filter(op == "~" & lhs %in% c("cy2", "cy3", "cy4")) %>% 
  mutate(pvalue=gtools::stars.pval(pvalue),
         ci = paste0("[", round(ci.lower, 3), "-", round(ci.upper, 3), "]")) %>% 
  select(-c(z, ci.lower, ci.upper))

# 4. Save and export ------------------------------------------------------

save(df_study1_wide, m1, m1_con, gofdt1, gofdt1_con,
     comp1, comp2, comp3, comp4,
     file = here("output/models/clpm_models.RData"),
     compress = TRUE)
