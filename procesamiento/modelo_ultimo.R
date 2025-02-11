# 0. Identification ---------------------------------------------------

# Title: Data analysis for research paper on Contact and Classism
# Institution: Centro de Estudios de Conflicto y Cohesión Social (COES)
# Responsible: Researcher

# Date: February 11, 2025

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


# CLPM d13, d15 y d18 sin controles -------------------------------------------

#CLPM requiere de al menos 2 olas para ser identificado. 
#RICLPM requiere de al menos 3 olas de datos, pero se recomiendan 4-5 olas.

df_study1_wide <- df_study1_wide %>% 
  select(idencuesta, freq_cont_lc2, freq_cont_lc3, freq_cont_lc4, 
         freq_neg_cont_lc2, freq_neg_cont_lc3, freq_neg_cont_lc4,
         aci2, aci3, aci4, sex, ess) %>% 
  rename(freq_cont_lc1 = freq_cont_lc2,
         freq_cont_lc2 = freq_cont_lc3,
         freq_cont_lc3 = freq_cont_lc4, 
         freq_neg_cont_lc1 = freq_neg_cont_lc2, 
         freq_neg_cont_lc2 = freq_neg_cont_lc3, 
         freq_neg_cont_lc3 = freq_neg_cont_lc4,
         aci1 = aci2, 
         aci2 = aci3,
         aci3 = aci4)

within  <- '
# Componentes within
   cx1 =~ 1*freq_cont_lc1
   cx2 =~ 1*freq_cont_lc2
   cx3 =~ 1*freq_cont_lc3

   cz1 =~ 1*freq_neg_cont_lc1
   cz2 =~ 1*freq_neg_cont_lc2
   cz3 =~ 1*freq_neg_cont_lc3

   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3

# Constreñir las varianzas del error de medicion a 0
   freq_cont_lc1 ~~ 0*freq_cont_lc1
   freq_cont_lc2 ~~ 0*freq_cont_lc2
   freq_cont_lc3 ~~ 0*freq_cont_lc3
   
   freq_neg_cont_lc1 ~~ 0*freq_neg_cont_lc1
   freq_neg_cont_lc2 ~~ 0*freq_neg_cont_lc2
   freq_neg_cont_lc3 ~~ 0*freq_neg_cont_lc3
 
   aci1 ~~ 0*aci1
   aci2 ~~ 0*aci2
   aci3 ~~ 0*aci3
           '

covarianzas  <- '
# Covarianza entre los componentes within t=1 con corr entre predictores
    cx1 ~~ cy1
    cz1 ~~ cy1
    cx1 ~~ cz1

# Covarianzas entre los residuos componente within con corr entre predictores
    cx2 ~~ cy2
    cx3 ~~ cy3
    
    cz2 ~~ cy2
    cz3 ~~ cy3

    cx2 ~~ cz2
    cx3 ~~ cz3

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
           '

## Model A: Autoregressive ----

a1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 

    cz2 ~ cz1 
    cz3 ~ cz2 

    cy2 ~ cy1
    cy3 ~ cy2
           '

a2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 

    cz2 ~ b*cz1 
    cz3 ~ b*cz2 

    cy2 ~ c*cy1
    cy3 ~ c*cy2
           '

## Model B: Forward ----

b1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 

    cz2 ~ cz1 
    cz3 ~ cz2 

    cy2 ~ cx1 + cz1 + cy1
    cy3 ~ cx2 + cz2 + cy2
           '
b2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 

    cz2 ~ b*cz1 
    cz3 ~ b*cz2 

    cy2 ~ d*cx1 + e*cz1 + c*cy1
    cy3 ~ d*cx2 + e*cz2 + c*cy2
           '

## Model C: Backward ----

c1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1
    cx3 ~ cx2 + cz2 + cy2

    cz2 ~ cz1 + cx1 + cy1
    cz3 ~ cz2 + cx2 + cy2

    cy2 ~ cy1
    cy3 ~ cy2
           '
c2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 + f*cz1 + g*cy1
    cx3 ~ a*cx2 + f*cz2 + g*cy2

    cz2 ~ b*cz1 + h*cx1 + i*cy1
    cz3 ~ b*cz2 + h*cx2 + i*cy2

    cy2 ~ c*cy1
    cy3 ~ c*cy2
           '

## Model D: Bidirectional ---- 

d1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1
    cx3 ~ cx2 + cz2 + cy2

    cz2 ~ cz1 + cx1 + cy1
    cz3 ~ cz2 + cx2 + cy2

    cy2 ~ cx1 + cz1 + cy1
    cy3 ~ cx2 + cz2 + cy2
           '

d2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 + f*cz1 + g*cy1
    cx3 ~ a*cx2 + f*cz2 + g*cy2

    cz2 ~ b*cz1 + h*cx1 + i*cy1
    cz3 ~ b*cz2 + h*cx2 + i*cy2

    cy2 ~ d*cx1 + e*cz1 + c*cy1
    cy3 ~ d*cx2 + e*cz2 + c*cy2
           '

## Estimation CLPM ----

# models
models <- c("a1","a2","b1","b2","c1","c2","d1","d2")

m4 <- list()
for (i in models){
  m4[[i]] <- lavaan(model = c(within,get(i),covarianzas),
                    data = df_study1_wide, 
                    estimator = "MLR", 
                    missing = "FIML",
                    meanstructure = T, 
                    int.ov.free = T)
}


gofdt4 <- list()
for (i in names(m4)){
  x <- fitMeasures(m4[[i]])[c("chisq.scaled", "df.scaled",
                              "pvalue.scaled", "cfi.scaled",
                              "tli.scaled", "rmsea.scaled",
                              "srmr_mplus", "aic",
                              "bic", "bic2",
                              "logl", "npar",
                              "scaling.factor.h0")]
  gofdt4[[i]] <- setNames(as.numeric(x),
                          c("X2","df",
                            "pvalue","CFI",
                            "TLI","RMSEA",
                            "SRMR","AIC",
                            "BIC","aBIC",
                            "LL","par",
                            "LLcorrectf"))}

gofdt4 <- data.table(m=names(gofdt4),dplyr::bind_rows(gofdt4))

gofdt4 <- gofdt4 %>% mutate(
  interpret_CFI = effectsize::interpret_cfi(CFI),
  interpret_RMSEA = effectsize::interpret_rmsea(RMSEA),
) %>% relocate(interpret_CFI, .after = CFI) %>%
  relocate(interpret_RMSEA, .after = RMSEA)

# comparacion modelos

# testear si efectos son iguales en el tiempo

comp17 <- gof.comp(data = gofdt4, 
                   pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# testear direccion relationes
comp18 <- gof.comp(data = gofdt4, 
                   pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                                c("b2","d2"), c("c2","d2"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))


summary(m4[["d1"]], fit.measures = T, ci = T, standardized = T)

summary(m4[["d2"]], fit.measures = T, ci = T, standardized = T)


# CLPM d13, d15 y d18 con controles -------------------------------------------

#CLPM requiere de al menos 2 olas para ser identificado. 
#RICLPM requiere de al menos 3 olas de datos, pero se recomiendan 4-5 olas.

within  <- '
# Componentes within
   cx1 =~ 1*freq_cont_lc1
   cx2 =~ 1*freq_cont_lc2
   cx3 =~ 1*freq_cont_lc3

   cz1 =~ 1*freq_neg_cont_lc1
   cz2 =~ 1*freq_neg_cont_lc2
   cz3 =~ 1*freq_neg_cont_lc3

   cy1 =~ 1*aci1
   cy2 =~ 1*aci2
   cy3 =~ 1*aci3

# Constreñir las varianzas del error de medicion a 0
   freq_cont_lc1 ~~ 0*freq_cont_lc1
   freq_cont_lc2 ~~ 0*freq_cont_lc2
   freq_cont_lc3 ~~ 0*freq_cont_lc3
   
   freq_neg_cont_lc1 ~~ 0*freq_neg_cont_lc1
   freq_neg_cont_lc2 ~~ 0*freq_neg_cont_lc2
   freq_neg_cont_lc3 ~~ 0*freq_neg_cont_lc3
 
   aci1 ~~ 0*aci1
   aci2 ~~ 0*aci2
   aci3 ~~ 0*aci3
           '

covarianzas  <- '
# Covarianza entre los componentes within t=1 con corr entre predictores
    cx1 ~~ cy1
    cz1 ~~ cy1
    cx1 ~~ cz1

# Covarianzas entre los residuos componente within con corr entre predictores
    cx2 ~~ cy2
    cx3 ~~ cy3
    
    cz2 ~~ cy2
    cz3 ~~ cy3

    cx2 ~~ cz2
    cx3 ~~ cz3

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
           '

## Model A: Autoregressive ----

a1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 

    cz2 ~ cz1 
    cz3 ~ cz2 

    cy2 ~ cy1
    cy3 ~ cy2
           '

a2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 

    cz2 ~ b*cz1 
    cz3 ~ b*cz2 

    cy2 ~ c*cy1
    cy3 ~ c*cy2
           '

## Model B: Forward ----

b1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 
    cx3 ~ cx2 

    cz2 ~ cz1 
    cz3 ~ cz2 

    cy2 ~ cx1 + cz1 + cy1 + sex + ess
    cy3 ~ cx2 + cz2 + cy2 + sex + ess
           '
b2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 
    cx3 ~ a*cx2 

    cz2 ~ b*cz1 
    cz3 ~ b*cz2 

    cy2 ~ d*cx1 + e*cz1 + c*cy1 + sexdep*sex + essdep*ess
    cy3 ~ d*cx2 + e*cz2 + c*cy2 + sexdep*sex + essdep*ess
           '

## Model C: Backward ----

c1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1 + sex + ess
    cx3 ~ cx2 + cz2 + cy2 + sex + ess

    cz2 ~ cz1 + cx1 + cy1 + sex + ess
    cz3 ~ cz2 + cx2 + cy2 + sex + ess

    cy2 ~ cy1
    cy3 ~ cy2
           '
c2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 + f*cz1 + g*cy1 + sexindepx*sex + sexindepx*ess
    cx3 ~ a*cx2 + f*cz2 + g*cy2 + sexindepx*sex + sexindepx*ess

    cz2 ~ b*cz1 + h*cx1 + i*cy1 + sexindepz*sex + sexindepz*ess
    cz3 ~ b*cz2 + h*cx2 + i*cy2 + sexindepz*sex + sexindepz*ess

    cy2 ~ c*cy1
    cy3 ~ c*cy2
           '

## Model D: Bidirectional ---- 

d1  <- '
# Estimar los efectos sin constreñir
    cx2 ~ cx1 + cz1 + cy1 + sex + ess
    cx3 ~ cx2 + cz2 + cy2 + sex + ess

    cz2 ~ cz1 + cx1 + cy1 + sex + ess
    cz3 ~ cz2 + cx2 + cy2 + sex + ess

    cy2 ~ cx1 + cz1 + cy1 + sex + ess
    cy3 ~ cx2 + cz2 + cy2 + sex + ess
           '

d2  <- '
# Estimar los efectos constreñidos
    cx2 ~ a*cx1 + f*cz1 + g*cy1 + sexindepx*sex + sexindepx*ess
    cx3 ~ a*cx2 + f*cz2 + g*cy2 + sexindepx*sex + sexindepx*ess

    cz2 ~ b*cz1 + h*cx1 + i*cy1 + sexindepz*sex + sexindepz*ess
    cz3 ~ b*cz2 + h*cx2 + i*cy2 + sexindepz*sex + sexindepz*ess

    cy2 ~ d*cx1 + e*cz1 + c*cy1 + sexdep*sex + essdep*ess
    cy3 ~ d*cx2 + e*cz2 + c*cy2 + sexdep*sex + essdep*ess
           '

## Estimation CLPM ----

# models
models <- c("a1","a2","b1","b2","c1","c2","d1","d2")

m4_con <- list()
for (i in models){
  m4_con[[i]] <- lavaan(model = c(within,get(i),covarianzas),
                       data = df_study1_wide, 
                       estimator = "MLR", 
                       missing = "FIML",
                       meanstructure = T, 
                       int.ov.free = T)
}

# fit measures
gofdt4_con <- list()
for (i in names(m4_con)){
  x <- fitMeasures(m4_con[[i]])[c("chisq.scaled", "df.scaled",
                                 "pvalue.scaled", "cfi.scaled",
                                 "tli.scaled", "rmsea.scaled",
                                 "srmr_mplus", "aic",
                                 "bic", "bic2",
                                 "logl", "npar",
                                 "scaling.factor.h0")]
  gofdt4_con[[i]] <- setNames(as.numeric(x),
                             c("X2","df",
                               "pvalue","CFI",
                               "TLI","RMSEA",
                               "SRMR","AIC",
                               "BIC","aBIC",
                               "LL","par",
                               "LLcorrectf"))}

gofdt4_con <- data.table(m=names(gofdt4_con),dplyr::bind_rows(gofdt4_con))


# comparacion modelos

# testear si efectos son iguales en el tiempo

comp19 <- gof.comp(data = gofdt4_con, 
                  pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))

# testear direccion relationes
comp20 <- gof.comp(data = gofdt4_con, 
                  pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                               c("b2","d2"), c("c2","d2"))) %>% 
  mutate(test_cfi_d = ifelse(CFI_D < 0.02, "meet invariance", "do not meet invariance"),
         test_rmsea_d = ifelse(RMSEA_D < 0.03, "meet invariance", "do not meet invariance"))


summary(m4_con[["d1"]], fit.measures = T, ci = T, standardized = T)

summary(m4_con[["d2"]], fit.measures = T, ci = T, standardized = T)

# Con o sin controles?

lavTestLRT(m4[["d2"]], m4_con[["d2"]]) # el modelo con conytoles no mejora significativamente el ajuste

performance::compare_performance(m4$d2, m4_con$d2)

data.table(parameterEstimates(m4_con[["d2"]])) %>% 
  filter(op == "~" & lhs %in% c("cy2")) %>% 
  mutate(pvalue=gtools::stars.pval(pvalue),
         ci = paste0("[", round(ci.lower, 3), "-", round(ci.upper, 3), "]")) %>% 
  select(-c(label, z, ci.lower, ci.upper))

