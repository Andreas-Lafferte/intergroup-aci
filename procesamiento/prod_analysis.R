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
               performance,
               data.table)

options(scipen=999)
options(survey.lonely.psu = "certainty")
rm(list = ls())

# 2. Data -----------------------------------------------------------------

load(file = here("input/data/proc/df_study1_wide.RData"))
load(file = here("input/data/proc/df_study1_long.RData"))

glimpse(df_study1_wide)

# 3. Analysis -------------------------------------------------------------

df_study1_wide <- df_study1_wide[,c(1:5,10:19)]

# 3.1 Descriptive ----

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
# Un 11% de la varianza de aci se asocia a diferencias entre 
# individuos, mientras que un 89% de su varianza corresponde al cambio en el tiempo en los mismos individuos
# Es posble usar un RI-CLPM ? ....


# 3.1 CLPM ----

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

## Model A: Independence hypothesis ----

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

## Model B: Spillover hypothesis ----

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

## Model C: Facilitation hypothesis ----

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



## Model D: Reciprocity hypothesis ---- 

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


# comparacion modelos
gof.comp  = function(data, pairs, measures = c("CFI", "TLI", "RMSEA", "SRMR", "AIC", 
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

# testear si efectos son iguales en el tiempo

comp1 <- gof.comp(data = gofdt, 
                  pairs = list(c("a2","a1"), c("b2","b1"), c("c2","c1"), c("d2","d1"))) %>% 
  mutate(decision = ifelse(CFI_D < 0.02 | RMSEA_D < 0.03, yes = "Accept", no = "Reject"))

# testear direccion relationes
comp2 <- gof.comp(data = gofdt, 
                  pairs = list(c("a2","b2"), c("a2","c2"), c("a2","d2"), 
                               c("b2","d2"), c("c2","d2"))) %>% 
  mutate(decision = ifelse(CFI_D < 0.02 | RMSEA_D < 0.03, yes = "Accept", no = "Reject"))

# conclusion, modelo d2 mejor modelo

summary(m[["d2"]], fit.measures = T, ci = T, standardized = T)
