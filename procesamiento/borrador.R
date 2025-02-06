

# 3.2 RI-CLPM ----

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
c1 <- ' # Sin constreñir
           cx2 ~ cx1 + cy1
           cx3 ~ cx2 + cy2 
           cx4 ~ cx3 + cy3
           cx5 ~ cx4 + cy4
           cy2 ~ cy1
           cy3 ~ cy2
           cy4 ~ cy3
           cy5 ~ cy4
          '

c2 <- ' # Constreñido
           cx2 ~ a*cx1 + b*cy1
           cx3 ~ a*cx2 + b*cy2 
           cx4 ~ a*cx3 + b*cy3
           cx5 ~ a*cx4 + b*cy4
           cy2 ~ d*cy1
           cy3 ~ d*cy2
           cy4 ~ d*cy3
           cy5 ~ d*cy4
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



m1 <- lavaan(model = c(within, a1, covarianzas),
             data = df_study1_wide, 
             estimator = "MLR", missing = "FIML",
             meanstructure = T, int.ov.free = T)

summary(m1, fit.measures = T, ci = T, standardized = T)




# efectos
param <- data.table(parameterEstimates(m[["d2"]]))

eff_d2 <- param %>% 
  filter(op == "~" & rhs %in% c("cx1","cz1","cy1") & lhs %in% c("cx2","cz2","cy2")) %>% 
  mutate(pvalue=gtools::stars.pval(pvalue),
         ci = paste0("[", round(ci.lower, 3), "-", round(ci.upper, 3), "]")) %>% 
  select(-c(label, z, ci.lower, ci.upper))

## autoregresivos
eff_d2[c(1,4,8),] %>% 
  kableExtra::kable(., "markdown")


## cross lagged
eff_d2 %>% 
  kableExtra::kable(., "markdown")

