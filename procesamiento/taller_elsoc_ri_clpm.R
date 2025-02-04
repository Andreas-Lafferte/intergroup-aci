library(data.table)
library(lavaan)

# Data
load(url("https://dataverse.harvard.edu/api/access/datafile/6160173"))
setDT(elsoc_long_2016_2021)
elsoc_long_2016_2021 <- elsoc_long_2016_2021[muestra==1,]
# Long to wide
elsoc <- dcast(elsoc_long_2016_2021, idencuesta ~ ola, 
               value.var = c("c08_02","c08_04","m0_edad","m01","m0_sexo"))
# Select variables
old <- c("idencuesta",paste0("c08_02_",1:5),paste0("c08_04_",1:5),
         "m0_edad_1","m01_1","m0_sexo_1")
new <- c("id",paste0("off",1:5),paste0("onl",1:5),"edad","educ","sexo")
setnames(elsoc,old,new)
elsoc <- elsoc[, ..new]
# Recode NAs
elsoc <- elsoc[, lapply(.SD, function(x) replace(x, 
                                                 which(x %in% c(-666,-888,-999)), NA))]
# Cleaning workspace
rm(elsoc_long_2016_2021,old,new)
elsoc <- as.data.frame(elsoc)

summarytools::descr(elsoc[-1], order="p", stats = "common", 
                    transpose = T, style = "rmarkdown", headings = F)

between <- ' 
# Crear los componentes between
   RI_x =~ 1*off1 + 1*off2 + 1*off3 + 1*off4 + 1*off5
   RI_y =~ 1*onl1 + 1*onl2 + 1*onl3 + 1*onl4 + 1*onl5
           '

within  <- '
# Crear los componentes within
   cx1 =~ 1*off1
   cx2 =~ 1*off2
   cx3 =~ 1*off3
   cx4 =~ 1*off4
   cx5 =~ 1*off5
   cy1 =~ 1*onl1
   cy2 =~ 1*onl2
   cy3 =~ 1*onl3
   cy4 =~ 1*onl4
   cy5 =~ 1*onl5
   
# constreñir las varianzas del error de medicion a 0
   off1 ~~ 0*off1
   off2 ~~ 0*off2
   off3 ~~ 0*off3
   off4 ~~ 0*off4
   off5 ~~ 0*off5
   onl1 ~~ 0*onl1
   onl2 ~~ 0*onl2
   onl3 ~~ 0*onl3
   onl4 ~~ 0*onl4
   onl5 ~~ 0*onl5
          '

efectos  <- '
# Estimar los efectos lagged
    cx2 ~ cx1 + cy1
    cx3 ~ cx2 + cy2
    cx4 ~ cx3 + cy3
    cx5 ~ cx4 + cy4
    cy2 ~ cx1 + cy1
    cy3 ~ cx2 + cy2
    cy4 ~ cx3 + cy3
    cy5 ~ cx4 + cy4
           '

covarianzas  <- '
# Covarianza entre los componentes within t=1
    cx1 ~~ cy1
   
# Covarianzas entre los residuos componente within
    cx2 ~~ cy2
    cx3 ~~ cy3
    cx4 ~~ cy4
    cx5 ~~ cy5
    
# Varianzas residuales del componente within
    cx1 ~~ cx1
    cy1 ~~ cy1 
    cx2 ~~ cx2
    cy2 ~~ cy2 
    cx3 ~~ cx3 
    cy3 ~~ cy3 
    cx4 ~~ cx4 
    cy4 ~~ cy4 
    cx5 ~~ cx5
    cy5 ~~ cy5    
    
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

fit_lavaan2 <- lavaan(model = c(between, within, efectos, covarianzas),
                      data = elsoc, 
                      estimator = "MLR", missing = "FIML",
                      meanstructure = T, int.ov.free = T)
summary(fit_lavaan2, fit.measures = T, ci = T, standardized = T)
