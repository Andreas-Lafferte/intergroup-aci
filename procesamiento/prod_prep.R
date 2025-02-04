#*******************************************************************************************************************
#
# 0. Identification ---------------------------------------------------
# Title: Data preparation for Research Paper
# Author: Andreas Laffert            
# Overview: Preparatrion of ELSOC 2018-2023      
# Date: 04-01-2025           
#
#******************************************************************************************************************

# 1. Packages  -----------------------------------------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse,
               here,
               magrittr,
               rio,
               sjmisc,
               sjlabelled,
               naniar,
               data.table)

options(scipen=999)
options(survey.lonely.psu = "certainty")


# 2. Data -----------------------------------------------------------------

load(url("https://dataverse.harvard.edu/api/access/datafile/10797987"))

glimpse(elsoc_long_2016_2023)

# 3. Processing --------------------------------------------------------------

elsoc_long_2016_2023 %>% 
  group_by(ola, muestra) %>% 
  count(d18) %>% 
  na.omit() %>% 
  print(n = nrow(.)) # var dep by wave and sample: d18 is only un refresh sample since 2018-2013


# 3.1 Select ----

db_proc <- elsoc_long_2016_2023 %>% 
  dplyr::select(idencuesta, ola, muestra, 
                ponderador_long_total,segmento, estrato,
                freq_cont_lc = d13, 
                positiv_contc_lc = d14, 
                freq_neg_cont_lc = d15, 
                aci = d18, 
                ess = d01_01, 
                sex = m0_sexo)

# 3.2 Filter ----
  
db_proc <- db_proc %>%  dplyr::filter(ola %in% c(3,4,6,7) & muestra == 2)
  
# 3.3 Recode and transform ----

frq(db_proc$freq_cont_lc)
frq(db_proc$positiv_contc_lc)
frq(db_proc$freq_neg_cont_lc)
frq(db_proc$aci)
frq(db_proc$ess)
frq(db_proc$sex)

## Identify NA's

db_proc <- db_proc %>% 
  mutate(
    across(.cols = c(freq_cont_lc, positiv_contc_lc, freq_neg_cont_lc, aci, ess, sex),
           .fns = ~ set_na(., na = c(-666,-777,-888,-999))
           
    )
  )

## Freq_cont_lc and freq_neg_cont_lc

labels1 <- c("Never" = 1, 
             "Hardly ever" = 2, 
             "Sometimes" = 3, 
             "Almost always" = 4, 
             "Always" = 5)

db_proc <- db_proc %>% 
  mutate(
    across(
      .cols = c(freq_cont_lc, freq_neg_cont_lc),
      .fns = ~ sjlabelled::set_labels(., labels = labels1)
      )
    )


## positiv_contc_lc

db_proc$positiv_contc_lc <- sjlabelled::set_labels(db_proc$positiv_contc_lc,
                                                   labels = c("Very unfriendly" = 1,
                                                              "Not very friendly" = 2,
                                                              "Neither friendly nor unfriendly" = 3,
                                                              "Fairly friendly" = 4,
                                                              "Very friendly" = 5))

## aci

db_proc$aci <- sjlabelled::set_labels(db_proc$aci,
                                      labels = c("Very negative" = 1,
                                                 "Negative" = 2,
                                                 "Neither negative nor positive" = 3,
                                                 "Positive" = 4,
                                                 "Very positive" = 5))

## sex

db_proc$sex <- car::recode(db_proc$sex, 
                           recodes = c("1='Male'; 2='Female'"), 
                           levels = c("Male", "Female"),
                           as.factor = T)


## wave 

db_proc <- db_proc %>% 
  rename(wave = ola) %>% 
  mutate(wave = case_when(wave == 3 ~ 1,
                          wave == 4 ~ 2,
                          wave == 6 ~ 3,
                          wave == 7 ~ 4,
                          TRUE ~ NA_real_))

# 3.4 Long to wide ----

df_study1_wide <- db_proc %>% 
  pivot_wider(id_cols = idencuesta,
              names_from = wave,
              values_from = c(freq_cont_lc, positiv_contc_lc, freq_neg_cont_lc, aci, ess, sex, ponderador_long_total,segmento, estrato),
              names_glue = "{.value}{wave}")


# 3.5 Missings values ----

colSums(is.na(df_study1_wide))

n_miss(df_study1_wide)

prop_miss(df_study1_wide)*100

miss_var_summary(df_study1_wide)

miss_var_table(df_study1_wide)

vis_miss(df_study1_wide) + theme(axis.text.x = element_text(angle=80))

# 3.6 Labels ----

vars1 <- c("freq_cont_lc1", "freq_cont_lc2", "freq_cont_lc3", "freq_cont_lc4")
labs1 <- c("Frequency of contact with lower class people W1",
          "Frequency of contact with lower class people W2",
          "Frequency of contact with lower class people W3",
          "Frequency of contact with lower class people W4")

for (i in seq_along(vars1)) {
  df_study1_wide[[vars1[i]]] <- set_label(df_study1_wide[[vars1[i]]], labs1[i])
}

vars2 <- c("positiv_contc_lc1", "positiv_contc_lc2", "positiv_contc_lc3", "positiv_contc_lc4")
labs2 <- c("Positive contact with lower class people W1",
           "Positive contact with lower class people W2",
           "Positive contact with lower class people W3",
           "Positive contact with lower class people W4")

for (i in seq_along(vars2)) {
  df_study1_wide[[vars2[i]]] <- set_label(df_study1_wide[[vars2[i]]], labs2[i])
}


vars3 <- c("freq_neg_cont_lc1", "freq_neg_cont_lc2", "freq_neg_cont_lc3", "freq_neg_cont_lc4")
labs3 <- c("Frequency of negative contact with lower class people W1",
           "Frequency of negative contact with lower class people W2",
           "Frequency of negative contact with lower class people W3",
           "Frequency of negative contact with lower class people W4")

for (i in seq_along(vars3)) {
  df_study1_wide[[vars3[i]]] <- set_label(df_study1_wide[[vars3[i]]], labs3[i])
}


vars4 <- c("aci1", "aci2", "aci3", "aci4")
labs4 <- c("Attitude towards lower class people W1",
           "Attitude towards lower class people W2",
           "Attitude towards lower class people W3",
           "Attitude towards lower class people W4")

for (i in seq_along(vars4)) {
  df_study1_wide[[vars4[i]]] <- set_label(df_study1_wide[[vars4[i]]], labs4[i])
}


vars5 <- c("ess1", "ess2", "ess3", "ess4")
labs5 <- c("Subjective Social Status W1",
           "Subjective Social Status W2",
           "Subjective Social Status W3",
           "Subjective Social Status W4")

for (i in seq_along(vars5)) {
  df_study1_wide[[vars5[i]]] <- set_label(df_study1_wide[[vars5[i]]], labs5[i])
}


vars6 <- c("sex1", "sex2", "sex3", "sex4")
labs6 <- c("Sex W1",
           "Sex W2",
           "Sex W3",
           "Sex W4")

for (i in seq_along(vars6)) {
  df_study1_wide[[vars6[i]]] <- set_label(df_study1_wide[[vars6[i]]], labs6[i])
}


# 4. Save and export  ----------------------------------------------------------------

df_study1_wide <- df_study1_wide %>% 
  select(idencuesta,
         starts_with(c("freq", "positiv", "aci")),
         sex = sex1,
         ess = ess1,
         starts_with(c("ponderador", "segmento", "estrato")))

save(df_study1_wide, file = here("input/data/proc/df_study1_wide.RData"))
save(db_proc, file = here("input/data/proc/df_study1_long.RData"))
