library(lubridate)
library(tidyr)
library(dplyr)
library(readxl)
library(ggplot2)
library(MuMIn)
library(DHARMa)


#load the demographic and clinical data
data <- read.csv("C:/Users/julio/OneDrive/Escritorio/maestria/tfm/data_gihsn.csv")
#cargo la informacion estandarizada de laboratorio
data_lab <- read.csv("C:/Users/julio/OneDrive/Escritorio/maestria/tfm/data_lab.csv")
#load the subgrouping and sequencing data
lab2_data <- read.csv("C:/Users/julio/OneDrive/Escritorio/maestria/tfm/lab2_data.csv")
#change the date format
data$Fecha.de.Recoleccion <- as.Date(data$Fecha.de.Registro.y.Recolección.de.Muestras, format="%Y-%m-%d") 
#join the demographic+clinical data with the standarized lab data
df_combined <- data %>%
  left_join(data_lab, by="Record.ID")
#I filter only to include the IDs starting with "C" or "Q", since those were actually enrolled
#Other IDs are failed enrollments or test IDs
df_combined <- df_combined %>%
  filter(grepl("^[CQ]", Record.ID))

#Select only the period to be analyzed
df_combined <- df_combined %>%
  filter(Fecha.de.Registro.y.Recolección.de.Muestras>= "2025-08-26"& Fecha.de.Registro.y.Recolección.de.Muestras<="2026-01-31")


#Create a new variable, where, based on the ID's first letter, collection site is assigned
df_combined <- df_combined %>%
  mutate(hospital= case_when(
    startsWith(Record.ID, "C")~"Hospital Nacional de Chimaltenango",
    startsWith(Record.ID, "Q")~"Hospital Nacional de Coatepeque"
  ))

#now I create a variable where I have the age in months
df_combined <- df_combined %>%
  mutate(edad_meses=case_when(
    Unidades.de.edad...1..Años...2..Meses...3..Días. == 3 ~ 0,
    Unidades.de.edad...1..Años...2..Meses...3..Días. == 2 ~ Edad..,
    TRUE~NA
  ))


######filter only for <5 years old
under_five <- df_combined %>%
  filter(edad_meses < 60)
under_five$date <- as.Date(under_five$Fecha.de.Registro.y.Recolección.de.Muestras, format="%Y-%m-%d")

table(under_five$Género)

#now the additional lab data
add_lab <- under_five %>%
  left_join(lab2_data, by="Record.ID")
####hospitalization length 
add_lab$Fecha.de.ingreso.al.hospital <- as.Date(add_lab$Fecha.de.ingreso.al.hospital, format="%Y-%m-%d")
add_lab$Fecha.de.alta...defunción. <- as.Date(add_lab$Fecha.de.alta...defunción., format="%Y-%m-%d")
add_lab <- add_lab %>%
  mutate(dias_hospitalizacion = as.numeric(difftime(Fecha.de.alta...defunción., Fecha.de.ingreso.al.hospital, units="days")))

#filtro a los viajeros en el tiempo

add_lab <- add_lab %>%
  filter(dias_hospitalizacion >= 0)

######asignar severidad segun las variables
####definicion 1: admision a UCI
#def2_ oxigeno sin ventilacion mecanica
#def3: uci/ventilacion mecanica/fallecido durante hosp
#def4: def3 + oxigeno 
add_lab <- add_lab %>%
  mutate(def1=if_else(Admisión.en.UCI=="Si", 1, 0))%>%
  mutate(def2=if_else(Oxigeno.suplementario.sin.ventilación.mecánica=="Si", 1, 0))%>%
  mutate(def3=if_else(Admisión.en.UCI=="Si"|Ventilación.mecánica=="Si"|Fallecido.durante.la.hospitalización.=="Si", 1, 0))%>%
  mutate(def4=if_else(Admisión.en.UCI=="Si"|Ventilación.mecánica=="Si"|Fallecido.durante.la.hospitalización.=="Si"|Oxigeno.suplementario.sin.ventilación.mecánica=="Si", 1, 0))

##asigno subgrupo de VSR
add_lab <- add_lab%>%
  mutate(subgrupo_vsr=case_when(
    Resultados.target.VSR.A=="Amplificó"~"A",
    Resultados.target.VSR.B=="Amplificó"~"B",
    TRUE~"No subagrupado"
  ))

#Only keep those with subgrouping info
con_info_sub <- add_lab %>%
  filter(Resultados.target.VSR.A=="Amplificó"|Resultados.target.VSR.B=="Amplificó")

con_info_ngs <- add_lab %>%
  filter(Envió.una.muestra.de.Virus.Sincitial.Respiratorio.a.la.base.de.datos.GISAID.EpiFlu..=="Si")


add_lab <- add_lab %>%
  mutate(
    otra_comorb = if_else(
      Enfermedad.cardiovascular == "Si" |
        EPOC == "Si" |
        Asma == "Si" |
        Enfermedad.neurológica.o.neuromuscular == "Si" |
        Obesidad == "Si" |
        Hemoglobinopatías == "Si",
      "Si", "No"
    )
  )

##########################################
add_lab$Resultado.de.influenza <- add_lab$El.paciente.tuvo.un.resultado.positivo.a.Influenza.
add_lab <- add_lab %>%
  mutate(
    Resultado.de.influenza = ifelse(
      Resultado.de.influenza == "Si",
      "Positivo",
      Resultado.de.influenza
    )
  ) 

add_lab <- add_lab %>%
  rowwise() %>%
  mutate(
    Virus_detectados = {
      
      vals <- c_across(starts_with("Resultado.de"))
      
      positivos <- names(
        dplyr::select(add_lab, starts_with("Resultado.de"))
      )[
        !is.na(vals) & trimws(tolower(vals)) == "positivo"
      ]
      
      positivos <- gsub("^Resultado\\.de\\.?", "", positivos)
      
      if(length(positivos) == 0){
        "Negativo"
      } else {
        paste(positivos, collapse = ", ")
      }
    }
  ) %>%
  ungroup()
    
###############################################################
    
add_lab <- add_lab %>%
  mutate(resultado_resp = case_when(
    Virus_detectados=="influenza"|Virus_detectados=="Adenovirus, influenza"|
      Virus_detectados=="Coronavirus.Humano, influenza"~"Influenza",
    Virus_detectados=="Coronavirus.Humano, Virus.sincitial.respiratorio, Adenovirus"|
      Virus_detectados=="Enterovirus.Rhinovirus, Metapneumovirus.humano, Virus.sincitial.respiratorio, Virus.parainfluenza"|
      Virus_detectados=="Enterovirus.Rhinovirus, Virus.sincitial.respiratorio"|Virus_detectados=="Metapneumovirus.humano, Virus.sincitial.respiratorio, Adenovirus"|
      Virus_detectados=="SARS.CoV.2, Virus.sincitial.respiratorio"|Virus_detectados=="Virus.sincitial.respiratorio, Adenovirus"|
      Virus_detectados=="Virus.sincitial.respiratorio, influenza"~"RSV coinfection",
    Virus_detectados=="Virus.sincitial.respiratorio"~"RSV",
    Virus_detectados=="Negativo"~"Negativo",
    TRUE~"ORV"
  ))
    
    
add_lab$resultado_resp <- as.factor(add_lab$resultado_resp)    
    
add_lab$resultado_resp<-factor(add_lab$resultado_resp, levels=c("Negativo", "Influenza", "ORV", "RSV", "RSV coinfection"))

levels(add_lab$resultado_resp)



####rsv positives
rsv_pos <- add_lab %>%
  filter(VSR=="Positivo")
###########rsv negatives
rsv_neg <- add_lab %>%
  filter(VSR=="Negativo")

add_lab_filter <- add_lab %>%
  filter(Resultado.de.Virus.sincitial.respiratorio!="No sabe")%>%
  filter(Congestión.nasal...Rinorrea!="No sabe")%>%
  filter(Sibilancias!="No sabe")%>%
  filter(Falta.de.aire.dificultad.para.respirar!="No sabe")%>%
  filter(Diarrea!="No sabe")

add_lab <- add_lab %>%
  filter(Resultado.de.Virus.sincitial.respiratorio!="No sabe")

##################optimizando
cont_vars <- c("edad_meses", "dias_hospitalizacion")

cat_vars <- c(
  "Género",
  "hospital",
  "Nació.prematuro....37.semanas.",
  "Desnutrición..solo.para...5.años.",
  "Fiebre...antecedentes.de.fiebre",
  "Náuseas.o.vómitos",
  "Oxigeno.suplementario.sin.ventilación.mecánica",
  "Admisión.en.UCI",
  "Ventilación.mecánica",
  "Fallecido.durante.la.hospitalización.",
  "otra_comorb"
)

cat_vars_filter <- c(
  "Congestión.nasal...Rinorrea",
  "Sibilancias",
  "Falta.de.aire.dificultad.para.respirar",
  "Diarrea"
)




library(dplyr)
library(purrr)

n <- length(rsv_pos$edad_meses)
mean_x <- mean(rsv_pos$edad_meses)
se <- sd(rsv_pos$edad_meses) / sqrt(n)

# 95% CI using t distribution
ci <- mean_x + qt(c(0.025, 0.975), df = n - 1) * se

ci


n <- length(rsv_neg$edad_meses)
mean_x <- mean(rsv_neg$edad_meses)
se <- sd(rsv_neg$edad_meses) / sqrt(n)

# 95% CI using t distribution
ci <- mean_x + qt(c(0.025, 0.975), df = n - 1) * se

ci

##########

n <- length(rsv_pos$dias_hospitalizacion)
mean_x <- mean(rsv_pos$dias_hospitalizacion)
se <- sd(rsv_pos$dias_hospitalizacion) / sqrt(n)

# 95% CI using t distribution
ci <- mean_x + qt(c(0.025, 0.975), df = n - 1) * se

ci


n <- length(rsv_neg$dias_hospitalizacion)
mean_x <- mean(rsv_neg$dias_hospitalizacion)
se <- sd(rsv_neg$dias_hospitalizacion) / sqrt(n)

# 95% CI using t distribution
ci <- mean_x + qt(c(0.025, 0.975), df = n - 1) * se

ci

##################################

rsv_a <- add_lab %>%
  filter(Resultados.target.VSR.A=="Amplificó")
rsv_b <- add_lab %>%
  filter(Resultados.target.VSR.B=="Amplificó")



n <- length(rsv_a$edad_meses)
mean_x <- mean(rsv_a$edad_meses)
se <- sd(rsv_a$edad_meses) / sqrt(n)

# 95% CI using t distribution
ci <- mean_x + qt(c(0.025, 0.975), df = n - 1) * se

ci


n <- length(rsv_b$edad_meses)
mean_x <- mean(rsv_b$edad_meses)
se <- sd(rsv_b$edad_meses) / sqrt(n)

# 95% CI using t distribution
ci <- mean_x + qt(c(0.025, 0.975), df = n - 1) * se

ci








n <- length(rsv_a$dias_hospitalizacion)
mean_x <- mean(rsv_a$dias_hospitalizacion)
se <- sd(rsv_a$dias_hospitalizacion) / sqrt(n)

# 95% CI using t distribution
ci <- mean_x + qt(c(0.025, 0.975), df = n - 1) * se

ci





n <- length(rsv_b$dias_hospitalizacion)
mean_x <- mean(rsv_b$dias_hospitalizacion)
se <- sd(rsv_b$dias_hospitalizacion) / sqrt(n)

# 95% CI using t distribution
ci <- mean_x + qt(c(0.025, 0.975), df = n - 1) * se

ci








bw <- round(2 * IQR(rsv_pos$edad_meses) / length(rsv_pos$edad_meses)^(1/3))

bw

breaks <- seq(
  floor(min(add_lab$edad_meses, na.rm = TRUE)),
  ceiling(max(add_lab$edad_meses, na.rm = TRUE)) + bw,
  by = bw
)

add_lab$bins <- cut(
  add_lab$edad_meses,
  breaks = breaks,
  include.lowest = TRUE,
  right = FALSE
)

rsv_summary <- add_lab %>%
  group_by(bins) %>%
  summarise(
    total = n(),
    positives = sum(Resultado.de.Virus.sincitial.respiratorio == "Positivo", na.rm = TRUE),
    proportion = positives / total
  )



a <- ggplot(rsv_summary, aes(x = bins)) +
  
  # counts
  geom_col(aes(y = total),
           fill = "grey80",
           color = "black") +
  
  # proportion line
  geom_line(aes(y = proportion * max(total),
                group = 1),
            linetype = "dashed") +
  
  geom_point(aes(y = proportion * max(total)),
             size = 3) +
  
  # primary y axis = counts
  scale_y_continuous(
    name = "Count",
    
    # secondary axis = proportion
    sec.axis = sec_axis(
      ~ . / max(rsv_summary$total),
      name = "RSV positivity proportion"
    )
  ) +
  
  xlab("Age in months") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.5)
  )

a
ggsave("graf_pos_edad.jpg",a, dpi=300, height = 4.5, width=8)



analyze_cont <- function(var){
  
  cat("\n====================\n")
  cat("Variable:", var, "\n")
  
  # Summary
  print(
    add_lab %>%
      group_by(Resultado.de.Virus.sincitial.respiratorio) %>%
      summarise(
        mean = mean(.data[[var]], na.rm = TRUE),
        sd   = sd(.data[[var]], na.rm = TRUE),
        median = median(.data[[var]], na.rm = TRUE),
        IQR = IQR(.data[[var]], na.rm = TRUE)
      )
  )
  
  # Shapiro
  print(by(add_lab[[var]],
           add_lab$Resultado.de.Virus.sincitial.respiratorio,
           shapiro.test))
  
  # Wilcoxon
  formula <- as.formula(paste(var, "~ Resultado.de.Virus.sincitial.respiratorio"))
  print(wilcox.test(formula, data = add_lab))
}

walk(cont_vars, analyze_cont)

analyze_cat <- function(var, data = add_lab){
  
  cat("\n====================\n")
  cat("Variable:", var, "\n")
  
  tab <- table(data[[var]], data$Resultado.de.Virus.sincitial.respiratorio)
  print(tab)
  
  test <- chisq.test(tab)
  
  if(any(test$expected < 5)){
    cat("→ Fisher test\n")
    print(fisher.test(tab))
  } else {
    cat("→ Chi-square test\n")
    print(test)
  }
}

walk(cat_vars, analyze_cat)

analyze_cat_filter <- function(var, data = add_lab_filter){
  
  cat("\n====================\n")
  cat("Variable:", var, "\n")
  
  tab <- table(data[[var]], data$Resultado.de.Virus.sincitial.respiratorio)
  print(tab)
  
  test <- chisq.test(tab)
  
  if(any(test$expected < 5)){
    cat("→ Fisher test\n")
    print(fisher.test(tab))
  } else {
    cat("→ Chi-square test\n")
    print(test)
  }
}

walk(cat_vars_filter, analyze_cat_filter)


subagrupados <- rsv_pos %>%
  filter(subgrupo_vsr!="No subagrupado")


cont_vars <- c("edad_meses", "dias_hospitalizacion")

cat_vars <- c(
  "Género",
  "hospital",
  "Nació.prematuro....37.semanas.",
  "Desnutrición..solo.para...5.años.",
  "Fiebre...antecedentes.de.fiebre",
  "Náuseas.o.vómitos",
  "Oxigeno.suplementario.sin.ventilación.mecánica",
  "Admisión.en.UCI",
  "Ventilación.mecánica",
  "Fallecido.durante.la.hospitalización.",
  "otra_comorb"
)

cat_vars_filter <- c(
  "Congestión.nasal...Rinorrea",
  "Sibilancias",
  "Falta.de.aire.dificultad.para.respirar",
  "Diarrea"
)



library(dplyr)
library(purrr)

analyze_cont <- function(var){
  
  cat("\n====================\n")
  cat("Variable:", var, "\n")
  
  # Summary
  print(
    subagrupados %>%
      group_by(subgrupo_vsr) %>%
      summarise(
        mean = mean(.data[[var]], na.rm = TRUE),
        sd   = sd(.data[[var]], na.rm = TRUE),
        median = median(.data[[var]], na.rm = TRUE),
        IQR = IQR(.data[[var]], na.rm = TRUE)
      )
  )
  
  # Shapiro
  print(by(subagrupados[[var]],
           subagrupados$subgrupo_vsr,
           shapiro.test))
  
  # Wilcoxon
  formula <- as.formula(paste(var, "~ subgrupo_vsr"))
  print(wilcox.test(formula, data = subagrupados))
}

walk(cont_vars, analyze_cont)


analyze_cat <- function(var, data = subagrupados){
  
  cat("\n====================\n")
  cat("Variable:", var, "\n")
  
  tab <- table(data[[var]], data$subgrupo_vsr)
  print(tab)
  
  test <- chisq.test(tab)
  
  if(any(test$expected < 5)){
    cat("→ Fisher test\n")
    print(fisher.test(tab))
  } else {
    cat("→ Chi-square test\n")
    print(test)
  }
}

walk(cat_vars, analyze_cat)




##############################################


a<-ggplot(add_lab, aes(Resultado.de.Virus.sincitial.respiratorio, edad_meses))+
  geom_boxplot(linewidth=1)+
  ylab("Age in months")+
  xlab("RSV test result")+
  scale_y_continuous(
    limits = c(0, 60),
    breaks = seq(0, 60, 10)
  )+
  theme_minimal(base_size=16)
a
ggsave("edad_resutlado_vsr.jpg", a, dpi=300, width = 9, height = 4.5)  

b<-ggplot(subagrupados, aes(subgrupo_vsr, dias_hospitalizacion))+
  geom_boxplot(linewidth=1)+
  ylab("Days of hospitalization")+
  xlab("RSV subgroup")+
  theme_minimal(base_size=16)+
scale_y_continuous(
  limits = c(0, 70),
  breaks = seq(0, 70, 10)
)
b
ggsave("dias_hosp_subgrupo.jpg", b, dpi=300, width = 9, height = 4.5) 



  df_plot <- add_lab %>%
  count(hospital, Resultado.de.Virus.sincitial.respiratorio) %>%
  group_by(hospital) %>%
  mutate(
    prop = n / sum(n),
    label = scales::percent(prop, accuracy = 0.1)
  )

b<- ggplot(df_plot,
           aes(x = hospital,
               y = n,
               fill = Resultado.de.Virus.sincitial.respiratorio)) +
  geom_col(color="black") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5))  +
  labs(y = "Percentage", x = "Hospital")+
  ylab("Number of enrolled")+
  xlab("Location")+ labs(fill = "RSV result")
b

ggsave("hosp_result.jpg", b, dpi=300, width = 9, height = 4.5)  



df_plot5 <- subagrupados %>%
  count(Fiebre...antecedentes.de.fiebre, subgrupo_vsr) %>%
  group_by(subgrupo_vsr) %>%
  mutate(
    prop = n / sum(n),
    label = scales::percent(prop, accuracy = 0.1)
  )

f<-ggplot(df_plot5,
          aes(x = subgrupo_vsr,
              y = prop,
              fill = Fiebre...antecedentes.de.fiebre)) +
  geom_col(color="black", linewidth=1) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5))+
  ylab("Proportion")+
  xlab("RSV subgroup")+ labs(fill = "Presence of fever")+
  scale_fill_manual(
    values = c("Si" = "#A6CEE3", "No" = "#D9D9D9"),  # optional colors
    labels = c("Si" = "Yes", "No" = "No")
  )+
  theme_minimal(base_size=16)
f
ggsave("fiebre_subgrupo.jpg", f, dpi=300, width = 9, height = 4.5) 

#rinorrea en pos y neg


df_plot6 <- add_lab %>%
  filter(Congestión.nasal...Rinorrea!="No sabe")%>%
  count(Congestión.nasal...Rinorrea, Resultado.de.Virus.sincitial.respiratorio) %>%
  group_by(Resultado.de.Virus.sincitial.respiratorio) %>%
  mutate(
    prop = n / sum(n),
    label = scales::percent(prop, accuracy = 0.1)
  )
g<-ggplot(df_plot6,
          aes(Resultado.de.Virus.sincitial.respiratorio,
              y = prop,
              fill = Congestión.nasal...Rinorrea)) +
  geom_col(color="black",linewidth = 1) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5))+
  ylab("Proportion")+
  xlab("RSV result")+ labs(fill = "Presence of \nrhinorrea")+
  scale_fill_manual(
    values = c("Si" = "#A6CEE3", "No" = "#D9D9D9"),  # optional colors
    labels = c("Si" = "Yes", "No" = "No")
  )+
  theme_minimal(base_size = 16)
g

ggsave("rinorrea_resultado_vsr.jpg", g, dpi=300, width = 9, height = 4.5) 

df_plot7 <- add_lab %>%
  filter(Sibilancias!="No sabe")%>%
  count(Sibilancias, Resultado.de.Virus.sincitial.respiratorio) %>%
  group_by(Resultado.de.Virus.sincitial.respiratorio) %>%
  mutate(
    prop = n / sum(n),
    label = scales::percent(prop, accuracy = 0.1)
  )



h<-ggplot(df_plot7,
          aes(Resultado.de.Virus.sincitial.respiratorio,
              y = prop,
              fill = Sibilancias)) +
  geom_col(color="black",linewidth = 1) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5))+
  ylab("Proportion")+
  xlab("RSV result")+ labs(fill = "Presence of \nwheezing")+
  scale_fill_manual(
    values = c("Si" = "#A6CEE3", "No" = "#D9D9D9"),  # optional colors
    labels = c("Si" = "Yes", "No" = "No")
  )+
  theme_minimal(base_size = 16)
h

ggsave("sibilancias_resultado_vsr.jpg", h, dpi=300, width = 9, height = 4.5) 





df_plot8 <- add_lab %>%
  filter(Falta.de.aire.dificultad.para.respirar!="No sabe")%>%
  count(Falta.de.aire.dificultad.para.respirar, Resultado.de.Virus.sincitial.respiratorio) %>%
  group_by(Resultado.de.Virus.sincitial.respiratorio) %>%
  mutate(
    prop = n / sum(n),
    label = scales::percent(prop, accuracy = 0.1)
  )


i<-ggplot(df_plot8,
          aes(Resultado.de.Virus.sincitial.respiratorio,
              y = prop,
              fill = Falta.de.aire.dificultad.para.respirar)) +
  geom_col(color="black",linewidth = 1) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5))+
  ylab("Proportion")+
  xlab("RSV result")+ labs(fill = "Presence of \nshortness of breath")+
  scale_fill_manual(
    values = c("Si" = "#A6CEE3", "No" = "#D9D9D9"),  # optional colors
    labels = c("Si" = "Yes", "No" = "No")
  )+
  theme_minimal(base_size = 16)
i
ggsave("dificultad_respirar_resultado_vsr.jpg", i, dpi=300, width = 9, height = 4.5) 

#modelaje previo a calcular variables nutricionales (varios no tienen info de peso/talla)

#creo rangos
add_lab <- add_lab %>%
  mutate(rango_meses=case_when(
    edad_meses < 3 ~ "0-2",
    edad_meses >= 3 & edad_meses <7~ "3-6",
    edad_meses >= 7 & edad_meses <13~ "7-12",
    edad_meses >= 12 & edad_meses <= 24 ~ "12-24",
    edad_meses > 24 ~ ">24",
    TRUE~NA
  ))
add_lab$rango_meses <- factor(add_lab$rango_meses, 
                                      levels = c("0-2", "3-6", "7-12", "12-24", ">24"))

add_lab$rango_meses <- factor(add_lab$rango_meses, 
                              levels = c(">24", "12-24", "7-12", "3-6", "0-2"))


#generacion de modelos candidatos para def 1
cand.models1 <- list()
cand.models1[["Completo"]] = glm(def1 ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial")
cand.models1[["solo sub"]] = glm(def1 ~ Resultado.de.Virus.sincitial.respiratorio, data = add_lab, na.action=na.fail, family = "binomial")
cand.models1[["sin genero"]] = glm(def1 ~ Resultado.de.Virus.sincitial.respiratorio +hospital + rango_meses+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial")
cand.models1[["solo edad"]] = glm(def1 ~rango_meses, data = add_lab, na.action=na.fail, family = "binomial")
cand.models1[["solo prem"]] = glm(def1 ~Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial")
cand.models1[["edad y prem"]] = glm(def1 ~Nació.prematuro....37.semanas. +rango_meses, data = add_lab, na.action=na.fail, family = "binomial")
cand.models1[["demo"]] = glm(def1 ~ Género+rango_meses+hospital, data = add_lab, na.action=na.fail, family = "binomial")
cand.models1[["nulo"]] = glm(def1 ~ 1, data = add_lab, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models1,rank="AICc")
dredge.out
subset(dredge.out,delta<4)

avg_model <- model.avg(dredge.out, subset = delta < 4)

summary(avg_model)

confint(avg_model)

# Extract coefficients
coef_est <- coef(avg_model)

# Extract CI
ci <- confint(avg_model)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)

results

#validation
simulationOutput <- simulateResiduals(fittedModel = glm(def1 ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial"), plot = T)

#################definicion2
cand.models2 <- list()
cand.models2[["Completo"]] = glm(def2 ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial")
cand.models2[["solo sub"]] = glm(def2 ~ Resultado.de.Virus.sincitial.respiratorio, data = add_lab, na.action=na.fail, family = "binomial")
cand.models2[["sin genero"]] = glm(def2 ~ Resultado.de.Virus.sincitial.respiratorio +hospital + rango_meses+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial")
cand.models2[["solo edad"]] = glm(def2 ~rango_meses, data = add_lab, na.action=na.fail, family = "binomial")
cand.models2[["solo prem"]] = glm(def2 ~Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial")
cand.models2[["edad y prem"]] = glm(def2 ~Nació.prematuro....37.semanas. +rango_meses, data = add_lab, na.action=na.fail, family = "binomial")
cand.models2[["demo"]] = glm(def2 ~ Género+rango_meses+hospital, data = add_lab, na.action=na.fail, family = "binomial")
cand.models2[["nulo"]] = glm(def2 ~ 1, data = add_lab, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models2,rank="AICc")
dredge.out
subset(dredge.out,delta<2)

avg_model <- model.avg(dredge.out, subset = delta < 2)

summary(avg_model)

confint(avg_model)

# Extract coefficients
coef_est <- coef(avg_model)

# Extract CI
ci <- confint(avg_model)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)


results

#validation
simulationOutput2 <- simulateResiduals(fittedModel = glm(def2 ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial"), plot = T)



#################definicion3
cand.models3 <- list()
cand.models3[["Completo"]] = glm(def4 ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial")
cand.models3[["solo sub"]] = glm(def4 ~ Resultado.de.Virus.sincitial.respiratorio, data = add_lab, na.action=na.fail, family = "binomial")
cand.models3[["sin genero"]] = glm(def4 ~ Resultado.de.Virus.sincitial.respiratorio +hospital + rango_meses+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial")
cand.models3[["solo edad"]] = glm(def4 ~rango_meses, data = add_lab, na.action=na.fail, family = "binomial")
cand.models3[["solo prem"]] = glm(def4 ~Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial")
cand.models3[["edad y prem"]] = glm(def4 ~Nació.prematuro....37.semanas. +rango_meses, data = add_lab, na.action=na.fail, family = "binomial")
cand.models3[["demo"]] = glm(def4 ~ Género+rango_meses+hospital, data = add_lab, na.action=na.fail, family = "binomial")
cand.models3[["nulo"]] = glm(def4 ~ 1, data = add_lab, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models3,rank="AICc")
dredge.out
subset(dredge.out,delta<4)



avg_model <- model.avg(dredge.out, subset = delta < 4)

summary(avg_model)

confint(avg_model)

# Extract coefficients
coef_est <- coef(avg_model)

# Extract CI
ci <- confint(avg_model)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)


results

#validation
simulationOutput3 <- simulateResiduals(fittedModel = glm(def4 ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "binomial"), plot = T)


ggplot(add_lab, aes(dias_hospitalizacion))+
  geom_bar()


##############################
cand.models4 <- list()
cand.models4[["Completo"]] = glm(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "poisson")
cand.models4[["solo sub"]] = glm(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio, data = add_lab, na.action=na.fail, family = "poisson")
cand.models4[["sin genero"]] = glm(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio +hospital + rango_meses+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "poisson")
cand.models4[["solo edad"]] = glm(dias_hospitalizacion ~rango_meses, data = add_lab, na.action=na.fail, family = "poisson")
cand.models4[["solo prem"]] = glm(dias_hospitalizacion ~Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "poisson")
cand.models4[["edad y prem"]] = glm(dias_hospitalizacion ~Nació.prematuro....37.semanas. +rango_meses, data = add_lab, na.action=na.fail, family = "poisson")
cand.models4[["demo"]] = glm(dias_hospitalizacion ~ Género+rango_meses+hospital, data = add_lab, na.action=na.fail, family = "poisson")
cand.models4[["nulo"]] = glm(dias_hospitalizacion ~ 1, data = add_lab, na.action=na.fail, family = "poisson")



dredge.out<-model.sel(cand.models4,rank="AICc")
dredge.out
subset(dredge.out,delta<4)



avg_model <- model.avg(dredge.out, subset = delta < 4)

summary(avg_model)

confint(avg_model)

# Extract coefficients
coef_est <- coef(avg_model)

# Extract CI
ci <- confint(avg_model)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)


results

simulationOutput4 <- simulateResiduals(fittedModel = glm(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "poisson"), plot = T)

Modelo1 <-glm(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail, family = "poisson")
#test sobredispersi?n
E1 <- resid(Modelo1, type = "pearson")
sum(E1^2)/(Modelo1$df.res)
#If φ ≈ 1 → good (no overdispersion)
#If φ > 1 → overdispersion
#If φ < 1 → underdispersion (less common)

#bondad de ajuste de test sobredispersi?n
residual.variance <- sum(resid(Modelo1, type="pearson")^2)
binomial.variance <- Modelo1$df.residual
1-pchisq(residual.variance, binomial.variance)

#Small p-value (e.g. < 0.05) → evidence of overdispersion





cand.models5 <- list()
cand.models5[["Completo"]] = glm.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail)
cand.models5[["solo sub"]] = glm.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio, data = add_lab, na.action=na.fail)
cand.models5[["sin genero"]] = glm.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio +hospital + rango_meses+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail)
cand.models5[["solo edad"]] = glm.nb(dias_hospitalizacion ~rango_meses, data = add_lab, na.action=na.fail)
cand.models5[["solo prem"]] = glm.nb(dias_hospitalizacion ~Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail)
cand.models5[["edad y prem"]] = glm.nb(dias_hospitalizacion ~Nació.prematuro....37.semanas. +rango_meses, data = add_lab, na.action=na.fail)
cand.models5[["demo"]] = glm.nb(dias_hospitalizacion ~ Género+rango_meses+hospital, data = add_lab, na.action=na.fail)
cand.models5[["nulo"]] = glm.nb(dias_hospitalizacion ~ 1, data = add_lab, na.action=na.fail)



dredge.out<-model.sel(cand.models5,rank="AICc")
dredge.out
subset(dredge.out,delta<2)


avg_model <- model.avg(dredge.out, subset = delta < 2)

summary(avg_model)

confint(avg_model)


# Extract coefficients
coef_est <- coef(avg_model)

# Extract CI
ci <- confint(avg_model)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)


results

simulationOutput4 <- simulateResiduals(fittedModel = glm.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail), plot = T)

Modelo2 <- glm.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = add_lab, na.action=na.fail)
#test sobredispersi?n
E1 <- resid(Modelo2, type = "pearson")
sum(E1^2)/(Modelo2$df.res)
#If φ ≈ 1 → good (no overdispersion)
#If φ > 1 → overdispersion
#If φ < 1 → underdispersion (less common)

#bondad de ajuste de test sobredispersi?n
residual.variance <- sum(resid(Modelo2, type="pearson")^2)
binomial.variance <- Modelo2$df.residual
1-pchisq(residual.variance, binomial.variance)






library(glmmTMB)

mod_quasi<-glm(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., family = quasipoisson , data = add_lab, na.action=na.fail)

simulationOutput4 <- simulateResiduals(fittedModel = mod_quasi, plot = T)

#test sobredispersi?n
E1 <- resid(mod_quasi, type = "pearson")
sum(E1^2)/(mod_quasi$df.res)
#If φ ≈ 1 → good (no overdispersion)
#If φ > 1 → overdispersion
#If φ < 1 → underdispersion (less common)

#bondad de ajuste de test sobredispersi?n
residual.variance <- sum(resid(mod_quasi, type="pearson")^2)
binomial.variance <- mod_quasi$df.residual
1-pchisq(residual.variance, binomial.variance)

library(DHARMa)

sim_res <- simulateResiduals(mod_quasi)
plot(sim_res)

testDispersion(sim_res)

add_lab <- add_lab %>%
  mutate(
    semana = isoweek(date),
    year = isoyear(date)
  )


cand.models6 <- list()
cand.models6[["Completo"]] = glmer(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas.+(1|semana),family="poisson", data = add_lab, na.action=na.fail)
cand.models6[["solo sub"]] = glmer(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio+(1|semana),family="poisson", data = add_lab, na.action=na.fail)
cand.models6[["sin genero"]] = glmer(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio +hospital + rango_meses+Nació.prematuro....37.semanas.+(1|semana),family="poisson", data = add_lab, na.action=na.fail)
cand.models6[["solo edad"]] = glmer(dias_hospitalizacion ~rango_meses+(1|semana), data = add_lab,family="poisson", na.action=na.fail)
cand.models6[["solo prem"]] = glmer(dias_hospitalizacion ~Nació.prematuro....37.semanas.+(1|semana),family="poisson", data = add_lab, na.action=na.fail)
cand.models6[["edad y prem"]] = glmer(dias_hospitalizacion ~Nació.prematuro....37.semanas. +rango_meses+(1|semana),family="poisson", data = add_lab, na.action=na.fail)
cand.models6[["demo"]] = glmer(dias_hospitalizacion ~ Género+rango_meses+hospital+(1|semana),family="poisson", data = add_lab, na.action=na.fail)
cand.models6[["nulo"]] = glmer(dias_hospitalizacion ~ 1+(1|semana),family="poisson", data = add_lab, na.action=na.fail)


dredge.out<-model.sel(cand.models6,rank="AICc")
dredge.out
subset(dredge.out,delta<4)

mod <- glmer(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas.+(1|semana),family="poisson", data = add_lab, na.action=na.fail)
summary(mod)


# Extract coefficients
coef_est <- fixef(mod)

# Extract CI
ci <- confint(mod)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)


results

simulationOutput4 <- simulateResiduals(fittedModel =  mod, plot = T)





cand.models7 <- list()
cand.models7[["Completo"]] = glmer.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas.+(1|semana), data = add_lab, na.action=na.fail)
cand.models7[["solo sub"]] = glmer.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio+(1|semana), data = add_lab, na.action=na.fail)
cand.models7[["sin genero"]] = glmer.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio +hospital + rango_meses+Nació.prematuro....37.semanas.+(1|semana), data = add_lab, na.action=na.fail)
cand.models7[["solo edad"]] = glmer.nb(dias_hospitalizacion ~rango_meses+(1|semana), data = add_lab, na.action=na.fail)
cand.models7[["solo prem"]] = glmer.nb(dias_hospitalizacion ~Nació.prematuro....37.semanas.+(1|semana), data = add_lab, na.action=na.fail)
cand.models7[["edad y prem"]] = glmer.nb(dias_hospitalizacion ~Nació.prematuro....37.semanas. +rango_meses+(1|semana), data = add_lab, na.action=na.fail)
cand.models7[["demo"]] = glmer.nb(dias_hospitalizacion ~ Género+rango_meses+hospital+(1|semana), data = add_lab, na.action=na.fail)
cand.models7[["nulo"]] = glmer.nb(dias_hospitalizacion ~ 1+(1|semana), data = add_lab, na.action=na.fail)


dredge.out<-model.sel(cand.models7,rank="AICc")
dredge.out
subset(dredge.out,delta<4)

mod <- glmer.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas.+(1|semana), data = add_lab, na.action=na.fail)
summary(mod)


# Extract coefficients
coef_est <- fixef(mod)

# Extract CI
ci <- confint(mod)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)


results

simulationOutput4 <- simulateResiduals(fittedModel =  mod, plot = T)


boxplot(add_lab$dias_hospitalizacion)

mod <- glmmTMB(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas.+(1|semana), data = add_lab, na.action=na.fail, family = nbinom2, ziformula = ~1)

simulationOutput4 <- simulateResiduals(fittedModel =  mod, plot = T)


menos_de_mes <- add_lab %>%
  filter(dias_hospitalizacion<30)


cand.models8 <- list()
cand.models8[["Completo"]] = glm.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = menos_de_mes, na.action=na.fail)
cand.models8[["solo sub"]] = glm.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio, data = menos_de_mes, na.action=na.fail)
cand.models8[["sin genero"]] = glm.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio +hospital + rango_meses+Nació.prematuro....37.semanas., data = menos_de_mes, na.action=na.fail)
cand.models8[["solo edad"]] = glm.nb(dias_hospitalizacion ~rango_meses, data = menos_de_mes, na.action=na.fail)
cand.models8[["solo prem"]] = glm.nb(dias_hospitalizacion ~Nació.prematuro....37.semanas., data = menos_de_mes, na.action=na.fail)
cand.models8[["edad y prem"]] = glm.nb(dias_hospitalizacion ~Nació.prematuro....37.semanas. +rango_meses, data = menos_de_mes, na.action=na.fail)
cand.models8[["demo"]] = glm.nb(dias_hospitalizacion ~ Género+rango_meses+hospital, data = menos_de_mes, na.action=na.fail)
cand.models8[["nulo"]] = glm.nb(dias_hospitalizacion ~ 1, data = menos_de_mes, na.action=na.fail)


dredge.out<-model.sel(cand.models8,rank="AICc")
dredge.out
subset(dredge.out,delta<4)



avg_model <- model.avg(dredge.out, subset = delta < 4)

summary(avg_model)

confint(avg_model)

mod<-glm.nb(dias_hospitalizacion ~ Resultado.de.Virus.sincitial.respiratorio + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = menos_de_mes, na.action=na.fail)
simulationOutput4 <- simulateResiduals(fittedModel =  mod, plot = T)






############
prop_df <- add_lab %>%
  group_by(rango_meses) %>%
  summarise(
    n_total = n(),
    n_uci_pos = sum(def1 == 1),
    prop_uci = n_uci_pos / n_total
  )

a<-ggplot(prop_df, aes(rango_meses, prop_uci, group=1))+
  geom_point(size=3)+
  geom_line(size=1, linetype = "dashed")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ylab("Porcentaje de pacientes ingresados a la UCI")+
  xlab("Rango de edad (meses)")+
  theme_minimal(base_size = 14)
a
ggsave("rango_edad_uci.jpg", a, dpi=300, width = 8, height = 4.5)



prop_oxy <- add_lab %>%
  group_by(rango_meses, hospital) %>%
  summarise(
    n_total = n(),
    n_oxy_pos = sum(def2 == 1),
    prop_oxy = n_oxy_pos / n_total
  )

b<-ggplot(prop_oxy, aes(rango_meses, prop_oxy, group=1))+
  geom_point(size=3)+
  geom_line(size=1, linetype = "dashed")+
  facet_wrap(~hospital)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),limits = c(0, 1))+
  ylab("Porcentaje de pacientes que recibieron oxígeno")+
  xlab("Rango de edad (meses)")+
  theme_bw(base_size = 14)
b

ggsave("rango_edad_oxy.jpg", b, dpi=300, width = 8, height = 4.5)



prop_def4 <- add_lab %>%
  group_by(rango_meses, hospital) %>%
  summarise(
    n_total = n(),
    n_def4_pos = sum(def4 == 1),
    prop_def4 = n_def4_pos / n_total
  )



d<-ggplot(prop_def4, aes(rango_meses, prop_def4, group=1))+
  geom_point(size=3)+
  geom_line(size=1, linetype = "dashed")+
  facet_wrap(~hospital)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),limits = c(0, 1))+
  ylab("Porcentaje de pacientes que cumplen con \nla definición de severidad def4")+
  xlab("Rango de edad (meses)")+
  theme_bw(base_size = 14)
d

ggsave("rango_edad_def4.jpg", d, dpi=300, width = 8, height = 4.5)
#########################################################################
subagrupados <- add_lab %>%
  filter(subgrupo_vsr!="No subagrupado")
cand.models11 <- list()
cand.models11[["Completo"]] = glm(def1 ~ subgrupo_vsr + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial")
cand.models11[["solo sub"]] = glm(def1 ~ subgrupo_vsr, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models11[["sin genero"]] = glm(def1 ~ subgrupo_vsr +hospital + rango_meses+Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial")
cand.models11[["solo edad"]] = glm(def1 ~rango_meses, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models11[["solo prem"]] = glm(def1 ~Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial")
cand.models11[["edad y prem"]] = glm(def1 ~Nació.prematuro....37.semanas. +rango_meses, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models11[["demo"]] = glm(def1 ~ Género+rango_meses+hospital, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models11[["nulo"]] = glm(def1 ~ 1, data = subagrupados, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models11,rank="AICc")
dredge.out
subset(dredge.out,delta<4)

avg_model <- model.avg(dredge.out, subset = delta < 4)

summary(avg_model)

confint(avg_model)

# Extract coefficients
coef_est <- coef(avg_model)

# Extract CI
ci <- confint(avg_model)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)

results

#validation
simulationOutput <- simulateResiduals(fittedModel = glm(def1 ~ subgrupo_vsr + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial"), plot = T)




prop_df <- subagrupados %>%
  group_by(rango_meses) %>%
  summarise(
    n_total = n(),
    n_uci_pos = sum(def1 == 1),
    prop_uci = n_uci_pos / n_total
  )

a<-ggplot(prop_df, aes(rango_meses, prop_uci, group=1))+
  geom_point(size=3)+
  geom_line(size=1, linetype = "dashed")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ylab("Porcentaje de pacientes ingresados a la UCI")+
  xlab("Rango de edad (meses)")+
  theme_minimal(base_size = 14)
a
ggsave("rango_edad_uci_subgrupo.jpg", a, dpi=300, width = 8, height = 4.5)


cand.models12 <- list()
cand.models12[["Completo"]] = glm(def2 ~ subgrupo_vsr + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial")
cand.models12[["solo sub"]] = glm(def2 ~ subgrupo_vsr, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models12[["sin genero"]] = glm(def2 ~ subgrupo_vsr +hospital + rango_meses+Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial")
cand.models12[["solo edad"]] = glm(def2 ~rango_meses, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models12[["solo prem"]] = glm(def2 ~Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial")
cand.models12[["edad y prem"]] = glm(def2 ~Nació.prematuro....37.semanas. +rango_meses, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models12[["demo"]] = glm(def2 ~ Género+rango_meses+hospital, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models12[["nulo"]] = glm(def2 ~ 1, data = subagrupados, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models12,rank="AICc")
dredge.out
subset(dredge.out,delta<2)

avg_model <- model.avg(dredge.out, subset = delta < 2)

summary(avg_model)

confint(avg_model)

# Extract coefficients
coef_est <- coef(avg_model)

# Extract CI
ci <- confint(avg_model)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)

results

#validation
simulationOutput <- simulateResiduals(fittedModel = glm(def1 ~ subgrupo_vsr + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial"), plot = T)




prop_df <- subagrupados %>%
  group_by(rango_meses, Género, hospital) %>%
  summarise(
    n_total = n(),
    n_oxy_pos = sum(def2 == 1),
    prop_oxy = n_oxy_pos / n_total
  )





a<-ggplot(prop_df, aes(rango_meses, prop_oxy, colour=Género, group=Género))+
  geom_point(size=3)+
  facet_wrap(~hospital)+
  geom_line(size=1, linetype = "solid")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ylab("Porcentaje de pacientes que reciben oxígeno suplementario")+
  xlab("Rango de edad (meses)")+
  theme_bw(base_size = 14)
a
ggsave("rango_edad_oxy_subgrupo_genero.jpg", a, dpi=300, width = 8, height = 4.5)








cand.models13 <- list()
cand.models13[["Completo"]] = glm(def4 ~ subgrupo_vsr + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial")
cand.models13[["solo sub"]] = glm(def4 ~ subgrupo_vsr, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models13[["sin genero"]] = glm(def4 ~ subgrupo_vsr +hospital + rango_meses+Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial")
cand.models13[["solo edad"]] = glm(def4 ~rango_meses, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models13[["solo prem"]] = glm(def4 ~Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial")
cand.models13[["edad y prem"]] = glm(def4 ~Nació.prematuro....37.semanas. +rango_meses, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models13[["demo"]] = glm(def4 ~ Género+rango_meses+hospital, data = subagrupados, na.action=na.fail, family = "binomial")
cand.models13[["nulo"]] = glm(def4 ~ 1, data = subagrupados, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models13,rank="AICc")
dredge.out
subset(dredge.out,delta<2)

avg_model <- model.avg(dredge.out, subset = delta < 2)

summary(avg_model)

confint(avg_model)

# Extract coefficients
coef_est <- coef(avg_model)

# Extract CI
ci <- confint(avg_model)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)

results

#validation
simulationOutput <- simulateResiduals(fittedModel = glm(def1 ~ subgrupo_vsr + Género + rango_meses+hospital+Nació.prematuro....37.semanas., data = subagrupados, na.action=na.fail, family = "binomial"), plot = T)




prop_df <- subagrupados %>%
  group_by(rango_meses, Género, hospital) %>%
  summarise(
    n_total = n(),
    n_oxy_pos = sum(def2 == 1),
    prop_oxy = n_oxy_pos / n_total
  )





a<-ggplot(prop_df, aes(rango_meses, prop_oxy, colour=Género, group=Género))+
  geom_point(size=3)+
  facet_wrap(~hospital)+
  geom_line(size=1, linetype = "solid")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ylab("Porcentaje de pacientes que reciben oxígeno suplementario")+
  xlab("Rango de edad (meses)")+
  theme_bw(base_size = 14)
a
ggsave("rango_edad_oxy_subgrupo_genero.jpg", a, dpi=300, width = 8, height = 4.5)
















































#####################################################################3
#Prior to modelling, i will calculate the lfa, wfa and wfl zscores for the children
#uso de anthro
library(anthro)
library(readr)
#usare edad en meses

#edad en dias
add_lab$Fecha.de.nacimiento.del.sujeto <- as.Date(add_lab$Fecha.de.nacimiento.del.sujeto, format="%Y-%m-%d")
add_lab <- add_lab %>%
  mutate(edad_dias = as.numeric(difftime(Fecha.de.ingreso.al.hospital, Fecha.de.nacimiento.del.sujeto, units="days")))
###modificar variable sexo
add_lab <- add_lab %>%
  mutate(sexo=case_when(
    Género=="Masculino"~"M",
    Género=="Femenino"~"F",
    TRUE~NA
  ))
zscores <- anthro_zscores(sex=add_lab$sexo, age=add_lab$edad_dias, weight = add_lab$Peso.en.Kilogramos..DEJAR.VACIO.SI.NO.SABE., lenhei = add_lab$Talla.en.centímetros..DEJAR.VACIO.SI.NO.SABE.)
ggplot(zscores, aes("x",zlen))+
  geom_jitter()+
  geom_hline(yintercept=3)+
  geom_hline(yintercept=-6)
  
ggplot(zscores, aes("x",zwei))+
  geom_jitter()+
  geom_hline(yintercept=3)+
  geom_hline(yintercept=-5)



ggplot(zscores, aes("x",zwfl))+
  geom_jitter()+
  geom_hline(yintercept=4)+
  geom_hline(yintercept=-5)



#there are some extreme values that fall out of what has been reported in literature, will consider only those
#with laz between -6 and +3 

#necesito agregar un id a las entradas del df original
library(dplyr)

add_lab <- add_lab %>%
  mutate(id = row_number())

zscores <- zscores %>%
  mutate(id = row_number())

add_lab_con_des <- add_lab %>%
  left_join(zscores, by="id")%>%
  filter(!is.na(zlen))


add_lab_con_des <- add_lab_con_des %>%
  filter(zlen >=-6 & zlen<=3)%>%
  filter(zwei>=-5&zlen<=3)%>%
  filter(zwfl>=-5&zwfl<=4)

#clasifico como wasting, underweight y stunting
add_lab_con_des <- add_lab_con_des %>%
  mutate(underweight = if_else(zwei< -2, "Underweight", "No underweight"))%>%
  mutate(stunting = if_else(zlen< -2, "Stunting", "No stunting"))%>%
  mutate(wasting = if_else(zwfl< -2, "Wasting", "No wasting"))
         
####################
#modelos
####################
library(dplyr)
library(readxl)
library(ggplot2)
library(MASS)
library(MuMIn)
library(DHARMa)
library(lme4)
library(nlme)
#voy a utilizar como variables predictoras: sexo, subgrupo de vsr, edad en años,
#peso-talla?, hospital, prematuridad



#verificacion de colinealidad
#Test de Colinealidad
source("C:/Users/julio/OneDrive/Escritorio/HighstatLibV10.R")
Var <- c("zwei","zlen", "zwfl", "edad_dias")
pairs(add_lab_con_des[,Var], lower.panel = panel.smooth,
      upper.panel = panel.cor,
      diag.panel = panel.hist)


Var1 <- c("zwei","zlen", "zwfl", "edad_dias") 
corvif(add_lab_con_des[,Var1])

Var2 <- c("zlen", "zwfl", "edad_dias") 
corvif(add_lab_con_des[,Var2])


subagrupados_con_des <- add_lab_con_des %>%
  filter(subgrupo_vsr!="No subagrupado")



#generacion de modelos candidatos para def 1
cand.models1 <- list()
cand.models1[["Completo"]] = glm(def1 ~ subgrupo_vsr + Género + rango_meses+hospital+Nació.prematuro....37.semanas.+wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models1[["solo sub"]] = glm(def1 ~ subgrupo_vsr, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models1[["sin genero"]] = glm(def1 ~ subgrupo_vsr +hospital + rango_meses+Nació.prematuro....37.semanas.+wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models1[["sin edad"]] = glm(def1~ subgrupo_vsr +hospital+ Género+Nació.prematuro....37.semanas.+wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models1[["solo edad"]] = glm(def1 ~rango_meses, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models1[["solo prem"]] = glm(def1 ~Nació.prematuro....37.semanas., data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models1[["edad y prem"]] = glm(def1 ~Nació.prematuro....37.semanas. +rango_meses, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models1[["solo nutri"]] = glm(def1 ~ wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models1[["solo genero"]] = glm(def1 ~ Género, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models1[["demo"]] = glm(def1 ~ Género+rango_meses+hospital, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models1[["nulo"]] = glm(def1 ~ 1, data = subagrupados_con_des, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models1,rank="AICc")
dredge.out
subset(dredge.out,delta<4)











#generacion de modelos candidatos para def 2
cand.models2 <- list()
cand.models2[["Completo"]] = glm(def2 ~ subgrupo_vsr + Género + rango_meses+hospital+Nació.prematuro....37.semanas.+wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models2[["solo sub"]] = glm(def2 ~ subgrupo_vsr, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models2[["sin genero"]] = glm(def2 ~ subgrupo_vsr +hospital + rango_meses+Nació.prematuro....37.semanas.+wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models2[["sin edad"]] = glm(def2~ subgrupo_vsr +hospital+ Género+Nació.prematuro....37.semanas.+wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models2[["solo edad"]] = glm(def2 ~rango_meses, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models2[["solo prem"]] = glm(def2 ~Nació.prematuro....37.semanas., data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models2[["edad y prem"]] = glm(def2 ~Nació.prematuro....37.semanas. +rango_meses, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models2[["solo nutri"]] = glm(def2 ~ wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models2[["solo genero"]] = glm(def2 ~ Género, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models2[["demo"]] = glm(def2 ~ Género+rango_meses+hospital, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models2[["nulo"]] = glm(def2 ~ 1, data = subagrupados_con_des, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models2,rank="AICc")
dredge.out
subset(dredge.out,delta<4)


avg_model <- model.avg(dredge.out, subset = delta < 4)

summary(avg_model)

confint(avg_model)

# Extract coefficients
coef_est <- coef(avg_model)

# Extract CI
ci <- confint(avg_model)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)

results


#generacion de modelos candidatos para def 4
cand.models3 <- list()
cand.models3[["Completo"]] = glm(def4 ~ subgrupo_vsr + Género + rango_meses+hospital+Nació.prematuro....37.semanas.+wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models3[["solo sub"]] = glm(def4 ~ subgrupo_vsr, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models3[["sin genero"]] = glm(def4 ~ subgrupo_vsr +hospital + rango_meses+Nació.prematuro....37.semanas.+wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models3[["sin edad"]] = glm(def4~ subgrupo_vsr +hospital+ Género+Nació.prematuro....37.semanas.+wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models3[["solo edad"]] = glm(def4 ~rango_meses, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models3[["solo prem"]] = glm(def4 ~Nació.prematuro....37.semanas., data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models3[["edad y prem"]] = glm(def4 ~Nació.prematuro....37.semanas. +rango_meses, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models3[["solo nutri"]] = glm(def4 ~ wasting+stunting, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models3[["solo genero"]] = glm(def4 ~ Género, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models3[["demo"]] = glm(def4 ~ Género+rango_meses+hospital, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
cand.models3[["nulo"]] = glm(def4 ~ 1, data = subagrupados_con_des, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models3,rank="AICc")
dredge.out
subset(dredge.out,delta<4)


summary(glm(def4 ~ Género+rango_meses+hospital, data = subagrupados_con_des, na.action=na.fail, family = "binomial"))

confint(glm(def4 ~ Género+rango_meses+hospital, data = subagrupados_con_des, na.action=na.fail, family = "binomial"))
mod <- glm(def4 ~ Género+rango_meses+hospital, data = subagrupados_con_des, na.action=na.fail, family = "binomial")
# Extract coefficients
coef_est <- coef(mod)

# Extract CI
ci <- confint(mod)

# Convert to OR
OR <- exp(coef_est)
OR_low <- exp(ci[,1])
OR_high <- exp(ci[,2])

# Combine
results <- data.frame(
  Variable = names(coef_est),
  Beta = coef_est,
  OR = OR,
  CI_low = OR_low,
  CI_high = OR_high
)

results






































cand.models3 <- list()
cand.models3[["Completo"]] = glm(def4 ~ VSR + Género + rango_meses+hospital+Nació.prematuro....37.semanas.+wasting+stunting, data = add_lab_con_des, na.action=na.fail, family = "binomial")
cand.models3[["solo sub"]] = glm(def4 ~ VSR, data = add_lab_con_des, na.action=na.fail, family = "binomial")
cand.models3[["sin genero"]] = glm(def4 ~ VSR +hospital + rango_meses+Nació.prematuro....37.semanas.+wasting+stunting, data = add_lab_con_des, na.action=na.fail, family = "binomial")
cand.models3[["sin edad"]] = glm(def4~ VSR +hospital+ Género+Nació.prematuro....37.semanas.+wasting+stunting, data = add_lab_con_des, na.action=na.fail, family = "binomial")
cand.models3[["solo edad"]] = glm(def4 ~rango_meses, data = add_lab_con_des, na.action=na.fail, family = "binomial")
cand.models3[["solo prem"]] = glm(def4 ~Nació.prematuro....37.semanas., data = add_lab_con_des, na.action=na.fail, family = "binomial")
cand.models3[["edad y prem"]] = glm(def4 ~Nació.prematuro....37.semanas. +rango_meses, data = add_lab_con_des, na.action=na.fail, family = "binomial")
cand.models3[["solo nutri"]] = glm(def4 ~ wasting+stunting, data = add_lab_con_des, na.action=na.fail, family = "binomial")
cand.models3[["solo genero"]] = glm(def4 ~ Género, data = add_lab_con_des, na.action=na.fail, family = "binomial")
cand.models3[["demo"]] = glm(def4 ~ Género+rango_meses+hospital, data = add_lab_con_des, na.action=na.fail, family = "binomial")
cand.models3[["nulo"]] = glm(def4 ~ 1, data = add_lab_con_des, na.action=na.fail, family = "binomial")


dredge.out<-model.sel(cand.models3,rank="AICc")
dredge.out
subset(dredge.out,delta<2)


mejor_modelo_def3 <- glm(def4 ~ Género+rango_meses+hospital, data = add_lab_con_des, na.action=na.fail, family = "binomial")
confint(mejor_modelo_def3)

ggplot(add_lab_con_des, aes(VSR))+
  geom_bar()+
  facet_wrap(~def4)

simulationOutput3 <- simulateResiduals(fittedModel = mejor_modelo_def3, plot = T)


################

subagrupados$def1 <- as.factor(subagrupados$def1)
#######solo coatepeque
sub_coat <- subagrupados_con_des %>%
  filter(hospital=="Hospital Nacional de Coatepeque")

#generacion de modelos candidatos para def 2

cand.models2 <- list()
cand.models2[["Completo"]] = glm(def2 ~ subgrupo_vsr + Género + rango_meses+Nació.prematuro....37.semanas.+wasting+stunting, data = sub_coat, na.action=na.fail, family = "binomial")
cand.models2[["solo sub"]] = glm(def2 ~ subgrupo_vsr, data = sub_coat, na.action=na.fail, family = "binomial")
cand.models2[["sin genero"]] = glm(def2 ~ subgrupo_vsr  + rango_meses+Nació.prematuro....37.semanas.+wasting+stunting, data = sub_coat, na.action=na.fail, family = "binomial")
cand.models2[["sin edad"]] = glm(def2 ~ subgrupo_vsr + Género+Nació.prematuro....37.semanas.+wasting+stunting, data = sub_coat, na.action=na.fail, family = "binomial")
cand.models2[["solo edad"]] = glm(def2 ~rango_meses , data = sub_coat, na.action=na.fail, family = "binomial")
cand.models2[["solo prem"]] = glm(def2 ~Nació.prematuro....37.semanas. , data = sub_coat, na.action=na.fail, family = "binomial")
cand.models2[["edad y prem"]] = glm(def2 ~Nació.prematuro....37.semanas. +rango_meses, data = sub_coat, na.action=na.fail, family = "binomial")
cand.models2[["solo nutri"]] = glm(def2 ~ wasting+stunting, data = sub_coat, na.action=na.fail, family = "binomial")
cand.models2[["solo genero"]] = glm(def2 ~ Género, data = sub_coat, na.action=na.fail, family = "binomial")
cand.models2[["demo"]] = glm(def2 ~ Género+rango_meses, data = sub_coat, na.action=na.fail, family = "binomial")
cand.models2[["nulo"]] = glm(def2 ~ 1, data = sub_coat, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models2,rank="AICc")
dredge.out
subset(dredge.out,delta<2)

mejor_modelo <- glm(def2 ~ subgrupo_vsr + Género + rango_meses+Nació.prematuro....37.semanas.+wasting+stunting, data = sub_coat, na.action=na.fail, family = "binomial")
confint(mejor_modelo)

#######solo chimaltenango
sub_chimal <- subagrupados_con_des %>%
  filter(hospital=="Hospital Nacional de Chimaltenango")

#generacion de modelos candidatos para def 2

cand.models2 <- list()
cand.models2[["Completo"]] = glm(def2 ~ subgrupo_vsr + Género + rango_meses+Nació.prematuro....37.semanas.+wasting+stunting, data = sub_chimal, na.action=na.fail, family = "binomial")
cand.models2[["solo sub"]] = glm(def2 ~ subgrupo_vsr, data = sub_chimal, na.action=na.fail, family = "binomial")
cand.models2[["sin genero"]] = glm(def2 ~ subgrupo_vsr  + rango_meses+Nació.prematuro....37.semanas.+wasting+stunting, data = sub_chimal, na.action=na.fail, family = "binomial")
cand.models2[["sin edad"]] = glm(def2 ~ subgrupo_vsr + Género+Nació.prematuro....37.semanas.+wasting+stunting, data = sub_chimal, na.action=na.fail, family = "binomial")
cand.models2[["solo edad"]] = glm(def2 ~rango_meses , data = sub_chimal, na.action=na.fail, family = "binomial")
cand.models2[["solo prem"]] = glm(def2 ~Nació.prematuro....37.semanas. , data = sub_chimal, na.action=na.fail, family = "binomial")
cand.models2[["edad y prem"]] = glm(def2 ~Nació.prematuro....37.semanas. +rango_meses, data = sub_chimal, na.action=na.fail, family = "binomial")
cand.models2[["solo nutri"]] = glm(def2 ~ wasting+stunting, data = sub_chimal, na.action=na.fail, family = "binomial")
cand.models2[["solo genero"]] = glm(def2 ~ Género, data = sub_chimal, na.action=na.fail, family = "binomial")
cand.models2[["demo"]] = glm(def2 ~ Género+rango_meses, data = sub_chimal, na.action=na.fail, family = "binomial")
cand.models2[["nulo"]] = glm(def2 ~ 1, data = sub_chimal, na.action=na.fail, family = "binomial")

dredge.out<-model.sel(cand.models2,rank="AICc")
dredge.out
subset(dredge.out,delta<2)

mejor_modelo <- glm(def2 ~ subgrupo_vsr + Género + rango_meses+Nació.prematuro....37.semanas.+wasting+stunting, data = sub_coat, na.action=na.fail, family = "binomial")
confint(mejor_modelo)



subagrupados_con_des <- subagrupados_con_des %>%
  mutate(congestion=case_when(
    Congestión.nasal...Rinorrea=="Si"~"1",
    Congestión.nasal...Rinorrea=="No"~"0",
    TRUE~NA
  ))%>%
  filter(!is.na(congestion))

subagrupados_con_des <- subagrupados_con_des %>%
  mutate(Sibilancias_F=case_when(
    Sibilancias=="Si"~"1",
    Sibilancias=="No"~"0",
    TRUE~NA
  ))%>%
  filter(!is.na(Sibilancias_F))

subagrupados_con_des <- subagrupados_con_des %>%
  mutate(Fiebre=case_when(
    Fiebre...antecedentes.de.fiebre=="Si"~"1",
    Fiebre...antecedentes.de.fiebre=="No"~"0",
    TRUE~NA
  ))%>%
  filter(!is.na(Fiebre))

########figures
#assign to each sample its epiweek and year
under_five <- under_five %>%
  mutate(
    semana = isoweek(date),
    year = isoyear(date)
  )

#group by week, year, hospital
weekly_counts<-under_five%>%
  group_by(semana, year, hospital)%>%
  summarise(
    total_samples=n(),
    rsv_pos = sum(Resultado.de.Virus.sincitial.respiratorio=="Positivo", na.rm=TRUE))




#epi curve combining both sites

under_five <- under_five %>%
  filter(Resultado.de.Virus.sincitial.respiratorio!="No sabe")

under_five$Resultado.de.Virus.sincitial.respiratorio <-
  recode(under_five$Resultado.de.Virus.sincitial.respiratorio,
         "Negativo" = "Negative",
         "Positivo" = "Positive")

i<-ggplot(under_five, aes(factor(paste(year, semana, sep = "-")), fill=Resultado.de.Virus.sincitial.respiratorio))+
  geom_bar(color="black")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 45))+
  xlab("Epidemiological week")+
  scale_fill_manual(values = c(
    "Negative" = "grey", 
    "Positive" = "steelblue"
  ))+
  labs(fill="RSV result", y="Number of enrolled")
i
ggsave("curva_hosp_ambos.jpg", i, dpi=300, height = 4.5, width = 8)

#plot the epi curve for both sites independently

h<-ggplot(under_five, aes(factor(paste(year, semana, sep = "-")), fill=Resultado.de.Virus.sincitial.respiratorio))+
  geom_bar(color="black")+
  xlab("Epidemiological week")+
  scale_fill_manual(values = c(
    "Negative" = "grey", 
    "Positive" = "steelblue"
  ))+
  labs(fill="RSV result", y="Number of enrolled")+
  facet_wrap(~hospital, ncol=1)+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
h
ggsave("curva_hosp.jpg", h, dpi=300, height = 4.5, width = 8)

weekly_pos <- under_five %>%
  group_by(hospital, year, semana) %>%
  summarise(
    total = n(),
    positive = sum(Resultado.de.Virus.sincitial.respiratorio == "Positive",
                   na.rm = TRUE),
    positivity = positive / total
  )




######heatmap de sintomas
#re-locate some of the variables
under_five <- under_five %>%
  relocate(Definición.de.caso.utilizada..Consulte.el.Anexo.2...1..Definición.de.caso.de.SARI..2..Definición.ampliada.de.caso.de.SARI..3..Definición.de.caso.modificada.del.ECDC, .after = last_col())


sintomas_long <- under_five %>%
  pivot_longer(
    cols = Fiebre...antecedentes.de.fiebre:Dolor.de.pecho,
    names_to = "sintoma",
    values_to = "presente"
  ) 
sintomas_long <- sintomas_long %>%
  mutate(presente = ifelse(presente == "Si", 1, 0))

heatmap_data <- sintomas_long %>%
  group_by(sintoma, Resultado.de.Virus.sincitial.respiratorio) %>%
  summarise(
    frecuencia = mean(presente, na.rm = TRUE)
  )

heatmap_data <- heatmap_data %>%
  filter(sintoma=="Congestión.nasal...Rinorrea"|sintoma=="Fiebre...antecedentes.de.fiebre"|
           sintoma=="Falta.de.aire.dificultad.para.respirar"|sintoma=="Diarrea"|sintoma=="Sibilancias"|
           sintoma=="Náuseas.o.vómitos")

ggplot(heatmap_data, aes(x = sintoma, y = 1, fill = frecuencia)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(
    title = "Frecuencia de síntomas en GIHSN",
    x = "Virus",
    y = "Síntoma",
    fill = "Frecuencia"
  ) +
  facet_wrap(~Resultado.de.Virus.sincitial.respiratorio, ncol=1)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=45))



