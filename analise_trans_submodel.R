library(tidyverse)

# the original study by marks
# https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30985-3/fulltext
# later analyzed by Marc et al.
df_bahv <- read_csv("./elife-69302-data1-v2-converted.csv")
# Elife

length(unique(df_bahv$case_id))

df_bahv %>%
  filter(!is.na(case_day0_VL)) -> df_inf

length(unique(df_inf$case_id))

df_bahv %>%
  select(case_id, case_sex) %>%
  distinct(case_id, case_sex) -> xx2

table(xx2$case_sex)
# 202 women

# this should have problems
glm1 <- glm(infected_contact ~ case_day0_VL, family = "binomial", data = df_inf)
summary(glm1)
