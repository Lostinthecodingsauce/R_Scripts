canopus <- read_csv("canopus_summary.tsv.csv")

canopus_binary <- canopus[2:ncol(canopus)]
canopus_binary[canopus_binary > 0.75] <- 1
canopus_binary[canopus_binary < 0.75] <- 0
canopus_binary <- cbind(canopus$name, canopus_binary)


write_csv(canopus_binary, "test.csv")


matric <- matrix(norm(5*6, mean = 20, sd = 2),5,6)
'>5'


Irina_yes <- Irina_fixed%>%
  dplyr::select(-c(newcolumn))

Irina_important <- Irina_yes%>%
  dplyr::filter(RF_Group1_Top500 =="important")

Irina_others <-Irina_yes%>%
  dplyr::filter(!RF_Group1_Top500 =="important")


Irina_fixed$no <- sapply(Irina_fixed$Count, function(x) x-1000)


bind_rows(Irina_others, Irina_important)