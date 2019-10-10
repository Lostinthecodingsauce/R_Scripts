#PO_filt, etc. Taken directly from "Making Pie Charts.R script because im lazy.
#Fucking. Deal. With. It.
PO_filt_05 <- PO_filt%>%
  group_by(Organism, Order)%>%
    dplyr::filter(RA, sum(RA[RA>0.02]))%>%
    ungroup(PO_filt_05)%>%
  arrange(Order)

HR_filt_05 <- HR_filt%>%
  group_by(Organism, Order)%>%
  dplyr::filter(RA, sum(RA[RA>0.02]))%>%
  ungroup(HR_filt_05)%>%
  arrange(Order)

PO_unf_05 <- PO_unf%>%
  group_by(Organism, Order)%>%
  dplyr::filter(RA, sum(RA[RA>0.02]))%>%
  ungroup(PO_unf_05)%>%
  arrange(Order)

HR_unf_05 <- HR_unf%>%
  group_by(Organism, Order)%>%
  dplyr::filter(RA, sum(RA[RA>0.02]))%>%
  ungroup(HR_unf_05)%>%
  arrange(Order)

PO_CCA_05 <- PO_CCA%>%
  group_by(Organism, Order)%>%
  dplyr::filter(RA, sum(RA[RA>0.02]))%>%
  ungroup(PO_CCA_05)%>%
  arrange(Order)

HR_CCA_05 <- HR_CCA%>%
  group_by(Organism, Order)%>%
  dplyr::filter(RA, sum(RA[RA>0.02]))%>%
  ungroup(HR_CCA_05)%>%
  arrange(Order)