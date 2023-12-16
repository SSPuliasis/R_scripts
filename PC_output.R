library(tidyr)
library(readxl)

PC_digest <- read_excel("AT1G67120_PC.xlsx")
View(PC_digest)
attach(PC_digest)

filled_in <- fill(PC_digest, `Position of cleavage site`:`Peptide mass [Da]`, .direction=c("down"))
colnames(AT1G67120_PC_table) <- c("cleavage_position", "enzyme", "sequence", "peptide_size", "mol_weight")

write.csv(filled_in, "AT1G67120_PC_table.csv")

detach(PC_digest)
