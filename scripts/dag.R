library(dagitty)

dag <-dagitty('dag{
  "Population structure" -> "Host genetics" -> "Cirrhosis" -> "Immune cell counts"
  "Population structure" -> "Immune cell counts"
  "Host genetics" -> "Immune cell counts"
  Cirrhosis <- Sex -> "Immune cell counts"
  Cirrhosis <- Age -> "Immune cell counts"
  Cirrhosis <- "IFNL4 SNP" -> "Immune cell counts"
  "Host genetics" -> "Viral load" <- "Sex"
  "Host genetics" -> "IFNL4 SNP" -> "Viral load" -> "Immune cell counts"
  Cirrhosis <- "Viral load" <- "Previous treatment" -> "Immune cell counts"
  Cirrhosis <- "Previous treatment" -> "Immune cell counts"
}')

adjustmentSets( dag, exposure="Host genetics", outcome="Immune cell counts", type = "minimal", effect="direct")
plot(dag)



