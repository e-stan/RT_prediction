options(java.parameters = "-Xmx4096m")
library(rcdk)
setwd("D:/covid/ID/RetentionTimes/")

data = read.csv("MoNA_smiles_of_interest.csv")

mols <- parse.smiles(as.character(data[,2]))

descNames = unique(unlist(sapply(get.desc.categories(), get.desc.names)))
first = TRUE

for (i in 1:length(mols))
{
  tryCatch({
    if (first){
      df = eval.desc(mols[i],descNames)
      df = cbind(df,data[i,1])
      first = FALSE
    }
    else{
      df = rbind(df,cbind(eval.desc(mols[i],descNames),data[i,1]))#,data[i,3]))
    }
  }
  ,error = function(e) {
    
  })
}

write.csv(df,"molecular_properties_mona_of_interest.csv")
