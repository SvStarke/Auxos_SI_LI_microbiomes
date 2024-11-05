library(stringr)
library(data.table)

mero_domains <- readLines("/mnt/nuuk/Resources/MEROPS/domain.sql", encoding = "UTF-8") 
MEROlines <- grep("MER[0-9]{7}", mero_domains)

dfs <- lapply(MEROlines, function(x) {
  print(x)
  values_part <- str_extract(mero_domains[x], "(?<=VALUES \\().*(?=\\))")
  records <- str_split(values_part, "\\),\\(")[[1]]
  clean_records <- gsub("^\\(|\\)$", "", records)
  clean_records <- gsub("NULL", "NA", clean_records)
  records_list <- str_split(clean_records, ",(?=(?:[^']*'[^']*')*[^']*$)")
  records_list <- lapply(records_list, function(x) gsub("^\\'|\\'$","",x))
  as.data.table(do.call(rbind, records_list))
})
dfsBound <- rbindlist(dfs)
dfsBound <- dfsBound[V7 != "inhibitor"]

Merocolnames <- mero_domains[26:42]
Merocolnames <- str_extract(Merocolnames,  "(?<=`)[^`]+(?=`)")

setnames(dfsBound, Merocolnames)

fwrite(dfsBound, "Data/domain2.csv.gz")
