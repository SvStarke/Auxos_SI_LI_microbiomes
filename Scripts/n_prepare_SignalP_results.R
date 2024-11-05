library(data.table)

sgp <- dir("/mnt/nuuk/2021/HRGM/signalP/", full.names = T)

sgp <- lapply(sgp, function(x) fread(x))

sgp <- rbindlist(sgp)

sgp <- sgp[Prediction != "OTHER"]

sgp[, gene := gsub(" .*$","",`# ID`)]

fwrite(sgp, "Data/signalP_predictions.tsv.gz", sep = "\t", quote = FALSE)
rm(sgp)
