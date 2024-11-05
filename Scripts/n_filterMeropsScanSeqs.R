merof <- readAAStringSet("/mnt/nuuk/Resources/MEROPS/merops_scan.lib")
mero_fam <- fread("Data/domain2.csv.gz")

names(merof) <- substr(names(merof), 1,10)

merof <- merof[names(merof) %in% mero_fam$mernum]
writeXStringSet(merof, "/mnt/nuuk/Resources/MEROPS/merops_scan_filt.lib")

fwrite(data.table(mernum = names(merof),
                  len = width(merof)), "Data/merops_scan_filt_length.csv")
