library(data.table)

### analysis of groups and isolation source in metadata 
meta <- fread("Data/Metadata_Reimagine.txt")

# Extract sampling location
meta[, location := NA_character_]
meta[grepl("[0-9]+J", `Sample Name`), location := "Jejunum"] # regulärer ausdruck: Sample name beginnt mit mindestens einer Zahl "[0-9]+"; gefolgt von einem großen "J" 
meta[grepl("[0-9]+D", `Sample Name`), location := "Duodenum"]
meta[grepl("[0-9]+I", `Sample Name`), location := "Ileum"]
meta[grepl("[0-9]+S", `Sample Name`), location := "Stool"]

# extract patient ID
meta[, PatNo := NA_character_]
meta[, PatNo := gsub("(I|J|D|S)-DBE$","",`Sample Name`)] # Regulärer Ausdruck: Ab einem J,I,S oder D bis zum Ende des Strings wird alles abgeschnitten


# Einteilung in Gruppen
meta[, Group1 := "Duodenum" %in% location & "Colon" %in% location , by = PatNo]
#Correction of value --> too much in Group 1 !!!!

meta[, Group2 := "Duodenum" %in% location & "Ileum" %in% location & "Jejunum" %in% location, by = PatNo]
meta[, Group3 := "Duodenum" %in% location & "Ileum" %in% location & "Jejunum" %in% location & "Stool" %in% location, by = PatNo]
meta[Group3 == TRUE,]
