dRaw <- read.csv("../../Data/Daten_für_Grafik_Version_2_3.csv", header=TRUE, sep=';', na.strings='')
dRaw <- dRaw[!is.na(dRaw$diameter), ]

dRaw$phenotypes <- substring(dRaw$phenotypes, 2)  # removes leading ";"
dRaw$phenotypes[is.na(dRaw$phenotypes == "")] = 'NA'
dRaw$phenotypes <- factor(dRaw$phenotypes)
new <- sapply(strsplit(levels(dRaw$phenotypes), ';'), function(x) {if(length(x) > 1) {paste(x[-1], collapse='; ')} else {x}})
names(new) = levels(dRaw$phenotypes)
dRaw$phenotypes <- revalue(dRaw$phenotypes, new)

dRaw$checkin <- factor(dRaw$checkin)

dRaw$hDiscr <- factor(round(dRaw$h))


dRaw %>% select(organism, checkin, atb_acronim, phenotypes) %>% filter(phenotypes %in% paste(c("AG", "BL", "FA", "MLS", "RIF", "SXT"), "non-wild-type")) %>% distinct %>% group_by(organism, atb_acronim, phenotypes) %>% summarise(count=n()) %>% data.frame %>% write.table(file='unspecifiedNonWildType.csv', sep=';', row.names=FALSE)


dRaw %>% select(organism, checkin, atb_acronim, phenotypes) %>% filter(phenotypes %in% paste(c("AG", "BL", "FA", "MLS", "OXD", "QL", "RIF", "SXT", "TET"), "wild-type")) %>% distinct %>% group_by(organism, atb_acronim, phenotypes) %>% summarise(count=n()) %>% data.frame %>% write.table(file='nWildType.csv', sep=';', row.names=FALSE)
