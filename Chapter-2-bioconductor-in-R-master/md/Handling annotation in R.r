source("http://www.bioconductor.org/biocLite.R")
biocLite("hgu133a.db")
library(hgu133a.db)
# `hgu133a.db` is the chip annotation package for Affymetrix Human Genome U133 Set, 
# which contains mappings between an Affymetrix's identifiers and accessions. 


mymap <- hgu133aENTREZID 
mymap
# `hgu133aENTREZID` provides mappings of manufacturer identifiers to a vector of entrez gene identifiers.
# mymap shows a map between ENtrez ids and probes, the mapped probes.

mapped_probes <- mappedkeys(mymap)
mapped_probes
# this shows a list of mapped probes.
# 'mappedkeys' function gets the probe identifiers that are mapped to a entrez gene id as an object and later the selected IDs 
# can be extracted.
           
myentrez <- as.list(mymap[mapped_probes[1:5]])
myentrez
# myentrez shows the first five entrez ids (mapped_probes) extracted in the mapped_probes list
# another way to extract entrez ids:

mylength = 3          # for entire list, 'mylength = length(myentrez)'
for(i in 1: mylength){
myentrez[[i]]
}
# be sure that `mylength` is less than 'length(myentrez)'

ls("package:hgu133a.db")
# check the 'hgu133a.db' package
# 'ls("package:<package_name>") is used to have an overview of all possible mappings in the package.

