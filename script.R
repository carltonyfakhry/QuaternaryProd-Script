require(QuaternaryProd)
require(fdrtool)

source('./utils.R')

args = commandArgs(TRUE);

ents.all.file      = args[[1]] ## KB ent file
rels.all.file      = args[[2]] ## KB rel file
evidence.file      = args[[3]] ## up or down regulated genes
out.file           = args[[4]] ## output

L = processData(ents.all.file, rels.all.file, evidence.file)
Tab = generateCREtable(L$ents, L$rels, L$evidence)

write.table(Tab, out.file, quote = F, sep = '\t', row.names = F, col.names = T)

