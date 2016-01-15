#################################################################################
## This file contains formating and utility files for running the the 
## QuaternaryProd package. It works with stringDB as well as CRE KB.
##
##
## Last Modified: 1/11/2016
#################################################################################


## This function cleans up the KB.
##
## Arguments:
##
## ents.file     -  entry files from the knowledge base
##
## rels.file     -  relation files from the knowledge base
##
## L             -  the returned list of pruned entries and relations
##
## Note: It is assumed that the entity file is a tab seperated text file with the
## following columns:
##
## 1 - uid. This is a unique identifier for the entity (character).
## 2 - name. This is the entity name (character).
## 3 - id. This should be the entrez id for mRNAs and an string (arbitrary) for others.
## 4 - type. This is the entity type. Must be one of mRNA, Protein or Compound.
##
## The relation file is assumed to have the following columns.
##
## 1 - uid. This is a unique identifier for the relation (character).
## 2 - srcuid. This is the source uid as deteremined by the entity file.
## 3 - trguid. This is the target uid as deteremined by the entity file.
## 4 - type. This is the type of regulation. One of "increase", "decrease", "conflict", or ""
## 5 - pmid. This is the pmid of the paper reporting the relation (can be "").
## 6 - nls. This is a text describing the relation (can be "").

processKB = function(ents.file, rels.file){
  ## Generating the one-level network form the knowledge base
  ## Note: There are ' and # characters in KB that mess up the tables. The following will clean them up.
  ents = read.table(ents.file, header = T, stringsAsFactors = F, strip.white=T, sep = '\t', quote = NULL,
                    comment.char = '')
  colnames(ents) = c('uid', 'name', 'id', 'type')
  rels = read.table(rels.file, header = T, stringsAsFactors = F, strip.white=T, sep = '\t', quote = NULL,
                    comment.char = '')
  colnames(rels) = c('uid', 'srcuid', 'trguid', 'type', 'pmids', 'nls')

  ents = data.frame(apply(ents, 2, function(x) gsub('\"', '', x)), stringsAsFactors = F)
  ents = data.frame(apply(ents, 2, function(x) gsub('\'', 'p', x)), stringsAsFactors = F)
  ents = data.frame(apply(ents, 2, function(x) gsub('#', '_', x)), stringsAsFactors = F)
  ents = data.frame(apply(ents, 2, function(x) gsub(' ', '', x)), stringsAsFactors = F)

  rels = data.frame(apply(rels, 2, function(x) gsub('\"', '', x)), stringsAsFactors = F)
  rels = data.frame(apply(rels, 2, function(x) gsub('\'', 'p', x)), stringsAsFactors = F)
  rels = data.frame(apply(rels, 2, function(x) gsub('#', '_', x)), stringsAsFactors = F)
  rels = data.frame(apply(rels, 2, function(x) gsub(' ', '', x)), stringsAsFactors = F)

  rownames(rels) = 1:nrow(rels)

  ## Convert uids to integers
  uid.orig = ents$uid
  uid.new = seq(1, length(ents$uid))
  id.map = data.frame(uid.orig = uid.orig, uid.new = uid.new, stringsAsFactors = F)
  ents$uid = uid.new
  rels$srcuid = uid.new[match(rels$srcuid, uid.orig)]
  rels$trguid = uid.new[match(rels$trguid, uid.orig)]

  L = list(ents = ents, rels = rels, id.map = id.map)

  L
}

## genOneLevel:
##
## Given a knowledge-base, this function generates a one-level causal network from the protein/chemical
## /mRNA entries. The hypothesis layer consists of protiens and chemicals while the evidence layer
## consistes of mRNAs only. The i-->i and cyles are removed from the causal network as well.
##
## Arguments:
##
## ents.file     -  entry files from the knowledge base
##
## rels.file     -  relation files from the knowledge base
##
## L                 -  the returned list of pruned entries and relations
genOneLevel = function(ents.file, rels.file){
  ## pre-process KB
  print('processing the network...')
  L = processKB(ents.file, rels.file)

  ents.all = L$ents
  rels.all = L$rels

  ## Protein, Compound or mRNA entries
  ents = unique(ents.all[which(ents.all$type %in% c('mRNA', 'Protein', 'Compound')),])

  ## (unique) Protein/Compound entries
  ents.pc = unique(ents[which(ents$type %in% c('Protein', 'Compound')),])
  ## (unique) mRNA entries
  ents.mRNA = unique(ents[which(ents$type == 'mRNA'),])

  ## src has to be protiens or a compunds
  rels = rels.all[which(rels.all$srcuid %in% ents.pc$uid),]
  rownames(rels) = 1:nrow(rels)
  ## target has to be mRNA
  rels = rels[which(rels$trguid %in% ents.mRNA$uid),]
  rownames(rels) = 1:nrow(rels)

  ## Only consider entries that are participating in valid relations
  ents = ents[which(ents$uid %in% c(rels$srcuid, rels$trguid)), ]
  rownames(ents) = 1:nrow(ents)

  ## Identify duplicated rows
  dup.rels = rels[duplicated(rels[,c(2,3)]), ]
  ## Take one example from each
  dup.rels.uniq = dup.rels[!duplicated(dup.rels[,c(2,3)]), ]

  for(i in 1:nrow(dup.rels.uniq)){
    ind = which(rels$srcuid == dup.rels.uniq$srcuid[i] & rels$trguid == dup.rels.uniq$trguid[i])
    if(all(c("increase","decrease") %in% unique(rels$type[ind]))){
      rels$type[ind] = 'conflict'
    }else if("increase" %in% unique(rels$type[ind])){
      rels$type[ind] = 'increase'
    }else if("decrease" %in% unique(rels$type[ind])){
      rels$type[ind] = 'decrease'
    }else{
      rels$type[ind] = 'conflict'
    }
  }
  rels$type[which(!(rels$type %in% c('increase', 'decrease', 'conflict')))] = 'conflict'

  ## Remove the duplicates
  rels = rels[!duplicated(rels[,c(2,3)]),]
  rownames(rels) = 1:nrow(rels)

  ## Filter regulators that have less than 5 genes underneath
  u.hyps = unique(rels$srcuid)
  num.targs = unlist(lapply(u.hyps, function(x) length(which(rels$srcuid == x))))
  keep.ind  = which(num.targs >= 5)
  rels = rels[which(rels$srcuid %in% u.hyps[keep.ind]), ]
  rownames(rels) = 1:nrow(rels)
  ## Update ents
  ents = ents[which(ents$uid %in% c(rels$srcuid, rels$trguid)), ]
  rownames(ents) = 1:nrow(ents)

  print('processed network dimenssions:')
  print(paste('ents:', dim(ents)[1]))
  print(paste('rels:', dim(rels)[1]))

  L = list(ents = ents, rels = rels)

  L
}


## Given a list of gene uids and the gene expression file, this function reads the values
## of the genes from the gene expression values (if they exist)
##
## uids       -   uids of the list of genes
##
## evidence   -   gene expression data
getGeneVals = function(uids, evidence){
  val = rep(0, length(uids))
  ind = match(uids, evidence$uid)
  ind2 = na.omit(ind)
  if(length(ind2) > 0){
    val[which(!is.na(ind))] = evidence$val[ind2]
  }
  val
}


## processData
## This function prepares the CRE input format data
##
## Arguments
##
## ents.all.file     -       entries
##
## rels.all.file     -       relationships
##
## evidence.file     -       evidence file
processData = function(ents.all.file, rels.all.file, evidence.file){
  
  ## Generating the one level-network form the knowledge base
  print('pre-processing...')
  L = genOneLevel(ents.all.file, rels.all.file)
  ents = L$ents
  rels = L$rels
  ents.mRNA = ents[which(ents$type == 'mRNA'),]
  
  
  ## Read in the evidence 
  evidence = read.table(evidence.file, header = T, sep = '\t', stringsAsFactors = F)
  if(ncol(evidence) == 3){
    pval.ind = grep('pval|p.val|p-val|P-val', colnames(evidence), ignore.case = T)
    fc.ind = grep('fc|FC|fold', colnames(evidence), ignore.case = T)
    id.ind = grep('id|entr', colnames(evidence), ignore.case = T)
    if(length(id.ind) == 0 | length(fc.ind) == 0 | length(pval.ind) == 0){
      print('Please make sure the expression files column names are labled as entrez, fc, pvalue')
      quit(save = "no", status = 1, runLast = FALSE)
    }
    
    p.thresh = 0.05
    isLog = grep('log', colnames(evidence)[fc.ind], ignore.case = T)
    if(length(isLog) > 0){
      fc.thresh = log2(1.3)
    }else{
      fc.thresh = 1.3
    }
    
    #fc.thresh = 0
    evidence = evidence[which(abs(evidence[,fc.ind]) >= fc.thresh & 
                                evidence[,pval.ind] <= p.thresh ),c(id.ind, fc.ind)]
    evidence[,2] = ifelse(evidence[,2] > 0, 1, -1)
    
    evidence = data.frame(id = suppressWarnings(as.numeric(evidence[,1])), val = as.numeric(evidence[,2]),
                          stringsAsFactors = F)
    evidence = evidence[!is.na(evidence[,1]),]
  }else if(ncol(evidence) == 2){
    id.ind = grep('id|entr|Entr', colnames(evidence))
    val.ind = grep('val', colnames(evidence))
    if(length(id.ind) == 0 | length(val.ind) == 0){
      print('Please provide a valid expression data file:')
      print('1: a file containing columns: entrez, fc, pvalue')
      print('2: a file containing columns: entrez, val, whre val is 1, -1, or 0')
      quit(save = "no", status = 1, runLast = FALSE)
    }
    colnames(evidence) = c('id', 'val')
  }else{
    print('Please provide a valid expression data file:')
    print('1: a file containing columns: entrez, fc, pvalue')
    print('2: a file containing columns: entrez, val, whre val is 1, -1, or 0')
    quit(save = "no", status = 1, runLast = FALSE)
  }
  evidence = evidence[!duplicated(evidence[,1]), ]
  #evidence = read.table(evidence.file, header = T, stringsAsFactors = F, strip.white=T, sep = '\t')
  n.e1 = nrow(evidence)
  ## Make sure evidence id is in the ents. Same gene with multiple ids may exist
  ## in the original ents. 
  evidence = evidence[which(evidence[,1] %in% ents.mRNA$id),]
  n.e2 = nrow(evidence)
  print(paste((n.e1-n.e2), "evidence removed!"))
  
  ## Change id back to uid
  evidence.tmp = merge(evidence, ents.mRNA, by.x = 1, by.y = 3)
  if(nrow(evidence.tmp) > nrow(evidence)){
    print('Warning')
    print('Entrez ids in evidence file mapped to multiple uids in KB')
    print(paste('total number of duplicates:' , (nrow(evidence.tmp) - nrow(evidence))))
    print('selecting one uid at random')
  }
  evidence = data.frame(uid = evidence.tmp$uid, val = evidence.tmp$val, stringsAsFactors = F)
  evidence = evidence[!duplicated(evidence),]
  
  print('dim of input file:')
  print(dim(evidence))
  
  L = list(ents = ents, rels = rels, evidence = evidence)
  
  L
}  


## This function generates the prior pvalues based on the CRE analysis
##
## Arguments
##
## ents              -       entries in the KB
##
## rels              -       relationships in the KB
##
## evidence          -       up or down regulated genes. See getDEGs.
generateCREtable = function(ents, rels, evidence){
  ents.mRNA = ents[which(ents$type == 'mRNA'), ]
  ## For each hypothesis, identify the children and non-children and thier evidence values
  u.hyps = unique(rels$srcuid)
  child.uid = lapply(u.hyps, function(x) rels$trguid[which(rels$srcuid == x)])
  child.sgn = lapply(u.hyps, function(x) ifelse(rels$type[which(rels$srcuid == x)] == 'increase',
                                                1, ifelse(rels$type[which(rels$srcuid == x)] == 'decrease', -1, 0)))

  child.val = lapply(child.uid, function(x) getGeneVals(x, evidence))

  non.child.uid = lapply(child.uid, function(x) unique(ents.mRNA$uid[which(!(ents.mRNA$uid %in% x))]))
  non.child.val = lapply(non.child.uid, function(x) getGeneVals(x, evidence))

  ## Get the data slices corresponding to each hypothesis
  child.id = lapply(child.uid, function(x) as.numeric(ents.mRNA$id[match(x,ents.mRNA$uid)])) ## to get the id


  print('computing pvalues')
  results = data.frame(matrix(0, nrow  = 2 * length(u.hyps), ncol = 17), stringsAsFactors = F)
  colnames(results) = c('uid', 'name', 'regulation', 'correct.pred', 'incorrect.pred', 'score',
                        'total.reachable', 'significant.reachable', 'total.ambiguous', 'significant.ambiguous',
                        'unknow', 'quaternary.pval', 'ternary.pval', 'enrichment.pval',
                        'quaternary.qval', 'ternary.qval', 'enrichment.qval')

  print(paste('Total number of hypothesis to consider:', length(u.hyps)))
  for(p.s in 1:length(u.hyps)){
    #cat('.')
    results[(2*(p.s-1)+1), 1] = u.hyps[p.s]
    results[(2*p.s), 1]       = u.hyps[p.s]
    results[(2*(p.s-1)+1), 2] = ents$name[which(ents$uid == u.hyps[p.s])]
    results[(2*p.s), 2]       = ents$name[which(ents$uid == u.hyps[p.s])]
    results[(2*(p.s-1)+1), 3] = 'up'
    results[(2*p.s), 3]       = 'down'

    npp = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == 1))
    npm = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == -1))
    npz = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == 0))

    nmp = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == 1))
    nmm = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == -1))
    nmz = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == 0))

    nrp = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == 1))
    nrm = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == -1))
    nrz = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == 0))

    nzp = length(which(non.child.val[[p.s]] == 1))
    nzm = length(which(non.child.val[[p.s]] == -1))
    nzz = length(which(non.child.val[[p.s]] == 0))

    pvalq = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Quaternary')
    pvalt = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Ternary')
    pvale = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Enrichment')
    ##pvalf = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Fisher')


    qPlus  = npp + npm + npz
    qMinus = nmp + nmm + nmz
    qR     = nrp + nrm + nrz
    qZero  = nzp + nzm + nzz

    results[(2*(p.s-1)+1), 4]  = npp + nmm
    results[(2*(p.s-1)+1), 5]  = npm + nmp
    results[(2*(p.s-1)+1), 6]  = npp + nmm - (npm + nmp)
    results[(2*(p.s-1)+1), 7]  = qPlus + qMinus + qR
    results[(2*(p.s-1)+1), 8]  = npp + npm + nmp + nmm + nrp + nrm
    results[(2*(p.s-1)+1), 9]  = qR
    results[(2*(p.s-1)+1), 10] = nrp + nrm
    results[(2*(p.s-1)+1), 11] = qZero
    results[(2*(p.s-1)+1), 12] = pvalq$pval.up
    results[(2*(p.s-1)+1), 13] = pvalt$pval.up
    results[(2*(p.s-1)+1), 14] = pvale$pval.up

    results[(2*p.s), 4]  = nmp + npm
    results[(2*p.s), 5]  = npp + nmm
    results[(2*p.s), 6]  = nmp + npm - (npp + nmm)
    results[(2*p.s), 7]  = qPlus + qMinus + qR
    results[(2*p.s), 8]  = npp + npm + nmp + nmm + nrp + nrm
    results[(2*p.s), 9]  = qR
    results[(2*p.s), 10] = nrp + nrm
    results[(2*p.s), 11] = qZero
    results[(2*p.s), 12] = pvalq$pval.down
    results[(2*p.s), 13] = pvalt$pval.down
    results[(2*p.s), 14] = pvale$pval.down
  }


  results[,15] = fdrtool(results[,12], "pvalue", F, F, F, "fndr")$q
  results[,16] = fdrtool(results[,13], "pvalue", F, F, F, "fndr")$q
  results[,17] = fdrtool(results[,14], "pvalue", F, F, F, "fndr")$q

  #results = results[order(as.numeric(results$quaternary.pval)), ]
  #rownames(results) = 1:nrow(results)

  results
}

## This function runs the CRE based on the version specified
runCRE = function(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = ''){

  if(method == 'Quaternary'){
    qR     = nrp + nrm + nrz
    qZero  = nzp + nzm + nzz
    nPlus  = npp + nmp + nrp + nzp
    nMinus = npm + nmm + nrm + nzm
    nZero  = npz + nmz + nrz + nzz

    ## Assume up-regulated
    qPlus  = npp + npm + npz
    qMinus = nmp + nmm + nmz
    score  = npp + nmm + nrp + nrm - (npm + nmp)
    pval.up   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                          q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)

    ## Assume down-regulated
    qPlus  = nmp + nmm + nmz
    qMinus = npp + npm + npz
    score  = nmp + npm + nrp + nrm - (npp + nmm)
    pval.down   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                            q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
  }else if(method == 'Ternary'){
    qR     = 0
    qZero  = nzp + nzm + nzz
    nPlus  = npp + nmp + nzp
    nMinus = npm + nmm + nzm
    nZero  = npz + nmz + nzz

    ## Assume up-regulated
    qPlus  = npp + npm + npz
    qMinus = nmp + nmm + nmz
    score  = npp + nmm - (npm + nmp)
    pval.up   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                          q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)

    ## Assume down-regulated
    qPlus  = nmp + nmm + nmz
    qMinus = npp + npm + npz
    score  = nmp + npm - (npp + nmm)
    pval.down   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                            q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)

  }else if(method == 'Enrichment'){
    nrp    = npp + nmp + nrp
    nrm    = npm + nmm + nrm
    nrz    = npz + nmz + nrz

    qPlus  = 0
    qMinus = 0
    qR     = nrp + nrm + nrz
    qZero  = nzp + nzm + nzz

    nPlus  = nrp + nzp
    nMinus = nrm + nzm
    nZero  = nrz + nzz

    score  = nrp + nrm

    pval.up   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                          q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    pval.down = pval.up
  }else if(method == 'Fisher'){
    nrp    = npp + nmp + nrp
    nrm    = npm + nmm + nrm
    nrz    = npz + nmz + nrz

    M = matrix(0, nrow = 2, ncol = 2)
    M[1,1] = nrp + nrm
    M[1,2] = nrz
    M[2,1] = nzp + nzm
    M[2,2] = nzz

    pval.up = fisher.test(M, alternative = 'greater')$p.val
    pval.down = pval.up

  }

  pval = list(pval.up = pval.up, pval.down = pval.down)
}

