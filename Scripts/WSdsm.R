# WSdsm.R  (WordSpace / Distributional Semantic Model)
#
# Script constitué par un ensemble de fonctions destinées à faciliter l'usage 
# de la bibliothèque R 'wordspace' (Stefan Evert), à partir d'un corpus enregistré sous CWB.

# Un premier groupe de fonctions est destiné à créer un DSM et à calculer, à partir de ce DSM,
# les champs des lemmes choisis, avec visualisation par analyse factorielle des correspondances (AFC),
# et à extraire les mots-clés d'un ensemble à partir des valences lexicales généralisées pondérées.

# Le second groupe permet d'appliquer les mêmes procédures sur un corpus
# découpé en tranches : l'objectif est l'analyse de l'évolution d'un champ sémantique,
# la visualisation est conçue pour faire ressortir les éléments liés plus particulièrement
# à telle ou telle période. Les deux groupes doivent être employés de manière complémentaire.

# version pré-alpha 0.3    AG novembre 2015 - mars 2017. GPL3

# TODO : autres méthodes d'examen des évolutions.

#########################################################################################
# premier groupe : analyses globales > champs sémantiques
#########################################################################################

corpus2scan <- function(corp, dis=3, posA="QLF|SUB|VBE", posB="QLF|SUB|VBE", objetA= "lemma", objetB = "lemma", D=0, F="", attr="", val="", destination , flag=TRUE ) {

# Création d'un fichier-somme du scan complet d'un corpus
# ou d'une partie de corpus,
# résultant de l'application de 'cwb-scan-corpus' à une fenêtre
# de la largeur choisie (de part et d'autre du pivot).
#
# Double contrainte : taille de mémoire et temps d'exécution.
# le programme scanne 2 colonnes et décompte toutes les paires identiques ;
# on prend les colonnes successivement pour balayer toute la fenêtre choisie
# et on enregistre au fur et à mesure sur le DD ;
# après quoi, on récupère les fichiers un par un et on les concatène.
# Les affichages pendant l'exécution sont très approximatifs, il s'agit seulement 
# de faire patienter !

if (flag==TRUE){
  t1 <- Sys.time()
}

if (destination == "") {
  stop(" Indiquer une destination pour le scan ", call.=FALSE)
}
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(scipen=999)  # supprimer la notation scientifique (pb avec cqp)

efftt <- size(corpus(corp))
effpart <- 0
if (F==""){
  F <- efftt
}
if (D!=0 | F!="") {
  effpart <- F-D
}
if (attr!="") {
  def.scorp <- paste('[lemma=".*" %cd]', "::match.", attr, "=\"", val, "\"", sep="")
  CRP <- corpus(corp)
  crp <- subcorpus(CRP, def.scorp)
  effpart <- size(crp)
}

# boucle : scans par colonne
for (i in 0:(dis*2)) {
  if (i==dis){
    next()
    }
  # création des paramètres pour la ligne de commande / paramètre -b excessif ??
  params <- paste("-b 200000000 -q -s ",D,  sep="")
  if (F != efftt){
    params <- paste(params, " -e ",F, sep="")
  }
#  if (reg != ""){
#    params <- paste(params, " -r '",reg,"' ", sep="")
#  }
  params <- paste(params, "  ",corp, " ", objetA,"+",dis," '?pos+",dis,"=/", posA, "/' pos+",dis," ",objetB,"+",i," '?pos+",i,"=/", posB, "/' pos+",i,  sep="")
  if (attr != "" & val != ""){
    params <- paste(params,"  '?", attr, "=/", val,"/'", sep="")
  }
  
  sortie <- paste("/tmp/xyzxyz",i,".tsv", sep="")
  # exécution (sortie sur disque automatique)
  system2(command="cwb-scan-corpus", args=params, stdout=sortie) 
  cat("scan =",i, "sur", dis*2, "\n")
  gc()
}

# rassemblement en un seul fichier (sur disque)
commd <- paste("cat /tmp/xyzxyz* > ", destination, sep="")
system(command=commd)
commd2 <- paste("rm /tmp/xyzxyz*") # nettoyage des fichiers provisoires
system(command=commd2)

# création et enregistrement d'un fichier d'infos sur le scan
destination2 <- paste(destination, "_params", sep="")
parametres <- c("corpus","eff.total","eff.actuel","distance","posA","posB","objetA","objetB","D","F","attr","val")
valeurs <- c(corp,efftt,effpart,dis,posA,posB,objetA,objetB,D,F,attr,val)
infos <- cbind(parametres,valeurs)
write.table(infos,file=destination2, quote=FALSE,sep="\t",row.names=FALSE)

if (flag==TRUE) {
 t2 <- Sys.time()
 td <- difftime(t1,t2)
 cat("\n","Temps écoulé :", round(as.numeric(td),2), units(td), "\n")
}
}

#################################################################################

scan2dsm <- function(scan, seuil= 9, coef="simple-ll", transf="log", nproj="", flag=TRUE) {

# récupération sur le DD d'un fichier issu de corpus2scan()
# + paramètres
# nettoyage, scoring, création d'un objet WS exploitable

library(wordspace, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)

if (flag==TRUE) {
  t1 <- Sys.time()
}
gc()

scanp <- paste(scan, "_params", sep="")
params.tripl <- read.table(scanp, header=FALSE, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)

tripl <- read.table(scan, header=FALSE, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)
tripl <- tripl[, c(2,4,1)]
names(tripl) <- c("target", "feature", "eff") # esthétique !
# création de l'objet
triplobj <- dsm(target=tripl$target, feature=tripl$feature, score=tripl$eff, N=as.numeric(params.tripl[4,2]), raw.freq=TRUE, sort=TRUE)
rm(tripl) # nettoyage
gc()
# élagage
triplobj <- subset(triplobj, nnzero > seuil, nnzero > seuil, recursive=TRUE)
# scoring (filtrage des cooccurrents significatifs)
triplobjS <- dsm.score(triplobj, score= coef, transform=transf, normalize=TRUE)
# réduction des dimensions de la matrice
if (nproj != "") {
  triplobjS <- dsm.projection(triplobjS, method="rsvd", n=nproj, oversampling=4)
}

# enregistrement des infos ($globals) > dsm documenté !
triplobjS$globals$corpus <- params.tripl[2,2]
triplobjS$globals$nblignes <- length(triplobjS$rows$term)
triplobjS$globals$nbcols <- length(triplobjS$cols$term)
triplobjS$globals$posA <- params.tripl[6,2]
triplobjS$globals$posB <- params.tripl[7,2]
triplobjS$globals$objetA <- params.tripl[8,2]
triplobjS$globals$objetB <- params.tripl[9,2]
triplobjS$globals$dis <- params.tripl[5,2]
triplobjS$globals$effactuel <- as.numeric(params.tripl[4,2])
triplobjS$globals$D <- as.numeric(params.tripl[10,2])
triplobjS$globals$F <- as.numeric(params.tripl[11,2])
if (triplobjS$globals$F==Inf) triplobjS$globals$F <- triplobjS$globals$N-1
triplobjS$globals$attr <- params.tripl[12,2]
triplobjS$globals$val <- params.tripl[13,2]
triplobjS$globals$effcorpus <- params.tripl[3,2]
triplobjS$globals$seuil <- seuil
triplobjS$globals$coef <- coef
triplobjS$globals$transf <- transf
triplobjS$globals$nproj <- nproj

if (flag==T) {
  t2 <- Sys.time()
  td <- difftime(t1,t2)
  cat("\n","Temps écoulé :", round(as.numeric(td),2), units(td), "\n\n")
}
cat("lignes   : ", length(triplobjS$rows$term), "\n")
cat("colonnes : ", length(triplobjS$cols$term), "\n")
return(triplobjS)
}

###############################################################################

corpus2dsm <- function(corp, dis=5, posA="QLF|SUB|VBE", posB="QLF|SUB|VBE", objetA= "lemma", objetB = "lemma",D=0, F="", attr="", val="", destination, seuil= 9, coef="simple-ll",transf="log", nproj=""){

# regroupement de l'ensemble des opérations
# on part d'un corpus, d'une largeur de fenêtre
# et d'un choix des POS (pivot et cooc) ;
# on obtient un objet DSM prêt à l'emploi.

t1 <- Sys.time()
options(warn=-1)

if (destination == "") {
  stop(" Indiquer une destination pour le scan ", call.=FALSE)
}
# 1. scan
corpus2scan(corp=corp, dis=dis, posA=posA, posB=posB, objetA=objetA, objetB = objetB ,D=D, F=F, attr=attr, val=val, destination=destination, flag=FALSE)
cat("\n","Traitements...","\n")
gc() # nettoyage
# 2. construction de l'objet
res <- scan2dsm(scan=destination, seuil=seuil, coef=coef,transf=transf, nproj=nproj, flag=FALSE)

t2 <- Sys.time()
td <- difftime(t1,t2)
  cat("\n","Temps écoulé :", round(as.numeric(td),2), units(td), "\n")

res
}

################################################################################

dsm2af <- function(dsm, lm, nppv=40, cex=.9, decal=TRUE) {

# Recherche, dans un dsm donné, des p.p.voisins d'un lemme,
# et représentation par AFC de la matrice des distances.
# Le graphique fournit quelque chose d'analogue au Wortfeld
# au sens de Jost Trier. Les points sont répartis selon
# leurs distances réciproques : les divers 'nuages' correspondent
# aux sous-ensembles du champ.
# Création d'un objet contenant tous les éléments intermédiaires.

opar <- par(mar=par("mar"))
on.exit(par(opar))
library(wordspace, quietly=TRUE, warn.conflicts=FALSE)
library(ade4, quietly=TRUE, warn.conflicts=FALSE)
library(circular, quietly=TRUE, warn.conflicts=FALSE)
#library(MASS, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1) 

t1 <- Sys.time()

# recherche des p.p.voisins
vec.ppvoisins <- nearest.neighbours(M=dsm, term=lm, n=nppv)
ppv.names <- names(vec.ppvoisins)
val.ppvoisins <- cbind(as.character(ppv.names), as.numeric(vec.ppvoisins))
row.names(val.ppvoisins) <- NULL
mat.ppvoisins <- nearest.neighbours(M=dsm, term=lm, n=nppv, skip.missing=TRUE, dist.matrix=TRUE)

res <- list(NULL)
res[[1]] <- mat.ppvoisins
res[[2]] <- val.ppvoisins

# AFC sur le tableau des colonnes (la matrice est symétrique : on utilise les noms de ligne)
af.mat.ppvoisins <- dudi.coa(mat.ppvoisins, scannf=FALSE)
af.util <- af.mat.ppvoisins$co

# éviter les recouvrements d'étiquettes (appel à la fonction lisible())
if (decal==TRUE){
  ymax <- max(af.util[,2])
  ymin <- min(af.util[,2])
  Tbon <- lisible(af.util[,1],af.util[,2],lab=row.names(af.util),mn=ymin, mx=ymax,cex=(cex+.1))
  af.util[,1] <- Tbon[,1]
  af.util[,2] <- Tbon[,2]
}

res[[3]] <- af.util
names(res) <- c("matrice_distances", "vecteur_ppvoisins", "coordonnees")

# affichage de l'AF
par(mar=c(0.5,0.5,1.7,0.5))
plot(af.util, type="n", asp=1, axes=FALSE, frame.plot=TRUE)
text(af.util[1,], labels=row.names(af.util[1,]), cex=(cex+.2), col="red", font=2)
af.util <- af.util[-1,]
text(af.util, labels=row.names(af.util), cex=cex, col="blue")
# affichage d'un titre
nbr <- length(dsm$rows[,1])
nbc <- length(dsm$cols[,1])
nm.obj <- deparse(substitute(dsm))
mn <- paste("DSM d'origine : ",nm.obj,"  (matrice de ", nbc , " sur ", nbr ,"). ",nppv, " éléments.", sep = "")
title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
titranal <- paste("STRUCTURE GLOBALE DU CHAMP SÉMANTIQUE  de  *",lm,"*", sep="")
mtext(titranal, 3, line=0,cex=.8, font=1, adj=0)

#write.matrix(val.ppvoisins)

for (i in 1:nppv){
  cat(names(res$vecteur_ppvoisins)[i], "\n")
}

class(res) <- "NPPV"

t2 <- Sys.time()
td <- difftime(t1,t2)
cat("\n","Temps écoulé :", round(as.numeric(td),2), units(td), "\n")

res
}

##################################################################################

dsm2carte <- function(dsm, seuil= "", mincoo=3, stopw="",nseg=50, decal=TRUE, cex=.8) {

# 1. calcul des mots-clés enn fonction de la valence lexicale pondérée
# 2. représentation factorielle de l'ensemble (ACP sur indices de cooccurrence)
#
# calcul d'une liste de lemmes, évaluée à partir de 2 paramètres :
# seuil = nb minimal de cooc dans chaque case du tableau (calcul pour chaque ligne du nombre de cases > seuil)
# mincoo = nb minimal de cases > 0 dans chaque ligne (tri des lignes en fonction du nbe de cases retenues)
# stopw = fichier de mots-outils ou assimilés, un mot par ligne

t1 <- Sys.time()

gc()
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
library(ade4, quietly=TRUE, warn.conflicts=FALSE)
library(circular, quietly=TRUE, warn.conflicts=FALSE)
library(wordspace, quietly=TRUE, warn.conflicts=FALSE)
library(MASS, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
options(scipen=999)  # supprimer la notation scientifique (pb avec cqp)
if (!inherits(dsm, "dsm"))   stop("en entrée : un objet de classe dsm")

corp <- dsm$globals$corpus # nom du corpus
attr <- dsm$globals$attr
val <- dsm$globals$val
D <- dsm$globals$D
F <- dsm$globals$F

dsmm <- dsm$M              # matrice des effectifs de coocs bruts
cat("cooc.freq.max = ", max(dsmm), "\n")

# calculs (= tris en fonction des paramètres choisis)
S1 <- apply(dsmm, 1, function(x) length(x[x>seuil])) # nbe par ligne de cases > seuil
cat("nb.somme.cooc > seuil = ",length(S1[S1>0]),"\n")
S2 <- S1[S1>mincoo]                                  # tri des lignes à somme > mincoo
rm(S1)
gc()
S2 <- as.data.frame(cbind(names(S2),S2, stringsAsFactors=FALSE))
names(S2) <- c("names", "valbr")
S2$valbr <- as.numeric(as.character(S2$valbr))
S2$names <- as.character(S2$names)

# application d'une pondération aux valeurs brutes (par les fréquences totales)
# on calcule ces fréquences dans l'ensemble considéré, corpus ou sous-corpus
CRP <- corpus(corp)
if (dsm$globals$attr=="" & dsm$globals$N==dsm$globals$effactuel) {
  crp <- subcorpus(CRP, '[lemma=".*" & (pos="SUB"|pos="VBE"|pos="QLF")]')
}
else {
  def.scorp <- paste("abc:[lemma=\".*\" & (pos=\"SUB\"|pos=\"VBE\"|pos=\"QLF\") & _.", attr, "=\"", val, "\"]::abc >=",D," & abc <=",F, sep="")
  crp <- subcorpus(CRP, def.scorp)
}

Clist <- cqp_flist(crp, "match", "lemma")
Clist2 <- Clist[1:length(Clist)]
rm(Clist)
gc()
Clist <- as.data.frame(cbind(names(Clist2),Clist2, stringsAsFactors=FALSE))
names(Clist) <- c("names", "freqtt")
Clist$freqtt <- as.numeric(as.character(Clist$freqtt))
Clist$names <- as.character(Clist$names)
S3 <- merge(S2, Clist, by.x="names", by.y="names", all.x=TRUE, all.y=FALSE)
rm(S2)
gc()
S3$valence <- S3$valbr / S3$freqtt
S3 <- S3[,c(1,2,4,6)]
S3 <- S3[rev(order(S3[,4])),]
gc()

# tri des stop-words (par défaut : QLF et VBE latins inutiles)
if (stopw==""){
stopw <- c("--","ago","aio","alius","audio","debeo","dico1","dico2","facio","fio","habeo","inquio","ipse1","loquor","meus","multus","nihil","nolo","noster","nullus","omnis","pono","possum","quidam","sequor","sum","suus","talis","tantus","totus","tuus","uenio","uester","uolo2","hic2","hic1","iste","ille","diuersus","inquantus","alter","ceterus","quisque","ullus") 
}
else {
  stopw2 <- read.csv2(stopw, header=FALSE, stringsAsFactors=FALSE, quote="", fill=TRUE)
  stopw <- stopw2[,1]
}
lsttri <- setdiff(S3[,1],stopw)
S3 <- S3[(S3[,1]%in%lsttri),]
dsm2 <- subset(dsm, subset=(term %in% S3[,1]), select=(term %in% S3[,1]))
dsm3 <- as.matrix(dsm2$M)
dsm4 <- as.matrix(dsm2$S)

cat("nb.kw =  ", ncol(dsm3), "\n\n", sep="")
if (ncol(dsm3) < 20) {
  cat("Moins de 20 lemmes retenus : baissez les paramètres !", "\n\n", sep="")
  res <- list(S3,dsm3,dsm4)
  class(res) <- "carte"
  names(res) <- c("valences","mat.brute","mat.coeff")
  write.matrix(S3)
  t2 <- Sys.time()
  td <- difftime(t1,t2)
  cat("\n","Temps écoulé :", round(as.numeric(td),2), units(td), "\n")
  return(res)
}
if (ncol(dsm3) > 250) {
  cat("Plus de 250 lemmes retenus : augmentez les paramètres !", "\n\n", sep="")
  res <- list(S3,dsm3,dsm4)
  class(res) <- "carte"
  names(res) <- c("valences","mat.brute","mat.coeff")
  write.matrix(S3)
  t2 <- Sys.time()
  td <- difftime(t1,t2)
  cat("\n","Temps écoulé :", round(as.numeric(td),2), units(td), "\n")
  return(res)
}

# ACP sur le tableau des colonnes (la matrice est (à peu près) symétrique)
af.mat.coeff <- dudi.pca(dsm4, scannf=FALSE, nf=2) # on peut utiliser la transposée...
af.util <- af.mat.coeff$co

# éviter les recouvrements d'étiquettes (appel à la fonction lisible())
if (decal==TRUE){
  ymax <- max(af.util[,2])
  ymin <- min(af.util[,2])
  Tbon <- lisible(af.util[,1],af.util[,2],lab=row.names(af.util),mn=ymin, mx=ymax,cex=(cex+.1))
  af.util[,1] <- Tbon[,1]
  af.util[,2] <- Tbon[,2]
}
rm(Tbon)
gc()

# calcul des distances les plus importantes
distab <- data.frame(NULL)
cpt <- 1
nbcol <- ncol(dsm4)
for (i in 1:(nbcol-1)){
  for (j in (i+1):nbcol){
    distab[cpt,1] <- rownames(dsm4)[i]
    distab[cpt,2] <- colnames(dsm4)[j]
    distab[cpt,3] <- dsm4[i,j]
    cpt <- cpt+1
  }
}
distab.tr <- distab[order(distab[,3],decreasing=TRUE),]
distab <- distab.tr[1:nseg,]
# coordonnées d'affichage
R1r <- match(distab[,1], rownames(af.util))
R2r <- match(distab[,2], rownames(af.util))
distab[,4] <- af.util[R1r, 1]
distab[,5] <- af.util[R1r, 2]
distab[,6] <- af.util[R2r, 1]
distab[,7] <- af.util[R2r, 2]

# affichage de l'AF
par(mar=c(0.5,0.5,1.7,0.5))
plot(af.util, type="n", asp=1, axes=FALSE, frame.plot=TRUE)
nb.kw <- nrow(af.util)
text(af.util, labels=row.names(af.util), cex=cex, col="#005500")
segments(distab[,4], distab[,5], distab[,6], distab[,7], lwd=1, col="grey")

# affichage d'un titre
nbr <- length(dsm$rows[,1])
nbc <- length(dsm$cols[,1])
nm.obj <- deparse(substitute(dsm))
mn <- paste("DSM d'origine : ",nm.obj,"  (matrice de ", nbc , " sur ", nbr ,")    effectif : ",dsm$globals$N, " tokens", sep = "")
title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
if (dsm$globals$effactuel==dsm$globals$effcorpus & dsm$globals$attr=="") {
  titranal <- paste("CARTE SÉMANTIQUE DU CORPUS  *",corp, "*     ",nb.kw, " mots-clés", sep="")
}
else if (dsm$globals$attr!="") {
  titranal <- paste("CARTE DU SOUS-CORPUS  *",corp,"* attribut = ",dsm$globals$attr,"  valeur = ",dsm$globals$val,"  ",nb.kw, " mots-clés", sep="")
}
else {
  titranal <- paste("CARTE DU SOUS-CORPUS  *",corp,"* D = ",dsm$globals$D," F = ",dsm$globals$F,"  ",nb.kw, " mots-clés", sep="")
} 
mtext(titranal, 3, line=0,cex=.8, font=1, adj=0)

# création d'une liste en sortie
res <- list(S3,dsm3,dsm4,af.util,distab)
class(res) <- "carte"
names(res) <- c("valences","mat.brute","mat.coeff","acp","distab")
write.matrix(S3)

t2 <- Sys.time()
td <- difftime(t1,t2)
  cat("\n","Temps écoulé :", round(as.numeric(td),2), units(td), "\n")

return(res)

}


################################################################################
#  Second groupe : analyses par tranches > évolutions
################################################################################


corpus2scanm <- function(corp, dis=5, posA="QLF|SUB|VBE", posB="QLF|SUB|VBE", objetA="lemma", objetB ="lemma", attr="", val="", destination, trnch=5, flag=TRUE){

# construction d'une série de scans
# correspondant aux tranches successives d'un corpus
# enregistrés sur DD

if (flag==TRUE) {
  t1 <- Sys.time()
}
if (destination == "") {
  stop(" Indiquer une destination pour le scan ", call.=FALSE)
}
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)

# découpage en tranches égales
efftt <- size(corpus(corp))
bornes <- seq(1, efftt, length.out=trnch+1)

# boucle : scans des tranches
for (i in 1:trnch) {
  D <- bornes[i]
  F <- bornes[i+1]
  destin <- paste(destination, "_", i, sep="")
  corpus2scan(corp=corp, dis=dis, posA=posA, posB=posB, objetA=objetA, objetB=objetB, D=D, F=F, attr=attr, val=val, destination=destin, flag=FALSE)
  cat("Tranche ",i, " sur ",trnch, "terminée","\n\n")
}

if (flag==TRUE) {
t2 <- Sys.time()
td <- difftime(t1,t2)
cat("\n","Temps écoulé :", round(as.numeric(td),2), units(td), "\n")
}
}

###############################################################################

scanm2dsmm <- function(scan, seuil= 5, coef="simple-ll", nproj="", trnch, flag=TRUE) {

# récupération d'une série de scans ;
# constrution d'une série correspondante de DSM regroupés dans un objet list.

# boucle : récup des scans > liste de DSM
res <- list(NULL)
for (i in 1:trnch) {
  scanm <- paste(scan, "_", i, sep="")
  res[[i]] <- scan2dsm(scan=scanm, seuil=seuil, coef=coef, nproj=nproj, flag=FALSE)
}

res
}

#############################################################################

corpus2dsmm <- function(corp, dis=5, posA="QLF|SUB|VBE", posB="QLF|SUB|VBE", objetA= "lemma", objetB = "lemma", trnch=5, attr="", val="", destination, seuil= 5, coef="simple-ll", nproj=""){

# regroupement de deux scripts
# permettant d'effectuer à la suite un scan par tranches
# et la création d'un dsm multiple correspondant.

t1 <- Sys.time()
if (destination == "") {
  stop(" Indiquer une destination pour le scan ", call.=FALSE)
}

# 1. scans
corpus2scanm(corp=corp, dis=dis, posA=posA, posB=posB, objetA=objetA, objetB = objetB , trnch= trnch, attr=attr, val=val, destination=destination, flag=FALSE)
cat("\n","Traitements...","\n")
gc() # nettoyage
# 2. construction d'une série d'objets
res <- scanm2dsmm(scan=destination, seuil=seuil, coef=coef, nproj=nproj, flag=FALSE, trnch=trnch)

t2 <- Sys.time()
td <- difftime(t1,t2)
cat("\n","Temps écoulé :", round(as.numeric(td),2), units(td), "\n")

res
}

#################################################################################

dsmm2af <- function(dsmm, lm, nppv, xax=1, yax=1, cex=.9, decal=TRUE) {

# Récupération des p.p.voisins d'un lemme dans une suite de dsmm ;
# construction d'une matrice des distances ;
# visualisation par AFC.
# Strictement complémentaire de dsm2af() : 
# cette visualisation est seulement destinée à éclaircir
# les évolutions, - pas les sous-ensembles.
# Difficulté : éliminer les outliers qui bloquent l'AF.
# TODO : afficher ces outliers en éléments supplémentaires.

t1 <- Sys.time()
library(wordspace, quietly=TRUE, warn.conflicts=FALSE)
library(ade4, quietly=TRUE, warn.conflicts=FALSE)
library(circular, quietly=TRUE, warn.conflicts=FALSE)
#library(MASS, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)

cooc.tt <- as.vector(NULL)
trnch <- length(dsmm)

# premier passage : récupérer les lemmes pour chaque tranche
for (i in 1:trnch) {
  cooc.raw   <- nearest.neighbours(dsmm[[i]], lm, n=nppv)
  cooc.nam <- names(cooc.raw)
  cooc.tt <- c(cooc.tt, cooc.nam)
}

# établir la liste complète (vecteur char)
cooc.tt <- unique(cooc.tt)
nb.cooc <- length(cooc.tt)
piv <- rep(lm, times=nb.cooc)

# rechercher les distances pour tous les lemmes dans toutes les tranches > matrice
distmat <- matrix(ncol = trnch, nrow = nb.cooc, 0)
vec.ncol <- as.vector(NULL)
for (i in 1:trnch) {
  distmat[,i] <- pair.distances(piv, cooc.tt, method="cosine", dsmm[[i]], convert=FALSE)
  vec.ncol <- c(vec.ncol, paste("PER_", i, sep=""))
}
colnames(distmat) <- vec.ncol # on met les noms
rownames(distmat) <- cooc.tt
# calcul de la variance de chaque ligne (lemmes)
diff0.var <- apply(distmat, 1, var, na.rm=TRUE)
distmat0 <- cbind(distmat, diff0.var*100)

res <- list(NULL)    # création d'une liste pour les sorties
res[[1]] <- distmat0  # matrice brute

# nettoyer les lignes incomplètes (contenant au moins une valeur Inf)
li.sum <- apply(distmat, 1, sum)
li.sum[is.finite(li.sum)] <- TRUE
li.sum[is.infinite(li.sum)] <- FALSE
li.sum <- as.logical(li.sum)
distmat <- distmat[li.sum,]

# nettoyage des lignes dont la variance dépasse 2 écarts-types
diff.var <- apply(distmat, 1, sd, na.rm=TRUE)
var.mean <- mean(diff.var, na.rm=TRUE)
var.sd <- sd(diff.var, na.rm=TRUE)
li.rm <- (diff.var < (var.mean + (2*var.sd)) & diff.var > (var.mean - (2*var.sd)))
li.rm[is.na(li.rm)] <- FALSE
distmat <- distmat[li.rm,]

# réorganisation (moyennes réciproques)
li.m <- rep(0, nrow(distmat))
for (j in 1:ncol(distmat)){
   for (i in 1:nrow(distmat)){
      li.m[i] <- li.m[i]+(distmat[i,j]*j)
   }
   }
li.m <- li.m / rowSums(distmat)
distmat <- distmat[rev(sort.list(li.m)),]
res[[2]] <- distmat # matrice nettoyée et réorganisée

dis.mean <- apply(distmat,1,mean)
dis.sd <- apply(distmat,1,sd)
dis.util <- cbind(row.names(distmat),dis.mean,dis.sd)
colnames(dis.util) <- c("lemmes", "coeff.moy.", "sd")
cat("\n")
write.matrix(dis.util)
res[[3]] <- dis.util[,1:2]

# lissage (avec lowess()) : indispensable !
nb <- nrow(distmat)
dist.colnames <- colnames(distmat)
dist.rownames <- rownames(distmat)
ls.coeff2 <- matrix(0, nrow=nb, ncol=trnch)
for (i in 1:nb) {
  ls.coeff2[i,] <- lowess(distmat[i,])$y
}
distmat <- ls.coeff2
colnames(distmat) <- dist.colnames
rownames(distmat) <- dist.rownames
distmat[distmat<0] <- 0

# calcul de la variance par ligne et par colonne
diff.var <- sort(apply(distmat, 1, var, na.rm=TRUE))
diff2.var <- apply(distmat, 2, var, na.rm=TRUE)
res[[4]] <- diff.var
res[[5]] <- diff2.var

# analyse factorielle (AFC)
af.distmat <- dudi.coa(distmat, scannf=FALSE)
af.util.co <- af.distmat$co
colnames(af.util.co) <- c("axe1", "axe2")
af.util.li <- af.distmat$li
colnames(af.util.li) <- c("axe1", "axe2")
af.util.tt <- rbind(af.util.co, af.util.li)
res[[6]] <- af.util.tt
names(res) <- c("matrice_brute", "matrice_nettoyee", "vecteur_ppvoisins", "variances_lignes", "variances_colonnes", "coordonnees")
co.nm <- colnames(distmat)
li.nm <- rownames(distmat)
tt.nm <- c(co.nm, li.nm)

# éviter les recouvrements d'étiquettes (appel à la fonction lisible())
if (decal==TRUE){
  ymax <- max(af.util.tt[,2])
  ymin <- min(af.util.tt[,2])
  Tbon <- lisible(af.util.tt[,1],af.util.tt[,2],lab=row.names(af.util.tt),mn=ymin, mx=ymax,cex=(cex+.1))
  af.util.tt[,1] <- Tbon[,1]
  af.util.tt[,2] <- Tbon[,2]
}
af.util.tt[,1] <- af.util.tt[,1]*xax # contrôle de l'orientation des axes
af.util.tt[,2] <- af.util.tt[,2]*yax
# distinguer les lignes et colonnes
af.util.co <- af.util.tt[(1:trnch),]
af.util.li <- af.util.tt[((trnch+1):(length(af.util.tt[,2]))),]

# affichage de l'AF
par(mar=c(0.5,0.5,1.7,0.5))
if (asp==1){
plot(af.util.tt, asp=1, type="n", axes=FALSE, frame.plot=TRUE) # cadre
}
else {
plot(af.util.tt, type="n", axes=FALSE, frame.plot=TRUE) # cadre
}

#lines(af.util.co, col="grey", lwd=3)                           # trace
text(af.util.co, labels=co.nm, cex=cex, col="red", font=2)     # points-colonnes
text(af.util.li, labels=li.nm, cex=cex, col="blue")            # points-lignes
nm.obj <- deparse(substitute(dsmm))
nbcoocr <- length(af.util.li[,1])
mn <- paste("DSM multiple d'origine : ",nm.obj,"  (", trnch," tranches).   Lemme : ", lm, ".   ",nppv, " > ", nbcoocr,  " éléments.", sep = "")
title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)    # titre
titranal= "ÉVOLUTION DU CHAMP SÉMANTIQUE (sémantique distributionnelle)"
mtext(titranal, 3, line=0,cex=.8, font=1, adj=0)
spllines(af.util.co[,1], af.util.co[,2], col="red") # trace

class(res) <- "NPPVM"
return(res)
}

##########################################################################
##########################################################################

#####################
# suppression des recouvrements
# partant du centre, on écarte les points qui provoquent recouvrement,
# toujours vers l'extérieur (selon le quadrant), alternativement horizontalement
# et verticalement, de manière à éviter la déformation du nuage,
# en pondérant l'alternance par la proximité angulaire avec l'axe 1 ou 2
# peut durer de quelques secondes à quelques minutes !!!
#####################

lisible <- function (x, y, lab, mn, mx, cex=.2){

#on constitue le tab(leau de )dep(art)
library(circular, quietly=TRUE, warn.conflicts=FALSE)
eps <- 0.0000000001
tabdep <- as.data.frame(cbind(x,y,lab))
names(tabdep) <- c("x","y","lab")
row.names(tabdep) <- seq(1,nrow(tabdep))
tabdep$x <- as.numeric(as.character(tabdep[,1]))
tabdep$y <- as.numeric(as.character(tabdep[,2]))
tabdep$lab <- as.character(tabdep$lab)
htlet <- (mx-mn)/(30/cex)
lglet <- htlet*.5
H <- lglet/2
indx <- as.numeric(row.names(tabdep))
d2 <- (tabdep$x^2)+(tabdep$y^2)
drt <- tabdep$x + (H*nchar(tabdep$lab))
gau <- tabdep$x - (H*nchar(tabdep$lab))
angl <- deg(atan(tabdep$y/tabdep$x))/.9

tabdep <- as.data.frame(cbind(tabdep,indx,d2,drt,gau,angl))
tt <- length(x)
tabfin <- tabpro <- tabdep

# problème : points aux mêmes coordonnées
tabpro <- tabpro[sort.list(tabpro$d2),]
for (i in 2:nrow(tabpro)) {
  if (signif(tabpro[i,5],8) == signif(tabpro[i-1,5],8)) {
    tabpro[i,1] <- tabpro[i,1] + (tabpro[i,1]/10000)
  }
}
tabpro$d2 <- (tabpro$x^2)+(tabpro$y^2)
rn <- (runif(tt*100))*100

for (i in 1:tt){

# on trie et on évacue la première ligne >> tableau final
  tabpro <- tabpro[sort.list(tabpro$d2),]
  cnt <- (tabpro[1,])
  tabfin[i,] <- cnt
  tabpro <- tabpro[-1,]

# il faut repousser tout ce qui peut recouvrir le point actif  (cnt)
# constitution du rub(an) formé de tous les points à écarter
  if (nrow(tabpro)==0) next
  cnt[1] <- as.numeric(as.character(cnt[1]))-(eps*sign(as.numeric(as.character(cnt[1]))))
  cnt[2] <- as.numeric(as.character(cnt[2]))-(eps*sign(as.numeric(as.character(cnt[2]))))
  ruban <- tabpro[(abs(as.numeric(tabpro$y)-as.numeric(as.character(cnt[2])))< htlet),]

  if (nrow(ruban) == 0) next
  rubg <- ruban[(ruban$x < as.numeric(as.character(cnt[1])) & ruban$drt > as.numeric(as.character(cnt[7]))),]
  rubd <- ruban[(ruban$x > as.numeric(as.character(cnt[1])) & ruban$gau < as.numeric(as.character(cnt[6]))),]
  rub <- rbind(rubg,rubd)
  rub <- unique(rub)
  if (nrow(rub) == 0) next
  n <- nrow(rub)
  r <- 1

# on écarte tous les points du rub(an) alternativement horizontalement et verticalement, vers l'extérieur 
# du quadrant en combinant la valeur de l'angle et un nombre aléatoire (!)

for (j in 1:n){
    if (rub[j,1]>0 & rub[j,2]>0 & rub[j,8]<rn[r]) tabpro[(tabpro[,4]==rub[j,4]),1] <- cnt[6]+(H*nchar(rub[j,3]))
    if (rub[j,1]>0 & rub[j,2]>0 & rub[j,8]>=rn[r]) tabpro[(tabpro[,4]==rub[j,4]),2] <- cnt[2]+(htlet)
    
    if (rub[j,1]>0 & rub[j,2]<0 & abs(rub[j,8])<rn[r]) tabpro[(tabpro[,4]==rub[j,4]),1] <- cnt[6]+(H*nchar(rub[j,3]))
    if (rub[j,1]>0 & rub[j,2]<0 & abs(rub[j,8])>=rn[r]) tabpro[(tabpro[,4]==rub[j,4]),2] <- cnt[2]-(htlet)
    
    if (rub[j,1]<0 & rub[j,2]<0 & rub[j,8]<rn[r]) tabpro[(tabpro[,4]==rub[j,4]),1] <- cnt[7]-(H*nchar(rub[j,3]))
    if (rub[j,1]<0 & rub[j,2]<0 & rub[j,8]>=rn[r]) tabpro[(tabpro[,4]==rub[j,4]),2] <- cnt[2]-(htlet)
    
    if (rub[j,1]<0 & rub[j,2]>0 & abs(rub[j,8])<rn[r]) tabpro[(tabpro[,4]==rub[j,4]),1] <- cnt[7]-(H*nchar(rub[j,3]))
    if (rub[j,1]<0 & rub[j,2]>0 & abs(rub[j,8])>=rn[r]) tabpro[(tabpro[,4]==rub[j,4]),2] <- cnt[2]+(htlet)
    r <- r+1
    }
    
# on recalcule la position relative de tous les points restants
# de manière à être sûr d'attaquer le bon point au tour suivant
  tabpro$d2 <- (tabpro$x^2) + (tabpro$y^2)
  tabpro$drt <- tabpro$x + (H*nchar(tabpro$lab))
  tabpro$gau <- tabpro$x - (H*nchar(tabpro$lab))
  }

# on remet le tableau final dans l'ordre des lignes au départ (indx)
tabfin <- tabfin[sort.list(tabfin$indx),]
tabfin[,3] <- lab

return(tabfin)
}

###################################################################################

listevaleurs <- function(corp, attr) {

# utilitaire de listage des valeurs
# d'un attribut, avec calcul de l'effectif

t1 <- Sys.time()
gc()

library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)

efftt <- size(corpus(corp))
requ <- paste(corp,".",attr, sep="")

# liste des ids de l'attribut
idsattr <- unique(cqi_cpos2struc(requ, 0:(efftt-1)))
nb.idsattr <- length(idsattr)

# pour chaque id, les cpos-bornes et le nom
df.val <- data.frame(NULL)
for (i in 1:nb.idsattr) {
  df.val[i,1] <- idsattr[i]
  bornes <- cqi_struc2cpos(requ, idsattr[i])
  df.val[i,2] <- cqi_struc2str(requ, idsattr[i])
  df.val[i,3] <- bornes[2]-bornes[1]+1
}
names(df.val) <- c("id","nom","effectif")

# cumul des effectifs par valeur
prov <- df.val
prov[,2] <- as.factor(prov[,2])
df.valsum <- tapply(prov[,3],prov[,2],sum)

res <- list(df.val,df.valsum)
names(res) <- c("df.val","df.valsum")
cat("effectif total du corpus  ", efftt, "\n\n")
print(as.matrix(df.valsum))

t2 <- Sys.time()
td <- difftime(t1,t2)
  cat("\n","Temps écoulé :", round(as.numeric(td),2), units(td), "\n")

return(res)

}

