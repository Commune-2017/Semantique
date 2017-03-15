# cooc.R

# script "boîte à outils" d'analyse des cooccurrents 
# écrit pour traiter des corpus de textes anciens,
# spécialement de textes latins
# lemmatisés avec les paramètres OMNIA,
# mais commodément adaptable.

# alain guerreau - août 2012-mars 2016 domaine public ou GPL3...
# version pré-alpha !!  0.4
# guerreau@msh-paris.fr

##########################################################
# premier groupe : évolutions
##########################################################

freqlem <- function(corp, reg="", lm, trnch=20, attr="", val="", date) {

# recherche d'un lemme et affichage
# d'un histogramme simple (par tranches de mêmes effectifs)
# à utiliser pour une première observation globale
# A GUERREAU gpl3 janvier 2013

if (reg != "") {
  reg0 <- Sys.getenv("CORPUS_REGISTRY")
  Sys.setenv(CORPUS_REGISTRY=reg)
}
library(MASS, quietly=TRUE, warn.conflicts=FALSE)
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
LL <- "Lsub"
corpdate <- paste(corp, date, sep=".")

# vérifications élémentaires
if ((attr!="" & val=="") | (attr=="" & val!="")) stop ("Requête incomplète !", call.=FALSE)

provl <- paste(corp, ".lemma", sep="")
provp <- paste(corp, ".pos", sep="")
idprov <- cqi_regex2id(provl, lm)

if (length(idprov)==0){
  fl <- 1
  lm2 <- paste(lm, ".*", sep="")
  idprov2 <- cqi_regex2id(provl, lm2)

if (length(idprov2) == 0) stop(paste(lm, "n'existe pas !"), call. = FALSE)

strprov <- cqi_id2str(provl, idprov2)
strcposprov <- NULL
posprov <- NULL
freq <- NULL
cat("\n")
for (i in 1:length(idprov2)){
  strcposprov[i] <- cqi_id2cpos(provl, idprov2[i])
  posprov[i] <- cqi_cpos2str(provp, strcposprov[i])
  freq[i] <- cqi_id2freq(provl, idprov2[i])
  cat(strprov[i], "/", posprov[i], "/", freq[i], "\n")
}
cat("\n")
if (fl > 0) stop ("Choisissez parmi ces lemmes !", call.=FALSE)
}

# cas simple : corpus entier
if (attr=="" & val=="") {
quer1 <- paste ('[lemma="', lm, '" %cd]', sep ="")
cqi_query(corp, LL, quer1)
sc <- paste(corp, LL, sep = ":")
A1d <- cqi_dump_subcorpus(sc)
ef <- nrow(A1d)
cat("\n")
cat(lm, "effectif total : ", ef)

position <- classrang(corp = corp, lm = lm)
cat (paste("\n","classe = ",position[1]," sur ",position[2],"\n", sep=""))

flist <- cqi_fdist1(sc, "match", "word") # formes correspondant au lemme
sc2 <- paste(corp, ".word", sep="")
flistb <- data.frame(cqi_id2str(sc2, flist[,1]), flist[,2])
cat("\n", "\n")
write.matrix(flistb)

# découpage et affichage de l'histogramme
ln <- size(corpus(corp))
xl <- paste(corp, " / ", trnch, "  (effectif total = ", ln, ")", sep = "")
mn <- paste("lemme : ", lm, " (effectif du lemme = ", ef, ")", sep = "")
br <- seq(1, ln, length.out=trnch+1)
decal <- (ln/trnch)/2
br2 <- as.integer((br[1:(length(br))-1]) + decal)
listdate <- cqi_cpos2struc(corpdate, br2)
listdateb <- cqi_struc2str(corpdate, listdate)
opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(4.2,4.2,2,0.5))
H <- hist(A1d[,1], breaks = br, border = "grey", main = "", xlab = "", ylab = "Effectif par tranche")
axis(side=1, at=br2, labels=listdateb, line=2)
title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
mtext(xl, 3, line=0,cex=.8, font=1, adj=0)
}

# cas complexe : seulement une partie du corpus
# nécessite de créer un sous-corpus provisoire complet correspondant
# au critère de choix, pour pouvoir créer des tranches égales
# dans ce sous-corpus
if (attr!="" & val!="") {
cat("Le calcul de l'effectif d'une fraction du corpus total prend du temps... \n")
LL2 <- paste(LL, "2", sep="")
quer2 <- paste('[lemma="',lm,'" %cd]', "::match.", attr, "=\"", val, "\"", sep="")
cqi_query(corp, LL2, quer2)
sc2 <- paste(corp, LL2, sep = ":")
lsc <- cqi_subcorpus_size(sc2)
if (lsc==0) stop("Pas de réponse (valeur de l'attribut inexistante ?), reformulez la requête !", call.=FALSE)
A1d2 <- cqi_dump_subcorpus(sc2)
ef <- nrow(A1d2)
indexL <- as.integer(seq(1, ef, 1))
efI <- cbind(as.integer(A1d2[,1]), indexL)
LL3 <- paste(LL, "3", sep="")
quertt <- paste("[lemma=\".*\"]", "::match.", attr, "=\"", val, "\"", sep="")
cqi_query(corp, LL3, quertt)
sc3 <- paste(corp, LL3, sep = ":")
A2d3 <- cqi_dump_subcorpus(sc3)
ln2 <- nrow(A2d3)

# ici, une ruse diabolique !
# (nécessaire pour contourner la discontinuité des cpos)
indexT <- as.integer(seq(1, ln2, 1)) # index auxiliaire
corpI <- cbind(as.integer(A2d3[,1]), indexT)
choix <- corpI[,1] %in% efI[,1]
corpI2 <- cbind(corpI, choix)
indbon <- corpI2[(corpI2[,3]==TRUE),2]

# découpage et affichage de l'histogramme
mn <- paste("lemme : ", lm, " (effectif du lemme = ", ef, ")", sep = "")
sb <- paste(attr, " = ", val, " (effectif = ", ln2, ")", sep="")
br <- seq(1, ln2, length.out=trnch+1)
decal <- (ln2/trnch)/2
br2 <- as.integer((br[1:(length(br))-1]) + decal)
br3 <- corpI2[br2,1]
listdate <- cqi_cpos2struc(corpdate, br3)
listdateb <- cqi_struc2str(corpdate, listdate)
opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(4.2,4.2,2,0.5))
H <- hist(indbon, breaks = br, border = "grey", main = "", xlab = "", ylab = "Effectif par tranche")
axis(side=1, at=br2, labels=listdateb, line=2)
title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
mtext(sb, 3, line=0,cex=.8, font=1, adj=0)
}

# calcul et affichage d'une courbe de tendance (fonction lowess())
mov <- lowess(H$counts)
br2 <- (br[1:(length(br))-1]) + decal
lines(br2, mov$y, col="blue", lwd=1.5)
EFTR <<- H$counts
LB1 <<- br2

if (reg != "") {
  Sys.setenv(CORPUS_REGISTRY=reg0)
}


}

#####################################################################################

freqlem2 <- function(corp, reg="", lm, trnch=20, attr="", val="", date, D=0, F=Inf) {

# recherche d'un lemme et affichage
# d'un histogramme simple (par tranches de mêmes effectifs)
# à utiliser pour une première observation globale
# affichage *chronologique* (anamorphose des tranches)
# A GUERREAU gpl3 janvier 2013

if (reg != "") {
  reg0 <- Sys.getenv("CORPUS_REGISTRY")
  Sys.setenv(CORPUS_REGISTRY=reg)
}

library(MASS, quietly=TRUE, warn.conflicts=FALSE)
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
LL <- "Lsub"

# vérifications
at <- cqi_attributes(corp, "s")
if ((date %in% at) == FALSE) stop ('Indiquez *explicitement* le nom du paramètre "date" !', call.=FALSE)
corpdate <- paste(corp, date, sep=".")
if ((attr!="" & val=="") | (attr=="" & val!="")) stop ("Requête incomplète !", call.=FALSE)
ln <- size(corpus(corp))
if (F=="Inf") {
  F <- ln-1
}
lgcorp <- F - D # effectif du corpus traité

provl <- paste(corp, ".lemma", sep="")
provp <- paste(corp, ".pos", sep="")
idprov <- cqi_regex2id(provl, lm)

if (length(idprov)==0){
  fl <- 1
  lm2 <- paste(lm, ".*", sep="")
  idprov2 <- cqi_regex2id(provl, lm2)

if (length(idprov2) == 0) stop(paste(lm, "n'existe pas !"), call. = FALSE)

strprov <- cqi_id2str(provl, idprov2)
strcposprov <- NULL
posprov <- NULL
freq <- NULL
cat("\n")
for (i in 1:length(idprov2)){
  strcposprov[i] <- cqi_id2cpos(provl, idprov2[i])
  posprov[i] <- cqi_cpos2str(provp, strcposprov[i])
  freq[i] <- cqi_id2freq(provl, idprov2[i])
  cat(strprov[i], "/", posprov[i], "/", freq[i], "\n")
}
cat("\n")
if (fl > 0) stop ("Choisissez parmi ces lemmes !", call.=FALSE)
}

# cas simple
if (attr=="" & val=="") {
quer1 <- paste ('[lemma="', lm, '" %cd]', sep ="")
cqi_query(corp, LL, quer1)
sc <- paste(corp, LL, sep = ":")
A1d <- cqi_dump_subcorpus(sc)

# choix d'un sous-corpus en fonction des bornes D / F
cposLm <- A1d[,1]
cposLm <- cposLm[cposLm > D & cposLm < F]
cposLm <- as.integer(cposLm)
ef <- length(cposLm) # effectif du lemme

cat("\n")
cat(lm, " : effectif total dans le (sous)corpus considéré = ", ef)

position <- classrang(corp = corp, lm = lm)
cat (paste("\n","classe = ",position[1]," sur ",position[2],"\n", sep=""))

flist <- cqi_fdist1(sc, "match", "word") # formes correspondant au lemme
sc2 <- paste(corp, ".word", sep="")
lid <- cqi_cpos2id(sc2, cposLm)
forme <- cqi_id2str(sc2, lid)
flistb <- data.frame(lid,forme)
lifreq <- as.data.frame(aggregate(1:nrow(flistb), flistb, length))
lifreq <- lifreq[rev(order(lifreq[,3])),]
cat("\n", "\n")
write.matrix(lifreq[,2:3])

xl <- paste(corp, " / ", trnch, "  (effectif total dans le (sous)corpus considéré  = ", lgcorp, ")", sep = "")
mn <- paste("lemme : ", lm, " (effectif du lemme = ", ef, ")", sep = "")

br <- seq(D, F, length.out=trnch+1)
decal <- (lgcorp/trnch)/2
br2 <- as.integer((br[1:(length(br))-1]) + decal)
listdate2 <- cqi_cpos2struc(corpdate, br2)
listdateb2 <- cqi_struc2str(corpdate, listdate2)

# il faut essayer d'élininer tout ce qui n'est pas daté
listdate <- cqi_cpos2struc(corpdate, br)
listdateb <- cqi_struc2str(corpdate, listdate)
listdateb <- as.numeric(listdateb)
listdateb <- listdateb[listdateb < 9000]
if (nchar(tail(listdateb,1))==3) {
	listdateb <- listdateb[listdateb < 999]
	}
trnch2 <- length(listdateb)-1
listdateb2 <- listdateb2[1:trnch2]

opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(4.2,4.2,2,0.5))
H <- hist(cposLm, breaks = br, plot = FALSE)
efftrn <- H$counts
efftrn <- efftrn[1:trnch2]
EFTR <<- efftrn
LB1 <<- listdateb
LB2 <<- listdateb2
eff2 <- c(efftrn,1)
plot(listdateb, eff2, type="n", xlab="chronologie", ylab="effectif par tranche")
for (i in 1:trnch2) {
#i <- 1:trnch2
  rect(listdateb[i], 0, listdateb[i+1], efftrn[i], border="grey")
  }

# axis(side=1, at=br2, labels=listdateb, line=2)
title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
mtext(xl, 3, line=0,cex=.8, font=1, adj=0)
}

# cas complexe
if (attr!="" & val!="") {
cat("Le calcul de l'effectif d'une fraction du corpus total prend du temps... \n")
LL2 <- paste(LL, "2", sep="")
quer2 <- paste('[lemma="',lm,'" %cd]', "::match.", attr, "=\"", val, "\"", sep="")
cqi_query(corp, LL2, quer2)
sc2 <- paste(corp, LL2, sep = ":")
lsc <- cqi_subcorpus_size(sc2)
if (lsc==0) stop("Pas de réponse (valeur de l'attribut inexistante ?), reformulez la requête !", call.=FALSE)
A1d2 <- cqi_dump_subcorpus(sc2)

# choix d'un sous-corpus en fonction des bornes D / F
cposLm <- A1d2[,1]
cposLm <- cposLm[cposLm > D & cposLm < F]
cposLm <- as.integer(cposLm)
ef <- length(cposLm) # effectif du lemme
cat("\n")
cat(lm, " : effectif total dans le (sous)corpus considéré = ", ef)
cat("\n")

indexL <- as.integer(seq(1, ef, 1))
efI <- cbind(cposLm, indexL)
LL3 <- paste(LL, "3", sep="")
quertt <- paste("[lemma=\".*\"]", "::match.", attr, "=\"", val, "\"", sep="")
cqi_query(corp, LL3, quertt)
sc3 <- paste(corp, LL3, sep = ":")
A2d3 <- cqi_dump_subcorpus(sc3)

# restriction  D / F
cposSc <- A2d3[,1]
cposSc <- cposSc[cposSc > D & cposSc < F]
cposSc <- as.integer(cposSc)
ln2 <- length(cposSc)

# ici, une ruse diabolique !
indexT <- as.integer(seq(1, ln2, 1))
corpI <- cbind(as.integer(cposSc), indexT)
choix <- corpI[,1] %in% efI[,1]
corpI2 <- cbind(corpI, choix)
indbon <- corpI2[(corpI2[,3]==TRUE),2]

flist <- cqi_fdist1(sc2, "match", "word")
sc4 <- paste(corp, ".word", sep="")
lid <- cqi_cpos2id(sc4, cposLm)
forme <- cqi_id2str(sc4, lid)
flistb <- data.frame(lid,forme)
lifreq <- as.data.frame(aggregate(1:nrow(flistb), flistb, length))
lifreq <- lifreq[rev(order(lifreq[,3])),]
cat("\n", "\n")
write.matrix(lifreq[,2:3])

mn <- paste("lemme : ", lm, " (effectif du lemme = ", ef, ")", sep = "")
sb <- paste(attr, " = ", val, " (effectif = ", ln2, ")", sep="")

br <- seq(1, ln2, length.out=trnch+1)
decal <- (ln2/trnch)/2
br2 <- as.integer((br[1:(length(br))-1]) + decal)
br3 <- corpI2[br2,1]
listdate <- cqi_cpos2struc(corpdate, br3)
listdateb <- cqi_struc2str(corpdate, listdate)
listdateb2 <- listdateb[listdateb < 9000]
if (nchar(tail(listdateb,1))==3) {
	listdateb <- listdateb[listdateb < 999]
	}
trnch2 <- length(listdateb)-1
listdateb2 <- listdateb2[1:trnch2]

opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(4.2,4.2,2,0.5))
H <- hist(indbon, breaks = br,plot=FALSE)
efftrn <- H$counts
efftrn <- efftrn[1:trnch2]
EFTR <<- efftrn
LB1 <<- listdateb
LB2 <<- listdateb2
eff2 <- c(efftrn,1)
plot(listdateb, eff2, type="n", xlab="chronologie", ylab="effectif par tranche")
for (i in 1:trnch2) {
  rect(listdateb[i], 0, listdateb[i+1], efftrn[i], border="grey")
}

axis(side=1, at=br2, labels=listdateb, line=2)
title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
mtext(sb, 3, line=0,cex=.8, font=1, adj=0)
}

# courbe de tendance
mov <- lowess(efftrn)
lines(listdateb2, mov$y, col="blue", lwd=1.5)

if (reg != "") {
  Sys.setenv(CORPUS_REGISTRY=reg0)
}

}

#####################################################################################

freqfrm <- function(corp, frm, trnch=20, attr="", val="") {

# recherche d'une forme et affichage
# d'un histogramme simple (par tranches de mêmes effectifs)
# à utiliser pour une première observation globale

library(MASS, quietly=TRUE, warn.conflicts=FALSE)
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
LL <- "Lsub"

# vérifications
if ((attr!="" & val=="") | (attr=="" & val!="")) stop ("Requête incomplète !", call.=FALSE)

provl <- paste(corp, ".word", sep="")
provp <- paste(corp, ".pos", sep="")
provlm <- paste(corp, ".lemma", sep="")
idprov <- cqi_regex2id(provl, frm)

if (length(idprov)==0){
  fl <- 1
  lm2 <- paste(frm, ".*", sep="")
  idprov2 <- cqi_regex2id(provl, lm2)

if (length(idprov2) == 0) stop(paste(frm, "n'existe pas !"), call. = FALSE)

strprov <- cqi_id2str(provl, idprov2)
strcposprov <- NULL
posprov <- NULL
freq <- NULL
cat("\n")
for (i in 1:length(idprov2)){
  strcposprov[i] <- cqi_id2cpos(provl, idprov2[i])
  posprov[i] <- cqi_cpos2str(provp, strcposprov[i])
  freq[i] <- cqi_id2freq(provl, idprov2[i])
  cat(strprov[i], "/", posprov[i], "/", freq[i], "\n")
}
cat("\n")
if (fl > 0) stop ("Choisissez parmi ces lemmes !", call.=FALSE)
}

# cas simple
if (attr=="" & val=="") {
quer1 <- paste ("[word='", frm, "'%cd]", sep ="")
cqi_query(corp, LL, quer1)
sc <- paste(corp, LL, sep = ":")
A1d <- cqi_dump_subcorpus(sc)
ef <- nrow(A1d)
cat("\n")
cat(frm, "effectif total : ", ef)

stopos <- cqi_cpos2str(provp, A1d)
stolm <- cqi_cpos2str(provlm, A1d)
stodf <- as.data.frame(cbind(stopos, stolm))
stodf[,1] <- as.character(stodf[,1])
stodf[,2] <- as.character(stodf[,2])
stob <- as.data.frame(aggregate(1:nrow(stodf), stodf, length))
cat("\n", unique(stob[,1]), unique(stob[,2]))
flist <- cqi_fdist1(sc, "match", "word")
sc2 <- paste(corp, ".word", sep="")
flistb <- data.frame(cqi_id2str(sc2, flist[,1]), flist[,2])
names(flistb) <- c("forme", "effectif")
cat("\n", "\n")
write.matrix(flistb)

ln <- size(corpus(corp))
xl <- paste(corp, " / ", trnch, "  (effectif total = ", ln, ")", sep = "")
mn <- paste("forme : ", frm, " (effectif de la forme = ", ef, ")", sep = "")
br <- seq(1, ln, length.out=trnch+1)
decal <- (ln/trnch)/2
H <- hist(A1d[,1], breaks = br, border = "grey", main = mn, xlab = xl, ylab = "Effectif par tranche")
}

# cas complexe
if (attr!="" & val!="") {
cat("Le calcul de l'effectif d'une fraction du corpus total prend du temps... \n")
LL2 <- paste(LL, "2", sep="")
quer2 <- paste("[word=\"",frm,"\"%cd]", "::match.", attr, "=\"", val, "\"", sep="")
cqi_query(corp, LL2, quer2)
sc2 <- paste(corp, LL2, sep = ":")
lsc <- cqi_subcorpus_size(sc2)
if (lsc==0) stop("Pas de réponse (valeur de l'attribut inexistante ?), reformulez la requête !", call.=FALSE)
A1d2 <- cqi_dump_subcorpus(sc2)
ef <- nrow(A1d2)
cat("\n")
cat(frm, "effectif total : ", ef)
stopos <- cqi_cpos2str(provp, A1d2)
stolm <- cqi_cpos2str(provlm, A1d2)
stodf <- as.data.frame(cbind(stopos, stolm))
stodf[,1] <- as.character(stodf[,1])
stodf[,2] <- as.character(stodf[,2])
stob <- as.data.frame(aggregate(1:nrow(stodf), stodf, length))
cat("\n", unique(stob[,1]), unique(stob[,2]))


flist <- cqi_fdist1(sc2, "match", "word")
sc3 <- paste(corp, ".word", sep="")
flistb <- data.frame(cqi_id2str(sc3, flist[,1]), flist[,2])
names(flistb) <- c("forme", "effectif")
cat("\n", "\n")
write.matrix(flistb)

cat("\n")
indexL <- as.integer(seq(1, ef, 1))
efI <- cbind(as.integer(A1d2[,1]), indexL)
LL3 <- paste(LL, "3", sep="")
quertt <- paste("[word=\".*\"]", "::match.", attr, "=\"", val, "\"", sep="")
cqi_query(corp, LL3, quertt)
sc3 <- paste(corp, LL3, sep = ":")
A2d3 <- cqi_dump_subcorpus(sc3)
ln2 <- nrow(A2d3)
# ici, une ruse diabolique !
indexT <- as.integer(seq(1, ln2, 1))
corpI <- cbind(as.integer(A2d3[,1]), indexT)
choix <- corpI[,1] %in% efI[,1]
corpI2 <- cbind(corpI, choix)
indbon <- corpI2[(corpI2[,3]==TRUE),2]
mn <- paste("forme : ", frm, " (effectif de la forme = ", ef, ")", sep = "")
sb <- paste("index du corpus cible (effectif total = ", ln2, ")", sep="")
br <- seq(1, ln2, length.out=trnch+1)
decal <- (ln2/trnch)/2
H <- hist(indbon, breaks = br, border = "grey", main = mn, xlab = sb, ylab = "Effectif par tranche")
}

# courbe de tendance
mov <- lowess(H$counts)
br2 <- (br[1:(length(br))-1]) + decal
lines(br2, mov$y, col="blue", lwd=1.5)
EFTR <<- H$counts
LB1 <<- br2
}

#####################################################################################

freqexp <- function(corp, exp, trnch=25, attr="", lim=200) {

# recherche d'une expression (CQL) et affichage
# d'un histogramme simple (par tranches de mêmes effectifs)
# enregistrement de la liste des occurrences
# affichage réduit à 'lim' occurrences max

library(MASS, quietly=TRUE, warn.conflicts=FALSE)
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
if (exp=="") {
  stop("pas d'expression !")
}

# calculs et affichage de l'histogramme
LL <- "EXP"
provl <- paste(corp, ".word", sep="")
provp <- paste(corp, ".pos", sep="")
provlm <- paste(corp, ".lemma", sep="")
obj3 <- paste(corp, attr, sep=".")
cqi_query(corp, LL, exp)
sc <- paste(corp, LL, sep = ":")
A1d <- cqi_dump_subcorpus(sc)
res <- A1d[,1]
ef <- nrow(A1d)
cat("\n")
cat(exp, "effectif total : ", ef)
cat("\n")

ln <- size(corpus(corp))
xl <- paste(corp, " / ", trnch, "  (effectif total = ", ln, ")", sep = "")
mn <- paste("expression : ", exp, " (effectif de l'expression = ", ef, ")", sep = "")
br <- seq(1, ln, length.out=trnch+1)
decal <- (ln/trnch)/2
H <- hist(A1d[,1], breaks = br, border = "grey", main = mn, xlab = xl, ylab = "Effectif par tranche")
mov <- lowess(H$counts)
br2 <- (br[1:(length(br))-1]) + decal
lines(br2, mov$y, col="blue", lwd=1.5)
EFTR <<- H$counts
LB1 <<- br2

# recherche et affichage de l'attribut choisi
if (attr != "") {
res2 <- cqi_cpos2struc(obj3, res)
res3 <- cqi_struc2str(obj3, res2)
}
if (attr == "") {
res3 <- rep("", length(res))
}

d <- res - 6
vec <- rep(0, 20*length(d))
mat <- matrix(vec, nrow=length(d), ncol=20)
for (i in 1:20) {
  mat[,i] <- d + i
  }

ench2 <- function(x, provl){
  a <- cqi_cpos2str(provl, x)
  b <- paste(a[1:20], sep="", collapse=" ")
  return(b)
}

tabw <- apply(mat, 1, ench2, provl)
res4 <- cbind(as.matrix(tabw), as.matrix(res3))
if (nrow(res4) <= lim) {
  write.matrix(res4, file="", sep=" ")
}
if (nrow(res4) > lim){
  tirage <- sort(as.integer(runif(lim,1,nrow(res4))))
  res5 <- res4[tirage,]
  write.matrix(res5, file="", sep=" ")
}
return(res4)
}

##########################################################

freqcooc <- function(corp, lm1, lm2, dis=5, coeff="dice", trnch=10, D=0, F=Inf, date, selec=1) {

# juillet 2014 Alain GUERREAU -> ||GPL3||
# calculs des effectifs de lm1 et lm2 et de la cooc,
# division en tranches égales ;
# graphique lm1 lm2 coeff_cooc
# + courbes de tendance
# fenêtre = dis*2

library(MASS, quietly=TRUE, warn.conflicts=FALSE)
library(plyr, quietly=TRUE, warn.conflicts=FALSE)
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
gcinfo(verbose=FALSE)

# AGENDA éliminer dates 9999 / tri par attr + val


# vérifications préalables

if (dis < 1) {
  stop("si dis=0, pas de cooccurrent !", call.=FALSE)
}
vcoeff <- c("dice","poiss","daille","pmi","hyperg","ms")
if (coeff %in% vcoeff == FALSE) stop ("coefficient inconnu, choisir : dice/poiss/daille/pmi/hyperg/ms", call.=FALSE)
  
provl <- paste(corp, ".lemma", sep="")
a <- cqi_str2id(provl, lm1)
if (a < 0) stop(paste(lm1, "n'existe pas !"), call.=FALSE)
b <- cqi_str2id(provl, lm2)
if (b < 0) stop(paste(lm2, "n'existe pas !"), call.=FALSE)
# if ((attr!="" & val=="") | (attr=="" & val!="")) stop ("Requête incomplète !", call.=FALSE)
file <- lm
LL1 <- "Lsub1"
LL2 <- "Lsub2"
LL3 <- "Lsub3"
corpdate <- paste(corp, date, sep=".")

N <- size(corpus(corp))
Ntr <- N / trnch
br <- seq.int(1, N, length.out=trnch+1)
decal <- (N/trnch)/2
br2 <- as.integer((br[1:(length(br))-1]) + decal)

listdate <- cqi_cpos2struc(corpdate, br2)
listdateb <- cqi_struc2str(corpdate, listdate)

# cas simple
#if (attr=="" & val=="") {

# calculs des effectifs observés
quer1 <- paste ('[lemma="', lm1, '" %cd]', sep ="")
quer2 <- paste ('[lemma="', lm2, '" %cd]', sep ="")
quer3 <- paste ('([lemma="', lm1, '" %cd][]{0,',dis,'}[lemma="', lm2, '" %cd])|([lemma="', lm2, '" %cd][]{0,',dis,'}[lemma="', lm1, '" %cd])', sep ="")
sc1 <- paste(corp, LL1, sep = ":")
sc2 <- paste(corp, LL2, sep = ":")
sc3 <- paste(corp, LL3, sep = ":")

cqi_query(corp, LL1, quer1)
dump1 <- cqi_dump_subcorpus(sc1)[,1]
cqi_query(corp, LL2, quer2)
dump2 <- cqi_dump_subcorpus(sc2)[,1]
cqi_query(corp, LL3, quer3)
subsz <- cqi_subcorpus_size(sc3)
if (subsz == 0) stop ("Aucune cooccurrence !", call.=FALSE)
dump3 <- cqi_dump_subcorpus(sc3)[,1]

ef1 <- length(dump1)
ef2 <- length(dump2)
ef3 <- length(dump3)
ef311 <- as.numeric(ef1) * as.numeric(ef2) / N

cat("\n")
cat(lm1, ", effectif total : ", ef1, "\n")
cat(lm2, ", effectif total : ", ef2, "\n")
cat("cooc ",lm1, "/",lm2,  ", effectif total : ", ef3, " (indépendance =", ef311,")", "\n\n")
mn <- paste("corpus : ", corp, " (effectif :", N, ") /", trnch,"    ",lm1, " = ", ef1,"  ", lm2, " = ", ef2, "  cooc = ",ef3,sep = "")
# sb <- paste(lm1, " = ", ef1,"    ", lm2, " = ", ef2, "    cooc = ",ef3, sep="")

R1 <- hist(dump1, breaks=br, plot=FALSE)$counts
C1 <- hist(dump2, breaks=br, plot=FALSE)$counts
O11 <- hist(dump3, breaks=br, plot=FALSE)$counts 
RCO <- cbind(R1,C1,O11)
E11 <- R1*C1/Ntr
O11mE11 <- O11-E11
mindiff <- min(O11mE11)
if (mindiff < 0) cat ("tranche(s) sans cooccurrence ! \n\n") 

# calculs des coefficients
c.dice <- (O11*1e+10)/(C1+R1)^selec
c.poiss <- O11*(log(O11)-log(E11)-1)
c.daille <- log((O11^3)/E11)
c.pmi <- log(O11/E11)
c.hyperg <- -phyper(O11, C1, (N-C1), (2*dis*R1), lower.tail=FALSE, log.p=TRUE)
Or <- O11/R1
Oc <- O11/C1
OCR <- cbind(Or,Oc)
c.ms <- apply(OCR,1,min)

HDD <- cbind(RCO,c.dice,c.poiss,c.daille,c.pmi,c.hyperg,c.ms,O11, E11)
cat(N,Ntr,ef1,ef2,ef3,ef311,sum(E11),"\n")

# transformation en indices
R1 <- R1/mean(R1)*100
C1 <- C1/mean(C1)*100
O11 <- O11/mean(O11)*100
c.dice <- c.dice/mean(c.dice)*100
c.poiss <- c.poiss/mean(abs(c.poiss))*100
c.daille <- c.daille/mean(c.daille)*100
c.pmi <- c.pmi/mean(c.pmi)*100
c.hyperg <- c.hyperg/mean(c.hyperg)*100
c.ms <- c.ms/mean(c.ms)*100

# choix du coefficient
if (coeff=="dice") c.coeff <- c.dice
if (coeff=="poiss") c.coeff <- c.poiss
if (coeff=="daille") c.coeff <- c.daille
if (coeff=="pmi") c.coeff <- c.pmi
if (coeff=="hyperg") c.coeff <- c.hyperg
if (coeff=="ms") c.coeff <- c.ms

# préparation du graphique
mima <- range(c(R1,C1,O11))
opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(4.2,4.2,2,0.5))
plot(mima, main = "", xlab = "", ylab = "Indices (base 100)",xlim=c(1,N),ylim=mima, type="n")
axis(side=1, at=br2, labels=listdateb, line=2)
title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
# mtext(sb, 3, line=0,cex=.8, font=1, adj=0)

# affichage des valeurs observées
lines(br2,R1,type="o", col="blue", lty=3, pch=21)
lines(br2,C1, type="o", col="red", lty=3, pch=21)
lines(br2,O11, type="o", col= "green", lty=2, pch=19)

# affichage des tendances et de la légende
mov <- lowess(c.coeff)
lines(br2, mov$y, col="seagreen", lwd=4)
mov1 <- lowess(R1)
lines(br2, mov1$y, col="blue", lwd=1)
mov2 <- lowess(C1)
lines(br2, mov2$y, col="red", lwd=1)
legend("bottomright", legend=c(lm1,lm2,coeff), cex=.7, lty=c(1,1,1),lwd=c(1,1,4), col=c("blue","red","seagreen"))

return(HDD)
#}
}

#####################################################################################

stablem <- function(corp="PL", lm, trnch="", mod="R") {

# étude de la stabilité d'un lemme dans un corpus
# moyenne et variances cumulées
# NB les courbes aléatoires sont différentes à chaque essai : 
# les résultats de la fonction sample() sont différents à chaque éxécution
# >> bootstrap rapide
# AG 9/2014   GPL3

library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
t1 <- Sys.time()

LL <- "Lsub"

# vérification
provl <- paste(corp, ".lemma", sep="")
idprov <- cqi_regex2id(provl, lm)
if (length(idprov)==0) stop ("lemme inexistant dans ce corpus !", call.=FALSE)

# récupération des occurrences
quer1 <- paste ('[lemma="', lm, '"]', sep ="")
cqi_query(corp, LL, quer1)
sscorp <- paste(corp, LL, sep = ":")
sscorp.dump <- cqi_dump_subcorpus(sscorp)
sscorp.dump <- as.vector(sscorp.dump[,1])
sscorp.ef <- length(sscorp.dump)
corp.size <- size(corpus(corp))
sscorp.var <- NULL
sscorp.moy <- NULL
if (trnch=="") {
  trnch <- round(corp.size/12000) # par défaut, tranches de 12000 tokens
  }

# découpage
br <- round(seq(1, corp.size, length.out=trnch+1))
H <- hist(sscorp.dump, breaks=br, plot=FALSE)
val.trnch <- H$counts
rm(H)
gc()

# calcul de la variance cumulée
if (mod=="A") { # mise en ordre aléatoire (avec bootstrap rapide)
  cat("Les tirages aléatoires successifs demandent du temps ! \n")
  sscorp.vartab <- matrix(0, nrow=trnch, ncol=100)
  for (j in 1:100) {
  val.trnch<- sample(val.trnch)
  for (i in 1:trnch) {
  sscorp.vartab[i,j] <- var(val.trnch[1:i])
  }
  }
  sscorp.var <- apply(sscorp.vartab, 1, mean)

  sscorp.moytab <- matrix(0, nrow=trnch, ncol=100)
  for (j in 1:100) {
  val.trnch<- sample(val.trnch)
  for (i in 1:trnch) {
  sscorp.moytab[i,j] <- mean(val.trnch[1:i])
  }
  }
  sscorp.moy <- apply(sscorp.moytab, 1, mean)
 }

 if (mod=="R") {
for (i in 1:trnch) {
  sscorp.var[i] <- var(val.trnch[1:i])
}

for (i in 1:trnch) {
  sscorp.moy[i] <- mean(val.trnch[1:i])
}
}

sstitre <- "ordre réel (dans le corpus)"
if (mod=="A") {
  sstitre <- "ordre aléatoire (bootstrap rapide)"
  }

# affichages
layout(matrix(1:4, 2,2, byrow=TRUE))
titre2 <- paste("corpus : ",corp," ; lemme : ", lm, " ; effectif total = ", sscorp.ef, sep="")
H <- hist(sscorp.dump, breaks = br, freq=TRUE, border = "grey", main = "", xlab = "ordre du corpus", ylab = "Effectif par tranche")
title(main = titre2, line=1.5, cex.main=1, font.main=1, adj = 0)
mov <- lowess(H$counts)
L <<- mov$y
decal <- (corp.size/trnch)/2
br2 <- (br[1:(length(br)-1)]) + decal
B <<- br2
lines(br2, mov$y, col="green", type="l", lwd=1.5)

#dev.new()
titranal <- "évolution de la moyenne cumulée"
plot(x=(1:trnch),y=sscorp.moy, type="b", xlab = "nombre de tranches cumulées considérées", ylab="moyenne", col="blue", sub=sstitre)
title(main = titre2, line=1.5, cex.main=1, font.main=1, adj = 0)
mtext(titranal, 3, line=0.3,cex=1, font=1, adj=0)
#dev.new()
titranal <- "effectifs par tranche"
hist(val.trnch, xlab="histogramme des fréquences par tranche", ylab="effectifs", main="", sub=sstitre)
title(main = titre2, line=1.5, cex.main=1, font.main=1, adj = 0)
mtext(titranal, 3, line=0.3,cex=1, font=1, adj=0)
#
titranal <- "évolution de la variance cumulée"
plot(x=(1:trnch),y=sscorp.var, type="b", xlab = "nombre de tranches cumulées considérées", ylab="variance", col="red", sub=sstitre)
title(main = titre2, line=1.5, cex.main=1, font.main=1, adj = 0)
mtext(titranal, 3, line=0.3,cex=1, font=1, adj=0)

stabvar <- as.list(NULL)
stabvar$val <- val.trnch
stabvar$cum <- sscorp.var

t2 <- Sys.time()
td <- difftime(t1,t2, units="secs")
cat("\n","temps écoulé :",round(td), "secondes \n")

return (stabvar)
}


##########################################################
# deuxième groupe : cooccurrences
##########################################################

coocA <- function(corp, reg="", piv, poscooc="QLF SUB VBE", dis=5, coeff="dice", selec=1, nbs=30, D=0, F=Inf, attr="", val="") {

# novembre 2012 - juillet 2014  Alain GUERREAU -> ||GPL3||

# ===== COOCCURRENCES DE PREMIER ORDRE ======
#        = cooccurrents ordinaires

# choix du corpus, du pivot, des pos des cooccurrents, de la largeur de la fenêtre,
# du type de coefficient, du nombre de cooccurrents retenus
# sélection possible d'une tranche, d'une paire attribut/valeur.

if (reg != "") {
  reg0 <- Sys.getenv("CORPUS_REGISTRY")
  Sys.setenv(CORPUS_REGISTRY=reg)
}

library(MASS, quietly=TRUE, warn.conflicts=FALSE)
library(plyr, quietly=TRUE, warn.conflicts=FALSE)
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
options(stringsAsFactors=FALSE)
gcinfo(verbose=FALSE)
gc()
t1 <- Sys.time()

#
# vérifications de base
#

# du corpus
lst.corp <- cqi_list_corpora()
if (corp %in% lst.corp == FALSE) stop (paste(corp, ": aucun corpus de ce nom !"), call.=FALSE)

# du pivot
corp.lm <- paste(corp, ".lemma", sep="")
corp.ps <- paste(corp, ".pos", sep="")
id.piv <- cqi_str2id(corp.lm, piv)

if (id.piv < 0) {  # si lemme inconnu, recherche simple par joker affixé
  id.prov <- cqi_regex2id(corp.lm, piv) # id provisoire
  if (length(id.prov)==0){
    piv2 <- paste(piv, ".*", sep="")    # on ajoute un joker
    id.prov <- cqi_regex2id(corp.lm, piv2) # vecteur des id correspondant à piv+joker
    }
  if (length(id.prov)==0) stop(paste(piv, ": lemme inconnu ! utilisez la fonction flou() !"),call.=FALSE) 
  str.prov <- cqi_id2str(corp.lm, id.prov)
  cpos.prov <- NULL
  pos.prov <- NULL
  freq.prov <- NULL
  for (i in 1:length(id.prov)){
    cpos.prov[i] <- cqi_id2cpos(corp.lm, id.prov[i])
    pos.prov[i] <- cqi_cpos2str(corp.ps, cpos.prov[i])
    freq.prov[i] <- cqi_id2freq(corp.lm, id.prov[i])
    cat(str.prov[i], "/", pos.prov[i], "/", freq.prov[i], "\n")
    }
  if (length(id.prov) > 0) stop ("Choisissez parmi ces lemmes !", call.=FALSE)
  }

# du coefficient
lst.coeff <- c("dice","poiss","daille","pmi","hyperg","mins")
if (coeff %in% lst.coeff == FALSE) stop ("coefficient inconnu, choisir : daille/dice/hyperg/mins/pmi/poiss", call.=FALSE)

# de la borne supérieure
CORP <- corpus(corp) # corpus complet
CORP.ln <- size(CORP)
lgcorp <- CORP.ln
if (F == Inf) F <- CORP.ln
if (F > CORP.ln) F <- CORP.ln
lgcorp <- F-D

# de la paire attribut/valeur (validité de la paire : AGENDUM !!!)
if ((attr!="" & val=="") | (attr=="" & val!="")) stop ("Requête incomplète !", call.=FALSE)
# cat ("\n Vérifications simples : OK ! \n")

#
# établissement de la liste des cpos du pivot
#
LL <- "Lsub" # nom générique de sous-corpus
A2 <- ""

if (attr=="" & val=="") {
quer1 <- paste ("[lemma='", piv, "']", sep ="")
cqi_query(corp, LL, quer1)
sscorp <- paste(corp, LL, sep = ":")
if (cqi_subcorpus_size(sscorp)==0) stop ("Pas de réponse, vérifiez la requête !", call.=FALSE)
dump.pivot <- cqi_dump_subcorpus(sscorp)[,1]
cqi_drop_subcorpus(sscorp)
# CORPu.ln <- F-D
}
if (attr!="" & val!="") {
# cat("La recherche du pivot dans une fraction du corpus total prend du temps (beaucoup)... \n")
LL2 <- "Lsub2"
quer2 <- paste("[lemma='", piv, "' %cd]", "::match.", attr, "=\"", val, "\"", sep="")
cqi_query(corp, LL2, quer2)
sscorp2 <- paste(corp, LL2, sep = ":")
if (cqi_subcorpus_size(sscorp2)==0) stop ("Pas de réponse, vérifiez la requête !", call.=F)
dump.pivot <- cqi_dump_subcorpus(sscorp2)[,1]
cqi_drop_subcorpus(sscorp2)
}
gc()
cpos.pivot <- dump.pivot[dump.pivot>D & dump.pivot <F] # tranche simple
eff.pivot <- length(cpos.pivot)
if (eff.pivot == 0) stop("Aucune occurrence du pivot dans la sélection !")
pos.piv <- cqi_cpos2str(corp.ps ,cpos.pivot[1])

#
# établissement de la liste des POS retenus (lst.pos.r)
# acrobatie nécessaire pour le cas des pos subdivisés ("radical"+subdivision)
# e.g. VER = VER.cond + VER.futu + VER.impe + VER.impf .....
# > il faut passer des "c.hoisis" aux "r.éels"
#

lst.pos.tt <- cqp_flist(CORP, "pos") # liste complète des POS du corpus
lst.pos.tt <- sort(names(lst.pos.tt))
pos.coo <- unlist(strsplit(poscooc, " ")) # liste des POS choisis
nbpos.c <- length(pos.coo) # nombre des POS choisis
lst.pos.ri <- NULL
lst.pos.r <- NULL
lst.pos.c <- NULL

for (i in 1:nbpos.c) {  # on traite "radical" par "radical"
  lpp <- grep(pos.coo[i], lst.pos.tt) # vecteur des POS correspondant à un "radical" (indices)
  lpp2 <- rep(pos.coo[i], length(lpp))# liste parallèle des POS choisis
  lst.pos.ri <- c(lst.pos.ri, lpp) # constitution par étape du vecteur (indices) des POS r.éels
  lst.pos.c <- c(lst.pos.c, lpp2) # idem pour les POS c.hoisis (vecteur strings)
}

lst.pos.r <- lst.pos.tt[lst.pos.ri] # passage à un vecteur de strings
corr.pos.rc <- cbind(as.vector(lst.pos.r),as.vector(lst.pos.c)) # tableau de correspondance réels/choisis (strings)
nb.pos.r <- 1
if (length(lst.pos.r) > 1) {     # mises en ordre
corr.pos.rc <- corr.pos.rc[order(corr.pos.rc[,1]),]
lst.pos.r <- corr.pos.rc[,1]
nb.pos.r <- length(lst.pos.r)
}

#
# construction de la table des cooccurrents (ids)
# 2 sorties : * tableau lexical (une ligne par occ du pivot) * vecteur cumulatif
#

ids.GD <- as.vector(NULL)
ps.GD <-as.vector(NULL)
ids.tablex <- matrix(0, nrow=eff.pivot, ncol=(dis*2))
n <- 0
for (i in dis:1) {  # à gauche, traitement par colonne
  n <- n + 1
  cpos.G <- cpos.pivot - i
  ids.G <- cqi_cpos2id(corp.lm, cpos.G)  
  ids.tablex[,n] <- ids.G
  ps.G <- cqi_cpos2str(corp.ps, cpos.G)  
  ids.GD <- c(ids.GD, ids.G)
  ps.GD <- c(ps.GD, ps.G)
  }   
  for (i in 1:dis) {  # à droite, traitement par colonne
  n <- n + 1
  cpos.D <- cpos.pivot + i
  ids.D <- cqi_cpos2id(corp.lm, cpos.D)
  ids.tablex[,n] <- ids.D
  ps.D <- cqi_cpos2str(corp.ps, cpos.D)
  ids.GD <- c(ids.GD, ids.D)
  ps.GD <- c(ps.GD, ps.D)
  }
cooc.GD <- cbind(ids.GD,ps.GD) # ids et pos correspondants, 2 colonnes

#
# tri des ids en fonction des POS sélectionnés
#

cooc.GD.pos <- matrix(FALSE, nrow=eff.pivot, ncol=(dis*2)) # on joue sur le fait qu'une matrice est un vecteur...
for (i in 1:nb.pos.r) {
  logx <- ps.GD == lst.pos.r[i] # vecteur ajoute TRUE (=1) pour chaque POS sélectionné
  cooc.GD.pos <- cooc.GD.pos + logx # 
}
cooc.GD.pos <- as.logical(cooc.GD.pos) # vecteur de sélection
cooc.GD.tr <- cooc.GD[cooc.GD.pos,] # on ne garde que ce qui est TRUE
for (i in 1:nrow(corr.pos.rc)) {
  asub <- corr.pos.rc[i,1]
  cooc.GD.tr[(cooc.GD.tr[,2]==asub),2] <- corr.pos.rc[i,2] # table à 2 colonnes : ids/POS (noms choisis) (triée)
}
# TAB.IDS <<- ids.tablex
# TAB.TF <<- cooc.GD.pos
#
# calcul des fréquences par lemme et troncature
#

cooc.GD.tr <- as.data.frame(cooc.GD.tr)
cooc.final <- as.data.frame(aggregate(1:nrow(cooc.GD.tr), cooc.GD.tr, length)) # on fait la somme par lemme > 3e colonne
cooc.final <- cooc.final[rev(order(cooc.final[,3])),]
cooc.final <- cooc.final[cooc.final[,3] > 2,] # troncature : au moins 3 cooccurrents

# compléments et mise en forme du tableau
id.n <- as.character(cooc.final[,1])
id.n <- as.numeric(id.n)
cooc.final[,1] <- id.n
cooc.final[,4] <- cqi_id2str(corp.lm, id.n) # on ajoute le lemme (string) > 4e colonne
cooc.final <- cooc.final[(cooc.final[,4]!=piv),]
cooc.final <- cooc.final[(cooc.final[,4]!="--"),]
cooc.final <- cooc.final[(cooc.final[,4]!="-"),]
cooc.final <- cooc.final[(cooc.final[,4]!="unknown"),]
 
sup <- rep(1, nrow(cooc.final)) # nettoyage des lemmes parasites restants (?)
cooc.final[,5] <- sup 
cooc.final <- cooc.final[order(cooc.final[,4]),]

if (nrow(cooc.final)<2) stop(paste(piv, "/", eff.pivot,"/", "pas de cooccurrent significatif !"), call.=FALSE) 

for (i in 1:(nrow(cooc.final)-1)) {
  if (cooc.final[i,4] == cooc.final[i+1,4]){
      cooc.final[i,5] <- 0
      cooc.final[i+1,5] <- 0
  }
}
cooc.final <- cooc.final[(cooc.final[,5]==1),]

#
# calcul des fréquences *totales* des lemmes retenus dans le corpus considéré
# difficulté : les sélections de sous-corpus (par cpos ou/et par attribut)
#

if (attr=="" & val=="" & D==0 & F==CORP.ln) { #    cas élémentaire : corpus complet
  cooc.final[,5] <- cqi_id2freq(corp.lm, cooc.final[,1])
  eff.tt.final <- CORP.ln
  cat("Effectif du corpus considéré :",eff.tt.final, "\n")
  }

if (attr=="" & val=="" & (D!=0 | F!=CORP.ln)) { #    simple troncature par cpos
  options(scipen=999)  # supprimer la notation scientifique (pb avec cqp)
  quer5 <- paste("abc:[lemma=\".*\" %cd]::abc >=",D," & abc <=",F, sep="")
  LL3 <- "Lsub3"
  cqi_query(corp, LL3, quer5) # création du sous-corpus ad hoc
  sscorp3 <- paste(corp, LL3, sep=":")
  eff.tt.final <- cqi_subcorpus_size(sscorp3)
  cat(eff.tt.final, "\n")  
  lst.freq <- cqi_fdist1(sscorp3, "match", "lemma") # liste de fréquences du ss-corpus sélectionné
  lst.str <- cqi_id2str(corp.lm, lst.freq[,1]) # on passe aux strings
  lst.freq <- as.data.frame(cbind(lst.freq, lst.str))
  lst.freq <- lst.freq[(lst.freq[,1] %in% cooc.final[,1]),] # on ne retient que les lemmes qui cooccurrent
  lst.freq <- lst.freq[(order(lst.freq[,3])),] # ordre alpha des lemmes
  cooc.final[,5] <- as.numeric(as.character(lst.freq[,2])) # on ajoute les fréquences globales au tableau
  }

if (attr!="" & val!="") { # sélections
  options(scipen=999)  # supprimer la notation scientifique (pb avec cqp)
  quer5 <- paste("abc:[lemma=\".*\" %cd & _.", attr, "=\"", val, "\"]::abc >=",D," & abc <=",F, sep="")
  LL3 <- "Lsub3"
  cqi_query(corp, LL3, quer5) # création du sous-corpus ad hoc
  lgcorp <- cqi_subcorpus_size(paste(corp,LL3, sep=":"))
  sscorp3 <- paste(corp, LL3, sep=":")
  eff.tt.final <- cqi_subcorpus_size(sscorp3)
  cat(eff.tt.final, "\n")
  lst.freq <- cqi_fdist1(sscorp3, "match", "lemma") # liste de fréquences du ss-corpus sélectionné
  lst.str <- cqi_id2str(corp.lm, as.numeric(lst.freq[,1])) # on passe aux strings
  lst.freq <- as.data.frame(cbind(lst.freq, lst.str))
  lst.freq <- lst.freq[(lst.freq[,1] %in% cooc.final[,1]),] # on ne retient que les lemmes qui cooccurrent
  lst.freq <- lst.freq[(order(lst.freq[,3])),] # ordre alpha des lemmes
  cooc.final[,5] <- as.numeric(as.character(lst.freq[,2])) # on ajoute les fréquences globales au tableau
  A2 <- cqi_dump_subcorpus(sscorp3)
  }

#
# filtrage des cooccurrents "significatifs" par divers coefficients-filtres
#

eff.piv <- rep(eff.pivot, nrow(cooc.final))
cooc.final[,6] <- as.numeric(eff.piv)
E11 <- as.numeric(cooc.final[,5])*as.numeric(cooc.final[,6])/as.numeric(eff.tt.final)
cooc.final[,7] <- E11
nams.tab <- c("ids.GD","ps.GD", "O11", "strings.cooc", "C1", "R1", "E11")
names(cooc.final) <- nams.tab

# coefficients de Dice-Soerensen, de Poisson-Stirling, de Daille, hypergéométrique, PMI, smin
# (Stefan EVERT, _The Statistics of Word Cooccurrences_, 2005, pp. 75-94 // www.collocations.de)

attach(as.data.frame(cooc.final),warn.conflicts=FALSE) # manoeuvre risquée...
c.dice <- (O11*100000)/(C1+R1)^selec
c.poiss <- O11*(log(O11)-log(E11)-1)
c.daille <- log((O11^3)/E11)
c.pmi <- log(O11/E11)
c.hyperg <- as.numeric(as.character(round(-phyper(O11, C1, (eff.tt.final-C1), (2*dis*R1), lower.tail=FALSE, log.p=TRUE), digits=10)))
Or <- O11/R1
Oc <- O11/C1
OCR <- cbind(Or,Oc)
c.ms <- apply(OCR,1,min)
coeff.tt <- cbind(c.dice,c.poiss,c.daille,c.pmi,c.hyperg,c.ms)

detach()
cooc.final <- cbind(cooc.final,coeff.tt)
cat(paste(piv, pos.piv, eff.pivot, sep=" / "))
cat("\n")

#
# mise en forme de l'affichage
#

position <- classrang(corp = corp, lm = piv)
cat (paste("classe = ",position[1]," sur ",position[2]," dans le corpus entier","\n\n", sep=""))

if (coeff=="dice") numcol <- 8 
if (coeff=="poiss") numcol <- 9 
if (coeff=="daille") numcol <- 10
if (coeff=="pmi") numcol <- 11
if (coeff=="hyperg") numcol <- 12
if (coeff=="mins") numcol <- 13
resultat <- cooc.final[rev(order(cooc.final[,numcol])),]
if (nrow(resultat) > 9) {
  cat ("somme 1 à 10 : ", sum(resultat[1:10,3]),"\n")
}
if (nrow(resultat) > 19) {
cat ("somme 1 à 20 : ", sum(resultat[1:20,3]),"\n")
}
if (nrow(resultat) > 29) {
cat ("somme 1 à 30 : ", sum(resultat[1:30,3]),"\n")
}
if (nrow(resultat) > 49) {
cat ("somme 1 à 50 : ", sum(resultat[1:50,3]),"\n")
}
cat("\n")
nb.lignes <- min(nbs, nrow(resultat))
resultat <- resultat[1:nb.lignes,]
ord.res <- order(resultat[,2], -resultat[,numcol])
resultat <- resultat[ord.res,]

cooc.affich.p <- as.data.frame(cbind(resultat[,2],resultat[,4],resultat[,3],resultat[,5],resultat[,numcol]))
cooc.affich <- cooc.affich.p[1,]
nm.affich <- c("POS", "lemme", "eff.cooc", "eff.lemme", coeff)
cooc.affich <- rbind(nm.affich,rep(" ", 5), cooc.affich)

for (i in 2:nrow(cooc.affich.p)) {
  if (cooc.affich.p[i,1]==cooc.affich.p[i-1,1]) cooc.affich <- rbind(cooc.affich, cooc.affich.p[i,])
  if (cooc.affich.p[i,1]!=cooc.affich.p[i-1,1]) cooc.affich <- rbind(cooc.affich, rep(" ",5))  
  }

pp <- ncol(cooc.affich)
cat(format(t(cooc.affich)), file = "", sep = c(rep(" ", pp - 1L), "\n"))

# création de la liste finale / résultat
# bizarreries : noms de variables d'une version antérieure, conservés pour correspondance avec les programmes suivants

cooc0 <- data.frame(id=cooc.final[,1],pos=cooc.final[,2],cooc=cooc.final[,3], lem=cooc.final[,4],eff.lem=cooc.final[,5],eff.pivot=cooc.final[,6],Dice=cooc.final[,8],PoissonS=cooc.final[,9],Daille=cooc.final[,10],HyperG=cooc.final[,12],PMI=cooc.final[,11],cum="")
cooc0 <- cooc0[rev(order(cooc0[,7])),]
cooc1 <- data.frame(id=resultat[,1],pos=resultat[,2],cooc=resultat[,3], lem=resultat[,4],eff.lem=resultat[,5],eff.pivot=resultat[,6],Dice=resultat[,8],PoissonS=resultat[,9],Daille=resultat[,10],HyperG=resultat[,12],PMI=resultat[,11],cum="")
DF <- paste(D,F, sep="//")
result <- list(call=match.call(),cooc0=cooc0,cooc1=cooc1,A2d3=A2,lgcorp=lgcorp,lemme=piv,cposLm=cpos.pivot,corpus=corp,effectif=eff.pivot,debfin=DF,pospiv=pos.piv,poscooc=poscooc,dist=dis,attr=attr,val=val,nb=nbs,NN=CORP.ln,tlex=ids.tablex)

t2 <- Sys.time()
td <- difftime(t1,t2, units="secs")
cat("\n","temps écoulé :",round(td), "secondes \n")

if (reg != "") {
  Sys.setenv(CORPUS_REGISTRY=reg0)
}

class(result) <- "COOCA"
return (result)
}

#################################################################
#################################################################

coocB <- function(tc, cex=.8, marg=4, prox="dice", selec=1){

# septembre 2012 AG -> ||domaine public||
# cooccurrents de deuxième ordre
# création de 2 lexicogrammes :
# > cooc dans le contexte uniquement
# > cooc dans le corpus entier


if (!inherits(tc, "COOCA"))
   stop("en entrée : un objet produit par coocA()", call.=F)
library(ade4, quietly=TRUE, warn.conflicts=FALSE)
library(circular, quietly=TRUE, warn.conflicts=FALSE)
library(MASS, quietly=TRUE, warn.conflicts=FALSE)
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
library(plyr, quietly=TRUE, warn.conflicts=FALSE)
options(warn = -1)
t1 <- Sys.time()
graphics.off()

# mises en place générales
#attach(t, warn.conflicts=FALSE)
corp <- tc$corpus
GDbon <- tc$cooc0
lnbs <- (tc$nb^2 - tc$nb)/2
ef <- tc$effectif
lico <- tc$cooc1[,4]
fr <- tc$cooc1[,5]
w1 <- as.numeric(tc$cooc1[,7])
attr <- tc$attr
val <- tc$val
poscooc <- tc$poscooc
pospiv <- tc$pospiv
nb <- tc$nb
dist <- tc$dist
effectif <- tc$effectif
cposLm <- tc$cposLm
lgcorp <- tc$lgcorp
tlex <- tc$tlex
txtitr <- paste("corpus : ", tc$corpus, "   bornes : ", tc$debfin, "   stock cible :  ", tc$NN,  "   lemme : ", tc$lemme, "   fréqu. : ", tc$effectif, "   fenêtre = +/-", tc$dist, "   nb cooccurrents =  ", tc$nb, sep="")
clr <- c("green3", "blue2", "red2", "magenta", "red", "darkcyan", "darkorange", "yellow4", "darkorchid", "darksalmon", "steelblue4", "sienna2", "ivory4")
df <- list(corpus=tc$corpus, debfin=tc$debfin, lgcorp=tc$lgcorp, lemme=tc$lemme, effectif=tc$effectif, nb=tc$nb, dis=tc$dist, cposLm=tc$cposLm, tlex=tc$tlex, cooc1=tc$cooc1, NN=tc$NN, A2d3 =tc$A2d3)

####################################################
# COOCCURRENCES DE DEUXIÈME ORDRE **T**
# = dénombrement des paires cooccurrentes
# dans le corpus entier

# -> effectifs dans un tableau simple
# calcul de la distance de Dice pour toutes les paires

c <- corpus(corp)
Clm <- paste(corp, ".lemma", sep="")

larg <- length(lico)
vdis1 <- vdis2 <- vdis3 <- vdis4 <- vdis5 <- vdis6 <- 0
n <- 0

###### !

# cas simple
if (attr=="" & val=="") {
for (i in 1:(larg-1)) {
  prem <- i + 1
  cat("*")
  for (j in prem:larg) {
    n <- n + 1
    a <- lico[i]
    b <- lico[j]
#    dem1 <- paste("[lemma = '" , a , "'][]{0," , marg , "}[lemma = '" , b , "']", sep="")
#    dem2 <- paste("[lemma = '" , b , "'][]{0," , marg , "}[lemma = '" , a , "']", sep="")
#    N1 <- subcorpus(c, dem1)
#    N2 <- subcorpus(c, dem2)
#    eff <- size(N1) + size(N2)
# correction Krzysztof Nowak 2/2015
dem <- paste("MU(meet[lemma = '" , a , "'][lemma = '" , b , "']"," -", marg+1, " " , marg+1,")", sep="")
N <- subcorpus(c,dem)
eff <- size(N)

    vdis3[n] <- eff
    vdis1[n] <- lico[i]
    vdis2[n] <- lico[j]
    vdis4[n] <- fr[i]
    vdis5[n] <- fr[j]
    }
  }
}

# ici une possible vectorisation
# effet pratique quasi nul...
#tbin <- cbind(vdis1, vdis2)
#calcooc <- function(x, marg=marg) {
#  a <- x[1]
#  b <- x[2]
#  dem1 <- paste("[lemma = '" , a , "'][]{0," , marg , "}[lemma = '" , b , "']",  sep="")
#  dem2 <- paste("[lemma = '" , b , "'][]{0," , marg , "}[lemma = '" , a , "']",  sep="")
#  N1 <- subcorpus(c, dem1)
#  N2 <- subcorpus(c, dem2)
#  eff <- size(N1) + size(N2)
#  return(eff)
#}
#vdis3 <- apply(tbin, 1, calcooc, marg)

# cas complexe
if (attr!="" & val!="") {
for (i in 1:(larg-1)) {
  prem <- i + 1
  cat("*")
  for (j in prem:larg) {
    n <- n + 1
#    cpt <- cpt + 1
    a <- lico[i]
    b <- lico[j]
    dem1 <- paste("[lemma = '" , a , "'][]{0," , marg , "}[lemma = '" , b , "']", "::match.", attr, "=\"", val, "\"", sep="")
    dem2 <- paste("[lemma = '" , b , "'][]{0," , marg , "}[lemma = '" , a , "']", "::match.", attr, "=\"", val, "\"", sep="")
    N1 <- subcorpus(c, dem1)
    N2 <- subcorpus(c, dem2)
    eff <- size(N1) + size(N2)
    vdis3[n] <- eff
    vdis1[n] <- lico[i]
    vdis2[n] <- lico[j]
    vdis4[n] <- fr[i]
    vdis5[n] <- fr[j]
    }
  }
}

vdis6 <- vdis3 * 1000 / (vdis4 + vdis5)^selec
vdis6 <- sprintf("%06.3f", vdis6)
vdis <- cbind(vdis1, vdis2, vdis3, vdis4, vdis5, vdis6)
vdis <- as.data.frame(vdis, stringsAsFactors=FALSE)
names(vdis) <- c("lm1", "lm2", "cooc", "ef_lm1", "ef_lm2", "Dice")

# création de la matrice des distances
cat("\n")
dice <- vdis[,6]
bid <- rep(0,larg^2)
tabdisLm <- matrix(bid, nrow=larg, ncol=larg)
n <- 0
for (i in 1:(larg-1)) {
  prem <- i + 1
  for (j in prem:larg) {
    n <- n + 1
    tabdisLm[i,j] <- dice[n]
    tabdisLm[j,i] <- dice[n]
    }
  }
tabdisLm <- as.numeric(tabdisLm)
tabdisLm <- matrix(tabdisLm, nrow=larg, ncol=larg)
rownames(tabdisLm) <- lico
colnames(tabdisLm) <- lico

vdis <- vdis[rev(order(vdis[,6])),]
df$coocT <- vdis
df$distmatT <- tabdisLm
tabdisLm <- log1p(tabdisLm ^ .2)

# création d'un tableau de contingence des contextes
# une ligne par contexte, nbs colonnes (=coocccurrents retenus)
# le 1 signifie que le lemme est présent dans le contexte considéré
idlem <- tc$cooc1[,c(1,4)]
idlem <- idlem[(order(idlem[,2])),]
idutils <- idlem[,1]
tlexb <- tlex %in% idutils
tlexb <- matrix(tlexb, ncol=ncol(tlex), nrow=nrow(tlex))
tlexS <- apply(tlexb, 1, sum)
tlexu <- tlex[(tlexS>1),]
tlexcl <- alply(tlexu, 1, match, idutils)
tlexcl <- t(as.data.frame(tlexcl))

cretab <- function(x) {
  lig0 <- rep(0,50)
  lig1 <- x[!is.na(x)]
  lig0[lig1] <- 1
  return (lig0)
}

tabsup <- alply(tlexcl, 1, cretab)
tabsup <- t(as.data.frame(tabsup))
tabsup <- tabsup[,(1:nb)]
colnames(tabsup) <- idlem[,2]
# TABSUP <<- tabsup

# analyse en composantes principales
#  fréquences totales (dans le corpus)
tab.pca <- dudi.pca(tabdisLm, scannf=FALSE, nf=4)
contrib <- inertia.dudi(tab.pca, row.inertia=TRUE, col.inertia=TRUE)
df$ACPt <- tab.pca
df$CONTRt <- contrib
contrabs <- apply(contrib$col.cum, 1, sum)
df$CUMTOTt <- sort(contrabs)
#(hisden(contrabs))
T0 <- tab.pca$co
ymax <- max(T0[,2])
ymin <- min(T0[,2])
Tbon <- lisible(T0[,1],T0[,2],lab=row.names(T0),mn=ymin, mx=ymax,cex=(cex+.4))
T0[,1] <- Tbon[,1]
T0[,2] <- Tbon[,2]
T0[,5] <- tc$cooc1[,2]
poscoo <- unlist(strsplit(poscooc, " "))
nbpos <- length(poscoo)
for (i in 1:nbpos) {
  T0[(T0[,5]==poscoo[i]),6] <- clr[i]
}

# calcul des segments
if (prox=="spearman") {
  cortabdist <- cor(tabdisLm, method="spearman")}
if (prox=="dice") {
  cortabdist <- tabdisLm }
cortabdist <- 1 - cortabdist
rempl <- as.vector(rep(1, lnbs*3))
discor <- matrix(rempl, nrow=lnbs, ncol=3)
n <- 0
for (i in 1:(nb-1)) {
  prem <- i + 1
  for (j in prem:nb) {
    n <- n + 1
    discor[n,3] <- cortabdist[i,j]
    discor[n,1] <- rownames(cortabdist)[i]
    discor[n,2] <- colnames(cortabdist)[j]
    }
  }
  
discor <- discor[(order(discor[,3])),]
discor <- as.data.frame(discor)
R1r <- match(discor[,1], rownames(T0))
R2r <- match(discor[,2], rownames(T0))
discor[,4] <- T0[R1r, 1]
discor[,5] <- T0[R1r, 2]
discor[,6] <- T0[R2r, 1]
discor[,7] <- T0[R2r, 2]
discor1 <- discor[(1:12),]
discor2 <- discor[13:nb,]

# affichage des points-colonnes (lemmes)
# graphics.off()
dev.new()
titranal <- "ACP sur tableau des distances de Dice entre co-cooccurrents (corpus cible entier) + distances de Dice"
opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(0.5,0.5,1.7,0.5))
(plot(T0[,1], T0[,2], type="n", xlab="", ylab="", axes=FALSE, frame.plot=TRUE))
axis(side=1, labels=FALSE, tick=0)
axis(side=2, labels=FALSE, tick=0)
title(main = txtitr, line=1, cex.main=.8, font.main=1, adj = 0)
mtext(titranal, 3, line=0,cex=.8, font=1, adj=0)
segments(discor1[,4], discor1[,5], discor1[,6], discor1[,7], lwd=1.5, col="grey")
segments(discor2[,4], discor2[,5], discor2[,6], discor2[,7], lwd=.4, col="grey")

text(T0[,1], T0[,2], labels=row.names(T0), col=T0[,6], cex=cex)

# calcul des points-lignes
ttabsup <- t(tabsup)
suplignes <- supcol(tab.pca, ttabsup)


########################################################
# COOCCURRENCES DE DEUXIÈME ORDRE **C**
# = dénombrement des paires cooccurrentes
# calculs dans le contexte uniquement

if (ef < 200) {
  df$lgcorp <- lgcorp
  class(df) <- "COOCB"
  return(df)
}

# valeurs fixes
larg <- dist * 2
npair <- larg * (larg-1) / 2
long1 <- effectif
long2 <- effectif * npair
n <- 0
BINOM1 <- NULL
BINOM2 <- NULL
TBINOM <- matrix(seq(0, length(cposLm)*larg), nrow=length(cposLm), ncol=larg)

# création du tableau des ids des cooccurrents dans la fenêtre
for (i in dist:1) {
  cposG <- cposLm - i
  Gid <- cqi_cpos2id(Clm, cposG)
  n <- n+1
  TBINOM[,n] <- Gid
  }
  for (i in 1:dist) {
  cposD <- cposLm + i
  Did <- cqi_cpos2id(Clm, cposD)
  n <- n+1
  TBINOM[,n] <- Did
  }

# création de la suite des paires de cooccurrents
for (i in 1:(larg-1)) {
  cat (".", sep = "")
  prem <- i + 1
  for(j in prem:larg) {
    BINOM1 <- c(BINOM1,TBINOM[,i])
    BINOM2 <- c(BINOM2,TBINOM[,j])
  }
}
BINOM<- cbind(BINOM1, BINOM2)

# premier nettoyage
BINOM <- BINOM[(BINOM[,1]!=BINOM[,2]),]
BINOM <- BINOM[(BINOM[,1]%in%idutils),]
BINOM <- BINOM[(BINOM[,2]%in%idutils),]

if (nrow(BINOM) == 0) {
  cat("Pas de co-cooccurrents pour cette requête...")
  return (df)
  }
  
# permutations B-A -> A-B (si B < A)
# de manière à compter ensemble A-B et B-A
BINOMa <- BINOM[(BINOM[,1]<BINOM[,2]),]
BINOMb <- BINOM[(BINOM[,1]>BINOM[,2]),]
BINOMa <- as.data.frame(BINOMa)
BINOMb <- as.data.frame(BINOMb)
BINOMb[,3] <- BINOMb[,2]
BINOMb[,4] <- BINOMb[,1]
BINOMb <- BINOMb[,3:4]
names(BINOMb) <- c("BINOM1", "BINOM2")
BINOM  <- rbind(BINOMa, BINOMb)

# comptage des fréquences des paires
BINOM <- as.data.frame(BINOM)
liste <- as.data.frame(aggregate(1:nrow(BINOM), BINOM, length))
lnmax <- (nb^2 - nb)/2
df$densite <- nrow(liste)/lnmax

# compléments au tableau général et mise en forme
# liste <- liste[liste[,3]>4,]
liste[,4] <- cqi_id2str(Clm, liste[,1])
liste[,5] <- cqi_id2str(Clm, liste[,2])
i1 <- match(liste[,1], GDbon[,1])
liste[,6] <- GDbon[i1,5]
i2 <- match(liste[,2], GDbon[,1])
liste[,7] <- GDbon[i2,5]
Dice <- (liste[,3] * 1e+10) / (liste[,6] + liste[,7])^2
N <- sqrt(sum(Dice^2))
Dice <- Dice / N
liste[,8] <- Dice
liste <- liste[rev(order(liste[,8])),]

liste <- liste[(order(liste[,4])),]
nam <- c("id1", "id2", "F.cooc", "lem1", "lem2", "F.lm1", "F.lm2", "Dice")
names(liste) <- nam
or <- seq(1, nrow(liste))
row.names(liste) <- or

df$cooc2 <- liste
cooc2b <- liste[,c(3,4,5)]
cooc2b <- cooc2b[rev(order(cooc2b[,1])),]
cat("\n")
write.matrix(cooc2b[1:30,])

# CRÉATION DE LA MATRICE DE DISTANCES ENTRE LES COOCCURRENTS =====
#  = distances entre les "co-cooccurrents
li1 <- liste[,c(1,4)]
li2 <- liste[,c(2,5)]
nli <- c("ids", "lem")
names(li1) <- nli
names(li2) <- nli
listep1 <- rbind(li1, li2)
listep1 <- unique(listep1)
listep1 <- listep1[order(listep1[,2]),]
row.names(listep1) <- seq(1, nrow(listep1))
xn <- rep(0, (nrow(listep1)^2))
tabdisLm <- matrix(xn, nrow=nrow(listep1), ncol=nrow(listep1))
colnames(tabdisLm) <- listep1[,2]
rownames(tabdisLm) <- listep1[,2]

for (i in 1:nrow(liste)) {
  tabdisLm[(liste[i,4]),(liste[i,5])] <- liste[i,8]
  tabdisLm[(liste[i,5]),(liste[i,4])] <- liste[i,8]
  }
df$distmatc <- tabdisLm
tabdisLm <- log1p(tabdisLm ^ .2)

tabdisbrut <- matrix(xn, nrow=nrow(listep1), ncol=nrow(listep1))
colnames(tabdisbrut) <- listep1[,2]
rownames(tabdisbrut) <- listep1[,2]
for (i in 1:nrow(liste)) {
  tabdisbrut[(liste[i,4]),(liste[i,5])] <- liste[i,3]
  tabdisbrut[(liste[i,5]),(liste[i,4])] <- liste[i,3]
  }
df$distmatbrut <- tabdisbrut

w <- match(colnames(tabdisLm), tc$cooc1[,4])
w1 <- tc$cooc1[w,7]
w1 <- as.numeric(w1)
df$rowweights <- w1

# analyse en composantes principales
# fréquences dans le contexte uniquement
tab.pca <- dudi.pca(tabdisLm, scannf=FALSE, nf=4)
contrib <- inertia.dudi(tab.pca, row.inertia=TRUE, col.inertia=TRUE)
df$ACPc <- tab.pca
df$CONTRc <- contrib
contrabs <- apply(contrib$col.cum, 1, sum)
df$CUMTOTc <- sort(contrabs)
T0 <- tab.pca$co
ymax <- max(T0[,2])
ymin <- min(T0[,2])
Tbon <- lisible(T0[,1],T0[,2],lab=row.names(T0),mn=ymin, mx=ymax,cex=(cex+.4))
T0[,1] <- Tbon[,1]
T0[,2] <- Tbon[,2]
Tr <- match(row.names(tabdisLm), tc$cooc1[,4])
T0[,5] <- tc$cooc1[Tr,2]
poscoo <- unlist(strsplit(poscooc, " "))
nbpos <- length(poscoo)
for (i in 1:nbpos) {
  T0[(T0[,5]==poscoo[i]),6] <- clr[i]
}

# calcul des segments
if (prox=="spearman") {
  cortabdist <- cor(tabdisLm, method="spearman") }
if (prox=="dice") {
  cortabdist <- tabdisLm }

cortabdist <- 1 - cortabdist
rempl <- as.vector(rep(1, lnbs*3))
discor <- matrix(rempl, nrow=lnbs, ncol=3)
n <- 0
nb <- nrow(cortabdist)
for (i in 1:(nb-1)) {
  prem <- i + 1
  for (j in prem:nb) {
    n <- n + 1
    discor[n,3] <- cortabdist[i,j]
    discor[n,1] <- rownames(cortabdist)[i]
    discor[n,2] <- colnames(cortabdist)[j]
    }
  }
  
discor <- discor[(order(discor[,3])),]
discor <- as.data.frame(discor)
R1r <- match(discor[,1], rownames(T0))
R2r <- match(discor[,2], rownames(T0))
discor[,4] <- T0[R1r, 1]
discor[,5] <- T0[R1r, 2]
discor[,6] <- T0[R2r, 1]
discor[,7] <- T0[R2r, 2]
discor1 <- discor[(1:12),]
discor2 <- discor[13:nb,]

# affichage
dev.new()
titranal <- "ACP sur tableau des distances de Dice entre co-cooccurrents (contexte) + distances de Dice"
opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(0.5,0.5,1.7,0.5))
(plot(T0[,1], T0[,2], type="n", xlab="", ylab="", axes=FALSE, frame.plot=TRUE))
axis(side=1, labels=FALSE, tick=0)
axis(side=2, labels=FALSE, tick=0)
title(main = txtitr, line=1, cex.main=.8, font.main=1, adj = 0)
mtext(titranal, 3, line=0,cex=.8, font=1, adj=0)
segments(discor1[,4], discor1[,5], discor1[,6], discor1[,7], lwd=1.5, col="grey")
segments(discor2[,4], discor2[,5], discor2[,6], discor2[,7], lwd=.4, col="grey")

text(T0[,1], T0[,2], labels=row.names(T0), col=T0[,6], cex=cex)

# calcul des points-lignes
tabsup <- tabsup[,(listep1[,2])]
ttabsup <- t(tabsup)
df$tabsup2 <- tabsup
suplignes <- supcol(tab.pca, ttabsup)
df$tabsuplignes2 <- suplignes


t2 <- Sys.time()
td <- difftime(t1,t2, units="secs")
cat("\n","temps écoulé :",round(td), "secondes \n")

return(df)
}
################################################################################
# suppression des recouvrements dans un graphique factoriel.
# partant du centre, on écarte les points qui provoquent recouvrement,
# toujours vers l'extérieur (selon le quadrant), alternativement horizontalement
# et verticalement, de manière à éviter la déformation du nuage,
# en pondérant l'alternance par la proximité angulaire avec l'axe 1 ou 2 ;
# peut durer de quelques secondes à quelques minutes !!!
################################################################################

lisible <- function (x, y, lab,mn,mx, cex=.2){

# on constitue le tab(leau de )dep(art)
# library(circular, quietly=TRUE, warn.conflicts=FALSE)
eps <- 0.0000000001
tabdep <- as.data.frame(cbind(x,y,lab))
names(tabdep) <- c("x","y","lab")
row.names(tabdep) <- seq(1,nrow(tabdep))
tabdep$x <- as.numeric(as.character(tabdep[,1]))
tabdep$y <- as.numeric(as.character(tabdep[,2]))
tabdep$lab <- as.character(tabdep$lab)
htlet <- (mx-mn)/(30/cex) # choix empirique...
lglet <- htlet*.5         # choix empirique...
H <- lglet/2
indx <- as.numeric(row.names(tabdep))
d2 <- (tabdep$x^2)+(tabdep$y^2)
drt <- tabdep$x + (H*nchar(tabdep$lab))
gau <- tabdep$x - (H*nchar(tabdep$lab))
angl <- deg(atan(tabdep$y/tabdep$x))/.9

tabdep <- as.data.frame(cbind(tabdep,indx,d2,drt,gau,angl))
tt <- length(x)
# sens <- 0
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

# on écarte tous les points du rub(an) alternativement horizontalement et verticalement, vers l'extérieur du quadrant
# en combinant la valeur de l'angle et un nombre aléatoire (!)

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


##########################################################################
##########################################################################

coocE <- function(tc, trnch=3, cex=.8) {

# septembre 2012 AG -> ||GPL3||
# table des lemmes cooccurrents divisée en tranches chronologiques
# de même taille, classées par centres de gravité,
# graphe de l'évolution des fréquences
# + graphe Bertin-Cibois

library(ade4, quietly=TRUE, warn.conflicts=FALSE)
if (!inherits(tc, "COOCA"))
   stop("en entrée : un objet produit par coocA()", call.=F)
attach(tc, warn.conflicts=FALSE)
# tabdisLm <- distmatc
idutils <- cooc1[,1]
ef <- length(cposLm)


#### moyennes réciproques (sur colonnes)

mrecip1 <- function(x) {
li.m <- rep(0, nrow(x))
for (j in 1:ncol(x)){
   for (i in 1:nrow(x)){
      li.m[i] <- li.m[i]+(x[i,j]*j)
   }
   }
li.m <- li.m / rowSums(x)
x <- x[sort.list(li.m),]
return(x)
}

### affichage méthode Bertin-Cibois
### (assocplot() modifié)

repfac2 <- function (x, col = c("red", "white"), border= c("red", "blue"), main = NULL,
    xlab = NULL, ylab = NULL)  {
    opar <- par
    par(mar = c(2.5,4.8,2.0,0.1))
    n <- sum(x)
    f <- x[, NCOL(x):1]
    e <- outer(rowSums(f), colSums(f), "*")/n
    d <- as.matrix(f - e)/sqrt(e)
    e <- sqrt(e)
    x.w <- apply(e, 1, max)
    y.h <- apply(d, 2, max) - apply(d, 2, min)
    ylim <- c(0, sum(x.w))
    xlim <- c(0, sum(y.h))
    plot.new()
    plot.window(xlim, ylim, log = "")
    x.r <- cumsum(x.w)
    x.m <- (c(0, x.r[-NROW(f)]) + x.r)/2
    y.u <- cumsum(y.h)
    y.m <- y.u - apply(pmax(d, 0), 2, max)
    x.m <- rev(x.m)
    y.m <- rev(y.m)
    z <- expand.grid(x.m, y.m)
    rect(z[, 2], z[, 1] - e/2, z[, 2] + d, z[, 1] + e/2,col = col[1 + (d < 0)], border= border[1+(d<0)], lwd=2)
    txtitr <- paste("corpus : ", corpus, "   bornes : ", debfin, "   lemme : ", lemme, "   fréqu. : ", effectif, "   fenêtre = +/-", dist, "   nb cooccurrents =  ", nb, sep="")
    titranal <- "Graphe de Bertin-Cibois de la répartion des cooccurrents par tranches chronologiques (écarts à l'indépendance)"
    axis(2, at = x.m, labels = rownames(f), tick = FALSE, cex.axis=cex, las=2)
    axis(1, at = y.m, labels = colnames(f), tick = FALSE, cex.axis=1.0, las=1)
    abline(v = y.m, lty = 2)
    ndn <- names(dimnames(f))
    if (length(ndn) == 2) {
        if (is.null(xlab))
            xlab <- ndn[1]
        if (is.null(ylab))
            ylab <- ndn[2]
    }
    title(main = txtitr, line=1, xlab = xlab, ylab = ylab, cex.main=.8, font.main=1, adj=0)
    mtext(text=titranal, 3, line=0, cex=.8, font=1, adj=0)
    par <- opar
}
########################################################################

# création d'un tableau lexical sélectionné et trié ======

tablexT <- matrix(1, nrow=length(cposLm), ncol=nb)
for (i in 1:nb) {
  tabpro <- (tlex %in% idutils[i])
  tabpro <- matrix(tabpro, nrow=length(cposLm), ncol=dist*2)
  stab <- apply(tabpro, 1, sum)
  tablexT[,i] <- stab
  }
colnames(tablexT) <- cooc1[,4]
rownames(tablexT) <- cposLm
T <- tablexT

if (ef < 100) {
  cat("\n" , "ATTENTION ! moins de 100 occurrences du lemme choisi", "\n")
  cat("   Résultats à interpréter avec une extrême prudence !")
}

# découpage et traitement du tableau
# (on découpe le tableau en n groupes comportant le même nombre de lignes)
T <- tablexT
n <- trnch
eff <- nrow(T)
rt <- seq(1, eff, 1)
brn <- quantile(rt, probs=seq(0,1,1/n), type=1)
L <- list(NULL)
for (i in 1:n) {
  L[[i]] <- T[(brn[i]:brn[i+1]),]
}

rp <- (rep(0, (n*ncol(tablexT)) ))
TabC <- matrix(rp, ncol=n, nrow=ncol(tablexT))
for (i in 1:n) {
  S <- apply(L[[i]], 2, sum)
  TabC[,i] <- S
  }
rownames(TabC) <- cooc1[,4]
cnam <- NULL
for (i in 1:trnch) {
  rn <- paste("P", i, sep="")
  cnam <- c(cnam, rn)
}
colnames(TabC) <- cnam
txtitr <- paste("corpus : ", corpus, "   bornes : ", debfin, "   lemme : ", lemme, "   fréqu. : ", effectif, sep="")

# premier graphique : chronologie générale
if (attr=="" & val=="") {
opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(4.0,3.4,1.7,0.5))
trnch2 <- 20
cposbrn <- cposLm[brn]
br <- seq(1, NN, length.out=trnch2+1)
decal <- (NN/trnch2)/2
#mn <- paste("lemme : ", lemmme, " (effectif du lemme = ", effectif, ")", sep = "")
sb <- paste("index du corpus cible (effectif total = ", NN, ")", sep="")
H <- hist(cposLm, breaks = br, border = "grey", main = "", xlab = sb, ylab = "Effectif par tranche", cex.lab=.8, cex.axis=.7, mgp=c(2,1,0))
}

if (attr!="" & val!="") {
opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(4.0,3.4,1.7,0.5))
trnch2 <- 20
cposbrn <- cposLm[brn]
ef <- length(cposLm)
indexL <- as.integer(seq(1, ef, 1))
efI <- cbind(as.integer(cposLm, indexL))
# ici, une ruse diabolique !
ln2 <- nrow(A2d3)
indexT <- as.integer(seq(1, ln2, 1))
corpI <- cbind(as.integer(A2d3[,1]), indexT)
choix <- corpI[,1] %in% efI[,1]
corpI2 <- cbind(corpI, choix)
indbon <- corpI2[(corpI2[,3]==TRUE),2]
cposbrn <- indbon[brn]
#mn <- paste("lemme : ", lm, " (effectif du lemme = ", ef, ")", sep = "")
sb <- paste("index du corpus cible (effectif total = ", NN, ")", sep="")
br <- seq(1, ln2, length.out=trnch2+1)
decal <- (ln2/trnch2)/2
H <- hist(indbon, breaks = br, border = "grey", main = "", xlab = sb, ylab = "Effectif par tranche", cex.lab=.8, cex.axis=.7, mgp=c(2,1,0))
}

mov <- lowess(H$counts)
br2 <- (br[1:(length(br))-1]) + decal
lines(br2, mov$y, col="blue", lwd=1.5)
for (i in 2:(trnch)){
  points(cposbrn[i] ,0, pch="|",cex=3, col="turquoise")
}
title(main = txtitr, line=0, cex.main=.8, font.main=1, adj = 0)

# mise en ordre du tableau et calcul des écarts à l'indépendance
# affichage des écarts ordonnés
TAB <- mrecip1(TabC)
dev.new()
repfac2(TAB)
chi2 <- chisq.test(TAB)
dif <- round(chi2$observed-chi2$expected)
scT <- apply(TAB, 2, sum)
TAB <- rbind(TAB, scT)
srT <- apply(TAB, 1, sum)
TAB <- cbind(TAB, srT)

# affichages et sauvegardes
sink(file=lemme, append=TRUE)
cat("\n", date(), "\n\n")
print(TAB)
cat("\n\n")
print(dif)
sink(file=NULL)
cat("\n\n")
result <- cbind(row.names(TAB), TAB)
write.matrix(result, file="")
res <- list(observe=TAB, difference=dif)
class(res) <- "COOCC"

return(res)
}

##############################################################

coocC <- function(x, y) {

# mise en ordre de deux listes de cooccurrents
# pour comparaison simplifiée

if ( (!inherits(x, "COOCA")) | (!inherits(y, "COOCA")) )
   stop("en entrée : deux objets produits par coocA()", call.=F)

gau <- x$cooc1[,c(3,4,7)]
dro <- y$cooc1[,c(3,4,7)]
gd <- merge(gau, dro, by.x=2, by.y=2, all.x=TRUE, all.y=TRUE, incomparables= 0)
gd[is.na(gd)] <- 0
gd <- gd[order(gd[,5]),]
gdA <- gd[(gd[,4]==0),]
gdA <- gdA[(order(gdA[,3])),]
gd <- gd[(gd[,4]!=0),]
gd <- gd[rev(order(gd[,3])),]
gdB <- gd[(gd[,3]==0),]
gdB <- gdB[rev(order(gdB[,5])),]
gd <- gd[(gd[,3]!=0),]
rA <- rank(gd[,3])
rB <- rank(gd[,5])
rdiff <- rB - rA
gd <- gd[order(rdiff),]
gdf <- rbind(gdA,gd,gdB)
return(gdf)
}

###################################################################################

coocK2 <- function(corp, attr, lm1, lm2, dis=8, lim=200) {

# affichage d'une ligne de texte brut (20 tokens)
# pour une paire de lemmes (a - b  ou b - a)
# avec indication d'un attribut

library(plyr, quietly=TRUE, warn.conflicts=FALSE)
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
library(MASS, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)

obj  <- paste(corp, ".lemma", sep="")
obj2 <- paste(corp, ".word", sep="")
obj3 <- paste(corp, attr, sep=".")
a <- cqi_str2id(obj, lm1)
b <- cqi_str2id(obj, lm2)
if (a < 0) stop(paste(lm1, "n'existe pas !"), call.=FALSE)
if (b < 0) stop(paste(lm2, "n'existe pas !"), call.=FALSE)
listattr <- cqi_attributes(corp, "s")
if ((attr != "") & ((attr %in% listattr)==FALSE)) {
  stop(paste(attr, "n'est pas un attribut de ce corpus !"), call.=FALSE)
}
LL1 <- "Lsub1"
LL2 <- "Lsub2"
dem1 <- paste("[lemma = '" , lm1 , "'][]{0," , dis , "}[lemma = '" , lm2 , "']", sep="")
dem2 <- paste("[lemma = '" , lm2 , "'][]{0," , dis , "}[lemma = '" , lm1 , "']", sep="")

cqi_query(corp, LL1, dem1)
cqi_query(corp, LL2, dem2)
sc1 <- paste(corp, LL1, sep = ":")
sc2 <- paste(corp, LL2, sep = ":")
if (cqi_subcorpus_size(sc1) > 0) {
A1 <- cqi_dump_subcorpus(sc1)}
else A1 <- NULL
if (cqi_subcorpus_size(sc2) > 0) {
A2 <- cqi_dump_subcorpus(sc2)}
else A2 <- NULL
AA <- rbind(A1, A2)
if (is.null(AA)) stop("Aucun résultat !", call. = FALSE)
if (nrow(AA)> 1) AA <- AA[order(AA[,1]),]
res <- AA[,1]
res <- unique(res)
cat("\n","nombre d'occurrences :", length(res), "\n\n")


if (attr != "") {
res2 <- cqi_cpos2struc(obj3, res)
res3 <- cqi_struc2str(obj3, res2)
}
if (attr == "") {
res3 <- rep("", length(res))
}

d <- res - 6
vec <- rep(0, 20*length(d))
mat <- matrix(vec, nrow=length(d), ncol=20)
for (i in 1:20) {
  mat[,i] <- d + i
  }

ench2 <- function(x, obj2){
  a <- cqi_cpos2str(obj2, x)
  b <- paste(a[1:20], sep="", collapse=" ")
  return(b)
}

tabw <- apply(mat, 1, ench2, obj2)
res4 <- cbind(as.matrix(tabw), as.matrix(res3))
if (nrow(res4) <= lim) {
  write.matrix(res4, file="", sep=" ")
}
if (nrow(res4) > lim){
  tirage <- sort(as.integer(runif(lim,1,nrow(res4))))
  res5 <- res4[tirage,]
  write.matrix(res5, file="", sep=" ")
}
return(res4)
}

#############################################################################

coocK3 <- function(corp, attr, lm1, lm2, lm3, dis=8, lim=200) {

# affichage d'une ligne de texte brut (24 tokens)
# pour un triplet de lemmes (dans un ordre quelconque)
# avec indication d'un attribut

library(plyr, quietly=TRUE, warn.conflicts=FALSE)
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
library(MASS, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)

obj  <- paste(corp, ".lemma", sep="")
obj2 <- paste(corp, ".word", sep="")
obj3 <- paste(corp, attr, sep=".")
a <- cqi_str2id(obj, lm1)
b <- cqi_str2id(obj, lm2)
c <- cqi_str2id(obj, lm3)
if (a < 0) stop(paste(lm1, "n'existe pas !"), call.=FALSE)
if (b < 0) stop(paste(lm2, "n'existe pas !"), call.=FALSE)
if (c < 0) stop(paste(lm3, "n'existe pas !"), call.=FALSE)
listattr <- cqi_attributes(corp, "s")
if ((attr != "") & ((attr %in% listattr)==FALSE)) {
  stop(paste(attr, "n'est pas un attribut de ce corpus !"), call.=FALSE)
}

LL1 <- "Lsub1"
LL2 <- "Lsub2"
LL3 <- "Lsub3"
LL4 <- paste(LL1, LL2, sep="")
LL5 <- paste(LL1, LL3, sep="")
LL6 <- paste(LL2, LL3, sep="")
dem1 <- paste("[lemma = '" , lm1 , "'][]{0," , dis , "}[lemma = '" , lm2 , "'][]{0," , dis , "}[lemma = '" , lm3 , "']", sep="")
dem2 <- paste("[lemma = '" , lm1 , "'][]{0," , dis , "}[lemma = '" , lm3 , "'][]{0," , dis , "}[lemma = '" , lm2 , "']", sep="")
dem3 <- paste("[lemma = '" , lm2 , "'][]{0," , dis , "}[lemma = '" , lm1 , "'][]{0," , dis , "}[lemma = '" , lm3 , "']", sep="")
dem4 <- paste("[lemma = '" , lm2 , "'][]{0," , dis , "}[lemma = '" , lm3 , "'][]{0," , dis , "}[lemma = '" , lm1 , "']", sep="")
dem5 <- paste("[lemma = '" , lm3 , "'][]{0," , dis , "}[lemma = '" , lm1 , "'][]{0," , dis , "}[lemma = '" , lm2 , "']", sep="")
dem6 <- paste("[lemma = '" , lm3 , "'][]{0," , dis , "}[lemma = '" , lm2 , "'][]{0," , dis , "}[lemma = '" , lm1 , "']", sep="")

cqi_query(corp, LL1, dem1)
cqi_query(corp, LL2, dem2)
cqi_query(corp, LL3, dem3)
cqi_query(corp, LL4, dem4)
cqi_query(corp, LL5, dem5)
cqi_query(corp, LL6, dem6)

sc1 <- paste(corp, LL1, sep = ":")
sc2 <- paste(corp, LL2, sep = ":")
sc3 <- paste(corp, LL3, sep = ":")
sc4 <- paste(corp, LL4, sep = ":")
sc5 <- paste(corp, LL5, sep = ":")
sc6 <- paste(corp, LL6, sep = ":")
if (cqi_subcorpus_size(sc1) > 0) {
A1 <- cqi_dump_subcorpus(sc1)}
else A1 <- NULL
if (cqi_subcorpus_size(sc2) > 0) {
A2 <- cqi_dump_subcorpus(sc2)}
else A2 <- NULL
if (cqi_subcorpus_size(sc3) > 0) {
A3 <- cqi_dump_subcorpus(sc3)}
else A3 <- NULL
if (cqi_subcorpus_size(sc4) > 0) {
A4 <- cqi_dump_subcorpus(sc4)}
else A4 <- NULL
if (cqi_subcorpus_size(sc5) > 0) {
A5 <- cqi_dump_subcorpus(sc5)}
else A5 <- NULL
if (cqi_subcorpus_size(sc6) > 0) {
A6 <- cqi_dump_subcorpus(sc6)}
else A6 <- NULL
AA <- rbind(A1, A2, A3, A4, A5, A6)
if (is.null(AA)) stop("Aucun résultat !", call. = FALSE)
if (nrow(AA) > 1)   AA <- AA[order(AA[,1]),]
res <- AA[,1]
res <- unique(res)
cat("\n","nombre d'occurrences :", length(res), "\n\n")

if (attr != "") {
res2 <- cqi_cpos2struc(obj3, res)
res3 <- cqi_struc2str(obj3, res2)
}
if (attr == "") {
res3 <- rep("", length(res))
}

d <- res - 6
vec <- rep(0, 20*length(d))
mat <- matrix(vec, nrow=length(d), ncol=20)
for (i in 1:20) {
  mat[,i] <- d + i
  }

ench2 <- function(x, obj2){
  a <- cqi_cpos2str(obj2, x)
  b <- paste(a[1:20], sep="", collapse=" ")
  return(b)
}

tabw <- apply(mat, 1, ench2, obj2)
res4 <- cbind(as.matrix(tabw), as.matrix(res3))
if (nrow(res4) <= lim) {
  write.matrix(res4, file="", sep=" ")
}
if (nrow(res4) > lim){
  tirage <- sort(as.integer(runif(lim,1,nrow(res4))))
  res5 <- res4[tirage,]
  write.matrix(res5, file="", sep=" ")
}
return(res4)
}

#########################################################################
# troisième groupe : fonctions utilitaires
#########################################################################

flou <- function(corp, type, larg=0, mot) {

# fonction de recherche floue
# choix du corpus, du type (lemma ou word)
# et de la distance tolérée
# correction mars 2016 :
# liste globale par appel à fonction externe
# beaucoup plus efficace

library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
atr <- cqi_attributes(corp, "p")
if (!(type %in% atr)) stop ('paramètre "type" inutilisable !', call.=FALSE)

params1 <- paste("-f -P ",type," ", corp,  sep="")
rang1 <-  system2(command="cwb-lexdecode", args=params1, stdout=TRUE) 
rang1b <- substr(rang1, 9, 25) # les effectifs sont sur les 7 premiers car+1

fuzzy <- agrep(mot, rang1b, max=larg, ignore.case=TRUE, value=TRUE)
# print (fuzzy)
if (length(fuzzy)==0) stop ('pas de réponse, modifiez le paramètre "larg" !', call.=FALSE)

# affichages différents pour word et lemma : 2 routines distinctes

if (type == "word") {
  Leff <- Llm <- Lpos <- NULL
  cql <- paste(corp,".", "lemma", sep="")
  cqp <- paste(corp,".", "pos", sep="")
  cq <- paste(corp,".",type, sep="")
  for (i in 1:length(fuzzy)){
    Nid <- cqi_str2id(cq, fuzzy[i])
    Ncpos <- cqi_id2cpos(cq, Nid)
    Leff[i] <- as.numeric(length(Ncpos))
    Lpos[i] <- cqi_cpos2str(cqp, Ncpos[1])
    Llm[i] <- cqi_cpos2str(cql, Ncpos[1])
  }
  Ltt <- cbind(fuzzy, Leff, Lpos, Llm)
  Ltt <- Ltt[order(Ltt[,1]),]
  return (Ltt)
}

if (type == "lemma") {
  Leff <- Lpos <- NULL
  cqp <- paste(corp,".", "pos", sep="")
  cq <- paste(corp,".",type, sep="")
  for (i in 1:length(fuzzy)){
    Nid <- cqi_str2id(cq, fuzzy[i])
    Ncpos <- cqi_id2cpos(cq, Nid)
    Leff[i] <- as.numeric(length(Ncpos))
    Lpos[i] <- cqi_cpos2str(cqp, Ncpos[1])
  }
  Ltt <- cbind(fuzzy, Leff, Lpos)
  Ltt <- Ltt[order(Ltt[,1]),]
  return (Ltt)
}
}

####################################################################

classrang <- function(corp, lm, c=2) {

# calcul de la classe, méthode Leipzig-Quasthoff
# a = effectif du premier
# b = effectif à classer
# formule : a = c^x.b 
# on utilise une commande externe bien plus efficace
# pour obtenir les effectifs.
# on utilise un log en base 2 par défaut,
# en augmentant cette valeur, on réduit le nb de classes

params1 <- paste("-f -P lemma ",corp,  sep="")
rang1 <-  system2(command="cwb-lexdecode", args=params1, stdout=TRUE) 
rang1b <- substr(rang1, 1, 7) # les effectifs sont sur les 7 premiers car.
rang1c <- rev(sort(as.numeric(rang1b)))

paramsN <- paste("-f -P lemma -p '",lm,"' ",corp, sep="")
rangN <-   system2(command="cwb-lexdecode", args=paramsN, stdout=TRUE) 
rangN1 <- as.numeric(substr(rangN, 1,7))

x <- (log(rang1c[1]/rangN1[1]))/log(c)
xmax <- (log(rang1c[1]/1))/log(c)
x <- round(x, 0)
xmax <- round(xmax, 0)
res <- c(x, xmax) # on renvoie un vecteur

return (res)
}

####################################################################"

VB0 <- function(corp, attr, num=0){

# création des objets nécessaires
# aux analyses globales du vocabulaire.
# si l'attribut est numérique,
# mettre num=1

library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
gcinfo(verbose=FALSE)

cqi_query(corp, "Tt", '[lemma=".*"]')
quer <- paste(corp, "Tt", sep=":")
dumprov <- cqi_dump_subcorpus(quer)
cps <- dumprov[,1]
rm(dumprov)
gc()

cplem <- paste(corp, "lemma", sep=".")
cpfrm <- paste(corp, "word", sep=".")
cppos <- paste(corp, "pos", sep=".")
cpattr <- paste(corp, attr, sep=".")

tabcp <- as.data.frame(cps)
tabcp[,2] <- cqi_cpos2id(cplem, cps)
tabcp[,3] <- cqi_cpos2id(cpfrm, cps)
tabcp[,4] <- cqi_cpos2id(cppos, cps)
tabcp[,5] <- cqi_cpos2str(cplem, cps)
tabcp[,6] <- cqi_cpos2str(cpfrm, cps)

iddt <- cqi_cpos2struc(cpattr, cps)
tabcp[,7] <- cqi_struc2str(cpattr, iddt)
if (num == 1) {
  tabcp[,7] <- as.numeric(tabcp[,7])
}
names(tabcp) <- c("cpos", "idlem", "idfrm", "idpos", "lem", "frm", "attr")

cat("\n", "Effectif (tokens) : ", nrow(tabcp), "\n")
cat("Nombre de lemmes : ", length(unique(tabcp[,2])), "\n")
cat("Nombre de formes : ", length(unique(tabcp[,3])), "\n\n")

class(tabcp) <- "VBASE"
return(tabcp)

}
##############################################

VBZ <- function(vbase, D="", F="", val="", hz="") {

# en entrée un objet créé par VBase()
# troncature éventuelle par les dates
# et/ou un autre attribut
# possibilité de tirage au hasard
# création des objets >> zipfR
# le paramètre hz est utilisé pour générer un tirage au hasard
# selon une distribution de poisson ; pour obtenir un effectif
# donné, procéder par essais successifs
# en gros, la proportion souhaitée (essayer entre 0.1 et 1)

library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
library(zipfR, quietly=TRUE, warn.conflicts=FALSE)
options(warn=-1)
gc()

if (!inherits(vbase, "VBASE")) stop("en entrée : un objet produit par VB0()", call.=F)
vbase <- unclass(vbase)
vbase <- as.data.frame(vbase)

if (D != "") {
  D <- as.numeric(D)
  vbase <- vbase[(vbase[,1]>=D),]
}
if (F != "") {
  F <- as.numeric(F)
  vbase <-vbase[(vbase[,1]<=F),]
}
if (val != "") {
  vbase <- vbase[(vbase[,7]==val),]
}

cposd <- vbase[1,1]
cposf <- vbase[(nrow(vbase)),1]

if (hz != "") {
  hz <- as.numeric(hz)
  rg <- cposf - cposd
  hzd <- round(rpois(rg, hz))
  vbase <- vbase[(hzd>0),]
}

Ltabtfl <- vec2tfl(vbase[,5])
Ltabvgc <- vec2vgc(vbase[,5])
Ltabspc <- vec2spc(vbase[,5])
Ftabtfl <- vec2tfl(vbase[,6])
Ftabvgc <- vec2vgc(vbase[,6])
Ftabspc <- vec2spc(vbase[,6])

cat("\n", "Nombre de tokens retenus : ", nrow(vbase))

vb1 <- list(cposd = cposd, cposf=cposf, Ltfl=Ltabtfl, Lvgc=Ltabvgc, Lspc=Ltabspc, Ftfl=Ftabtfl, Fvgc=Ftabvgc, Fspc=Ftabspc)
gc()
return(vb1)
}

#####################################################

#### fonction de dessin de courbes splines ouvertes

"spllines" <-
function(x, y = NULL, n = 20, type="l", col=par("col"), lty=par("lty"), ...  ){

   xy <- cbind(x,y)
   nP <- length(xy[,1])
   if((nP < 3) || (is.na(any(xy[,1]))))
       return(list(x=numeric(0), y=numeric(0)))

  ## else :

   nP <- length(xy[,1])
   z <- n*(nP-1)

   i <- 1:nP
   Sx <- splinefun(i, xy[,1], method = "natural")
   Sy <- splinefun(i, xy[,2], method = "natural")
   ti <- seq(1, nP, length = z)
   opspl <- cbind(x = Sx(ti), y = Sy(ti))

   lines(opspl, type=type, col=col, lty=lty, ...  )
}
