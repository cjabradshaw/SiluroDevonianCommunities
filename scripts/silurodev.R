## Siluro-Devonian vertebrate/invertebrate/plant geological determinants of presence
## Corey Bradshaw & Richard Cloutier
## May 2024

rm(list = ls())

library(performance)
library(sjPlot)
library(lme4)
library(sjPlot)
library(ggplot2)
library(stringr)
library(vegan)
library(parallel)
library(vcdExtra)
library(plyr)
library(dplyr)

cl <- detectCores() - 1

# source files
source("new_lmer_AIC_tables3.R")
source("r.squared.R")

# functions
AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

# import data
dat <- read.table("sildevcomp2.txt", sep="\t", header=T, as.is=T)
colnames(dat)
head(dat)

###############
## vertebrates 
###############
vertcolnames <- colnames(dat[,c(34,35,37,39,41,43,45,47,49,51)])
vertcolnames
vertcols <- which(colnames(dat) %in% vertcolnames)
dat$vertcomb <- ifelse(apply(dat[,vertcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

vertNA <- ifelse(apply(is.na(dat[,vertcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(vertNA, MARGIN=2, sum)

# conglomerate
vertcong <- xtabs(dat$vertcomb ~ dat$Conglomerate)
vertcong

# sandstone
vertsand <- xtabs(dat$vertcomb ~ dat$Sandstone)
vertsand

# siltstone
vertsilt <- xtabs(dat$vertcomb ~ dat$Siltstone)
vertsilt

# shale
vertshal <- xtabs(dat$vertcomb ~ dat$Shale)
vertshal

# mudstone
vertmud <- xtabs(dat$vertcomb ~ dat$Mudstone)
vertmud

# limestone
vertlime <- xtabs(dat$vertcomb ~ dat$Limestone)
vertlime

vert.summ <- as.data.frame(rbind(vertcong,vertsand,vertsilt,vertshal,vertmud,vertlime))
stone.labs <- c("conglomerate","sandstone","siltstone","shale","mudstone","limestone")
rownames(vert.summ) <- stone.labs
colnames(vert.summ) <- c("abs","pres")
vert.summ
vert.summ$ppres <- vert.summ$pres / apply(vert.summ, MARGIN=1, sum)
vert.summ

ggplot(vert.summ, aes(x=rownames(vert.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")


#################
## invertebrates 
#################
invcolnames <- colnames(dat[,c(59:78)])
invcolnames
invcols <- which(colnames(dat) %in% invcolnames)
dat$invcomb <- ifelse(apply(dat[,invcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

invNA <- ifelse(apply(is.na(dat[,invcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(invNA, MARGIN=2, sum)

# conglomerate
invcong <- xtabs(dat$invcomb ~ dat$Conglomerate)
invcong

# sandstone
invsand <- xtabs(dat$invcomb ~ dat$Sandstone)
invsand

# siltstone
invsilt <- xtabs(dat$invcomb ~ dat$Siltstone)
invsilt

# shale
invshal <- xtabs(dat$invcomb ~ dat$Shale)
invshal

# mudstone
invmud <- xtabs(dat$invcomb ~ dat$Mudstone)
invmud

# limestone
invlime <- xtabs(dat$invcomb ~ dat$Limestone)
invlime

inv.summ <- as.data.frame(rbind(invcong,invsand,invsilt,invshal,invmud,invlime))
rownames(inv.summ) <- stone.labs
colnames(inv.summ) <- c("abs","pres")
inv.summ
inv.summ$ppres <- inv.summ$pres / apply(inv.summ, MARGIN=1, sum)
inv.summ

ggplot(inv.summ, aes(x=rownames(inv.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")


##########
## plants 
##########
plntcolnames <- colnames(dat[,c(80:87)])
plntcolnames
plntcols <- which(colnames(dat) %in% plntcolnames)
dat$plntcomb <- ifelse(apply(dat[,plntcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

plntNA <- ifelse(apply(is.na(dat[,plntcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(plntNA, MARGIN=2, sum)

# conglomerate
plntcong <- xtabs(dat$plntcomb ~ dat$Conglomerate)
plntcong

# sandstone
plntsand <- xtabs(dat$plntcomb ~ dat$Sandstone)
plntsand

# siltstone
plntsilt <- xtabs(dat$plntcomb ~ dat$Siltstone)
plntsilt

# shale
plntshal <- xtabs(dat$plntcomb ~ dat$Shale)
plntshal

# mudstone
plntmud <- xtabs(dat$plntcomb ~ dat$Mudstone)
plntmud

# limestone
plntlime <- xtabs(dat$plntcomb ~ dat$Limestone)
plntlime

plnt.summ <- as.data.frame(rbind(plntcong,plntsand,plntsilt,plntshal,plntmud,plntlime))
rownames(plnt.summ) <- stone.labs
colnames(plnt.summ) <- c("abs","pres")
plnt.summ
plnt.summ$ppres <- plnt.summ$pres / apply(plnt.summ, MARGIN=1, sum)
plnt.summ

ggplot(plnt.summ, aes(x=rownames(plnt.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

## combined data
vert.summ.pres <- data.frame("geol"=rownames(vert.summ), "pres"="pres", "freq"=vert.summ$pres)
vert.summ.abs <- data.frame("geol"=rownames(vert.summ), "pres"="abs", "freq"=vert.summ$abs)
vert.summ2 <- (rbind(vert.summ.pres, vert.summ.abs))
vert.summ2$tax <- "vert"

inv.summ.pres <- data.frame("geol"=rownames(inv.summ), "pres"="pres", "freq"=inv.summ$pres)
inv.summ.abs <- data.frame("geol"=rownames(inv.summ), "pres"="abs", "freq"=inv.summ$abs)
inv.summ2 <- (rbind(inv.summ.pres, inv.summ.abs))
inv.summ2$tax <- "inv"

plnt.summ.pres <- data.frame("geol"=rownames(plnt.summ), "pres"="pres", "freq"=plnt.summ$pres)
plnt.summ.abs <- data.frame("geol"=rownames(plnt.summ), "pres"="abs", "freq"=plnt.summ$abs)
plnt.summ2 <- (rbind(plnt.summ.pres, plnt.summ.abs))
plnt.summ2$tax <- "plnt"

comb.summ <- rbind(vert.summ2, inv.summ2, plnt.summ2)
comb.exp <- expand.dft(comb.summ, freq="freq")
comb.tab <- table(comb.exp$tax, comb.exp$geol, comb.exp$pres)
plot(comb.tab)

# mosaic
comb.struct <- structable(data=comb.tab)
mosaic(comb.struct)
mosaic(comb.struct, type="expected")
assoc(comb.struct, compress=F)

# logistic weighted models
# prepare data
comb.exp$pres1 <- as.integer(ifelse(comb.exp$pres == "pres", 1, 0))
comb.exp$geol <- as.factor(comb.exp$geol)
comb.exp$tax <- as.factor(comb.exp$tax)
str(comb.exp)
table(comb.exp$pres1)

# model set
#m1 <- "pres1 ~ geol + tax + geol*tax"
m2 <- "pres1 ~ geol + tax"
m3 <- "pres1 ~ geol"
m4 <- "pres1 ~ tax"
m5 <- "pres1 ~ 1"

## model vector
mod.vec <- c(m2,m3,m4,m5)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=comb.exp, family=binomial, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=T),]
summary.table

ERgeoltax <- round(sumtable[2,5]/sumtable[3,5], 3)
ERgeoltax


## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), data=comb.exp, na.action=na.omit)

check_model(fit)
plot_model(fit, show.values=T, vline.color = "purple")


## split by major period(Silurian, lower, mid, upper Devonian)
sil.ind <- grep("silur", dat$Age, ignore.case=T)
dat[sil.ind,]$Age
dat2 <- dat[-sil.ind,]
earldev.ind <- grep("lower", dat2$Age, ignore.case=T)
dat2[earldev.ind,]$Age
dat3 <- dat2[-earldev.ind,]
middev.ind <- grep("mid", dat3$Age, ignore.case=T)
dat3[middev.ind,]$Age
dat4 <- dat3[-middev.ind,]
latedev.ind <- grep("upp", dat4$Age, ignore.case=T)
dat4[latedev.ind,]$Age


############
## Silurian
############
Sil.dat <- dat[sil.ind,]
dim(Sil.dat)

###############
## vertebrates 
###############
vertcolnames <- colnames(Sil.dat[,c(34,35,37,39,41,43,45,47,49,51)])
vertcolnames
vertcols <- which(colnames(Sil.dat) %in% vertcolnames)
Sil.dat$vertcomb <- ifelse(apply(Sil.dat[,vertcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

vertNA <- ifelse(apply(is.na(Sil.dat[,vertcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(vertNA, MARGIN=2, sum)

# conglomerate
vertcong <- xtabs(Sil.dat$vertcomb ~ Sil.dat$Conglomerate)
vertcong

# sandstone
vertsand <- xtabs(Sil.dat$vertcomb ~ Sil.dat$Sandstone)
vertsand

# siltstone
vertsilt <- xtabs(Sil.dat$vertcomb ~ Sil.dat$Siltstone)
vertsilt

# shale
vertshal <- xtabs(Sil.dat$vertcomb ~ Sil.dat$Shale)
vertshal

# mudstone
vertmud <- xtabs(Sil.dat$vertcomb ~ Sil.dat$Mudstone)
vertmud

# limestone
vertlime <- xtabs(Sil.dat$vertcomb ~ Sil.dat$Limestone)
vertlime

Sil.vert.summ <- as.data.frame(rbind(vertcong,vertsand,vertsilt,vertshal,vertmud,vertlime))
rownames(Sil.vert.summ) <- stone.labs
colnames(Sil.vert.summ) <- c("abs","pres")
Sil.vert.summ
Sil.vert.summ$ppres <- Sil.vert.summ$pres / apply(Sil.vert.summ, MARGIN=1, sum)
Sil.vert.summ

ggplot(Sil.vert.summ, aes(x=rownames(Sil.vert.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

#################
## invertebrates 
#################
invcolnames <- colnames(Sil.dat[,c(59:78)])
invcolnames
invcols <- which(colnames(Sil.dat) %in% invcolnames)
Sil.dat$invcomb <- ifelse(apply(Sil.dat[,invcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

invNA <- ifelse(apply(is.na(Sil.dat[,invcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(invNA, MARGIN=2, sum)

# conglomerate
invcong <- xtabs(Sil.dat$invcomb ~ Sil.dat$Conglomerate)
invcong

# sandstone
invsand <- xtabs(Sil.dat$invcomb ~ Sil.dat$Sandstone)
invsand

# siltstone
invsilt <- xtabs(Sil.dat$invcomb ~ Sil.dat$Siltstone)
invsilt

# shale
invshal <- xtabs(Sil.dat$invcomb ~ Sil.dat$Shale)
invshal

# mudstone
invmud <- xtabs(Sil.dat$invcomb ~ Sil.dat$Mudstone)
invmud

# limestone
invlime <- xtabs(Sil.dat$invcomb ~ Sil.dat$Limestone)
invlime

Sil.inv.summ <- as.data.frame(rbind(invcong,invsand,invsilt,invshal,invmud,invlime))
rownames(Sil.inv.summ) <- stone.labs
colnames(Sil.inv.summ) <- c("abs","pres")
Sil.inv.summ
Sil.inv.summ$ppres <- Sil.inv.summ$pres / apply(Sil.inv.summ, MARGIN=1, sum)
Sil.inv.summ

ggplot(Sil.inv.summ, aes(x=rownames(Sil.inv.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

##########
## plants 
##########
plntcolnames <- colnames(Sil.dat[,c(80:87)])
plntcolnames
plntcols <- which(colnames(Sil.dat) %in% plntcolnames)
Sil.dat$plntcomb <- ifelse(apply(Sil.dat[,plntcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

plntNA <- ifelse(apply(is.na(Sil.dat[,plntcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(plntNA, MARGIN=2, sum)

# conglomerate
plntcong <- xtabs(Sil.dat$plntcomb ~ Sil.dat$Conglomerate)
plntcong

# sandstone
plntsand <- xtabs(Sil.dat$plntcomb ~ Sil.dat$Sandstone)
plntsand

# siltstone
plntsilt <- xtabs(Sil.dat$plntcomb ~ Sil.dat$Siltstone)
plntsilt

# shale
plntshal <- xtabs(Sil.dat$plntcomb ~ Sil.dat$Shale)
plntshal

# mudstone
plntmud <- xtabs(Sil.dat$plntcomb ~ Sil.dat$Mudstone)
plntmud

# limestone
plntlime <- xtabs(Sil.dat$plntcomb ~ Sil.dat$Limestone)
plntlime

Sil.plnt.summ <- as.data.frame(rbind(plntcong,plntsand,plntsilt,plntshal,plntmud,plntlime))
rownames(Sil.plnt.summ) <- stone.labs
colnames(Sil.plnt.summ) <- c("abs","pres")
Sil.plnt.summ
Sil.plnt.summ$ppres <- Sil.plnt.summ$pres / apply(Sil.plnt.summ, MARGIN=1, sum)
Sil.plnt.summ

ggplot(Sil.plnt.summ, aes(x=rownames(Sil.plnt.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

## combined data
Sil.vert.summ.pres <- data.frame("geol"=rownames(Sil.vert.summ), "pres"="pres", "freq"=Sil.vert.summ$pres)
Sil.vert.summ.abs <- data.frame("geol"=rownames(Sil.vert.summ), "pres"="abs", "freq"=Sil.vert.summ$abs)
Sil.vert.summ2 <- (rbind(Sil.vert.summ.pres, Sil.vert.summ.abs))
Sil.vert.summ2$tax <- "vert"

Sil.inv.summ.pres <- data.frame("geol"=rownames(Sil.inv.summ), "pres"="pres", "freq"=Sil.inv.summ$pres)
Sil.inv.summ.abs <- data.frame("geol"=rownames(Sil.inv.summ), "pres"="abs", "freq"=Sil.inv.summ$abs)
Sil.inv.summ2 <- (rbind(Sil.inv.summ.pres, Sil.inv.summ.abs))
Sil.inv.summ2$tax <- "inv"

Sil.plnt.summ.pres <- data.frame("geol"=rownames(Sil.plnt.summ), "pres"="pres", "freq"=Sil.plnt.summ$pres)
Sil.plnt.summ.abs <- data.frame("geol"=rownames(Sil.plnt.summ), "pres"="abs", "freq"=Sil.plnt.summ$abs)
Sil.plnt.summ2 <- (rbind(Sil.plnt.summ.pres, Sil.plnt.summ.abs))
Sil.plnt.summ2$tax <- "plnt"

Sil.comb.summ <- rbind(Sil.vert.summ2, Sil.inv.summ2, Sil.plnt.summ2)
Sil.comb.exp <- expand.dft(Sil.comb.summ, freq="freq")
Sil.comb.tab <- table(Sil.comb.exp$tax, Sil.comb.exp$geol, Sil.comb.exp$pres)
plot(Sil.comb.tab)

# mosaic
Sil.comb.struct <- structable(data=Sil.comb.tab)
mosaic(Sil.comb.struct)
mosaic(Sil.comb.struct, type="expected")
assoc(Sil.comb.struct, compress=F)

# logistic weighted models
# prepare data
Sil.comb.exp$pres1 <- as.integer(ifelse(Sil.comb.exp$pres == "pres", 1, 0))
Sil.comb.exp$geol <- as.factor(Sil.comb.exp$geol)
Sil.comb.exp$tax <- as.factor(Sil.comb.exp$tax)
str(Sil.comb.exp)
table(Sil.comb.exp$pres1)

## model vector
mod.vec <- c(m2,m3,m4,m5)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=Sil.comb.exp, family=binomial, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

Sil.sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(Sil.sumtable) <- mod.vec
Sil.summary.table <- Sil.sumtable[order(Sil.sumtable[,5],decreasing=T),]
Sil.summary.table

ERgeoltax <- round(Sil.sumtable[2,5]/Sil.sumtable[3,5], 3)
ERgeoltax

## saturated residual diagnostic
i <- 2
Sil.fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), data=Sil.comb.exp, na.action=na.omit)

check_model(Sil.fit)
plot_model(Sil.fit, show.values=T, vline.color = "purple")


########################
## Lower/Early Devonian
########################
EarlyDev.dat <- dat2[earldev.ind,]
dim(EarlyDev.dat)

###############
## vertebrates 
###############
vertcolnames <- colnames(EarlyDev.dat[,c(34,35,37,39,41,43,45,47,49,51)])
vertcolnames
vertcols <- which(colnames(EarlyDev.dat) %in% vertcolnames)
EarlyDev.dat$vertcomb <- ifelse(apply(EarlyDev.dat[,vertcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

vertNA <- ifelse(apply(is.na(EarlyDev.dat[,vertcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(vertNA, MARGIN=2, sum)

# conglomerate
vertcong <- xtabs(EarlyDev.dat$vertcomb ~ EarlyDev.dat$Conglomerate)
vertcong

# sandstone
vertsand <- xtabs(EarlyDev.dat$vertcomb ~ EarlyDev.dat$Sandstone)
vertsand

# siltstone
vertsilt <- xtabs(EarlyDev.dat$vertcomb ~ EarlyDev.dat$Siltstone)
vertsilt

# shale
vertshal <- xtabs(EarlyDev.dat$vertcomb ~ EarlyDev.dat$Shale)
vertshal

# mudstone
vertmud <- xtabs(EarlyDev.dat$vertcomb ~ EarlyDev.dat$Mudstone)
vertmud

# limestone
vertlime <- xtabs(EarlyDev.dat$vertcomb ~ EarlyDev.dat$Limestone)
vertlime

EarlyDev.vert.summ <- as.data.frame(rbind(vertcong,vertsand,vertsilt,vertshal,vertmud,vertlime))
rownames(EarlyDev.vert.summ) <- stone.labs
colnames(EarlyDev.vert.summ) <- c("abs","pres")
EarlyDev.vert.summ
EarlyDev.vert.summ$ppres <- EarlyDev.vert.summ$pres / apply(EarlyDev.vert.summ, MARGIN=1, sum)
EarlyDev.vert.summ

ggplot(EarlyDev.vert.summ, aes(x=rownames(EarlyDev.vert.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

#################
## invertebrates 
#################
invcolnames <- colnames(EarlyDev.dat[,c(59:78)])
invcolnames
invcols <- which(colnames(EarlyDev.dat) %in% invcolnames)
EarlyDev.dat$invcomb <- ifelse(apply(EarlyDev.dat[,invcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

invNA <- ifelse(apply(is.na(EarlyDev.dat[,invcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(invNA, MARGIN=2, sum)

# conglomerate
invcong <- xtabs(EarlyDev.dat$invcomb ~ EarlyDev.dat$Conglomerate)
invcong

# sandstone
invsand <- xtabs(EarlyDev.dat$invcomb ~ EarlyDev.dat$Sandstone)
invsand

# siltstone
invsilt <- xtabs(EarlyDev.dat$invcomb ~ EarlyDev.dat$Siltstone)
invsilt

# shale
invshal <- xtabs(EarlyDev.dat$invcomb ~ EarlyDev.dat$Shale)
invshal

# mudstone
invmud <- xtabs(EarlyDev.dat$invcomb ~ EarlyDev.dat$Mudstone)
invmud

# limestone
invlime <- xtabs(EarlyDev.dat$invcomb ~ EarlyDev.dat$Limestone)
invlime

EarlyDev.inv.summ <- as.data.frame(rbind(invcong,invsand,invsilt,invshal,invmud,invlime))
rownames(EarlyDev.inv.summ) <- stone.labs
colnames(EarlyDev.inv.summ) <- c("abs","pres")
EarlyDev.inv.summ
EarlyDev.inv.summ$ppres <- EarlyDev.inv.summ$pres / apply(EarlyDev.inv.summ, MARGIN=1, sum)
EarlyDev.inv.summ

ggplot(EarlyDev.inv.summ, aes(x=rownames(EarlyDev.inv.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

##########
## plants 
##########
plntcolnames <- colnames(EarlyDev.dat[,c(80:87)])
plntcolnames
plntcols <- which(colnames(EarlyDev.dat) %in% plntcolnames)
EarlyDev.dat$plntcomb <- ifelse(apply(EarlyDev.dat[,plntcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

plntNA <- ifelse(apply(is.na(EarlyDev.dat[,plntcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(plntNA, MARGIN=2, sum)

# conglomerate
plntcong <- xtabs(EarlyDev.dat$plntcomb ~ EarlyDev.dat$Conglomerate)
plntcong

# sandstone
plntsand <- xtabs(EarlyDev.dat$plntcomb ~ EarlyDev.dat$Sandstone)
plntsand

# siltstone
plntsilt <- xtabs(EarlyDev.dat$plntcomb ~ EarlyDev.dat$Siltstone)
plntsilt

# shale
plntshal <- xtabs(EarlyDev.dat$plntcomb ~ EarlyDev.dat$Shale)
plntshal

# mudstone
plntmud <- xtabs(EarlyDev.dat$plntcomb ~ EarlyDev.dat$Mudstone)
plntmud

# limestone
plntlime <- xtabs(EarlyDev.dat$plntcomb ~ EarlyDev.dat$Limestone)
plntlime

EarlyDev.plnt.summ <- as.data.frame(rbind(plntcong,plntsand,plntsilt,plntshal,plntmud,plntlime))
rownames(EarlyDev.plnt.summ) <- stone.labs
colnames(EarlyDev.plnt.summ) <- c("abs","pres")
EarlyDev.plnt.summ
EarlyDev.plnt.summ$ppres <- EarlyDev.plnt.summ$pres / apply(EarlyDev.plnt.summ, MARGIN=1, sum)
EarlyDev.plnt.summ

ggplot(EarlyDev.plnt.summ, aes(x=rownames(EarlyDev.plnt.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

## combined data
EarlyDev.vert.summ.pres <- data.frame("geol"=rownames(EarlyDev.vert.summ), "pres"="pres", "freq"=EarlyDev.vert.summ$pres)
EarlyDev.vert.summ.abs <- data.frame("geol"=rownames(EarlyDev.vert.summ), "pres"="abs", "freq"=EarlyDev.vert.summ$abs)
EarlyDev.vert.summ2 <- (rbind(EarlyDev.vert.summ.pres, EarlyDev.vert.summ.abs))
EarlyDev.vert.summ2$tax <- "vert"

EarlyDev.inv.summ.pres <- data.frame("geol"=rownames(EarlyDev.inv.summ), "pres"="pres", "freq"=EarlyDev.inv.summ$pres)
EarlyDev.inv.summ.abs <- data.frame("geol"=rownames(EarlyDev.inv.summ), "pres"="abs", "freq"=EarlyDev.inv.summ$abs)
EarlyDev.inv.summ2 <- (rbind(EarlyDev.inv.summ.pres, EarlyDev.inv.summ.abs))
EarlyDev.inv.summ2$tax <- "inv"

EarlyDev.plnt.summ.pres <- data.frame("geol"=rownames(EarlyDev.plnt.summ), "pres"="pres", "freq"=EarlyDev.plnt.summ$pres)
EarlyDev.plnt.summ.abs <- data.frame("geol"=rownames(EarlyDev.plnt.summ), "pres"="abs", "freq"=EarlyDev.plnt.summ$abs)
EarlyDev.plnt.summ2 <- (rbind(EarlyDev.plnt.summ.pres, EarlyDev.plnt.summ.abs))
EarlyDev.plnt.summ2$tax <- "plnt"

EarlyDev.comb.summ <- rbind(EarlyDev.vert.summ2, EarlyDev.inv.summ2, EarlyDev.plnt.summ2)
EarlyDev.comb.exp <- expand.dft(EarlyDev.comb.summ, freq="freq")
EarlyDev.comb.tab <- table(EarlyDev.comb.exp$tax, EarlyDev.comb.exp$geol, EarlyDev.comb.exp$pres)
plot(EarlyDev.comb.tab)

# mosaic
EarlyDev.comb.struct <- structable(data=EarlyDev.comb.tab)
mosaic(EarlyDev.comb.struct)
mosaic(EarlyDev.comb.struct, type="expected")
assoc(EarlyDev.comb.struct, compress=F)

# logistic weighted models
# prepare data
EarlyDev.comb.exp$pres1 <- as.integer(ifelse(EarlyDev.comb.exp$pres == "pres", 1, 0))
EarlyDev.comb.exp$geol <- as.factor(EarlyDev.comb.exp$geol)
EarlyDev.comb.exp$tax <- as.factor(EarlyDev.comb.exp$tax)
str(EarlyDev.comb.exp)
table(EarlyDev.comb.exp$pres1)

## model vector
mod.vec <- c(m2,m3,m4,m5)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=EarlyDev.comb.exp, family=binomial, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

EarlyDev.sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(EarlyDev.sumtable) <- mod.vec
EarlyDev.summary.table <- EarlyDev.sumtable[order(EarlyDev.sumtable[,5],decreasing=T),]
EarlyDev.summary.table

ERgeoltax <- round(EarlyDev.sumtable[2,5]/EarlyDev.sumtable[3,5], 3)
ERgeoltax

## saturated residual diagnostic
i <- 2
EarlyDev.fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), data=EarlyDev.comb.exp, na.action=na.omit)

check_model(EarlyDev.fit)
plot_model(EarlyDev.fit, show.values=T, vline.color = "purple")


###############
# Mid Devonian
###############
MidDev.dat <- dat3[middev.ind,]
dim(MidDev.dat)

###############
## vertebrates 
###############
vertcolnames <- colnames(MidDev.dat[,c(34,35,37,39,41,43,45,47,49,51)])
vertcolnames
vertcols <- which(colnames(MidDev.dat) %in% vertcolnames)
MidDev.dat$vertcomb <- ifelse(apply(MidDev.dat[,vertcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

vertNA <- ifelse(apply(is.na(MidDev.dat[,vertcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(vertNA, MARGIN=2, sum)

# conglomerate
vertcong <- xtabs(MidDev.dat$vertcomb ~ MidDev.dat$Conglomerate)
vertcong

# sandstone
vertsand <- xtabs(MidDev.dat$vertcomb ~ MidDev.dat$Sandstone)
vertsand

# siltstone
vertsilt <- xtabs(MidDev.dat$vertcomb ~ MidDev.dat$Siltstone)
vertsilt

# shale
vertshal <- xtabs(MidDev.dat$vertcomb ~ MidDev.dat$Shale)
vertshal

# mudstone
vertmud <- xtabs(MidDev.dat$vertcomb ~ MidDev.dat$Mudstone)
vertmud

# limestone
vertlime <- xtabs(MidDev.dat$vertcomb ~ MidDev.dat$Limestone)
vertlime

MidDev.vert.summ <- as.data.frame(rbind(vertcong,vertsand,vertsilt,vertshal,vertmud,vertlime))
rownames(MidDev.vert.summ) <- stone.labs
colnames(MidDev.vert.summ) <- c("abs","pres")
MidDev.vert.summ
MidDev.vert.summ$ppres <- MidDev.vert.summ$pres / apply(MidDev.vert.summ, MARGIN=1, sum)
MidDev.vert.summ

ggplot(MidDev.vert.summ, aes(x=rownames(MidDev.vert.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

#################
## invertebrates 
#################
invcolnames <- colnames(MidDev.dat[,c(59:78)])
invcolnames
invcols <- which(colnames(MidDev.dat) %in% invcolnames)
MidDev.dat$invcomb <- ifelse(apply(MidDev.dat[,invcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

invNA <- ifelse(apply(is.na(MidDev.dat[,invcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(invNA, MARGIN=2, sum)

# conglomerate
invcong <- xtabs(MidDev.dat$invcomb ~ MidDev.dat$Conglomerate)
invcong

# sandstone
invsand <- xtabs(MidDev.dat$invcomb ~ MidDev.dat$Sandstone)
invsand

# siltstone
invsilt <- xtabs(MidDev.dat$invcomb ~ MidDev.dat$Siltstone)
invsilt

# shale
invshal <- xtabs(MidDev.dat$invcomb ~ MidDev.dat$Shale)
invshal

# mudstone
invmud <- xtabs(MidDev.dat$invcomb ~ MidDev.dat$Mudstone)
invmud

# limestone
invlime <- xtabs(MidDev.dat$invcomb ~ MidDev.dat$Limestone)
invlime

MidDev.inv.summ <- as.data.frame(rbind(invcong,invsand,invsilt,invshal,invmud,invlime))
rownames(MidDev.inv.summ) <- stone.labs
colnames(MidDev.inv.summ) <- c("abs","pres")
MidDev.inv.summ
MidDev.inv.summ$ppres <- MidDev.inv.summ$pres / apply(MidDev.inv.summ, MARGIN=1, sum)
MidDev.inv.summ

ggplot(MidDev.inv.summ, aes(x=rownames(MidDev.inv.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

##########
## plants 
##########
plntcolnames <- colnames(MidDev.dat[,c(80:87)])
plntcolnames
plntcols <- which(colnames(MidDev.dat) %in% plntcolnames)
MidDev.dat$plntcomb <- ifelse(apply(MidDev.dat[,plntcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

plntNA <- ifelse(apply(is.na(MidDev.dat[,plntcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(plntNA, MARGIN=2, sum)

# conglomerate
plntcong <- xtabs(MidDev.dat$plntcomb ~ MidDev.dat$Conglomerate)
plntcong

# sandstone
plntsand <- xtabs(MidDev.dat$plntcomb ~ MidDev.dat$Sandstone)
plntsand

# siltstone
plntsilt <- xtabs(MidDev.dat$plntcomb ~ MidDev.dat$Siltstone)
plntsilt

# shale
plntshal <- xtabs(MidDev.dat$plntcomb ~ MidDev.dat$Shale)
plntshal

# mudstone
plntmud <- xtabs(MidDev.dat$plntcomb ~ MidDev.dat$Mudstone)
plntmud

# limestone
plntlime <- xtabs(MidDev.dat$plntcomb ~ MidDev.dat$Limestone)
plntlime

MidDev.plnt.summ <- as.data.frame(rbind(plntcong,plntsand,plntsilt,plntshal,plntmud,plntlime))
rownames(MidDev.plnt.summ) <- stone.labs
colnames(MidDev.plnt.summ) <- c("abs","pres")
MidDev.plnt.summ
MidDev.plnt.summ$ppres <- MidDev.plnt.summ$pres / apply(MidDev.plnt.summ, MARGIN=1, sum)
MidDev.plnt.summ

ggplot(MidDev.plnt.summ, aes(x=rownames(MidDev.plnt.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

## combined data
MidDev.vert.summ.pres <- data.frame("geol"=rownames(MidDev.vert.summ), "pres"="pres", "freq"=MidDev.vert.summ$pres)
MidDev.vert.summ.abs <- data.frame("geol"=rownames(MidDev.vert.summ), "pres"="abs", "freq"=MidDev.vert.summ$abs)
MidDev.vert.summ2 <- (rbind(MidDev.vert.summ.pres, MidDev.vert.summ.abs))
MidDev.vert.summ2$tax <- "vert"

MidDev.inv.summ.pres <- data.frame("geol"=rownames(MidDev.inv.summ), "pres"="pres", "freq"=MidDev.inv.summ$pres)
MidDev.inv.summ.abs <- data.frame("geol"=rownames(MidDev.inv.summ), "pres"="abs", "freq"=MidDev.inv.summ$abs)
MidDev.inv.summ2 <- (rbind(MidDev.inv.summ.pres, MidDev.inv.summ.abs))
MidDev.inv.summ2$tax <- "inv"

MidDev.plnt.summ.pres <- data.frame("geol"=rownames(MidDev.plnt.summ), "pres"="pres", "freq"=MidDev.plnt.summ$pres)
MidDev.plnt.summ.abs <- data.frame("geol"=rownames(MidDev.plnt.summ), "pres"="abs", "freq"=MidDev.plnt.summ$abs)
MidDev.plnt.summ2 <- (rbind(MidDev.plnt.summ.pres, MidDev.plnt.summ.abs))
MidDev.plnt.summ2$tax <- "plnt"

MidDev.comb.summ <- rbind(MidDev.vert.summ2, MidDev.inv.summ2, MidDev.plnt.summ2)
MidDev.comb.exp <- expand.dft(MidDev.comb.summ, freq="freq")
MidDev.comb.tab <- table(MidDev.comb.exp$tax, MidDev.comb.exp$geol, MidDev.comb.exp$pres)
plot(MidDev.comb.tab)

# mosaic
MidDev.comb.struct <- structable(data=MidDev.comb.tab)
mosaic(MidDev.comb.struct)
mosaic(MidDev.comb.struct, type="expected")
assoc(MidDev.comb.struct, compress=F)

# logistic weighted models
# prepare data
MidDev.comb.exp$pres1 <- as.integer(ifelse(MidDev.comb.exp$pres == "pres", 1, 0))
MidDev.comb.exp$geol <- as.factor(MidDev.comb.exp$geol)
MidDev.comb.exp$tax <- as.factor(MidDev.comb.exp$tax)
str(MidDev.comb.exp)
table(MidDev.comb.exp$pres1)

## model vector
mod.vec <- c(m2,m3,m4,m5)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=MidDev.comb.exp, family=binomial, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

MidDev.sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(MidDev.sumtable) <- mod.vec
MidDev.summary.table <- MidDev.sumtable[order(MidDev.sumtable[,5],decreasing=T),]
MidDev.summary.table

ERgeoltax <- round(MidDev.sumtable[2,5]/MidDev.sumtable[3,5], 3)
ERgeoltax

## saturated residual diagnostic
i <- 2
MidDev.fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), data=MidDev.comb.exp, na.action=na.omit)

check_model(MidDev.fit)
plot_model(MidDev.fit, show.values=T, vline.color = "purple")


################
# Late Devonian
################
LateDev.dat <- dat4[latedev.ind,]
dim(LateDev.dat)

###############
## vertebrates 
###############
vertcolnames <- colnames(LateDev.dat[,c(34,35,37,39,41,43,45,47,49,51)])
vertcolnames
vertcols <- which(colnames(LateDev.dat) %in% vertcolnames)
LateDev.dat$vertcomb <- ifelse(apply(LateDev.dat[,vertcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

vertNA <- ifelse(apply(is.na(LateDev.dat[,vertcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(vertNA, MARGIN=2, sum)

# conglomerate
vertcong <- xtabs(LateDev.dat$vertcomb ~ LateDev.dat$Conglomerate)
vertcong

# sandstone
vertsand <- xtabs(LateDev.dat$vertcomb ~ LateDev.dat$Sandstone)
vertsand

# siltstone
vertsilt <- xtabs(LateDev.dat$vertcomb ~ LateDev.dat$Siltstone)
vertsilt

# shale
vertshal <- xtabs(LateDev.dat$vertcomb ~ LateDev.dat$Shale)
vertshal

# mudstone
vertmud <- xtabs(LateDev.dat$vertcomb ~ LateDev.dat$Mudstone)
vertmud

# limestone
vertlime <- xtabs(LateDev.dat$vertcomb ~ LateDev.dat$Limestone)
vertlime

LateDev.vert.summ <- as.data.frame(rbind(vertcong,vertsand,vertsilt,vertshal,vertmud,vertlime))
rownames(LateDev.vert.summ) <- stone.labs
colnames(LateDev.vert.summ) <- c("abs","pres")
LateDev.vert.summ
LateDev.vert.summ$ppres <- LateDev.vert.summ$pres / apply(LateDev.vert.summ, MARGIN=1, sum)
LateDev.vert.summ

ggplot(LateDev.vert.summ, aes(x=rownames(LateDev.vert.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

#################
## invertebrates 
#################
invcolnames <- colnames(LateDev.dat[,c(59:78)])
invcolnames
invcols <- which(colnames(LateDev.dat) %in% invcolnames)
LateDev.dat$invcomb <- ifelse(apply(LateDev.dat[,invcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

invNA <- ifelse(apply(is.na(LateDev.dat[,invcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(invNA, MARGIN=2, sum)

# conglomerate
invcong <- xtabs(LateDev.dat$invcomb ~ LateDev.dat$Conglomerate)
invcong

# sandstone
invsand <- xtabs(LateDev.dat$invcomb ~ LateDev.dat$Sandstone)
invsand

# siltstone
invsilt <- xtabs(LateDev.dat$invcomb ~ LateDev.dat$Siltstone)
invsilt

# shale
invshal <- xtabs(LateDev.dat$invcomb ~ LateDev.dat$Shale)
invshal

# mudstone
invmud <- xtabs(LateDev.dat$invcomb ~ LateDev.dat$Mudstone)
invmud

# limestone
invlime <- xtabs(LateDev.dat$invcomb ~ LateDev.dat$Limestone)
invlime

LateDev.inv.summ <- as.data.frame(rbind(invcong,invsand,invsilt,invshal,invmud,invlime))
rownames(LateDev.inv.summ) <- stone.labs
colnames(LateDev.inv.summ) <- c("abs","pres")
LateDev.inv.summ
LateDev.inv.summ$ppres <- LateDev.inv.summ$pres / apply(LateDev.inv.summ, MARGIN=1, sum)
LateDev.inv.summ

ggplot(LateDev.inv.summ, aes(x=rownames(LateDev.inv.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

##########
## plants 
##########
plntcolnames <- colnames(LateDev.dat[,c(80:87)])
plntcolnames
plntcols <- which(colnames(LateDev.dat) %in% plntcolnames)
LateDev.dat$plntcomb <- ifelse(apply(LateDev.dat[,plntcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

plntNA <- ifelse(apply(is.na(LateDev.dat[,plntcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(plntNA, MARGIN=2, sum)

# conglomerate
plntcong <- xtabs(LateDev.dat$plntcomb ~ LateDev.dat$Conglomerate)
plntcong

# sandstone
plntsand <- xtabs(LateDev.dat$plntcomb ~ LateDev.dat$Sandstone)
plntsand

# siltstone
plntsilt <- xtabs(LateDev.dat$plntcomb ~ LateDev.dat$Siltstone)
plntsilt

# shale
plntshal <- xtabs(LateDev.dat$plntcomb ~ LateDev.dat$Shale)
plntshal

# mudstone
plntmud <- xtabs(LateDev.dat$plntcomb ~ LateDev.dat$Mudstone)
plntmud

# limestone
plntlime <- xtabs(LateDev.dat$plntcomb ~ LateDev.dat$Limestone)
plntlime

LateDev.plnt.summ <- as.data.frame(rbind(plntcong,plntsand,plntsilt,plntshal,plntmud,plntlime))
rownames(LateDev.plnt.summ) <- stone.labs
colnames(LateDev.plnt.summ) <- c("abs","pres")
LateDev.plnt.summ
LateDev.plnt.summ$ppres <- LateDev.plnt.summ$pres / apply(LateDev.plnt.summ, MARGIN=1, sum)
LateDev.plnt.summ

ggplot(LateDev.plnt.summ, aes(x=rownames(LateDev.plnt.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

## combined data
LateDev.vert.summ.pres <- data.frame("geol"=rownames(LateDev.vert.summ), "pres"="pres", "freq"=LateDev.vert.summ$pres)
LateDev.vert.summ.abs <- data.frame("geol"=rownames(LateDev.vert.summ), "pres"="abs", "freq"=LateDev.vert.summ$abs)
LateDev.vert.summ2 <- (rbind(LateDev.vert.summ.pres, LateDev.vert.summ.abs))
LateDev.vert.summ2$tax <- "vert"

LateDev.inv.summ.pres <- data.frame("geol"=rownames(LateDev.inv.summ), "pres"="pres", "freq"=LateDev.inv.summ$pres)
LateDev.inv.summ.abs <- data.frame("geol"=rownames(LateDev.inv.summ), "pres"="abs", "freq"=LateDev.inv.summ$abs)
LateDev.inv.summ2 <- (rbind(LateDev.inv.summ.pres, LateDev.inv.summ.abs))
LateDev.inv.summ2$tax <- "inv"

LateDev.plnt.summ.pres <- data.frame("geol"=rownames(LateDev.plnt.summ), "pres"="pres", "freq"=LateDev.plnt.summ$pres)
LateDev.plnt.summ.abs <- data.frame("geol"=rownames(LateDev.plnt.summ), "pres"="abs", "freq"=LateDev.plnt.summ$abs)
LateDev.plnt.summ2 <- (rbind(LateDev.plnt.summ.pres, LateDev.plnt.summ.abs))
LateDev.plnt.summ2$tax <- "plnt"

LateDev.comb.summ <- rbind(LateDev.vert.summ2, LateDev.inv.summ2, LateDev.plnt.summ2)
LateDev.comb.exp <- expand.dft(LateDev.comb.summ, freq="freq")
LateDev.comb.tab <- table(LateDev.comb.exp$tax, LateDev.comb.exp$geol, LateDev.comb.exp$pres)
plot(LateDev.comb.tab)

# mosaic
LateDev.comb.struct <- structable(data=LateDev.comb.tab)
mosaic(LateDev.comb.struct)
mosaic(LateDev.comb.struct, type="expected")
assoc(LateDev.comb.struct, compress=F)

# logistic weighted models
# prepare data
LateDev.comb.exp$pres1 <- as.integer(ifelse(LateDev.comb.exp$pres == "pres", 1, 0))
LateDev.comb.exp$geol <- as.factor(LateDev.comb.exp$geol)
LateDev.comb.exp$tax <- as.factor(LateDev.comb.exp$tax)
str(LateDev.comb.exp)
table(LateDev.comb.exp$pres1)

## model vector
mod.vec <- c(m2,m3,m4,m5)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=LateDev.comb.exp, family=binomial, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

LateDev.sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(LateDev.sumtable) <- mod.vec
LateDev.summary.table <- LateDev.sumtable[order(LateDev.sumtable[,5],decreasing=T),]
LateDev.summary.table

ERgeoltax <- round(LateDev.sumtable[2,5]/LateDev.sumtable[3,5], 3)
ERgeoltax

## saturated residual diagnostic
i <- 2
LateDev.fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), data=LateDev.comb.exp, na.action=na.omit)

check_model(LateDev.fit)
plot_model(LateDev.fit, show.values=T, vline.color = "purple")


###############################
## split by palaeo-environment
###############################
mar.ind <- grep("marin", dat$Environnement_interpretation, ignore.case=T)
mar.ind <- mar.ind[-18]
dat[mar.ind,]$Environnement_interpretation
dat2 <- dat[-mar.ind,]

alluv.ind <- grep("alluv", dat2$Environnement_interpretation, ignore.case=T)
dat2[alluv.ind,]$Environnement_interpretation

delt.ind <- grep("delt", dat2$Environnement_interpretation, ignore.case=T)
dat2[delt.ind,]$Environnement_interpretation

alluvdelt.ind <- sort(c(alluv.ind, delt.ind))
dat2[alluvdelt.ind,]$Environnement_interpretation
dat3 <- dat2[-alluvdelt.ind,]

fresh.ind <- grep("fresh", dat3$Environnement_interpretation, ignore.case=T)
dat3[fresh.ind,]$Environnement_interpretation

estuar.ind <- grep("estuar", dat3$Environnement_interpretation, ignore.case=T)
freshestuar.ind <- sort(c(fresh.ind, estuar.ind))
dat3[freshestuar.ind,]$Environnement_interpretation


##########
## marine
##########
mar.dat <- dat[mar.ind,]
dim(mar.dat)

###############
## vertebrates 
###############
vertcolnames <- colnames(mar.dat[,c(34,35,37,39,41,43,45,47,49,51)])
vertcolnames
vertcols <- which(colnames(mar.dat) %in% vertcolnames)
mar.dat$vertcomb <- ifelse(apply(mar.dat[,vertcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

vertNA <- ifelse(apply(is.na(mar.dat[,vertcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(vertNA, MARGIN=2, sum)

# conglomerate
vertcong <- xtabs(mar.dat$vertcomb ~ mar.dat$Conglomerate)
vertcong

# sandstone
vertsand <- xtabs(mar.dat$vertcomb ~ mar.dat$Sandstone)
vertsand

# siltstone
vertsilt <- xtabs(mar.dat$vertcomb ~ mar.dat$Siltstone)
vertsilt

# shale
vertshal <- xtabs(mar.dat$vertcomb ~ mar.dat$Shale)
vertshal

# mudstone
vertmud <- xtabs(mar.dat$vertcomb ~ mar.dat$Mudstone)
vertmud

# limestone
vertlime <- xtabs(mar.dat$vertcomb ~ mar.dat$Limestone)
vertlime

mar.vert.summ <- as.data.frame(rbind(vertcong,vertsand,vertsilt,vertshal,vertmud,vertlime))
rownames(mar.vert.summ) <- stone.labs
colnames(mar.vert.summ) <- c("abs","pres")
mar.vert.summ
mar.vert.summ$ppres <- mar.vert.summ$pres / apply(mar.vert.summ, MARGIN=1, sum)
mar.vert.summ

ggplot(mar.vert.summ, aes(x=rownames(mar.vert.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

#################
## invertebrates 
#################
invcolnames <- colnames(mar.dat[,c(59:78)])
invcolnames
invcols <- which(colnames(mar.dat) %in% invcolnames)
mar.dat$invcomb <- ifelse(apply(mar.dat[,invcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

invNA <- ifelse(apply(is.na(mar.dat[,invcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(invNA, MARGIN=2, sum)

# conglomerate
invcong <- xtabs(mar.dat$invcomb ~ mar.dat$Conglomerate)
invcong

# sandstone
invsand <- xtabs(mar.dat$invcomb ~ mar.dat$Sandstone)
invsand

# siltstone
invsilt <- xtabs(mar.dat$invcomb ~ mar.dat$Siltstone)
invsilt

# shale
invshal <- xtabs(mar.dat$invcomb ~ mar.dat$Shale)
invshal

# mudstone
invmud <- xtabs(mar.dat$invcomb ~ mar.dat$Mudstone)
invmud

# limestone
invlime <- xtabs(mar.dat$invcomb ~ mar.dat$Limestone)
invlime

mar.inv.summ <- as.data.frame(rbind(invcong,invsand,invsilt,invshal,invmud,invlime))
rownames(mar.inv.summ) <- stone.labs
colnames(mar.inv.summ) <- c("abs","pres")
mar.inv.summ
mar.inv.summ$ppres <- mar.inv.summ$pres / apply(mar.inv.summ, MARGIN=1, sum)
mar.inv.summ

ggplot(mar.inv.summ, aes(x=rownames(mar.inv.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

##########
## plants 
##########
plntcolnames <- colnames(mar.dat[,c(80:87)])
plntcolnames
plntcols <- which(colnames(mar.dat) %in% plntcolnames)
mar.dat$plntcomb <- ifelse(apply(mar.dat[,plntcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

plntNA <- ifelse(apply(is.na(mar.dat[,plntcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(plntNA, MARGIN=2, sum)

# conglomerate
plntcong <- xtabs(mar.dat$plntcomb ~ mar.dat$Conglomerate)
plntcong

# sandstone
plntsand <- xtabs(mar.dat$plntcomb ~ mar.dat$Sandstone)
plntsand

# siltstone
plntsilt <- xtabs(mar.dat$plntcomb ~ mar.dat$Siltstone)
plntsilt

# shale
plntshal <- xtabs(mar.dat$plntcomb ~ mar.dat$Shale)
plntshal

# mudstone
plntmud <- xtabs(mar.dat$plntcomb ~ mar.dat$Mudstone)
plntmud

# limestone
plntlime <- xtabs(mar.dat$plntcomb ~ mar.dat$Limestone)
plntlime

mar.plnt.summ <- as.data.frame(rbind(plntcong,plntsand,plntsilt,plntshal,plntmud,plntlime))
rownames(mar.plnt.summ) <- stone.labs
colnames(mar.plnt.summ) <- c("abs","pres")
mar.plnt.summ
mar.plnt.summ$ppres <- mar.plnt.summ$pres / apply(mar.plnt.summ, MARGIN=1, sum)
mar.plnt.summ

ggplot(mar.plnt.summ, aes(x=rownames(mar.plnt.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

## combined data
mar.vert.summ.pres <- data.frame("geol"=rownames(mar.vert.summ), "pres"="pres", "freq"=mar.vert.summ$pres)
mar.vert.summ.abs <- data.frame("geol"=rownames(mar.vert.summ), "pres"="abs", "freq"=mar.vert.summ$abs)
mar.vert.summ2 <- (rbind(mar.vert.summ.pres, mar.vert.summ.abs))
mar.vert.summ2$tax <- "vert"

mar.inv.summ.pres <- data.frame("geol"=rownames(mar.inv.summ), "pres"="pres", "freq"=mar.inv.summ$pres)
mar.inv.summ.abs <- data.frame("geol"=rownames(mar.inv.summ), "pres"="abs", "freq"=mar.inv.summ$abs)
mar.inv.summ2 <- (rbind(mar.inv.summ.pres, mar.inv.summ.abs))
mar.inv.summ2$tax <- "inv"

mar.plnt.summ.pres <- data.frame("geol"=rownames(mar.plnt.summ), "pres"="pres", "freq"=mar.plnt.summ$pres)
mar.plnt.summ.abs <- data.frame("geol"=rownames(mar.plnt.summ), "pres"="abs", "freq"=mar.plnt.summ$abs)
mar.plnt.summ2 <- (rbind(mar.plnt.summ.pres, mar.plnt.summ.abs))
mar.plnt.summ2$tax <- "plnt"

mar.comb.summ <- rbind(mar.vert.summ2, mar.inv.summ2, mar.plnt.summ2)
mar.comb.exp <- expand.dft(mar.comb.summ, freq="freq")
mar.comb.tab <- table(mar.comb.exp$tax, mar.comb.exp$geol, mar.comb.exp$pres)
plot(mar.comb.tab)

# mosaic
mar.comb.struct <- structable(data=mar.comb.tab)
mosaic(mar.comb.struct)
mosaic(mar.comb.struct, type="expected")
assoc(mar.comb.struct, compress=F)

# logistic weighted models
# prepare data
mar.comb.exp$pres1 <- as.integer(ifelse(mar.comb.exp$pres == "pres", 1, 0))
mar.comb.exp$geol <- as.factor(mar.comb.exp$geol)
mar.comb.exp$tax <- as.factor(mar.comb.exp$tax)
str(mar.comb.exp)
table(mar.comb.exp$pres1)

## model vector
mod.vec <- c(m2,m3,m4,m5)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=mar.comb.exp, family=binomial, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

mar.sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(mar.sumtable) <- mod.vec
mar.summary.table <- mar.sumtable[order(mar.sumtable[,5],decreasing=T),]
mar.summary.table

ERgeoltax <- round(mar.sumtable[2,5]/mar.sumtable[3,5], 3)
ERgeoltax

## saturated residual diagnostic
i <- 2
mar.fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), data=mar.comb.exp, na.action=na.omit)

check_model(mar.fit)
plot_model(mar.fit, show.values=T, vline.color = "purple")



####################
## alluvial/deltaic
####################
alluvdelt.dat <- dat[alluvdelt.ind,]
dim(alluvdelt.dat)

###############
## vertebrates 
###############
vertcolnames <- colnames(alluvdelt.dat[,c(34,35,37,39,41,43,45,47,49,51)])
vertcolnames
vertcols <- which(colnames(alluvdelt.dat) %in% vertcolnames)
alluvdelt.dat$vertcomb <- ifelse(apply(alluvdelt.dat[,vertcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

vertNA <- ifelse(apply(is.na(alluvdelt.dat[,vertcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(vertNA, MARGIN=2, sum)

# conglomerate
vertcong <- xtabs(alluvdelt.dat$vertcomb ~ alluvdelt.dat$Conglomerate)
vertcong

# sandstone
vertsand <- xtabs(alluvdelt.dat$vertcomb ~ alluvdelt.dat$Sandstone)
vertsand

# siltstone
vertsilt <- xtabs(alluvdelt.dat$vertcomb ~ alluvdelt.dat$Siltstone)
vertsilt

# shale
vertshal <- xtabs(alluvdelt.dat$vertcomb ~ alluvdelt.dat$Shale)
vertshal

# mudstone
vertmud <- xtabs(alluvdelt.dat$vertcomb ~ alluvdelt.dat$Mudstone)
vertmud

# limestone
vertlime <- xtabs(alluvdelt.dat$vertcomb ~ alluvdelt.dat$Limestone)
vertlime

alluvdelt.vert.summ <- as.data.frame(rbind(vertcong,vertsand,vertsilt,vertshal,vertmud,vertlime))
rownames(alluvdelt.vert.summ) <- stone.labs
colnames(alluvdelt.vert.summ) <- c("abs","pres")
alluvdelt.vert.summ
alluvdelt.vert.summ$ppres <- alluvdelt.vert.summ$pres / apply(alluvdelt.vert.summ, MARGIN=1, sum)
alluvdelt.vert.summ

ggplot(alluvdelt.vert.summ, aes(x=rownames(alluvdelt.vert.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

#################
## invertebrates 
#################
invcolnames <- colnames(alluvdelt.dat[,c(59:78)])
invcolnames
invcols <- which(colnames(alluvdelt.dat) %in% invcolnames)
alluvdelt.dat$invcomb <- ifelse(apply(alluvdelt.dat[,invcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

invNA <- ifelse(apply(is.na(alluvdelt.dat[,invcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(invNA, MARGIN=2, sum)

# conglomerate
invcong <- xtabs(alluvdelt.dat$invcomb ~ alluvdelt.dat$Conglomerate)
invcong

# sandstone
invsand <- xtabs(alluvdelt.dat$invcomb ~ alluvdelt.dat$Sandstone)
invsand

# siltstone
invsilt <- xtabs(alluvdelt.dat$invcomb ~ alluvdelt.dat$Siltstone)
invsilt

# shale
invshal <- xtabs(alluvdelt.dat$invcomb ~ alluvdelt.dat$Shale)
invshal

# mudstone
invmud <- xtabs(alluvdelt.dat$invcomb ~ alluvdelt.dat$Mudstone)
invmud

# limestone
invlime <- xtabs(alluvdelt.dat$invcomb ~ alluvdelt.dat$Limestone)
invlime

alluvdelt.inv.summ <- as.data.frame(rbind(invcong,invsand,invsilt,invshal,invmud,invlime))
rownames(alluvdelt.inv.summ) <- stone.labs
colnames(alluvdelt.inv.summ) <- c("abs","pres")
alluvdelt.inv.summ
alluvdelt.inv.summ$ppres <- alluvdelt.inv.summ$pres / apply(alluvdelt.inv.summ, MARGIN=1, sum)
alluvdelt.inv.summ

ggplot(alluvdelt.inv.summ, aes(x=rownames(alluvdelt.inv.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

##########
## plants 
##########
plntcolnames <- colnames(alluvdelt.dat[,c(80:87)])
plntcolnames
plntcols <- which(colnames(alluvdelt.dat) %in% plntcolnames)
alluvdelt.dat$plntcomb <- ifelse(apply(alluvdelt.dat[,plntcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

plntNA <- ifelse(apply(is.na(alluvdelt.dat[,plntcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(plntNA, MARGIN=2, sum)

# conglomerate
plntcong <- xtabs(alluvdelt.dat$plntcomb ~ alluvdelt.dat$Conglomerate)
plntcong

# sandstone
plntsand <- xtabs(alluvdelt.dat$plntcomb ~ alluvdelt.dat$Sandstone)
plntsand

# siltstone
plntsilt <- xtabs(alluvdelt.dat$plntcomb ~ alluvdelt.dat$Siltstone)
plntsilt

# shale
plntshal <- xtabs(alluvdelt.dat$plntcomb ~ alluvdelt.dat$Shale)
plntshal

# mudstone
plntmud <- xtabs(alluvdelt.dat$plntcomb ~ alluvdelt.dat$Mudstone)
plntmud

# limestone
plntlime <- xtabs(alluvdelt.dat$plntcomb ~ alluvdelt.dat$Limestone)
plntlime

alluvdelt.plnt.summ <- as.data.frame(rbind(plntcong,plntsand,plntsilt,plntshal,plntmud,plntlime))
rownames(alluvdelt.plnt.summ) <- stone.labs
colnames(alluvdelt.plnt.summ) <- c("abs","pres")
alluvdelt.plnt.summ
alluvdelt.plnt.summ$ppres <- alluvdelt.plnt.summ$pres / apply(alluvdelt.plnt.summ, MARGIN=1, sum)
alluvdelt.plnt.summ

ggplot(alluvdelt.plnt.summ, aes(x=rownames(alluvdelt.plnt.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

## combined data
alluvdelt.vert.summ.pres <- data.frame("geol"=rownames(alluvdelt.vert.summ), "pres"="pres", "freq"=alluvdelt.vert.summ$pres)
alluvdelt.vert.summ.abs <- data.frame("geol"=rownames(alluvdelt.vert.summ), "pres"="abs", "freq"=alluvdelt.vert.summ$abs)
alluvdelt.vert.summ2 <- (rbind(alluvdelt.vert.summ.pres, alluvdelt.vert.summ.abs))
alluvdelt.vert.summ2$tax <- "vert"

alluvdelt.inv.summ.pres <- data.frame("geol"=rownames(alluvdelt.inv.summ), "pres"="pres", "freq"=alluvdelt.inv.summ$pres)
alluvdelt.inv.summ.abs <- data.frame("geol"=rownames(alluvdelt.inv.summ), "pres"="abs", "freq"=alluvdelt.inv.summ$abs)
alluvdelt.inv.summ2 <- (rbind(alluvdelt.inv.summ.pres, alluvdelt.inv.summ.abs))
alluvdelt.inv.summ2$tax <- "inv"

alluvdelt.plnt.summ.pres <- data.frame("geol"=rownames(alluvdelt.plnt.summ), "pres"="pres", "freq"=alluvdelt.plnt.summ$pres)
alluvdelt.plnt.summ.abs <- data.frame("geol"=rownames(alluvdelt.plnt.summ), "pres"="abs", "freq"=alluvdelt.plnt.summ$abs)
alluvdelt.plnt.summ2 <- (rbind(alluvdelt.plnt.summ.pres, alluvdelt.plnt.summ.abs))
alluvdelt.plnt.summ2$tax <- "plnt"

alluvdelt.comb.summ <- rbind(alluvdelt.vert.summ2, alluvdelt.inv.summ2, alluvdelt.plnt.summ2)
alluvdelt.comb.exp <- expand.dft(alluvdelt.comb.summ, freq="freq")
alluvdelt.comb.tab <- table(alluvdelt.comb.exp$tax, alluvdelt.comb.exp$geol, alluvdelt.comb.exp$pres)
plot(alluvdelt.comb.tab)

# mosaic
alluvdelt.comb.struct <- structable(data=alluvdelt.comb.tab)
mosaic(alluvdelt.comb.struct)
mosaic(alluvdelt.comb.struct, type="expected")
assoc(alluvdelt.comb.struct, compress=F)

# logistic weighted models
# prepare data
alluvdelt.comb.exp$pres1 <- as.integer(ifelse(alluvdelt.comb.exp$pres == "pres", 1, 0))
alluvdelt.comb.exp$geol <- as.factor(alluvdelt.comb.exp$geol)
alluvdelt.comb.exp$tax <- as.factor(alluvdelt.comb.exp$tax)
str(alluvdelt.comb.exp)
table(alluvdelt.comb.exp$pres1)

## model vector
mod.vec <- c(m2,m3,m4,m5)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=alluvdelt.comb.exp, family=binomial, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

alluvdelt.sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(alluvdelt.sumtable) <- mod.vec
alluvdelt.summary.table <- alluvdelt.sumtable[order(alluvdelt.sumtable[,5],decreasing=T),]
alluvdelt.summary.table

ERgeoltax <- round(alluvdelt.sumtable[2,5]/alluvdelt.sumtable[3,5], 3)
ERgeoltax

## saturated residual diagnostic
i <- 2
alluvdelt.fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), data=alluvdelt.comb.exp, na.action=na.omit)

check_model(alluvdelt.fit)
plot_model(alluvdelt.fit, show.values=T, vline.color = "purple")


########################
## freshwater/estuarine
########################
freshestuar.dat <- dat[freshestuar.ind,]
dim(freshestuar.dat)

###############
## vertebrates 
###############
vertcolnames <- colnames(freshestuar.dat[,c(34,35,37,39,41,43,45,47,49,51)])
vertcolnames
vertcols <- which(colnames(freshestuar.dat) %in% vertcolnames)
freshestuar.dat$vertcomb <- ifelse(apply(freshestuar.dat[,vertcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

vertNA <- ifelse(apply(is.na(freshestuar.dat[,vertcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(vertNA, MARGIN=2, sum)

# conglomerate
vertcong <- xtabs(freshestuar.dat$vertcomb ~ freshestuar.dat$Conglomerate)
vertcong

# sandstone
vertsand <- xtabs(freshestuar.dat$vertcomb ~ freshestuar.dat$Sandstone)
vertsand

# siltstone
vertsilt <- xtabs(freshestuar.dat$vertcomb ~ freshestuar.dat$Siltstone)
vertsilt

# shale
vertshal <- xtabs(freshestuar.dat$vertcomb ~ freshestuar.dat$Shale)
vertshal

# mudstone
vertmud <- xtabs(freshestuar.dat$vertcomb ~ freshestuar.dat$Mudstone)
vertmud

# limestone
vertlime <- xtabs(freshestuar.dat$vertcomb ~ freshestuar.dat$Limestone)
vertlime

freshestuar.vert.summ <- as.data.frame(rbind(vertcong,vertsand,vertsilt,vertshal,vertmud,vertlime))
rownames(freshestuar.vert.summ) <- stone.labs
colnames(freshestuar.vert.summ) <- c("abs","pres")
freshestuar.vert.summ
freshestuar.vert.summ$ppres <- freshestuar.vert.summ$pres / apply(freshestuar.vert.summ, MARGIN=1, sum)
freshestuar.vert.summ

ggplot(freshestuar.vert.summ, aes(x=rownames(freshestuar.vert.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

#################
## invertebrates 
#################
invcolnames <- colnames(freshestuar.dat[,c(59:78)])
invcolnames
invcols <- which(colnames(freshestuar.dat) %in% invcolnames)
freshestuar.dat$invcomb <- ifelse(apply(freshestuar.dat[,invcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

invNA <- ifelse(apply(is.na(freshestuar.dat[,invcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(invNA, MARGIN=2, sum)

# conglomerate
invcong <- xtabs(freshestuar.dat$invcomb ~ freshestuar.dat$Conglomerate)
invcong

# sandstone
invsand <- xtabs(freshestuar.dat$invcomb ~ freshestuar.dat$Sandstone)
invsand

# siltstone
invsilt <- xtabs(freshestuar.dat$invcomb ~ freshestuar.dat$Siltstone)
invsilt

# shale
invshal <- xtabs(freshestuar.dat$invcomb ~ freshestuar.dat$Shale)
invshal

# mudstone
invmud <- xtabs(freshestuar.dat$invcomb ~ freshestuar.dat$Mudstone)
invmud

# limestone
invlime <- xtabs(freshestuar.dat$invcomb ~ freshestuar.dat$Limestone)
invlime

freshestuar.inv.summ <- as.data.frame(rbind(invcong,invsand,invsilt,invshal,invmud,invlime))
rownames(freshestuar.inv.summ) <- stone.labs
colnames(freshestuar.inv.summ) <- c("abs","pres")
freshestuar.inv.summ
freshestuar.inv.summ$ppres <- freshestuar.inv.summ$pres / apply(freshestuar.inv.summ, MARGIN=1, sum)
freshestuar.inv.summ

ggplot(freshestuar.inv.summ, aes(x=rownames(freshestuar.inv.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

##########
## plants 
##########
plntcolnames <- colnames(freshestuar.dat[,c(80:87)])
plntcolnames
plntcols <- which(colnames(freshestuar.dat) %in% plntcolnames)
freshestuar.dat$plntcomb <- ifelse(apply(freshestuar.dat[,plntcols], MARGIN=1, sum, na.rm=T) > 0, 1, 0)

plntNA <- ifelse(apply(is.na(freshestuar.dat[,plntcols]), MARGIN=1, is.na) == TRUE, 1, 0)
apply(plntNA, MARGIN=2, sum)

# conglomerate
plntcong <- xtabs(freshestuar.dat$plntcomb ~ freshestuar.dat$Conglomerate)
plntcong

# sandstone
plntsand <- xtabs(freshestuar.dat$plntcomb ~ freshestuar.dat$Sandstone)
plntsand

# siltstone
plntsilt <- xtabs(freshestuar.dat$plntcomb ~ freshestuar.dat$Siltstone)
plntsilt

# shale
plntshal <- xtabs(freshestuar.dat$plntcomb ~ freshestuar.dat$Shale)
plntshal

# mudstone
plntmud <- xtabs(freshestuar.dat$plntcomb ~ freshestuar.dat$Mudstone)
plntmud

# limestone
plntlime <- xtabs(freshestuar.dat$plntcomb ~ freshestuar.dat$Limestone)
plntlime

freshestuar.plnt.summ <- as.data.frame(rbind(plntcong,plntsand,plntsilt,plntshal,plntmud,plntlime))
rownames(freshestuar.plnt.summ) <- stone.labs
colnames(freshestuar.plnt.summ) <- c("abs","pres")
freshestuar.plnt.summ
freshestuar.plnt.summ$ppres <- freshestuar.plnt.summ$pres / apply(freshestuar.plnt.summ, MARGIN=1, sum)
freshestuar.plnt.summ

ggplot(freshestuar.plnt.summ, aes(x=rownames(freshestuar.plnt.summ), y=ppres)) +
  geom_bar(stat = "identity") +
  ylim(0,1) +
  ylab("% present") +
  xlab("geology")

## combined data
freshestuar.vert.summ.pres <- data.frame("geol"=rownames(freshestuar.vert.summ), "pres"="pres", "freq"=freshestuar.vert.summ$pres)
freshestuar.vert.summ.abs <- data.frame("geol"=rownames(freshestuar.vert.summ), "pres"="abs", "freq"=freshestuar.vert.summ$abs)
freshestuar.vert.summ2 <- (rbind(freshestuar.vert.summ.pres, freshestuar.vert.summ.abs))
freshestuar.vert.summ2$tax <- "vert"

freshestuar.inv.summ.pres <- data.frame("geol"=rownames(freshestuar.inv.summ), "pres"="pres", "freq"=freshestuar.inv.summ$pres)
freshestuar.inv.summ.abs <- data.frame("geol"=rownames(freshestuar.inv.summ), "pres"="abs", "freq"=freshestuar.inv.summ$abs)
freshestuar.inv.summ2 <- (rbind(freshestuar.inv.summ.pres, freshestuar.inv.summ.abs))
freshestuar.inv.summ2$tax <- "inv"

freshestuar.plnt.summ.pres <- data.frame("geol"=rownames(freshestuar.plnt.summ), "pres"="pres", "freq"=freshestuar.plnt.summ$pres)
freshestuar.plnt.summ.abs <- data.frame("geol"=rownames(freshestuar.plnt.summ), "pres"="abs", "freq"=freshestuar.plnt.summ$abs)
freshestuar.plnt.summ2 <- (rbind(freshestuar.plnt.summ.pres, freshestuar.plnt.summ.abs))
freshestuar.plnt.summ2$tax <- "plnt"

freshestuar.comb.summ <- rbind(freshestuar.vert.summ2, freshestuar.inv.summ2, freshestuar.plnt.summ2)
freshestuar.comb.exp <- expand.dft(freshestuar.comb.summ, freq="freq")
freshestuar.comb.tab <- table(freshestuar.comb.exp$tax, freshestuar.comb.exp$geol, freshestuar.comb.exp$pres)
plot(freshestuar.comb.tab)

# mosaic
freshestuar.comb.struct <- structable(data=freshestuar.comb.tab)
mosaic(freshestuar.comb.struct)
mosaic(freshestuar.comb.struct, type="expected")
assoc(freshestuar.comb.struct, compress=F)

# logistic weighted models
# prepare data
freshestuar.comb.exp$pres1 <- as.integer(ifelse(freshestuar.comb.exp$pres == "pres", 1, 0))
freshestuar.comb.exp$geol <- as.factor(freshestuar.comb.exp$geol)
freshestuar.comb.exp$tax <- as.factor(freshestuar.comb.exp$tax)
str(freshestuar.comb.exp)
table(freshestuar.comb.exp$pres1)

## model vector
mod.vec <- c(m2,m3,m4,m5)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=freshestuar.comb.exp, family=binomial, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

freshestuar.sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(freshestuar.sumtable) <- mod.vec
freshestuar.summary.table <- freshestuar.sumtable[order(freshestuar.sumtable[,5],decreasing=T),]
freshestuar.summary.table

ERgeoltax <- round(freshestuar.sumtable[2,5]/freshestuar.sumtable[3,5], 3)
ERgeoltax

## saturated residual diagnostic
i <- 2
freshestuar.fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), data=freshestuar.comb.exp, na.action=na.omit)

check_model(freshestuar.fit)
plot_model(freshestuar.fit, show.values=T, vline.color = "purple")


## permutation test
# taxa
vertcolnames
invcolnames
plntcolnames
taxa.sub <- c(vertcolnames, invcolnames, plntcolnames[-length(plntcolnames)])
dat[,taxa.sub]

# geology
str(dat)
geolID <- c("cng","snd","slt","shl","mud","lim")
geol.sub <- c(which(colnames(dat)=="Conglomerate"), which(colnames(dat)=="Sandstone"), 
              which(colnames(dat)=="Siltstone"), which(colnames(dat)=="Shale"),
              which(colnames(dat)=="Mudstone"), which(colnames(dat)=="Limestone"))
dat[,geol.sub]
head(geol.dat)
geolpc.sub <- c(which(colnames(dat)=="X._conglomerate"), which(colnames(dat)=="X._sandstone"), 
                which(colnames(dat)=="X._siltstone"), which(colnames(dat)=="X._shale"),
                which(colnames(dat)=="X._mudstone"), which(colnames(dat)=="X._limestone"))

table(dat$lDominant_lith)
dat$lDominant_lith <- revalue(dat$lDominant_lith, c("CON" = "cng",
                              "LIM" = "lim",
                              "MUD" = "mud",
                              "SAN" = "snd",
                              "SHA" = "shl",
                              "SHA LIM" = "shllim",
                              "SHA SAN" = "shlsnd",
                              "SHA SIL" = "shlslt",
                              "SIL" = "slt"))
table(dat$lDominant_lith)
domlith.sub <- which(colnames(dat)=="lDominant_lith")
geol.dat <- dat[, c(geol.sub, domlith.sub)]
head(geol.dat)

geolpc.dat <- as.data.frame(dat[, geolpc.sub])
str(geolpc.dat)
geolpc.dat$X._conglomerate <- as.numeric(geolpc.dat$X._conglomerate)
geolpc.dat$X._sandstone <- as.numeric(geolpc.dat$X._sandstone)
geolpc.dat$X._siltstone <- as.numeric(geolpc.dat$X._siltstone)
geolpc.dat$X._shale <- as.numeric(geolpc.dat$X._shale)
geolpc.dat$X._mudstone <- as.numeric(geolpc.dat$X._mudstone)
geolpc.dat$X._limestone <- as.numeric(geolpc.dat$X._limestone)
str(geolpc.dat)

geolpc.dom <- apply(geolpc.dat, MARGIN=1, max, na.rm=T)
geolpc.dom[which(geolpc.dom == -Inf)] <- NA
geolpc.dom[which(geolpc.dom == 0)] <- NA
geolpc.dom

geol.comb <- geoldom.cat <- rep(NA,length(geolpc.dom))
for (i in 1:length(geolpc.dom)) {
   geoldom.cat[i] <- geolpc$lDominant_lith[i]
   geol.comb[i] <- paste(geolID[which(geol.dat[i,]==1)], collapse="")
   geol.comb[i] <- ifelse(is.na(geoldom.cat[i])==F, geoldom.cat[i], geol.comb[i])
}
geol.comb[which(geol.comb=="")] <- NA
geol.comb
table(geol.comb)
length(table(geol.comb))

# period
## split by major period(Silurian, lower, mid, upper Devonian)
perID <- c("Silur","EarlDev","MidDev","LateDev")
percat <- rep(NA,dim(dat)[1])

sil.ind <- grep("silur", dat$Age, ignore.case=T)
dat[sil.ind,]$Age
percat[sil.ind] <- perID[1]

earldev.ind <- grep("lower", dat$Age, ignore.case=T)
dat[earldev.ind,]$Age
rem.prev1.sub <- which(earldev.ind %in% sil.ind == T)
which(earldev.ind[-rem.prev1.sub] %in% sil.ind == T)
percat[earldev.ind[-rem.prev1.sub]] <- perID[2]

middev.ind <- grep("mid", dat$Age, ignore.case=T)
dat[middev.ind,]$Age
sort(c(sil.ind, earldev.ind))
which(middev.ind %in% sort(c(sil.ind, earldev.ind)) == T)
percat[middev.ind] <- perID[3]

latedev.ind <- grep("Upp", dat$Age, ignore.case=T)
dat[latedev.ind,]$Age
sort(c(sil.ind, earldev.ind, middev.ind))
which(latedev.ind %in% sort(c(sil.ind, earldev.ind, middev.ind)) == T)
rem.prev2.sub <- which(latedev.ind %in% sort(c(sil.ind, earldev.ind, middev.ind)) == T)
dat[latedev.ind[-rem.prev2.sub],]$Age
percat[latedev.ind[-rem.prev2.sub]] <- perID[4]
percat

# combine taxa, geology, period
dat.reclass <- data.frame("ID"=dat$ID, "percat"=percat, "geolcat"=geol.comb, dat[,taxa.sub])
head(dat.reclass)
dat.reclass.red <- dat.reclass[-which(is.na(dat.reclass$percat)==T | is.na(dat.reclass$geolcat)==T),]
head(dat.reclass.red)
dim(dat.reclass.red)
table(dat.reclass.red$geolcat)

## run permutation linear models for community presence/absence data
dat.reclass.red.zero <- dat.reclass.red
dat.reclass.red.zero[is.na(dat.reclass.red.zero) == T] <- 0 # treat taxa NAs as absence (0)
tax.reclass <- dat.reclass.red.zero[,4:dim(dat.reclass.red.zero)[2]]
pergeo.reclass <- dat.reclass.red.zero[,2:3]


## permutation analysis of variance (Gower index)
## stochastic loop to select dominant geology class randomly
iter <- 1000
itdiv <- iter/10
permut <- 10000
st.time <- Sys.time()

fact.levels <- c("percatSilur", "percatLowDev", "percatMidDev", "percatUppDev",
                 "geolcatcng", "geolcatsnd", "geolcatslt", "geolcatshl", "geolcatmud", "geolcatlim")
factn <- length(fact.levels)

factlev.pmat <- factlev.R2mat <- matrix(data=NA, nrow=iter, ncol=factn)
colnames(factlev.R2mat) <- fact.levels
colnames(factlev.pmat) <- fact.levels
p.mat <- R2.mat <- matrix(data=NA, nrow=iter, ncol=3)

for (s in 1:iter) {
  
  # choose 'dominant' geo class randomly
  pergeo.reclass.it <- pergeo.reclass
  for (d in 1:dim(pergeo.reclass.it)[1]) {
    geo.split.vec <- as.vector(str_split(gsub("(.{3})", "\\1 ", pergeo.reclass.it$geolcat[d]), " ")[[1]])
    geosplit.vec <- geo.split.vec[-which(geo.split.vec == "")]
    pergeo.reclass.it$geolcat[d] <- sample(geosplit.vec,1)
  } # end d loop

  # resample dataset with replacement
  resamp.sub <- sort(sample(1:dim(pergeo.reclass.it)[1], round(dim(pergeo.reclass.it)[1]/1.5,0), replace=F))
  dat.resamp1 <- data.frame(pergeo.reclass.it, tax.reclass)
  dat.resamp <- dat.resamp1[resamp.sub,] 
  pergeo.resamp <- dat.resamp[,1:2]
  tax.resamp <- dat.resamp[,3:dim(dat.resamp)[2]]
  tax.resamp <- mutate_all(tax.resamp, function(x) as.integer(as.character(x)))
  
  # permutation analysis of variance (Cao index)
  div.iter <- adonis2(tax.resamp ~ percat*geolcat, data = pergeo.resamp,
                      permutations = permut, method="gower", binary=T, sqrt.dist=T, by="terms", parallel=cl)
  div.iter2 <- adonis2(tax.resamp ~ percat*geolcat, data = pergeo.resamp,
                      permutations = permut, method="gower", binary=T, sqrt.dist=T, by="onedf", parallel=cl)

  # store single-level p and R2 for factors
  for (n in 1:factn) {
    table.sub <- which(rownames(div.iter2[1]) == fact.levels[n])
    if (length(table.sub) > 0) {
      factlev.pmat[s, n] <- div.iter2$`Pr(>F)`[table.sub]
      factlev.R2mat[s, n] <- div.iter2$R2[table.sub]
    }
  }

  # save pr(F) results
  p.mat[s, ] <- c(div.iter$'Pr(>F)'[1], div.iter$'Pr(>F)'[2], div.iter$'Pr(>F)'[3])
  R2.mat[s, ] <- c(div.iter$R2[1], div.iter$R2[2], div.iter$R2[3])
  
  if (s %% itdiv==0) print(s) 
  
} # end s

end.time <- Sys.time()
proc.time <- end.time - st.time
proc.time

# set CI lims
# 0.67, 0.8, 0.95, 0.99
# i.e., 0.3333, 0.2, 0.05, 0.01
oneINx <- 20 # 3 # 5 # 20 # 100
CI.prob <- 1/oneINx
print(paste(round((1 - CI.prob)*100, 0), "% confidence interval", sep=""))
up.lim <- (1-CI.prob/2)
lo.lim <- 1 - up.lim

## period
# probabilities
exp(median(log(p.mat[, 1])))
print(c(exp(quantile(log(p.mat[, 1]), probs=lo.lim, na.rm=T)), exp(quantile(log(p.mat[, 1]), probs=up.lim, na.rm=T))))
hist(log(p.mat[, 1]), main="", xlab="log period prob")
abline(v=log(0.05), lty=1, lwd=3, col="red")
abline(v=(median(log(p.mat[, 1]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(p.mat[, 1]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(p.mat[, 1]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

# R2
exp(median(log(R2.mat[, 1])))
print(c(exp(quantile(log(R2.mat[, 1]), probs=lo.lim, na.rm=T)), exp(quantile(log(R2.mat[, 1]), probs=up.lim, na.rm=T))))
hist(log(R2.mat[, 1]), main="", xlab="log period R2")
abline(v=(median(log(R2.mat[, 1]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(R2.mat[, 1]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(R2.mat[, 1]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

## geology 
# probabilities
exp(median(log(p.mat[, 2])))
print(c(exp(quantile(log(p.mat[, 2]), probs=lo.lim, na.rm=T)), exp(quantile(log(p.mat[, 2]), probs=up.lim, na.rm=T))))
hist(log(p.mat[, 2]), main="", xlab="log geol prob")
abline(v=log(0.05), lty=1, lwd=3, col="red")
abline(v=(median(log(p.mat[, 2]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(p.mat[, 2]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(p.mat[, 2]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

# R2
exp(median(log(R2.mat[, 2])))
print(c(exp(quantile(log(R2.mat[, 2]), probs=lo.lim, na.rm=T)), exp(quantile(log(R2.mat[, 2]), probs=up.lim, na.rm=T))))
hist(log(R2.mat[, 2]), main="", xlab="log period R2")
abline(v=(median(log(R2.mat[, 2]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(R2.mat[, 2]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(R2.mat[, 2]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

## interaction 
# probabilities
exp(median(log(p.mat[, 3])))
print(c(exp(quantile(log(p.mat[, 3]), probs=lo.lim, na.rm=T)), exp(quantile(log(p.mat[, 3]), probs=up.lim, na.rm=T))))
hist(log(p.mat[, 3]), main="", xlab="log period*geology prob")
abline(v=log(0.05), lty=1, lwd=3, col="red")
abline(v=(median(log(p.mat[, 3]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(p.mat[, 3]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(p.mat[, 3]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

# R2
exp(median(log(R2.mat[, 3])))
print(c(exp(quantile(log(R2.mat[, 3]), probs=lo.lim, na.rm=T)), exp(quantile(log(R2.mat[, 3]), probs=up.lim, na.rm=T))))
hist(log(R2.mat[, 3]), main="", xlab="log period R2")
abline(v=(median(log(R2.mat[, 3]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(R2.mat[, 3]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(R2.mat[, 3]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

# single factor-level median probabilities & R2
fact.levels.sname <- c("Silur", "LowDev", "MidDev", "UppDev", "cng", "snd", "slt", "shl", "mud", "lim")
factlev.out <- data.frame(fact.levels.sname, apply(factlev.pmat, MARGIN=2, median, na.rm=T),
           apply(factlev.R2mat, MARGIN=2, median, na.rm=T))
colnames(factlev.out) <- c("factlev", "p", "R2")
rownames(factlev.out) <- NULL
factlev.sort <- factlev.out[order(factlev.out[,3],decreasing=T),]
factlev.sort

# percentages per taxon & epoch
# epoch
dim(dat.reclass)[2]
epoch.xtabs <- matrix(data=NA, ncol=4, nrow=dim(dat.reclass)[2] - 3)
for (x in 4:(length(colnames(dat.reclass[3:dim(dat.reclass)[2]]))+2)) {
  epoch.xtabs[x-3,] <- xtabs(dat.reclass[,x] ~ dat.reclass$percat)/sum(dat.reclass[,x], na.rm=T)
}
colnames(epoch.xtabs) <- c("LowDev", "MidDev", "Silur", "UppDev")
rownames(epoch.xtabs) <- colnames(dat.reclass[4:dim(dat.reclass)[2]])
epoch.xtabs

write.table(epoch.xtabs, "epochxtabs.csv", sep=",", row.names = T)

# geology
geolcols <- c(15,17,19,21,23,25)
geoltax.full <- data.frame(dat[,geolcols], dat[,vertcols], dat[,invcols], dat[,plntcols[-length(plntcols)]])
colnames(geoltax.full) <- c(colnames(dat[,geolcols], ), colnames(dat[,vertcols], ), colnames(dat[,invcols], ),
                            colnames(dat[,plntcols[-length(plntcols)]], ))
head(geoltax.full)

dim(geoltax.full)[2]
geol.xtabs <- matrix(data=NA, ncol=6, nrow=dim(geoltax.full)[2] - 6)
for (x in 7:(length(colnames(geoltax.full[7:dim(geoltax.full)[2]]))+6)) {
  geol.xtabs[x-6,1] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-6])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-6])))[2])
  geol.xtabs[x-6,2] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-5])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-5])))[2])
  geol.xtabs[x-6,3] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-4])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-4])))[2])
  geol.xtabs[x-6,4] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-3])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-3])))[2])
  geol.xtabs[x-6,5] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-2])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-2])))[2])
  geol.xtabs[x-6,6] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-1])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-1])))[2])
}

colnames(geol.xtabs) <- c("cng","snd","slt","shl","mud","lim")
rownames(geol.xtabs) <- colnames(geoltax.full[7:dim(geoltax.full)[2]])
geol.xtabs

write.table(geol.xtabs, "geolxtabs.csv", sep=",", row.names = T)


##############################################################################
## redo permutation analysis of variance, replacing epoch with age estimates
## permutation analysis of variance (Gower index)
## stochastic loop to select dominant geology class randomly
# combine taxa, geology, period
ageyr <- apply(dat[,c(which(colnames(dat)=="FAD"), which(colnames(dat)=="LAD"))],
               MARGIN=1, mean, na.rm=T)

dat.reclass2 <- data.frame("ID"=dat$ID, "ageyr"=ageyr, "geolcat"=geol.comb, dat[,taxa.sub])
head(dat.reclass2)
dat.reclass.red2 <- dat.reclass2[-which(is.na(dat.reclass2$ageyr)==T | is.na(dat.reclass2$geolcat)==T),]
head(dat.reclass.red2)
dim(dat.reclass.red2)
table(dat.reclass.red2$geolcat)

## run permutation linear models for community presence/absence data
dat.reclass.red.zero2 <- dat.reclass.red2
dat.reclass.red.zero2[is.na(dat.reclass.red.zero2) == T] <- 0 # treat taxa NAs as absence (0)
tax.reclass2 <- dat.reclass.red.zero2[,4:dim(dat.reclass.red.zero2)[2]]
pergeo.reclass2 <- dat.reclass.red.zero2[,2:3]

iter <- 1000
itdiv <- iter/10
permut <- 10000
st.time <- Sys.time()

fact.levels <- c("geolcatcng", "geolcatsnd", "geolcatslt", "geolcatshl", "geolcatmud", "geolcatlim")
factn <- length(fact.levels)

factlev.pmat <- factlev.R2mat <- matrix(data=NA, nrow=iter, ncol=factn)
colnames(factlev.R2mat) <- fact.levels
colnames(factlev.pmat) <- fact.levels
p.mat <- R2.mat <- matrix(data=NA, nrow=iter, ncol=3)

for (s in 1:iter) {
  
  # choose 'dominant' geo class randomly
  pergeo.reclass.it <- pergeo.reclass2
  for (d in 1:dim(pergeo.reclass.it)[1]) {
    geo.split.vec <- as.vector(str_split(gsub("(.{3})", "\\1 ", pergeo.reclass.it$geolcat[d]), " ")[[1]])
    geosplit.vec <- geo.split.vec[-which(geo.split.vec == "")]
    pergeo.reclass.it$geolcat[d] <- sample(geosplit.vec,1)
  } # end d loop
  
  # resample dataset with replacement
  resamp.sub <- sort(sample(1:dim(pergeo.reclass.it)[1], round(dim(pergeo.reclass.it)[1]/1.5,0), replace=F))
  dat.resamp1 <- data.frame(pergeo.reclass.it, tax.reclass2)
  dat.resamp <- dat.resamp1[resamp.sub,] 
  pergeo.resamp <- dat.resamp[,1:2]
  tax.resamp <- dat.resamp[,3:dim(dat.resamp)[2]]
  tax.resamp <- mutate_all(tax.resamp, function(x) as.integer(as.character(x)))
  
  # permutation analysis of variance (Cao index)
  div.iter <- adonis2(tax.resamp ~ ageyr*geolcat, data = pergeo.resamp,
                      permutations = permut, method="gower", binary=T, sqrt.dist=T, by="terms", parallel=cl)
  div.iter2 <- adonis2(tax.resamp ~ ageyr*geolcat, data = pergeo.resamp,
                       permutations = permut, method="gower", binary=T, sqrt.dist=T, by="onedf", parallel=cl)
  
  # store single-level p and R2 for factors
  for (n in 1:factn) {
    table.sub <- which(rownames(div.iter2[1]) == fact.levels[n])
    if (length(table.sub) > 0) {
      factlev.pmat[s, n] <- div.iter2$`Pr(>F)`[table.sub]
      factlev.R2mat[s, n] <- div.iter2$R2[table.sub]
    }
  }
  
  # save pr(F) results
  p.mat[s, ] <- c(div.iter$'Pr(>F)'[1], div.iter$'Pr(>F)'[2], div.iter$'Pr(>F)'[3])
  R2.mat[s, ] <- c(div.iter$R2[1], div.iter$R2[2], div.iter$R2[3])
  
  if (s %% itdiv==0) print(s) 
  
} # end s

end.time <- Sys.time()
proc.time <- end.time - st.time
proc.time

# set CI lims
# 0.67, 0.8, 0.95, 0.99
# i.e., 0.3333, 0.2, 0.05, 0.01
oneINx <- 20 # 3 # 5 # 20 # 100
CI.prob <- 1/oneINx
print(paste(round((1 - CI.prob)*100, 0), "% confidence interval", sep=""))
up.lim <- (1-CI.prob/2)
lo.lim <- 1 - up.lim

## age
# probabilities
exp(median(log(p.mat[, 1])))
print(c(exp(quantile(log(p.mat[, 1]), probs=lo.lim, na.rm=T)), exp(quantile(log(p.mat[, 1]), probs=up.lim, na.rm=T))))
hist(log(p.mat[, 1]), main="", xlab="log age prob")
abline(v=log(0.05), lty=2, lwd=2, col="red")
abline(v=log(0.05), lty=1, lwd=3, col="red")
abline(v=(median(log(p.mat[, 1]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(p.mat[, 1]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(p.mat[, 1]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

# R2
exp(median(log(R2.mat[, 1])))
print(c(exp(quantile(log(R2.mat[, 1]), probs=lo.lim, na.rm=T)), exp(quantile(log(R2.mat[, 1]), probs=up.lim, na.rm=T))))
hist(log(R2.mat[, 1]), main="", xlab="log period R2")
abline(v=(median(log(R2.mat[, 1]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(R2.mat[, 1]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(R2.mat[, 1]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

## geology 
# probabilities
exp(median(log(p.mat[, 2])))
print(c(exp(quantile(log(p.mat[, 2]), probs=lo.lim, na.rm=T)), exp(quantile(log(p.mat[, 2]), probs=up.lim, na.rm=T))))
hist(log(p.mat[, 2]), main="", xlab="log geol prob")
abline(v=log(0.05), lty=2, lwd=2, col="red")
abline(v=log(0.05), lty=1, lwd=3, col="red")
abline(v=(median(log(p.mat[, 2]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(p.mat[, 2]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(p.mat[, 2]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

# R2
exp(median(log(R2.mat[, 2])))
print(c(exp(quantile(log(R2.mat[, 2]), probs=lo.lim, na.rm=T)), exp(quantile(log(R2.mat[, 2]), probs=up.lim, na.rm=T))))
hist(log(R2.mat[, 2]), main="", xlab="log period R2")
abline(v=(median(log(R2.mat[, 2]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(R2.mat[, 2]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(R2.mat[, 2]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

## interaction 
# probabilities
exp(median(log(p.mat[, 3])))
print(c(exp(quantile(log(p.mat[, 3]), probs=lo.lim, na.rm=T)), exp(quantile(log(p.mat[, 3]), probs=up.lim, na.rm=T))))
hist(log(p.mat[, 3]), main="", xlab="log age*geology prob")
abline(v=log(0.05), lty=1, lwd=3, col="red")
abline(v=(median(log(p.mat[, 3]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(p.mat[, 3]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(p.mat[, 3]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

# R2
exp(median(log(R2.mat[, 3])))
print(c(exp(quantile(log(R2.mat[, 3]), probs=lo.lim, na.rm=T)), exp(quantile(log(R2.mat[, 3]), probs=up.lim, na.rm=T))))
hist(log(R2.mat[, 3]), main="", xlab="log period R2")
abline(v=(median(log(R2.mat[, 3]))), col="darkgreen", lwd=3, lty=1)
abline(v=quantile(log(R2.mat[, 3]), probs=lo.lim, na.rm=T), lwd=3, lty=3, col="purple")
abline(v=quantile(log(R2.mat[, 3]), probs=up.lim, na.rm=T), lwd=3, lty=3, col="purple")

# single factor-level median probabilities & R2
fact.levels.sname <- c("cng", "snd", "slt", "shl", "mud", "lim")
factlev.out <- data.frame(fact.levels.sname, apply(factlev.pmat, MARGIN=2, median, na.rm=T),
                          apply(factlev.R2mat, MARGIN=2, median, na.rm=T))
colnames(factlev.out) <- c("factlev", "p", "R2")
rownames(factlev.out) <- NULL
factlev.sort <- factlev.out[order(factlev.out[,3],decreasing=T),]
factlev.sort

# percentages per taxon & age & geology
# age
dim(dat.reclass2)[2]

# vertebrates
vertcols2 <- 4:13
vert.reclass2 <- dat.reclass2[,vertcols2]
grp.xtabs.names <- paste(colnames(vert.reclass2), ".ageXtabs", sep="")

for (g in 1:length(grp.xtabs.names)) {
  grp.it <- as.data.frame(xtabs(vert.reclass2[,g] ~ dat.reclass2[,2]))
  colnames(grp.it) <- c("ageyr","freq")
  grp.it$ageyr <- as.double(as.character(grp.it$ageyr))
  grp.it$freq <- ifelse(grp.it$freq == 0, "NA", grp.it$freq)
  plot(grp.it$ageyr, grp.it$freq, pch=19)
  assign(grp.xtabs.names[g], grp.it)
  write.table(get(grp.xtabs.names[g]), paste(grp.xtabs.names[g], ".csv", sep=""), sep=",", row.names = F)
}

# invertebrates
invcols2 <- 14:33
inv.reclass2 <- dat.reclass2[,invcols2]
grp.xtabs.names <- paste(colnames(inv.reclass2), ".ageXtabs", sep="")

for (g in 1:length(grp.xtabs.names)) {
  grp.it <- as.data.frame(xtabs(inv.reclass2[,g] ~ dat.reclass2[,2]))
  colnames(grp.it) <- c("ageyr","freq")
  grp.it$ageyr <- as.double(as.character(grp.it$ageyr))
  grp.it$freq <- ifelse(grp.it$freq == 0, "NA", grp.it$freq)
  plot(grp.it$ageyr, grp.it$freq, pch=19)
  assign(grp.xtabs.names[g], grp.it)
  write.table(get(grp.xtabs.names[g]), paste(grp.xtabs.names[g], ".csv", sep=""), sep=",", row.names = F)
}

# plants
plntcols2 <- 34:40
plnt.reclass2 <- dat.reclass2[,plntcols2]
grp.xtabs.names <- paste(colnames(plnt.reclass2), ".ageXtabs", sep="")

for (g in 1:length(grp.xtabs.names)) {
  grp.it <- as.data.frame(xtabs(plnt.reclass2[,g] ~ dat.reclass2[,2]))
  colnames(grp.it) <- c("ageyr","freq")
  grp.it$ageyr <- as.double(as.character(grp.it$ageyr))
  grp.it$freq <- ifelse(grp.it$freq == 0, "NA", grp.it$freq)
  plot(grp.it$ageyr, grp.it$freq, pch=19)
  assign(grp.xtabs.names[g], grp.it)
  write.table(get(grp.xtabs.names[g]), paste(grp.xtabs.names[g], ".csv", sep=""), sep=",", row.names = F)
}


# geology
geolcols <- c(15,17,19,21,23,25)
geoltax.full <- data.frame(dat[,geolcols], dat[,vertcols], dat[,invcols], dat[,plntcols[-length(plntcols)]])
colnames(geoltax.full) <- c(colnames(dat[,geolcols], ), colnames(dat[,vertcols], ), colnames(dat[,invcols], ),
                            colnames(dat[,plntcols[-length(plntcols)]], ))
head(geoltax.full)

dim(geoltax.full)[2]
geol.xtabs <- matrix(data=NA, ncol=6, nrow=dim(geoltax.full)[2] - 6)
for (x in 7:(length(colnames(geoltax.full[7:dim(geoltax.full)[2]]))+6)) {
  geol.xtabs[x-6,1] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-6])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-6])))[2])
  geol.xtabs[x-6,2] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-5])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-5])))[2])
  geol.xtabs[x-6,3] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-4])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-4])))[2])
  geol.xtabs[x-6,4] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-3])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-3])))[2])
  geol.xtabs[x-6,5] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-2])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-2])))[2])
  geol.xtabs[x-6,6] <- as.numeric((xtabs(geoltax.full[,x] ~ geoltax.full[,x-1])/
                                     sum(xtabs(geoltax.full[,x] ~ geoltax.full[,x-1])))[2])
}

colnames(geol.xtabs) <- c("cng","snd","slt","shl","mud","lim")
rownames(geol.xtabs) <- colnames(geoltax.full[7:dim(geoltax.full)[2]])
geol.xtabs

write.table(geol.xtabs, "geolxtabs.csv", sep=",", row.names = T)
