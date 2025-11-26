library("ape")
library("maps")
library("paleotree")
library("geoscale")
library("strap")
library("phytools")
library("TreeTools")
remotes::install_github("uyedaj/treeplyr")
library("treeplyr")
library("ontologyIndex")
library("shiny")
library("visNetwork")
library("ontoFAST")
library("tibble")
library("stringr")
library("phangorn")
library("rphenoscape")
library("rphenoscate")
library("treesurgeon")

# import head matrix 1
nex1 <- ReadCharacters(file = "cranial.NEX")
mat1 <- bind_cols(taxon = rownames(nex1), as_tibble(nex1)) 
mat1 <- setNames(mat1, c("taxon", colnames(nex1)))

# import trunk matrix 1
nex2 <- ReadCharacters(file = "trunk.NEX")
mat2 <- bind_cols(taxon = rownames(nex2), as_tibble(nex2))
mat2 <- setNames(mat2, c("taxon", colnames(nex2))) 

# import head matrix 2
nex3 <- read.nexus.data(file = "cranial2.NEX")
mat3 <- as_tibble(do.call(rbind, nex3))  
mat3 <- mat3 %>%
  mutate(taxon = names(nex3)) %>%       
  select(taxon, everything())          
mat32 <- mat3

# import trunk matrix 2
nex4 <- read.nexus.data(file = "trunk2.nex")
mat4 <- as_tibble(do.call(rbind, nex4))   
mat4 <- mat4 %>%
  mutate(taxon = names(nex4)) %>%         
  select(taxon, everything())            
mat42 <- mat4

#importing trees
tree1 <- read.nexus(file = "tree.tre")
plot.phylo(tree1, cex = 0.4)
tree2 <- read.nexus(file = "tree2.tree")
plot.phylo(tree2, cex = 0.4)
tree3 <- read.nexus(file = "tree3.tree")
plot.phylo(tree3, cex = 0.4)

#combine 'tree' and 'matrix' into 'treedata'
attributes(mat1)$row.names <- mat1$taxon
mat1 <- mat1[match(tree1$tip.label, rownames(mat1)),] 
td1 <- make.treedata(tree1, mat1, name_column = "taxon", as.is = TRUE) 
attributes(td1$dat)$row.names <- td1$phy$tip.label

attributes(mat2)$row.names <- mat2$taxon
mat2 <- mat2[match(tree1$tip.label, rownames(mat2)),] 
mat2 <- mat2[, !is.na(colnames(mat2)) & colnames(mat2) != ""] 
td2 <- make.treedata(tree1, mat2, name_column = "taxon", as.is = TRUE) 
attributes(td2$dat)$row.names <- td2$phy$tip.label   

#head2 linked to tree2'
mat3$taxon <- gsub("_", " ", mat3$taxon)
attributes(mat3)$row.names <- mat3$taxon
tree2$tip.label <- gsub("_", " ", tree2$tip.label)
mat3 <- mat3[match(tree2$tip.label, rownames(mat3)),] 
td3 <- make.treedata(tree2, mat3, name_column = "taxon", as.is = TRUE) 
attributes(td3$dat)$row.names <- td3$phy$tip.label  

#head2 linked to tree3'
mat32$taxon <- gsub("_", " ", mat32$taxon)
attributes(mat32)$row.names <- mat32$taxon
tree3$tip.label <- gsub("_", " ", tree3$tip.label)
mat32 <- mat32[match(tree3$tip.label, rownames(mat32)),] 
td4 <- make.treedata(tree3, mat32, name_column = "taxon", as.is = TRUE) 
attributes(td4$dat)$row.names <- td4$phy$tip.label   

#trunk2 linked to tree2'
mat4$taxon <- gsub("_", " ", mat4$taxon)
attributes(mat4)$row.names <- mat4$taxon
mat4 <- mat4[match(tree2$tip.label, rownames(mat4)),] 
td5 <- make.treedata(tree2, mat4, name_column = "taxon", as.is = TRUE) 
attributes(td5$dat)$row.names <- td5$phy$tip.label   

#trunk2 linked to tree3'
mat42$taxon <- gsub("_", " ", mat42$taxon)
attributes(mat42)$row.names <- mat42$taxon
mat42 <- mat42[match(tree3$tip.label, rownames(mat42)),] 
td6 <- make.treedata(tree3, mat42, name_column = "taxon", as.is = TRUE) 
attributes(td6$dat)$row.names <- td6$phy$tip.label

#time scale trees
data <- read.csv("tip.age.csv", header = FALSE)
rowname <- data$V1
First.inst <- data$V2
Last.inst <- data$V3
tip_ages <- structure(list(First.inst = First.inst, Last.inst = Last.inst), row.names = rowname,  class = "data.frame")
tip_ages <- tip_ages[match(tree1$tip.label, rownames(tip_ages)),]
timetrees1 <- tree1
timetrees_scaled1 <- tree1
timetrees1 <- timePaleoPhy(timetrees1, timeData = tip_ages[, 1:2], type = "equal", vartime = 5)
timetrees_scaled1$edge.length <- timetrees1$edge.length / max(node.depth.edgelength(timetrees1))
geoscalePhylo(timetrees1, units = "Period", boxes = "Period", cex.tip = 0.4, width = 1.0)

data2 <- read.csv("tip.age2.csv", header = FALSE)
rowname2 <- data2$V1
First.inst2 <- data2$V2
Last.inst2 <- data2$V3
tip_ages2 <- structure(list(First.inst = First.inst2, Last.inst = Last.inst2), row.names = rowname2,  class = "data.frame")
tree2$tip.label <- gsub("_", " ", tree2$tip.label)
match(tree2$tip.label, rownames(tip_ages2))
timetrees2 <- tree2
timetrees_scaled2 <- tree2
timetrees2 <- timePaleoPhy(timetrees2, timeData = tip_ages2[, 1:2], type = "equal", vartime = 5)
timetrees_scaled2$edge.length <- timetrees2$edge.length / max(node.depth.edgelength(timetrees2))
geoscalePhylo(timetrees2, units = "Period", boxes = "Period", cex.tip = 0.4, width = 1.0)

data3 <- read.csv("tip.age2.csv", header = FALSE)
rowname3 <- data2$V1
First.inst3 <- data2$V2
Last.inst3 <- data2$V3
tip_ages3 <- structure(list(First.inst = First.inst3, Last.inst = Last.inst3), row.names = rowname3,  class = "data.frame")
tree3$tip.label <- gsub("_", " ", tree3$tip.label)
match(tree3$tip.label, rownames(tip_ages3))
timetrees3 <- tree3
timetrees_scaled3 <- tree3
timetrees3 <- timePaleoPhy(timetrees3, timeData = tip_ages3[, 1:2], type = "equal", vartime = 5)
timetrees_scaled3$edge.length <- timetrees3$edge.length / max(node.depth.edgelength(timetrees3))
geoscalePhylo(timetrees3, units = "Period", boxes = "Period", cex.tip = 0.4, width = 1.0)

#define steps
D <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE, dimnames =list( c("D*", "D"), c("D*", "D")) )
E <- matrix(c(0, 2, 2, 0), 2, 2, byrow = TRUE, dimnames =list( c("E*", "E"), c("E*", "E")) )
L <- matrix(c(0, 3, 3, 0), 2, 2, byrow = TRUE, dimnames =list( c("L1", "LM"), c("L1", "LM")) )
C <- matrix(c(0, 4, 4, 4, 0, 4, 4, 4, 0), 3, 3, byrow = TRUE, dimnames =list( c("C1", "C2", "C3"), c("C1", "C2","C3")) )
R <- matrix(c(0, 5, 5, 0), 2, 2, byrow = TRUE, dimnames =list( c("R*", "R"), c("R*", "R")) )

step1 <- amaSMM(E, L)
step2 <- amaSMM(step1, C)
step3 <- amaSMM(step2, R)
step4 <- amaED(Qc = D, Qd = step3, type = "ql", diag = 0)

rownames(step4)

#convert character matrix into tip priors 
tp <- get_tip_priors(nex1)  
tp2 <- get_tip_priors(nex2)
tp3 <- get_tip_priors(nex3) 
tp4 <- get_tip_priors(nex3) 
tp5 <- get_tip_priors(nex4) 
tp6 <- get_tip_priors(nex4) 

tps <- list(tp, tp2, tp3, tp4, tp5, tp6)
col_lists <- list(D, E, L, C, R)
for (k in seq_along(tps)) {
  for (i in seq_along(col_lists)) {
    colnames(tps[[k]][[i]]) <- colnames(col_lists[[i]])
  }
}
tp  <- tps[[1]]
tp2 <- tps[[2]]
tp3 <- tps[[3]]
tp4 <- tps[[4]]
tp5 <- tps[[5]]
tp6 <- tps[[6]]

#amalgamate tip priors in the same order
tp11 <- amal_tip_priors(tp_1 = tp[[2]], tp_2 = tp[[3]], type = "SMM")
tp12 <- amal_tip_priors(tp_1 = tp[[5]], tp_2 = tp11, type = "SMM")
tp13 <- amal_tip_priors(tp_1 = tp[[4]], tp_2 = tp12, type = "SMM")
tp14 <- amal_tip_priors(tp_1 = tp[[1]], tp_2 = tp13, type = "ED")

tp21 <- amal_tip_priors(tp_1 = tp2[[2]], tp_2 = tp2[[3]], type = "SMM")
tp22 <- amal_tip_priors(tp_1 = tp2[[5]], tp_2 = tp21, type = "SMM")
tp23 <- amal_tip_priors(tp_1 = tp2[[4]], tp_2 = tp22, type = "SMM")
tp24 <- amal_tip_priors(tp_1 = tp2[[1]], tp_2 = tp23, type = "ED")

tp31 <- amal_tip_priors(tp_1 = tp3[[2]], tp_2 = tp3[[3]], type = "SMM")
tp32 <- amal_tip_priors(tp_1 = tp3[[5]], tp_2 = tp31, type = "SMM")
tp33 <- amal_tip_priors(tp_1 = tp3[[4]], tp_2 = tp32, type = "SMM")
tp34 <- amal_tip_priors(tp_1 = tp3[[1]], tp_2 = tp33, type = "ED")

tp41 <- amal_tip_priors(tp_1 = tp4[[2]], tp_2 = tp4[[3]], type = "SMM")
tp42 <- amal_tip_priors(tp_1 = tp4[[5]], tp_2 = tp41, type = "SMM")
tp43 <- amal_tip_priors(tp_1 = tp4[[4]], tp_2 = tp42, type = "SMM")
tp44 <- amal_tip_priors(tp_1 = tp4[[1]], tp_2 = tp43, type = "ED")

tp51 <- amal_tip_priors(tp_1 = tp5[[2]], tp_2 = tp5[[3]], type = "SMM")
tp52 <- amal_tip_priors(tp_1 = tp5[[5]], tp_2 = tp51, type = "SMM")
tp53 <- amal_tip_priors(tp_1 = tp5[[4]], tp_2 = tp52, type = "SMM")
tp54 <- amal_tip_priors(tp_1 = tp5[[1]], tp_2 = tp53, type = "ED")

tp61 <- amal_tip_priors(tp_1 = tp6[[2]], tp_2 = tp6[[3]], type = "SMM")
tp62 <- amal_tip_priors(tp_1 = tp6[[5]], tp_2 = tp61, type = "SMM")
tp63 <- amal_tip_priors(tp_1 = tp6[[4]], tp_2 = tp62, type = "SMM")
tp64 <- amal_tip_priors(tp_1 = tp6[[1]], tp_2 = tp63, type = "ED")

#fit model
fitRes <- fitMk(tree = timetrees1, x = tp14, model = step4)
fitRes2 <- fitMk(tree = timetrees1, x = tp24, model = step4)
rownames(tp34) <- gsub("_", " ", rownames(tp34))
fitRes3 <- fitMk(tree = timetrees2, x = tp34, model = step4)
rownames(tp44) <- gsub("_", " ", rownames(tp44))
fitRes4 <- fitMk(tree = timetrees3, x = tp44, model = step4)
rownames(tp54) <- gsub("_", " ", rownames(tp54))
fitRes5 <- fitMk(tree = timetrees2, x = tp54, model = step4)
rownames(tp64) <- gsub("_", " ", rownames(tp64))
fitRes6 <- fitMk(tree = timetrees3, x = tp64, model = step4)

#ancestral state estimation 
ancRes <- ancr(fitRes)   
ancRes2 <- ancr(fitRes2) 
ancRes3 <- ancr(fitRes3) 
ancRes4 <- ancr(fitRes4) 
ancRes5 <- ancr(fitRes5) 
ancRes6 <- ancr(fitRes6) 

library(randomcoloR)
cols <- distinctColorPalette(ncol(ancRes$ace))

#plot all states 
#head dataset1 timetree1
pdf("headtree1.pdf", width = 8.27, height = 11.69)
plot(ancRes, args.plotTree = list(fsize = 0.3, mar=c(0.1,1.1,2.1,1.1)), args.tiplabels=list(cex = 0.3), args.nodelabels=list(piecol=cols, cex = 0.4), legend = "topright")
title("cranial")
mtext(
  "D: Dermal skeleton; C1: Vascular canals state 1; C2: Vascular canals state 2; C3: Vascular canals state 3; R: Odontogenic resorption; E: Enamel; L1: single-layered enamel; LM: Multi-layered enamel; *: absent",
  side = 1,    
  line = 6,   
  cex = 0.5
)
dev.off()

#plot result
#enamel layer
pdf("headtree1 enamellayer.pdf", width = 8.27, height = 11.69)
Eanc <- as.matrix(data.frame(absent = rowSums(ancRes$ace[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(ancRes$ace[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(ancRes$ace[,c(5, 9, 13, 17, 21, 25)])))
Etip <- as.matrix(data.frame(absent = rowSums(tp14[,c(1, 2, 3, 6, 7, 10, 11, 14, 15)]), singlelayered = rowSums(tp14[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(tp14[,c(5, 9, 13, 17, 21, 25)])))
Etip_aligned <- Etip[tree1$tip.label, ]
geoscalePhylo(timetrees1, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 1.0)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "purple")
nodelabels(pie = Eanc, piecol = cols, cex = 0.2)
tiplabels(pie = Etip_aligned, piecol = cols, cex = 0.2)
title("Enamel (cranial dermal skeleton, GAP clade stem-Sarcop.)", cex.main = 0.8)
legend_labels <- c("absent", "single-layered", "multi-layered")
legend("topright", legend = legend_labels, pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#resorption
pdf("headtree1 resorption.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame(absent = rowSums(ancRes$ace[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(ancRes$ace[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip <- as.matrix(data.frame(absent = rowSums(tp14[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(tp14[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip_aligned <- Rtip[tree1$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees1, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange")
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Odontogenic resorption (cranial dermal skeleton, GAP clade stem-Sarcop.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#canal system
pdf("headtree1 canalsystem.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame("a dentinal system" = rowSums(ancRes$ace[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(ancRes$ace[,c(10:17)]), "an additional true pore canal system" = rowSums(ancRes$ace[,c(18:25)])))
Rtip <- as.matrix(data.frame("a dentinal system" = rowSums(tp14[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(tp14[,c(10:17)]), "an additional true pore canal system" = rowSums(tp14[,c(18:25)])))
Rtip_aligned <- Rtip[tree1$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees1, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "purple")
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Vascular canals (cranial dermal skeleton, GAP clade stem-Sarcop.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#cosmine
pdf("headtree1 cosmine.pdf", width = 8.27, height = 11.69)
Canc <- as.matrix(data.frame(absent = rowSums(ancRes$ace[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(ancRes$ace[,c(12, 13, 16, 17)]), "strict definition" = rowSums(ancRes$ace[,c(20, 21, 24, 25)])))
Ctip <- as.matrix(data.frame(absent = rowSums(tp14[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(tp14[,c(12, 13, 16, 17)]), "strict definition" = rowSums(tp14[,c(20, 21, 24, 25)])))
Ctip[Ctip >= 1] <- 1
Ctip_aligned <- Ctip[tree1$tip.label, ]
cols <- c("lightblue", "orange", "purple")
geoscalePhylo(timetrees1, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
nodelabels(pie = Canc, piecol = cols, cex = 0.2)
tiplabels(pie = Ctip_aligned, piecol = cols, cex = 0.2)
title("Cosmine (cranial dermal skeleton, GAP clade stem-Sarcop.)", cex.main = 0.8)
legend("topright", legend = colnames(Canc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#trunk dataset1 timetree1
#enamel layer
pdf("trunktree1 enamellayer.pdf", width = 8.27, height = 11.69)
Eanc <- as.matrix(data.frame(absent = rowSums(ancRes2$ace[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(ancRes2$ace[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(ancRes2$ace[,c(5, 9, 13, 17, 21, 25)])))
Etip <- as.matrix(data.frame(absent = rowSums(tp24[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(tp24[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(tp24[,c(5, 9, 13, 17, 21, 25)])))
Etip_aligned <- Etip[tree1$tip.label, ]
geoscalePhylo(timetrees1, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 1.0)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "mediumpurple1")
nodelabels(pie = Eanc, piecol = cols, cex = 0.2)
tiplabels(pie = Etip_aligned, piecol = cols, cex = 0.2)
title("Enamel (postcranial dermal skeleton, GAP clade stem-Sarcop.)", cex.main = 0.8)
legend_labels <- c("absent", "single-layered", "multi-layered")
legend("topright", legend = legend_labels, pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#canal system
pdf("trunktree1 canalsystem.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame("a dentinal system" = rowSums(ancRes2$ace[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(ancRes2$ace[,c(10:17)]), "an additional true pore canal system" = rowSums(ancRes2$ace[,c(18:25)])))
Rtip <- as.matrix(data.frame("a dentinal system" = rowSums(tp24[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(tp24[,c(10:17)]), "an additional true pore canal system" = rowSums(tp24[,c(18:25)])))
Rtip_aligned <- Rtip[tree1$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees1, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "purple")
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Vascular canals (postcranial dermal skeleton, GAP clade stem-Sarcop.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#resorption
pdf("trunktree1 resorption.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame(absent = rowSums(ancRes2$ace[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(ancRes2$ace[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip <- as.matrix(data.frame(absent = rowSums(tp24[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(tp24[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip_aligned <- Rtip[tree1$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees1, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange")
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Odontogenic resorption (postcranial dermal skeleton, GAP clade stem-Sarcop.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#cosmine
pdf("trunktree1 cosmine.pdf", width = 8.27, height = 11.69)
Canc <- as.matrix(data.frame(absent = rowSums(ancRes2$ace[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(ancRes2$ace[,c(12, 13, 16, 17)]), "strict definition" = rowSums(ancRes2$ace[,c(20, 21, 24, 25)])))
Ctip <- as.matrix(data.frame(absent = rowSums(tp24[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(tp24[,c(12, 13, 16, 17)]), "strict definition" = rowSums(tp24[,c(20, 21, 24, 25)])))
Ctip[Ctip >= 1] <- 1
Ctip_aligned <- Ctip[tree1$tip.label, ]
cols <- c("lightblue", "orange", "purple")
geoscalePhylo(timetrees1, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
nodelabels(pie = Canc, piecol = cols, cex = 0.2)
tiplabels(pie = Ctip_aligned, piecol = cols, cex = 0.2)
title("Cosmine (postcranial dermal skeleton, GAP clade stem-Sarcop.)", cex.main = 0.8)
legend("topright", legend = colnames(Canc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#head dataset2 timetree2
#enamel layer
pdf("head2tree2 enamellayer.pdf", width = 8.27, height = 11.69)
Eanc <- as.matrix(data.frame(absent = rowSums(ancRes3$ace[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(ancRes3$ace[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(ancRes3$ace[,c(5, 9, 13, 17, 21, 25)])))
Etip <- as.matrix(data.frame(absent = rowSums(tp34[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(tp34[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(tp34[,c(5, 9, 13, 17, 21, 25)])))
Etip_aligned <- Etip[tree2$tip.label, ]
geoscalePhylo(timetrees2, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 1.0)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "mediumpurple1")
nodelabels(pie = Eanc, piecol = cols, cex = 0.2)
tiplabels(pie = Etip_aligned, piecol = cols, cex = 0.2)
title("Enamel (cranial dermal skeleton, GAP clade stem-Actinop.)", cex.main = 0.8)
legend_labels <- c("absent", "single-layered", "multi-layered")
legend("topright", legend = legend_labels, pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#resorption
pdf("head2tree2 resorption.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame(absent = rowSums(ancRes3$ace[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(ancRes3$ace[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip <- as.matrix(data.frame(absent = rowSums(tp34[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(tp34[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip_aligned <- Rtip[tree2$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees2, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange")
par(lwd = 0.5)
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Odontogenic resorption (cranial dermal skeleton, GAP clade stem-Actinop.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#canal system
pdf("head2tree2 canalsystem.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame("a dentinal system" = rowSums(ancRes3$ace[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(ancRes3$ace[,c(10:17)]), "an additional true pore canal system" = rowSums(ancRes3$ace[,c(18:25)])))
Rtip <- as.matrix(data.frame("a dentinal system" = rowSums(tp34[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(tp34[,c(10:17)]), "an additional true pore canal system" = rowSums(tp34[,c(18:25)])))
Rtip_aligned <- Rtip[tree2$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees2, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "purple")
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Vascular canals (cranial dermal skeleton, GAP clade stem-Actinop.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#cosmine
pdf("head2tree2 cosmine.pdf", width = 8.27, height = 11.69)
Canc <- as.matrix(data.frame(absent = rowSums(ancRes3$ace[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(ancRes3$ace[,c(12, 13, 16, 17)]), "strict definition" = rowSums(ancRes3$ace[,c(20, 21, 24, 25)])))
Ctip <- as.matrix(data.frame(absent = rowSums(tp34[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(tp34[,c(12, 13, 16, 17)]), "strict definition" = rowSums(tp34[,c(20, 21, 24, 25)])))
Ctip[Ctip >= 1] <- 1
Ctip_aligned <- Ctip[tree2$tip.label, ]
cols <- c("lightblue", "orange", "purple")
geoscalePhylo(timetrees2, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
nodelabels(pie = Canc, piecol = cols, cex = 0.2)
tiplabels(pie = Ctip_aligned, piecol = cols, cex = 0.2)
title("Cosmine (cranial dermal skeleton, GAP clade stem-Actinop.)", cex.main = 0.8)
legend("topright", legend = colnames(Canc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#head dataset2 timetree3
#enamel layer
pdf("headtree3 enamellayer.pdf", width = 8.27, height = 11.69)
Eanc <- as.matrix(data.frame(absent = rowSums(ancRes4$ace[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(ancRes4$ace[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(ancRes4$ace[,c(5, 9, 13, 17, 21, 25)])))
Etip <- as.matrix(data.frame(absent = rowSums(tp44[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(tp44[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(tp44[,c(5, 9, 13, 17, 21, 25)])))
Etip_aligned <- Etip[tree3$tip.label, ]
geoscalePhylo(timetrees3, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 1.0)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "mediumpurple1")
nodelabels(pie = Eanc, piecol = cols, cex = 0.2)
tiplabels(pie = Etip_aligned, piecol = cols, cex = 0.2)
title("Enamel (cranial dermal skeleton, GAP clade stem-Ostei.)", cex.main = 0.8)
legend_labels <- c("absent", "single-layered", "multi-layered")
legend("topright", legend = legend_labels, pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#canal system
pdf("headtree3 canalsystem.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame("a dentinal system" = rowSums(ancRes4$ace[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(ancRes4$ace[,c(10:17)]), "an additional true pore canal system" = rowSums(ancRes4$ace[,c(18:25)])))
Rtip <- as.matrix(data.frame("a dentinal system" = rowSums(tp44[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(tp44[,c(10:17)]), "an additional true pore canal system" = rowSums(tp44[,c(18:25)])))
Rtip_aligned <- Rtip[tree3$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees3, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "purple")
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Vascular canals (cranial dermal skeleton, GAP clade stem-Ostei.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#resorption
pdf("head2tree3 resorption.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame(absent = rowSums(ancRes4$ace[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(ancRes4$ace[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip <- as.matrix(data.frame(absent = rowSums(tp44[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(tp44[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip_aligned <- Rtip[tree3$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees3, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange")
par(lwd = 0.5)
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Odontogenic resorption (cranial dermal skeleton, GAP clade stem-Ostei.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#cosmine
pdf("head2tree3 cosmine.pdf", width = 8.27, height = 11.69)
Canc <- as.matrix(data.frame(absent = rowSums(ancRes4$ace[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(ancRes4$ace[,c(12, 13, 16, 17)]), "strict definition" = rowSums(ancRes4$ace[,c(20, 21, 24, 25)])))
Ctip <- as.matrix(data.frame(absent = rowSums(tp44[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(tp44[,c(12, 13, 16, 17)]), "strict definition" = rowSums(tp44[,c(20, 21, 24, 25)])))
Ctip[Ctip >= 1] <- 1
Ctip_aligned <- Ctip[tree3$tip.label, ]
cols <- c("lightblue", "orange", "purple")
geoscalePhylo(timetrees3, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
nodelabels(pie = Canc, piecol = cols, cex = 0.2)
tiplabels(pie = Ctip_aligned, piecol = cols, cex = 0.2)
title("Cosmine (cranial dermal skeleton, GAP clade stem-Ostei.)", cex.main = 0.8)
legend("topright", legend = colnames(Canc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#trunk dataset2 timetree2
#enamel layer
pdf("trunk2tree2 enamellayer.pdf", width = 8.27, height = 11.69)
Eanc <- as.matrix(data.frame(absent = rowSums(ancRes5$ace[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(ancRes5$ace[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(ancRes5$ace[,c(5, 9, 13, 17, 21, 25)])))
Etip <- as.matrix(data.frame(absent = rowSums(tp54[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(tp54[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(tp54[,c(5, 9, 13, 17, 21, 25)])))
Etip_aligned <- Etip[tree2$tip.label, ]
geoscalePhylo(timetrees2, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 1.0)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "mediumpurple1")
nodelabels(pie = Eanc, piecol = cols, cex = 0.2)
tiplabels(pie = Etip_aligned, piecol = cols, cex = 0.2)
title("Enamel (postcranial dermal skeleton, GAP clade stem-Actinop.)", cex.main = 0.8)
legend_labels <- c("absent", "single-layered", "multi-layered")
legend("topright", legend = legend_labels, pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#canal system
pdf("trunk2tree2 canalsystem.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame("a dentinal system" = rowSums(ancRes5$ace[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(ancRes5$ace[,c(10:17)]), "an additional true pore canal system" = rowSums(ancRes5$ace[,c(18:25)])))
Rtip <- as.matrix(data.frame("a dentinal system" = rowSums(tp54[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(tp54[,c(10:17)]), "an additional true pore canal system" = rowSums(tp54[,c(18:25)])))
Rtip_aligned <- Rtip[tree2$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees2, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "purple")
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Vascular canals (postcranial dermal skeleton, GAP clade stem-Actinop.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#resorption
pdf("trunk2tree2 resorption.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame(absent = rowSums(ancRes5$ace[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(ancRes5$ace[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip <- as.matrix(data.frame(absent = rowSums(tp54[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(tp54[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip_aligned <- Rtip[tree2$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees2, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange")
par(lwd = 0.5)
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Odontogenic resorption (postcranial dermal skeleton, GAP clade stem-Actinop.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#cosmine
pdf("trunk2tree2 cosmine.pdf", width = 8.27, height = 11.69)
Canc <- as.matrix(data.frame(absent = rowSums(ancRes5$ace[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(ancRes5$ace[,c(12, 13, 16, 17)]), "strict definition" = rowSums(ancRes5$ace[,c(20, 21, 24, 25)])))
Ctip <- as.matrix(data.frame(absent = rowSums(tp54[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(tp54[,c(12, 13, 16, 17)]), "strict definition" = rowSums(tp54[,c(20, 21, 24, 25)])))
Ctip[Ctip >= 1] <- 1
Ctip_aligned <- Ctip[tree2$tip.label, ]
cols <- c("lightblue", "orange", "purple")
geoscalePhylo(timetrees2, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
nodelabels(pie = Canc, piecol = cols, cex = 0.2)
tiplabels(pie = Ctip_aligned, piecol = cols, cex = 0.2)
title("Cosmine (postcranial dermal skeleton, GAP clade stem-Actinop.)", cex.main = 0.8)
legend("topright", legend = colnames(Canc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#trunk dataset2 timetree3
#enamel layer
pdf("trunk2tree3 enamellayer.pdf", width = 8.27, height = 11.69)
Eanc <- as.matrix(data.frame(absent = rowSums(ancRes6$ace[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(ancRes6$ace[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(ancRes6$ace[,c(5, 9, 13, 17, 21, 25)])))
Etip <- as.matrix(data.frame(absent = rowSums(tp64[,c(1, 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23)]), singlelayered = rowSums(tp64[,c(4, 8, 12, 16, 20, 24)]), multilayered = rowSums(tp64[,c(5, 9, 13, 17, 21, 25)])))
Etip_aligned <- Etip[tree3$tip.label, ]
geoscalePhylo(timetrees3, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 1.0)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "mediumpurple1")
nodelabels(pie = Eanc, piecol = cols, cex = 0.2)
tiplabels(pie = Etip_aligned, piecol = cols, cex = 0.2)
title("Enamel (postcranial dermal skeleton, GAP clade stem-Ostei.)", cex.main = 0.8)
legend_labels <- c("absent", "single-layered", "multi-layered")
legend("topright", legend = legend_labels, pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#canal system
pdf("trunk2tree3 canalsystem.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame("a dentinal system" = rowSums(ancRes6$ace[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(ancRes6$ace[,c(10:17)]), "an additional true pore canal system" = rowSums(ancRes6$ace[,c(18:25)])))
Rtip <- as.matrix(data.frame("a dentinal system" = rowSums(tp64[,c(1:9)]), "an additional rudimentary pore canal system" = rowSums(tp64[,c(10:17)]), "an additional true pore canal system" = rowSums(tp64[,c(18:25)])))
Rtip_aligned <- Rtip[tree3$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees3, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
cols <- c("lightblue", "orange", "purple")
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Vascular canals (postcranial dermal skeleton, GAP clade stem-Ostei.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#resorption
pdf("trunk2tree3 resorption.pdf", width = 8.27, height = 11.69)
Ranc <- as.matrix(data.frame(absent = rowSums(ancRes6$ace[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(ancRes6$ace[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip <- as.matrix(data.frame(absent = rowSums(tp64[,c(1, 2, 3, 4, 5, 10, 11, 12, 13, 18, 19, 20, 21)]), present = rowSums(tp64[,c(6, 7, 8, 9, 14, 15, 16, 17, 22, 23, 24, 25)])))
Rtip_aligned <- Rtip[tree3$tip.label, ]
Rtip[Rtip >= 1] <- 1
geoscalePhylo(timetrees3, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
cols <- c("lightblue", "orange")
par(mar=c(0.1,1.1,2.1,1.1))
nodelabels(pie = Ranc, piecol = cols, cex = 0.2)
tiplabels(pie = Rtip_aligned, piecol = cols, cex = 0.2)
title("Odontogenic resorption (postcranial dermal skeleton, GAP clade stem-Ostei.)", cex.main = 0.8)
legend("topright", legend = colnames(Ranc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

#cosmine
pdf("trunk2tree3 cosmine.pdf", width = 8.27, height = 11.69)
Canc <- as.matrix(data.frame(absent = rowSums(ancRes6$ace[,c(1:11, 14, 15, 18, 19, 22, 23)]),  "broad definition" = rowSums(ancRes6$ace[,c(12, 13, 16, 17)]), "strict definition" = rowSums(ancRes6$ace[,c(20, 21, 24, 25)])))
Ctip <- as.matrix(data.frame(absent = rowSums(tp64[,c(1:11, 14, 15, 18, 19, 22, 23)]), "broad definition" = rowSums(tp64[,c(12, 13, 16, 17)]), "strict definition" = rowSums(tp64[,c(20, 21, 24, 25)])))
Ctip[Ctip >= 1] <- 1
Ctip_aligned <- Ctip[tree3$tip.label, ]
cols <- c("lightblue", "orange", "purple")
geoscalePhylo(timetrees3, units = "Period", boxes = "Period", cex.age = 0.5, cex.ts = 1, cex.tip = 0.5, label.offset = 1, width = 0.5)
par(mar=c(0.1,1.1,2.1,1.1))
nodelabels(pie = Canc, piecol = cols, cex = 0.2)
tiplabels(pie = Ctip_aligned, piecol = cols, cex = 0.2)
title("Cosmine (postcranial dermal skeleton, GAP clade stem-Ostei.)", cex.main = 0.8)
legend("topright", legend = colnames(Canc), pch = 21, pt.bg = cols, col = "black", pt.cex = 1.7, cex = 1, y.intersp = 1.5, inset = c(0.05, 0.05))
dev.off()

