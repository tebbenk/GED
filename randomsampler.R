library("stringr")
library("dplyr")
library("plyr")
library("data.table")

#Read in 1:1:1 orthologs and subset to the species you are interested in 
data <- read.delim("orthologs_pbpvpf.txt", header = 1)
data <- data[str_detect(data$Organism, "Plasmodium berghei ANKA|Plasmodium vivax P01|Plasmodium falciparum 3D7"), ]
data <- data[,c(1,3)]

data <- data[str_detect(data$Organism, "Plasmodium berghei ANKA|Plasmodium vivax P01|Plasmodium falciparum 3D7"), ]

#Read in Pf sc dataset 
pf_t <- read.delim("pf_t.txt", header = TRUE)

#Generate mixtures with increasing % trophs mixed with rings 
rings <- subset(pf_t, stage == "ring")
rings_50 <- rings %>% group_by(stage) %>% slice_sample(n = 50)
rings_125 <- rings %>% group_by(stage) %>% slice_sample(n = 125)
rings_250 <- rings %>% group_by(stage) %>% slice_sample(n = 250)
rings_375 <- rings %>% group_by(stage) %>% slice_sample(n = 375)
rings_450 <- rings %>% group_by(stage) %>% slice_sample(n = 450)
rings_500 <- rings %>% group_by(stage) %>% slice_sample(n = 500)

trophs <- subset(pf_t, stage == "troph")
trophs_50 <- trophs %>% group_by(stage) %>% slice_sample(n = 50)
trophs_125 <- trophs %>% group_by(stage) %>% slice_sample(n = 125)
trophs_250 <- trophs %>% group_by(stage) %>% slice_sample(n = 250)
trophs_375 <- trophs %>% group_by(stage) %>% slice_sample(n = 375)
trophs_450 <- trophs %>% group_by(stage) %>% slice_sample(n = 450)
trophs_500 <- trophs %>% group_by(stage) %>% slice_sample(n = 500)

rings_100percent <- rings_500
rings_100percent <- as.data.frame(t(rings_100percent))
colnames(rings_100percent) <- rings_100percent[1,]
rings_100percent <- rings_100percent[-1,]
rings_100percent[] <- sapply(rings_100percent, as.numeric)
rings_100percent$mix <- rowSums(rings_100percent)
mixA <- as.data.frame(rings_100percent$mix)
rownames(mixA) <- rownames(rings_100percent)

trophs10_rings90 <- rbind(trophs_50, rings_450)
trophs10_rings90 <- as.data.frame(t(trophs10_rings90))
trophs10_rings90 <- trophs10_rings90[-1,]
trophs10_rings90[] <- sapply(trophs10_rings90, as.numeric)
trophs10_rings90$mix <- rowSums(trophs10_rings90)
mixB <- as.data.frame(trophs10_rings90$mix)
rownames(mixB) <- rownames(trophs10_rings90)

trophs25_rings75 <- rbind(trophs_125, rings_375)
trophs25_rings75 <- as.data.frame(t(trophs25_rings75))
trophs25_rings75 <- trophs25_rings75[-1,]
trophs25_rings75[] <- sapply(trophs25_rings75, as.numeric)
trophs25_rings75$mix <- rowSums(trophs25_rings75)
mixC <- as.data.frame(trophs25_rings75$mix)
rownames(mixC) <- rownames(trophs25_rings75)

trophs50_rings50 <- rbind(trophs_250, rings_250)
trophs50_rings50 <- as.data.frame(t(trophs50_rings50))
trophs50_rings50 <- trophs50_rings50[-1,]
trophs50_rings50[] <- sapply(trophs50_rings50, as.numeric)
trophs50_rings50$mix <- rowSums(trophs50_rings50)
mixD <- as.data.frame(trophs50_rings50$mix)
rownames(mixD) <- rownames(trophs50_rings50)

trophs75_rings25 <- rbind(trophs_375, rings_125)
trophs75_rings25 <- as.data.frame(t(trophs75_rings25))
trophs75_rings25 <- trophs75_rings25[-1,]
trophs75_rings25[] <- sapply(trophs75_rings25, as.numeric)
trophs75_rings25$mix <- rowSums(trophs75_rings25)
mixE <- as.data.frame(trophs75_rings25$mix)
rownames(mixE) <- rownames(trophs75_rings25)

trophs90_rings10 <- rbind(trophs_450, rings_50)
trophs90_rings10 <- as.data.frame(t(trophs90_rings10))
trophs90_rings10 <- trophs90_rings10[-1,]
trophs90_rings10[] <- sapply(trophs90_rings10, as.numeric)
trophs90_rings10$mix <- rowSums(trophs90_rings10)
mixF <- as.data.frame(trophs90_rings10$mix)
rownames(mixF) <- rownames(trophs90_rings10)

trophs100percent <- trophs_500
trophs100percent <- as.data.frame(t(trophs100percent))
trophs100percent <- trophs100percent[-1,]
trophs100percent[] <- sapply(trophs100percent, as.numeric)
trophs100percent$mix <- rowSums(trophs100percent)
mixG <- as.data.frame(trophs100percent$mix)
rownames(mixG) <- rownames(trophs100percent)

mixes_A <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_A) <- c("100_rings", "trophs10_rings90", "trophs25_rings75", "trophs50_rings50", "trophs75_rings25", "trophs90_rings10","trophs100")

write.table(mixes_A, file = "mixA_cpm.txt", sep = "\t", quote = FALSE)

#Generate mixtures with increasing % schizonts mixed with rings
rings <- subset(pf_t, stage == "ring")
rings_50 <- rings %>% group_by(stage) %>% slice_sample(n = 50)
rings_125 <- rings %>% group_by(stage) %>% slice_sample(n = 125)
rings_250 <- rings %>% group_by(stage) %>% slice_sample(n = 250)
rings_375 <- rings %>% group_by(stage) %>% slice_sample(n = 375)
rings_450 <- rings %>% group_by(stage) %>% slice_sample(n = 450)
rings_500 <- rings %>% group_by(stage) %>% slice_sample(n = 500)

schiz <- subset(pf_t, stage == "schizont")
schiz_50 <- schiz %>% group_by(stage) %>% slice_sample(n = 50)
schiz_125 <- schiz %>% group_by(stage) %>% slice_sample(n = 125)
schiz_250 <- schiz %>% group_by(stage) %>% slice_sample(n = 250)
schiz_375 <- schiz %>% group_by(stage) %>% slice_sample(n = 375)
schiz_450 <- schiz %>% group_by(stage) %>% slice_sample(n = 450)
schiz_500 <- schiz %>% group_by(stage) %>% slice_sample(n = 500)

rings_100percent <- rings_500
rings_100percent <- as.data.frame(t(rings_100percent))
rings_100percent <- rings_100percent[-1,]
rings_100percent[] <- sapply(rings_100percent, as.numeric)
rings_100percent$mix <- rowSums(rings_100percent)
mixA <- as.data.frame(rings_100percent$mix)
rownames(mixA) <- rownames(rings_100percent)

schiz10_rings90 <- rbind(schiz_50, rings_450)
schiz10_rings90 <- as.data.frame(t(schiz10_rings90))
schiz10_rings90 <- schiz10_rings90[-1,]
schiz10_rings90[] <- sapply(schiz10_rings90, as.numeric)
schiz10_rings90$mix <- rowSums(schiz10_rings90)
mixB <- as.data.frame(schiz10_rings90$mix)
rownames(mixB) <- rownames(schiz10_rings90)

schiz25_rings75 <- rbind(schiz_125, rings_375)
schiz25_rings75 <- as.data.frame(t(schiz25_rings75))
schiz25_rings75 <- schiz25_rings75[-1,]
schiz25_rings75[] <- sapply(schiz25_rings75, as.numeric)
schiz25_rings75$mix <- rowSums(schiz25_rings75)
mixC <- as.data.frame(schiz25_rings75$mix)
rownames(mixC) <- rownames(schiz25_rings75)

schiz50_rings50 <- rbind(schiz_250, rings_250)
schiz50_rings50 <- as.data.frame(t(schiz50_rings50))
schiz50_rings50 <- schiz50_rings50[-1,]
schiz50_rings50[] <- sapply(schiz50_rings50, as.numeric)
schiz50_rings50$mix <- rowSums(schiz50_rings50)
mixD <- as.data.frame(schiz50_rings50$mix)
rownames(mixD) <- rownames(schiz50_rings50)

schiz75_rings25 <- rbind(schiz_375, rings_125)
schiz75_rings25 <- as.data.frame(t(schiz75_rings25))
schiz75_rings25 <- schiz75_rings25[-1,]
schiz75_rings25[] <- sapply(schiz75_rings25, as.numeric)
schiz75_rings25$mix <- rowSums(schiz75_rings25)
mixE <- as.data.frame(schiz75_rings25$mix)
rownames(mixE) <- rownames(schiz75_rings25)

schiz90_rings10 <- rbind(schiz_450, rings_50)
schiz90_rings10 <- as.data.frame(t(schiz90_rings10))
schiz90_rings10 <- schiz90_rings10[-1,]
schiz90_rings10[] <- sapply(schiz90_rings10, as.numeric)
schiz90_rings10$mix <- rowSums(schiz90_rings10)
mixF <- as.data.frame(schiz90_rings10$mix)
rownames(mixF) <- rownames(schiz90_rings10)

schiz100percent <- schiz_500
schiz100percent <- as.data.frame(t(schiz100percent))
schiz100percent <- schiz100percent[-1,]
schiz100percent[] <- sapply(schiz100percent, as.numeric)
schiz100percent$mix <- rowSums(schiz100percent)
mixG <- as.data.frame(schiz100percent$mix)
rownames(mixG) <- rownames(schiz100percent)

mixes_B <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_B) <- c("100_rings", "schiz10_rings90", "schiz25_rings75", "schiz50_rings50", "schiz75_rings25", "schiz90_rings10","schiz100")

write.table(mixes_B, file = "mixB_cpm.txt", sep = "\t", quote = FALSE)

##Generate mixtures with increasing % rings mixed with trophs
rings <- subset(pf_t, stage == "ring")
rings_50 <- rings %>% group_by(stage) %>% slice_sample(n = 50)
rings_125 <- rings %>% group_by(stage) %>% slice_sample(n = 125)
rings_250 <- rings %>% group_by(stage) %>% slice_sample(n = 250)
rings_375 <- rings %>% group_by(stage) %>% slice_sample(n = 375)
rings_450 <- rings %>% group_by(stage) %>% slice_sample(n = 450)
rings_500 <- rings %>% group_by(stage) %>% slice_sample(n = 500)

trophs <- subset(pf_t, stage == "troph")
trophs_50 <- trophs %>% group_by(stage) %>% slice_sample(n = 50)
trophs_125 <- trophs %>% group_by(stage) %>% slice_sample(n = 125)
trophs_250 <- trophs %>% group_by(stage) %>% slice_sample(n = 250)
trophs_375 <- trophs %>% group_by(stage) %>% slice_sample(n = 375)
trophs_450 <- trophs %>% group_by(stage) %>% slice_sample(n = 450)
trophs_500 <- trophs %>% group_by(stage) %>% slice_sample(n = 500)

trophs100percent <- trophs_500
trophs100percent <- as.data.frame(t(trophs100percent))
trophs100percent <- trophs100percent[-1,]
trophs100percent[] <- sapply(trophs100percent, as.numeric)
trophs100percent$mix <- rowSums(trophs100percent)
mixA <- as.data.frame(trophs100percent$mix)
rownames(mixA) <- rownames(trophs100percent)

rings10_trophs90 <- rbind(rings_50, trophs_450)
rings10_trophs90 <- as.data.frame(t(rings10_trophs90))
rings10_trophs90 <- rings10_trophs90[-1,]
rings10_trophs90[] <- sapply(rings10_trophs90, as.numeric)
rings10_trophs90$mix <- rowSums(rings10_trophs90)
mixB <- as.data.frame(rings10_trophs90$mix)
rownames(mixB) <- rownames(rings10_trophs90)

rings25_trophs75 <- rbind(rings_125, trophs_375)
rings25_trophs75 <- as.data.frame(t(rings25_trophs75))
rings25_trophs75 <- rings25_trophs75[-1,]
rings25_trophs75[] <- sapply(rings25_trophs75, as.numeric)
rings25_trophs75$mix <- rowSums(rings25_trophs75)
mixC <- as.data.frame(rings25_trophs75$mix)
rownames(mixC) <- rownames(rings25_trophs75)

rings50_trophs50 <- rbind(rings_250, trophs_250)
rings50_trophs50 <- as.data.frame(t(rings50_trophs50))
rings50_trophs50 <- rings50_trophs50[-1,]
rings50_trophs50[] <- sapply(rings50_trophs50, as.numeric)
rings50_trophs50$mix <- rowSums(rings50_trophs50)
mixD <- as.data.frame(rings50_trophs50$mix)
rownames(mixD) <- rownames(rings50_trophs50)

rings75_trophs25 <- rbind(rings_375, trophs_125)
rings75_trophs25 <- as.data.frame(t(rings75_trophs25))
rings75_trophs25 <- rings75_trophs25[-1,]
rings75_trophs25[] <- sapply(rings75_trophs25, as.numeric)
rings75_trophs25$mix <- rowSums(rings75_trophs25)
mixE <- as.data.frame(rings75_trophs25$mix)
rownames(mixE) <- rownames(rings75_trophs25)

rings90_trophs10 <- rbind(rings_450, trophs_50)
rings90_trophs10 <- as.data.frame(t(rings90_trophs10))
rings90_trophs10 <- rings90_trophs10[-1,]
rings90_trophs10[] <- sapply(rings90_trophs10, as.numeric)
rings90_trophs10$mix <- rowSums(rings90_trophs10)
mixF <- as.data.frame(rings90_trophs10$mix)
rownames(mixF) <- rownames(rings90_trophs10)

rings100percent <- rings_500
rings100percent <- as.data.frame(t(rings100percent))
rings100percent <- rings100percent[-1,]
rings100percent[] <- sapply(rings100percent, as.numeric)
rings100percent$mix <- rowSums(rings100percent)
mixG <- as.data.frame(rings100percent$mix)
rownames(mixG) <- rownames(rings100percent)

mixes_C <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_C) <- c("100_trophs", "rings10_trophs90", "rings25_trophs75", "rings50_trophs50", "rings75_trophs25", "rings90_trophs10","rings100")

write.table(mixes_C, file = "mixC_cpm.txt", sep = "\t", quote = FALSE)

##Generate mixtures with increasing % schizonts mixed with trophs 
schiz <- subset(pf_t, stage == "schizont")
schiz_50 <- schiz %>% group_by(stage) %>% slice_sample(n = 50)
schiz_125 <- schiz %>% group_by(stage) %>% slice_sample(n = 125)
schiz_250 <- schiz %>% group_by(stage) %>% slice_sample(n = 250)
schiz_375 <- schiz %>% group_by(stage) %>% slice_sample(n = 375)
schiz_450 <- schiz %>% group_by(stage) %>% slice_sample(n = 450)
schiz_500 <- schiz %>% group_by(stage) %>% slice_sample(n = 500)

trophs <- subset(pf_t, stage == "troph")
trophs_50 <- trophs %>% group_by(stage) %>% slice_sample(n = 50)
trophs_125 <- trophs %>% group_by(stage) %>% slice_sample(n = 125)
trophs_250 <- trophs %>% group_by(stage) %>% slice_sample(n = 250)
trophs_375 <- trophs %>% group_by(stage) %>% slice_sample(n = 375)
trophs_450 <- trophs %>% group_by(stage) %>% slice_sample(n = 450)
trophs_500 <- trophs %>% group_by(stage) %>% slice_sample(n = 500)

trophs100percent <- trophs_500
trophs100percent <- as.data.frame(t(trophs100percent))
trophs100percent <- trophs100percent[-1,]
trophs100percent[] <- sapply(trophs100percent, as.numeric)
trophs100percent$mix <- rowSums(trophs100percent)
mixA <- as.data.frame(trophs100percent$mix)
rownames(mixA) <- rownames(trophs100percent)

schiz10_trophs90 <- rbind(schiz_50, trophs_450)
schiz10_trophs90 <- as.data.frame(t(schiz10_trophs90))
schiz10_trophs90 <- schiz10_trophs90[-1,]
schiz10_trophs90[] <- sapply(schiz10_trophs90, as.numeric)
schiz10_trophs90$mix <- rowSums(schiz10_trophs90)
mixB <- as.data.frame(schiz10_trophs90$mix)
rownames(mixB) <- rownames(schiz10_trophs90)

schiz25_trophs75 <- rbind(schiz_125, trophs_375)
schiz25_trophs75 <- as.data.frame(t(schiz25_trophs75))
schiz25_trophs75 <- schiz25_trophs75[-1,]
schiz25_trophs75[] <- sapply(schiz25_trophs75, as.numeric)
schiz25_trophs75$mix <- rowSums(schiz25_trophs75)
mixC <- as.data.frame(schiz25_trophs75$mix)
rownames(mixC) <- rownames(schiz25_trophs75)

schiz50_trophs50 <- rbind(schiz_250, trophs_250)
schiz50_trophs50 <- as.data.frame(t(schiz50_trophs50))
schiz50_trophs50 <- schiz50_trophs50[-1,]
schiz50_trophs50[] <- sapply(schiz50_trophs50, as.numeric)
schiz50_trophs50$mix <- rowSums(schiz50_trophs50)
mixD <- as.data.frame(schiz50_trophs50$mix)
rownames(mixD) <- rownames(schiz50_trophs50)

schiz75_trophs25 <- rbind(schiz_375, trophs_125)
schiz75_trophs25 <- as.data.frame(t(schiz75_trophs25))
schiz75_trophs25 <- schiz75_trophs25[-1,]
schiz75_trophs25[] <- sapply(schiz75_trophs25, as.numeric)
schiz75_trophs25$mix <- rowSums(schiz75_trophs25)
mixE <- as.data.frame(schiz75_trophs25$mix)
rownames(mixE) <- rownames(schiz75_trophs25)

schiz90_trophs10 <- rbind(schiz_450, trophs_50)
schiz90_trophs10 <- as.data.frame(t(schiz90_trophs10))
schiz90_trophs10 <- schiz90_trophs10[-1,]
schiz90_trophs10[] <- sapply(schiz90_trophs10, as.numeric)
schiz90_trophs10$mix <- rowSums(schiz90_trophs10)
mixF <- as.data.frame(schiz90_trophs10$mix)
rownames(mixF) <- rownames(schiz90_trophs10)

schiz100percent <- schiz_500
schiz100percent <- as.data.frame(t(schiz100percent))
schiz100percent <- schiz100percent[-1,]
schiz100percent[] <- sapply(schiz100percent, as.numeric)
schiz100percent$mix <- rowSums(schiz100percent)
mixG <- as.data.frame(schiz100percent$mix)
rownames(mixG) <- rownames(schiz100percent)

mixes_D <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_D) <- c("100_trophs", "schiz10_trophs90", "schiz25_trophs75", "schiz50_trophs50", "schiz75_trophs25", "schiz90_trophs10","schiz100")

write.table(mixes_D, file = "mixD_cpm.txt", sep = "\t", quote = FALSE)

##Generate mixtures with increasing % rings mixed with schizonts 
schiz <- subset(pf_t, stage == "schizont")
schiz_50 <- schiz %>% group_by(stage) %>% slice_sample(n = 50)
schiz_125 <- schiz %>% group_by(stage) %>% slice_sample(n = 125)
schiz_250 <- schiz %>% group_by(stage) %>% slice_sample(n = 250)
schiz_375 <- schiz %>% group_by(stage) %>% slice_sample(n = 375)
schiz_450 <- schiz %>% group_by(stage) %>% slice_sample(n = 450)
schiz_500 <- schiz %>% group_by(stage) %>% slice_sample(n = 500)

rings <- subset(pf_t, stage == "ring")
rings_50 <- rings %>% group_by(stage) %>% slice_sample(n = 50)
rings_125 <- rings %>% group_by(stage) %>% slice_sample(n = 125)
rings_250 <- rings %>% group_by(stage) %>% slice_sample(n = 250)
rings_375 <- rings %>% group_by(stage) %>% slice_sample(n = 375)
rings_450 <- rings %>% group_by(stage) %>% slice_sample(n = 450)
rings_500 <- rings %>% group_by(stage) %>% slice_sample(n = 500)

schiz100percent <- schiz_500
schiz100percent <- as.data.frame(t(schiz100percent))
schiz100percent <- schiz100percent[-1,]
schiz100percent[] <- sapply(schiz100percent, as.numeric)
schiz100percent$mix <- rowSums(schiz100percent)
mixA <- as.data.frame(schiz100percent$mix)
#schiz100percent$avg <- rowMeans(schiz100percent)
rownames(mixA) <- rownames(schiz100percent)

rings10_schiz90 <- rbind(rings_50, schiz_450)
rings10_schiz90 <- as.data.frame(t(rings10_schiz90))
rings10_schiz90 <- rings10_schiz90[-1,]
rings10_schiz90[] <- sapply(rings10_schiz90, as.numeric)
rings10_schiz90$mix <- rowSums(rings10_schiz90)
mixB <- as.data.frame(rings10_schiz90$mix)
rownames(mixB) <- rownames(rings10_schiz90)

rings25_schiz75 <- rbind(rings_125, schiz_375)
rings25_schiz75 <- as.data.frame(t(rings25_schiz75))
rings25_schiz75 <- rings25_schiz75[-1,]
rings25_schiz75[] <- sapply(rings25_schiz75, as.numeric)
rings25_schiz75$mix <- rowSums(rings25_schiz75)
mixC <- as.data.frame(rings25_schiz75$mix)
rownames(mixC) <- rownames(rings25_schiz75)

rings50_schiz50 <- rbind(rings_250, schiz_250)
schiz50_rings50 <- as.data.frame(t(schiz50_rings50))
schiz50_rings50 <- schiz50_rings50[-1,]
schiz50_rings50[] <- sapply(schiz50_rings50, as.numeric)
schiz50_rings50$mix <- rowSums(schiz50_rings50)
mixD <- as.data.frame(schiz50_rings50$mix)
rownames(mixD) <- rownames(schiz50_rings50)

rings75_schiz25 <- rbind(rings_375, schiz_125)
rings75_schiz25 <- as.data.frame(t(rings75_schiz25))
rings75_schiz25 <- rings75_schiz25[-1,]
rings75_schiz25[] <- sapply(rings75_schiz25, as.numeric)
rings75_schiz25$mix <- rowSums(rings75_schiz25)
mixE <- as.data.frame(rings75_schiz25$mix)
rownames(mixE) <- rownames(rings75_schiz25)

rings90_schiz10 <- rbind(rings_450, schiz_50)
rings90_schiz10 <- as.data.frame(t(rings90_schiz10))
rings90_schiz10 <- rings90_schiz10[-1,]
rings90_schiz10[] <- sapply(rings90_schiz10, as.numeric)
rings90_schiz10$mix <- rowSums(rings90_schiz10)
mixF <- as.data.frame(rings90_schiz10$mix)
rownames(mixF) <- rownames(rings90_schiz10)

rings100percent <- rings_500
rings100percent <- as.data.frame(t(rings100percent))
rings100percent <- rings100percent[-1,]
rings100percent[] <- sapply(rings100percent, as.numeric)
rings100percent$mix <- rowSums(rings100percent)
mixG <- as.data.frame(rings100percent$mix)
rownames(mixG) <- rownames(rings100percent)

mixes_E <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_E) <- c("100schiz", "rings10_schiz90", "rings25_schiz75", "rings50_schiz50", "rings75_schiz25", "rings90_schiz10","rings100")

write.table(mixes_E, file = "mixE_cpm.txt", sep = "\t", quote = FALSE)

##Generate mixtures with increasing % trophs mixed with schizonts - F
trophs <- subset(pf_t, stage == "troph")
trophs_50 <- trophs %>% group_by(stage) %>% slice_sample(n = 50)
trophs_125 <- trophs %>% group_by(stage) %>% slice_sample(n = 125)
trophs_250 <- trophs %>% group_by(stage) %>% slice_sample(n = 250)
trophs_375 <- trophs %>% group_by(stage) %>% slice_sample(n = 375)
trophs_450 <- trophs %>% group_by(stage) %>% slice_sample(n = 450)
trophs_500 <- trophs %>% group_by(stage) %>% slice_sample(n = 500)

schiz <- subset(pf_t, stage == "schizont")
schiz_50 <- schiz %>% group_by(stage) %>% slice_sample(n = 50)
schiz_125 <- schiz %>% group_by(stage) %>% slice_sample(n = 125)
schiz_250 <- schiz %>% group_by(stage) %>% slice_sample(n = 250)
schiz_375 <- schiz %>% group_by(stage) %>% slice_sample(n = 375)
schiz_450 <- schiz %>% group_by(stage) %>% slice_sample(n = 450)
schiz_500 <- schiz %>% group_by(stage) %>% slice_sample(n = 500)

schiz100percent <- schiz_500
schiz100percent <- as.data.frame(t(schiz100percent))
schiz100percent <- schiz100percent[-1,]
schiz100percent[] <- sapply(schiz100percent, as.numeric)
schiz100percent$mix <- rowSums(schiz100percent)
mixA <- as.data.frame(schiz100percent$mix)
rownames(mixA) <- rownames(schiz100percent)

trophs10_schiz90 <- rbind(trophs_50, schiz_450)
trophs10_schiz90 <- as.data.frame(t(trophs10_schiz90))
trophs10_schiz90 <- trophs10_schiz90[-1,]
trophs10_schiz90[] <- sapply(trophs10_schiz90, as.numeric)
trophs10_schiz90$mix <- rowSums(trophs10_schiz90)
mixB <- as.data.frame(trophs10_schiz90$mix)
rownames(mixB) <- rownames(trophs10_schiz90)

trophs25_schiz75 <- rbind(trophs_125, schiz_375)
trophs25_schiz75 <- as.data.frame(t(trophs25_schiz75))
trophs25_schiz75 <- trophs25_schiz75[-1,]
trophs25_schiz75[] <- sapply(trophs25_schiz75, as.numeric)
trophs25_schiz75$mix <- rowSums(trophs25_schiz75)
mixC <- as.data.frame(trophs25_schiz75$mix)
rownames(mixC) <- rownames(trophs25_schiz75)

trophs50_schiz50 <- rbind(trophs_250, schiz_250)
trophs50_schiz50 <- as.data.frame(t(trophs50_schiz50))
trophs50_schiz50 <- trophs50_schiz50[-1,]
trophs50_schiz50[] <- sapply(trophs50_schiz50, as.numeric)
trophs50_schiz50$mix <- rowSums(trophs50_schiz50)
mixD <- as.data.frame(trophs50_schiz50$mix)
rownames(mixD) <- rownames(trophs50_schiz50)

trophs75_schiz25 <- rbind(trophs_375, schiz_125)
trophs75_schiz25 <- as.data.frame(t(trophs75_schiz25))
trophs75_schiz25 <- trophs75_schiz25[-1,]
trophs75_schiz25[] <- sapply(trophs75_schiz25, as.numeric)
trophs75_schiz25$mix <- rowSums(trophs75_schiz25)
mixE <- as.data.frame(trophs75_schiz25$mix)
rownames(mixE) <- rownames(trophs75_schiz25)

trophs90_schiz10 <- rbind(trophs_450, schiz_50)
trophs90_schiz10 <- as.data.frame(t(trophs90_schiz10))
trophs90_schiz10 <- trophs90_schiz10[-1,]
trophs90_schiz10[] <- sapply(trophs90_schiz10, as.numeric)
trophs90_schiz10$mix <- rowSums(trophs90_schiz10)
mixF <- as.data.frame(trophs90_schiz10$mix)
rownames(mixF) <- rownames(trophs90_schiz10)

trophs100percent <- trophs_500
trophs100percent <- as.data.frame(t(trophs100percent))
trophs100percent <- trophs100percent[-1,]
trophs100percent[] <- sapply(trophs100percent, as.numeric)
trophs100percent$mix <- rowSums(trophs100percent)
mixG <- as.data.frame(trophs100percent$mix)
rownames(mixG) <- rownames(trophs100percent)

mixes_F <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_F) <- c("100schiz", "trophs10_schiz90", "trophs25_schiz75", "trophs50_schiz50", "trophs75_schiz25", "trophs90_schiz10","trophs100")

write.table(mixes_F, file = "mixF_cpm.txt", sep = "\t", quote = FALSE)

######Using pure populations of P berghei
#Run each section 3 times and merge into 1 file so 3 replicates of random mixtures per panel 

pb <- read.csv("pb10xIDC_counts_purepops.csv", header = 1, check.names = FALSE)

#Increasing % trophs mixed with rings
rings <- subset(pb, stage == "ring")
rings_10 <- rings %>% group_by(stage) %>% slice_sample(n = 10)
rings_25 <- rings %>% group_by(stage) %>% slice_sample(n = 25)
rings_50 <- rings %>% group_by(stage) %>% slice_sample(n = 50)
rings_75 <- rings %>% group_by(stage) %>% slice_sample(n = 75)
rings_90 <- rings %>% group_by(stage) %>% slice_sample(n = 90)
rings_100 <- rings %>% group_by(stage) %>% slice_sample(n = 100)

trophs <- subset(pb, stage == "troph")
trophs_10 <- trophs %>% group_by(stage) %>% slice_sample(n = 10)
trophs_25 <- trophs %>% group_by(stage) %>% slice_sample(n = 25)
trophs_50 <- trophs %>% group_by(stage) %>% slice_sample(n = 50)
trophs_75 <- trophs %>% group_by(stage) %>% slice_sample(n = 75)
trophs_90 <- trophs %>% group_by(stage) %>% slice_sample(n = 90)
trophs_100 <- trophs %>% group_by(stage) %>% slice_sample(n = 100)

rings_100percent <- rings_100
rings_100percent <- as.data.frame(t(rings_100percent))
colnames(rings_100percent) <- rings_100percent[1,]
rings_100percent <- rings_100percent[-1,]
rings_100percent[] <- sapply(rings_100percent, as.numeric)
rings_100percent$mix <- rowSums(rings_100percent)
mixA <- as.data.frame(rings_100percent$mix)
rownames(mixA) <- rownames(rings_100percent)

trophs10_rings90 <- rbind(trophs_10, rings_90)
trophs10_rings90 <- as.data.frame(t(trophs10_rings90))
trophs10_rings90 <- trophs10_rings90[-1,]
trophs10_rings90[] <- sapply(trophs10_rings90, as.numeric)
trophs10_rings90$mix <- rowSums(trophs10_rings90)
mixB <- as.data.frame(trophs10_rings90$mix)
rownames(mixB) <- rownames(trophs10_rings90)

trophs25_rings75 <- rbind(trophs_25, rings_75)
trophs25_rings75 <- as.data.frame(t(trophs25_rings75))
trophs25_rings75 <- trophs25_rings75[-1,]
trophs25_rings75[] <- sapply(trophs25_rings75, as.numeric)
trophs25_rings75$mix <- rowSums(trophs25_rings75)
mixC <- as.data.frame(trophs25_rings75$mix)
rownames(mixC) <- rownames(trophs25_rings75)

trophs50_rings50 <- rbind(trophs_50, rings_50)
trophs50_rings50 <- as.data.frame(t(trophs50_rings50))
trophs50_rings50 <- trophs50_rings50[-1,]
trophs50_rings50[] <- sapply(trophs50_rings50, as.numeric)
trophs50_rings50$mix <- rowSums(trophs50_rings50)
mixD <- as.data.frame(trophs50_rings50$mix)
rownames(mixD) <- rownames(trophs50_rings50)

trophs75_rings25 <- rbind(trophs_75, rings_25)
trophs75_rings25 <- as.data.frame(t(trophs75_rings25))
trophs75_rings25 <- trophs75_rings25[-1,]
trophs75_rings25[] <- sapply(trophs75_rings25, as.numeric)
trophs75_rings25$mix <- rowSums(trophs75_rings25)
mixE <- as.data.frame(trophs75_rings25$mix)
rownames(mixE) <- rownames(trophs75_rings25)

trophs90_rings10 <- rbind(trophs_90, rings_10)
trophs90_rings10 <- as.data.frame(t(trophs90_rings10))
trophs90_rings10 <- trophs90_rings10[-1,]
trophs90_rings10[] <- sapply(trophs90_rings10, as.numeric)
trophs90_rings10$mix <- rowSums(trophs90_rings10)
mixF <- as.data.frame(trophs90_rings10$mix)
rownames(mixF) <- rownames(trophs90_rings10)

trophs100percent <- trophs_100
trophs100percent <- as.data.frame(t(trophs100percent))
trophs100percent <- trophs100percent[-1,]
trophs100percent[] <- sapply(trophs100percent, as.numeric)
trophs100percent$mix <- rowSums(trophs100percent)
mixG <- as.data.frame(trophs100percent$mix)
rownames(mixG) <- rownames(trophs100percent)

mixes_A_1 <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_A_1) <- c("100_rings", "trophs10_rings90", "trophs25_rings75", "trophs50_rings50", "trophs75_rings25", "trophs90_rings10","trophs100")

mixes_A_2 <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_A_2) <- c("100_rings_2", "trophs10_rings90_2", "trophs25_rings75_2", "trophs50_rings50_2", "trophs75_rings25_2", "trophs90_rings10_2","trophs100_2")

mixes_A_3 <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_A_3) <- c("100_rings_3", "trophs10_rings90_3", "trophs25_rings75_3", "trophs50_rings50_3", "trophs75_rings25_3", "trophs90_rings10_3","trophs100_3")

mixes_A <- cbind(mixes_A_1, mixes_A_2, mixes_A_3)
mixes_A$gene <- rownames(mixes_A)
mixes_A <- merge(mixes_A, orthologs, by = "gene")
mixes_A <- mixes_A[,-c(23:24)]
write.table(mixes_A, file = "mix_inc_troph_w_rings.txt", sep = "\t", quote = FALSE)

#Increasing % schizonts mixed with rings
rings <- subset(pb, stage == "ring")
rings_10 <- rings %>% group_by(stage) %>% slice_sample(n = 10)
rings_25 <- rings %>% group_by(stage) %>% slice_sample(n = 25)
rings_50 <- rings %>% group_by(stage) %>% slice_sample(n = 50)
rings_75 <- rings %>% group_by(stage) %>% slice_sample(n = 75)
rings_90 <- rings %>% group_by(stage) %>% slice_sample(n = 90)
rings_100 <- rings %>% group_by(stage) %>% slice_sample(n = 100)

schiz <- subset(pb, stage == "schizont")
schiz_10 <- schiz %>% group_by(stage) %>% slice_sample(n = 10)
schiz_25 <- schiz %>% group_by(stage) %>% slice_sample(n = 25)
schiz_50 <- schiz %>% group_by(stage) %>% slice_sample(n = 50)
schiz_75 <- schiz %>% group_by(stage) %>% slice_sample(n = 75)
schiz_90 <- schiz %>% group_by(stage) %>% slice_sample(n = 90)
schiz_100 <- schiz %>% group_by(stage) %>% slice_sample(n = 100)

rings_100percent <- rings_100
rings_100percent <- as.data.frame(t(rings_100percent))
colnames(rings_100percent) <- rings_100percent[1,]
rings_100percent <- rings_100percent[-1,]
rings_100percent[] <- sapply(rings_100percent, as.numeric)
rings_100percent$mix <- rowSums(rings_100percent)
mixA <- as.data.frame(rings_100percent$mix)
rownames(mixA) <- rownames(rings_100percent)

schiz10_rings90 <- rbind(schiz_10, rings_90)
schiz10_rings90 <- as.data.frame(t(schiz10_rings90))
colnames(schiz10_rings90) <- schiz10_rings90[1,]
schiz10_rings90 <- schiz10_rings90[-1,]
schiz10_rings90[] <- sapply(schiz10_rings90, as.numeric)
schiz10_rings90$mix <- rowSums(schiz10_rings90)
mixB <- as.data.frame(schiz10_rings90$mix)
rownames(mixB) <- rownames(schiz10_rings90)

schiz25_rings75 <- rbind(schiz_25, rings_75)
schiz25_rings75 <- as.data.frame(t(schiz25_rings75))
colnames(schiz25_rings75) <- schiz25_rings75[1,]
schiz25_rings75 <- schiz25_rings75[-1,]
schiz25_rings75[] <- sapply(schiz25_rings75, as.numeric)
schiz25_rings75$mix <- rowSums(schiz25_rings75)
mixC <- as.data.frame(schiz25_rings75$mix)
rownames(mixC) <- rownames(schiz25_rings75)

schiz50_rings50 <- rbind(schiz_50, rings_50)
schiz50_rings50 <- as.data.frame(t(schiz50_rings50))
colnames(schiz50_rings50) <- schiz50_rings50[1,]
schiz50_rings50 <- schiz50_rings50[-1,]
schiz50_rings50[] <- sapply(schiz50_rings50, as.numeric)
schiz50_rings50$mix <- rowSums(schiz50_rings50)
mixD <- as.data.frame(schiz50_rings50$mix)
rownames(mixD) <- rownames(schiz50_rings50)

schiz75_rings25 <- rbind(schiz_75, rings_25)
schiz75_rings25 <- as.data.frame(t(schiz75_rings25))
colnames(schiz75_rings25) <- schiz75_rings25[1,]
schiz75_rings25 <- schiz75_rings25[-1,]
schiz75_rings25[] <- sapply(schiz75_rings25, as.numeric)
schiz75_rings25$mix <- rowSums(schiz75_rings25)
mixE <- as.data.frame(schiz75_rings25$mix)
rownames(mixE) <- rownames(schiz75_rings25)

schiz90_rings10 <- rbind(schiz_90, rings_10)
schiz90_rings10 <- as.data.frame(t(schiz90_rings10))
colnames(schiz90_rings10) <- schiz90_rings10[1,]
schiz90_rings10 <- schiz90_rings10[-1,]
schiz90_rings10[] <- sapply(schiz90_rings10, as.numeric)
schiz90_rings10$mix <- rowSums(schiz90_rings10)
mixF <- as.data.frame(schiz90_rings10$mix)
rownames(mixF) <- rownames(schiz90_rings10)

schiz100percent <- schiz_100
schiz100percent <- as.data.frame(t(schiz100percent))
colnames(schiz100percent) <- schiz100percent[1,]
schiz100percent <- schiz100percent[-1,]
schiz100percent[] <- sapply(schiz100percent, as.numeric)
schiz100percent$mix <- rowSums(schiz100percent)
mixG <- as.data.frame(schiz100percent$mix)
rownames(mixG) <- rownames(schiz100percent)

mixes_B_1 <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_B_1) <- c("100_rings", "schiz10_rings90", "schiz25_rings75", "schiz50_rings50", "schiz75_rings25", "schiz90_rings10","schiz100")

mixes_B_2 <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_B_2) <- c("100_rings_2", "schiz10_rings90_2", "schiz25_rings75_2", "schiz50_rings50_2", "schiz75_rings25_2", "schiz90_rings10_2","schiz100_2")

mixes_B_3 <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_B_3) <- c("100_rings_3", "schiz10_rings90_3", "schiz25_rings75_3", "schiz50_rings50_3", "schiz75_rings25_3", "schiz90_rings10_3","schiz100_3")

mixes_B <- cbind(mixes_B_1, mixes_B_2, mixes_B_3)
mixes_B$gene <- rownames(mixes_B)
mixes_B <- merge(mixes_B, orthologs, by = "gene")
mixes_B <- mixes_B[,-c(23:24)]
write.table(mixes_B, file = "mixB_inc_schiz_w_rings.txt", sep = "\t", quote = FALSE)

#Increasing % schizonts mixed with trophs
schiz <- subset(pb, stage == "schizont")
schiz_10 <- schiz %>% group_by(stage) %>% slice_sample(n = 10)
schiz_25 <- schiz %>% group_by(stage) %>% slice_sample(n = 25)
schiz_50 <- schiz %>% group_by(stage) %>% slice_sample(n = 50)
schiz_75 <- schiz %>% group_by(stage) %>% slice_sample(n = 75)
schiz_90 <- schiz %>% group_by(stage) %>% slice_sample(n = 90)
schiz_100 <- schiz %>% group_by(stage) %>% slice_sample(n = 100)

trophs <- subset(pb, stage == "troph")
trophs_10 <- trophs %>% group_by(stage) %>% slice_sample(n = 10)
trophs_25 <- trophs %>% group_by(stage) %>% slice_sample(n = 25)
trophs_50 <- trophs %>% group_by(stage) %>% slice_sample(n = 50)
trophs_75 <- trophs %>% group_by(stage) %>% slice_sample(n = 75)
trophs_90 <- trophs %>% group_by(stage) %>% slice_sample(n = 90)
trophs_100 <- trophs %>% group_by(stage) %>% slice_sample(n = 100)

trophs100percent <- trophs_100
trophs100percent <- as.data.frame(t(trophs100percent))
trophs100percent <- trophs100percent[-1,]
trophs100percent[] <- sapply(trophs100percent, as.numeric)
trophs100percent$mix <- rowSums(trophs100percent)
mixA <- as.data.frame(trophs100percent$mix)
rownames(mixA) <- rownames(trophs100percent)

schiz10_trophs90 <- rbind(schiz_10, trophs_90)
schiz10_trophs90 <- as.data.frame(t(schiz10_trophs90))
schiz10_trophs90 <- schiz10_trophs90[-1,]
schiz10_trophs90[] <- sapply(schiz10_trophs90, as.numeric)
schiz10_trophs90$mix <- rowSums(schiz10_trophs90)
mixB <- as.data.frame(schiz10_trophs90$mix)
rownames(mixB) <- rownames(schiz10_trophs90)

schiz25_trophs75 <- rbind(schiz_25, trophs_75)
schiz25_trophs75 <- as.data.frame(t(schiz25_trophs75))
schiz25_trophs75 <- schiz25_trophs75[-1,]
schiz25_trophs75[] <- sapply(schiz25_trophs75, as.numeric)
schiz25_trophs75$mix <- rowSums(schiz25_trophs75)
mixC <- as.data.frame(schiz25_trophs75$mix)
rownames(mixC) <- rownames(schiz25_trophs75)

schiz50_trophs50 <- rbind(schiz_50, trophs_50)
schiz50_trophs50 <- as.data.frame(t(schiz50_trophs50))
schiz50_trophs50 <- schiz50_trophs50[-1,]
schiz50_trophs50[] <- sapply(schiz50_trophs50, as.numeric)
schiz50_trophs50$mix <- rowSums(schiz50_trophs50)
mixD <- as.data.frame(schiz50_trophs50$mix)
rownames(mixD) <- rownames(schiz50_trophs50)

schiz75_trophs25 <- rbind(schiz_75, trophs_25)
schiz75_trophs25 <- as.data.frame(t(schiz75_trophs25))
schiz75_trophs25 <- schiz75_trophs25[-1,]
schiz75_trophs25[] <- sapply(schiz75_trophs25, as.numeric)
schiz75_trophs25$mix <- rowSums(schiz75_trophs25)
mixE <- as.data.frame(schiz75_trophs25$mix)
rownames(mixE) <- rownames(schiz75_trophs25)

schiz90_trophs10 <- rbind(schiz_90, trophs_10)
schiz90_trophs10 <- as.data.frame(t(schiz90_trophs10))
schiz90_trophs10 <- schiz90_trophs10[-1,]
schiz90_trophs10[] <- sapply(schiz90_trophs10, as.numeric)
schiz90_trophs10$mix <- rowSums(schiz90_trophs10)
mixF <- as.data.frame(schiz90_trophs10$mix)
rownames(mixF) <- rownames(schiz90_trophs10)

schiz100percent <- schiz_100
schiz100percent <- as.data.frame(t(schiz100percent))
schiz100percent <- schiz100percent[-1,]
schiz100percent[] <- sapply(schiz100percent, as.numeric)
schiz100percent$mix <- rowSums(schiz100percent)
mixG <- as.data.frame(schiz100percent$mix)
rownames(mixG) <- rownames(schiz100percent)

mixes_C_1 <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_C_1) <- c("100_trophs", "schiz10_trophs90", "schiz25_trophs75", "schiz50_trophs50", "schiz75_trophs25", "schiz90_trophs10","schiz100")

mixes_C_2 <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_C_2) <- c("100_trophs_2", "schiz10_trophs90_2", "schiz25_trophs75_2", "schiz50_trophs50_2", "schiz75_trophs25_2", "schiz90_trophs10_2","schiz100_2")

mixes_C_3 <- cbind(mixA, mixB, mixC, mixD, mixE, mixF, mixG)
colnames(mixes_C_3) <- c("100_trophs_3", "schiz10_trophs90_3", "schiz25_trophs75_3", "schiz50_trophs50_3", "schiz75_trophs25_3", "schiz90_trophs10_3","schiz100_3")

mixes_C <- cbind(mixes_C_1, mixes_C_2, mixes_C_3)
mixes_C$gene <- rownames(mixes_C)
mixes_C <- merge(mixes_C, orthologs, by = "gene")
mixes_C <- mixes_C[,-c(23:24)]
write.table(mixes_C, file = "mix_inc_schiz_w_troph.txt", sep = "\t", quote = FALSE)

