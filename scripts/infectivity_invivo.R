library(dplyr)
library(lme4)
library(emmeans)
library(ggplot2)

dat <- read.csv("data/staph_invivo_eop_281123.csv",header=T)

dat <- na.omit(dat)

dat_m <- group_by(dat, phage, bact_time, block) %>%
  summarise(m_den = mean(lden))

ggplot(dat, aes(x = bact_time, y = lden, group = interaction(bact_time, phage))) +
  theme_bw() +
  MicrobioUoE::geom_pretty_boxplot(aes(col = phage, fill = phage)) +
  geom_point(position = position_dodge(0.75), shape = 21, fill = 'white', size = 2) +
  ylab("Phage density log10(PFU/mL)") +
  xlab("Phage") +
  theme(axis.text = element_text(size = 14, colour = "black"), axis.title = element_text(size = 14), legend.position = 'bottom', legend.text = element_text(size = 12), legend.title = element_text(size = 12))+
  palettetown::scale_color_poke(pokemon = "Roselia", spread = 11, labels = c("Day 0", "Day 0 (pus)", "Day (BAL)")) +
  palettetown::scale_fill_poke(pokemon = "Roselia", spread = 11, labels = c("Day 0", "Day 0 (pus)", "Day (BAL)")) +
  scale_x_discrete(labels = c("ISP","Day 6 (pre-dose)","Day 6 (post-dose)")) +
  labs(col = "Bacteria", fill = "Bacteria")

#definitely need to do a block correction

b1 <- filter(dat_m, block == 1)
b2 <- filter(dat_m, block == 2)

b1$b2 <- b2$m_den

#calculate percentage change between blocks for each phage
b1 <- mutate(b1, 
             change = b2 - m_den,
             pchange = (change/b2),
             pchange2 = abs(pchange))

b2dat <- filter(dat,block == 2)

b2dat$n[b2dat$bact_time == 'b1' & b2dat$phage == 'isp'] <- 0.4337660
b2dat$n[b2dat$bact_time == 'b1' & b2dat$phage == '3b'] <- 0.2984625
b2dat$n[b2dat$bact_time == 'b1' & b2dat$phage == '4b'] <- 0.2454854
b2dat$n[b2dat$bact_time == 'b2' & b2dat$phage == 'isp'] <- 0.4349424
b2dat$n[b2dat$bact_time == 'b2' & b2dat$phage == '3b'] <- 0.2738900
b2dat$n[b2dat$bact_time == 'b2' & b2dat$phage == '4b'] <- 0.2327821
b2dat$n[b2dat$bact_time == 'b5' & b2dat$phage == 'isp'] <- 0.4176180
b2dat$n[b2dat$bact_time == 'b5' & b2dat$phage == '3b'] <- 0.3452511
b2dat$n[b2dat$bact_time == 'b5' & b2dat$phage == '4b'] <- 0.3064925

b2dat <- mutate(b2dat,
                diff = lden * n,
                lden2 = lden + diff)

#tidy up datasets 

b2dat <- b2dat[,c(1,2,8,9,12)]
b1dat <- filter(dat, block == 1)
b1dat <- b1dat[,c(1,2,7:9)]

names(b2dat)[5] <- 'lden'

bdat_n <- merge(b1dat,b2dat,all=T)

bdat_n$phage2 <- bdat_n$phage
bdat_n$phage2[bdat_n$phage2 == 'isp'] <- '1isp'

ggplot(bdat_n, aes(x = phage2, y = lden, group = interaction(phage, bact_time))) +
  theme_bw() +
  MicrobioUoE::geom_pretty_boxplot(aes(col = bact_time, fill = bact_time)) +
  geom_point(position = position_dodge(0.75), shape = 21, fill = 'white', size = 2) +
  ylab("Estimated phage density log10(PFU/mL)") +
  xlab("Phage") +
  theme(axis.text = element_text(size = 14, colour = "black"), axis.title = element_text(size = 14), legend.position = 'bottom', legend.text = element_text(size = 12), legend.title = element_text(size = 12))+
  palettetown::scale_color_poke(pokemon = "Roselia", spread = 11, labels = c("Day 0", "Day 0 (pus)", "Day 6 (BAL)")) +
  palettetown::scale_fill_poke(pokemon = "Roselia", spread = 11, labels = c("Day 0", "Day 0 (pus)", "Day 6 (BAL)")) +
  scale_x_discrete(labels = c("ISP","Day 6 (pre-dose)","Day 6 (post-dose)")) +
  labs(col = "Bacteria", fill = "Bacteria")

ggplot(bdat_n, aes(x = bact_time, y = lden)) +
  theme_bw() +
  MicrobioUoE::geom_pretty_boxplot(col = 'grey', fill = 'grey')+
  geom_point(position = position_jitter(0.2), shape = 21, size = 2, aes(fill = phage2)) +
  ylab("Estimated phage density log10(PFU/mL)") +
  xlab("Bacteria sample") +
  theme(axis.text = element_text(size = 14, colour = "black"), axis.title = element_text(size = 14), legend.position = 'bottom', legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  palettetown::scale_fill_poke(pokemon = "Roselia", spread = 11, labels = c("ISP", "Day 6 (pleural; pre-dose)","Day 6 (pleural; post-dose)")) +
  scale_x_discrete(labels = c("Day 0\n(pleural fluid)","Day 0\n(pus)","Day 6\n(BAL)")) +
  labs(fill = "Phage")
