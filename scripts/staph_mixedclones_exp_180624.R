#bacteria and phage densities in single genotype and mixed staph from patient 3

library(tidyverse)
library(ggplot2)
library(emmeans)
library(lme4)

bact <- read.csv('data/staph_mixeddiversity_bactden_040624.csv',header=T)

bact <- na.omit(bact)

bact_m <- group_by(bact, treat, treat2, phage, time) %>%
  summarise(., m_den = mean(lden),
            se = sd(lden)/sqrt(length(lden)))

fac_labs <- c("(a) Isolate 1", "(b) Isolate 2", "(c) Isolate 3", "(d) Isolate 4", "(e) Mix of 24 Isolates")
names(fac_labs) <- c("clone1", "clone2", "clone3", "clone4", "mix")

ggplot() +
  theme_bw() +
  geom_point(data = bact, aes(x = time, y = lden, group = interaction(time, phage, treat), col = treat2), alpha = 0.5) +
  geom_line(data = bact, aes(x = time, y = lden, group = interaction(phage, treat, rep), linetype = phage, col = treat2), alpha = 0.2) +
  geom_point(data = bact_m, aes(x = time, y = m_den, group = interaction(time, phage, treat), col = treat2), size = 2, position = position_dodge(0.2)) +
  geom_line(data = bact_m, aes(x = time, y = m_den, group = interaction(phage, treat), linetype = phage, col = treat2), linewidth = 1, position = position_dodge(0.2)) +
  geom_errorbar(data = bact_m, aes(x = time, ymin = m_den-se, ymax = m_den+se, group = interaction(time, phage, treat), col = treat2), size = 1, width = 0.1, position = position_dodge(0.2)) +
  theme(axis.text.y = element_text(size = 13, colour = "black"), axis.title = element_text(size = 14), axis.text.x = element_text(size = 12, colour = "black"), title = element_text(size = 14),legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 12), strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0)) +
  ylab("Bacterial density log10(CFU/mL)") +
  xlab("Time (days)") +
  scale_x_discrete(labels = c("1","2","3","4","5","6")) +
  palettetown::scale_color_poke(pokemon = "charizard", spread = 5, labels = c('Isolate', 'Mixed-isolate')) +
  labs(linetype = "Phage", col = "Treatment") +
  facet_wrap(~treat, labeller = labeller(treat = fac_labs))+
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c('Absent', 'Present'))

#summary statistics of persister and recovered densities

bact_p <- filter(bact, phage == "Y" & dil_fact == 5 & time == 't6')
bact_p <- filter(bact_p, ! den == 0)
mean(bact_p$lden)
sd(bact_p$lden)/sqrt(length(bact_p$lden))

bact_r <- filter(bact, phage == "Y" & dil_fact == 40 & time == 't6')
mean(bact_r$lden)
sd(bact_r$lden)/sqrt(length(bact_r$lden))

bact_np <- filter(bact, phage == "N" & time == 't6')
mean(bact_np$lden)
sd(bact_np$lden)/sqrt(length(bact_np$lden))

#densities

bact$treat[bact$treat == 'clone1'] <- 'Isolate1'
bact$treat[bact$treat == 'clone2'] <- 'Isolate2'
bact$treat[bact$treat == 'clone3'] <- 'Isolate3'
bact$treat[bact$treat == 'clone4'] <- 'Isolate4'

bact$rep2 <- interaction(bact$treat, bact$rep)

m1 <- lmer(lden ~ treat * time * phage + (1|rep2), data = bact)
m2 <- lmer(lden ~ treat + time + phage + treat:time + treat:phage + time:phage + (1|rep2), data = bact)
anova(m1,m2)

m3 <- lmer(lden ~ treat + time + phage + treat:phage + time:phage + (1|rep2), data = bact) #treat:time

anova(m2,m3) #0.9816
m4 <- lmer(lden ~ treat + time + phage + treat:time + time:phage + (1|rep2), data = bact) #treat:phage
anova(m2,m4) #treat:phage very sig < 2.2e-16
m5 <- lmer(lden ~ treat + time + phage + treat:time + treat:phage + (1|rep2), data = bact) #time:phage
anova(m2,m5) #0.02696
#model cannot be simplified further. Best model is m3
emmeans::emmeans(m3, pairwise ~ treat|phage) #most treatments different to mix. Only clone 2 isn't, maybe bc this one has two reps recover and 1 persist (rest either extinct (clone 4) or only one recovered/1 persist)
emmeans::emmeans(m3, pairwise ~ time|phage) #effect driven by populations beginning to recover at day 6 generally

#Table of contrasts

mens2 <- data.frame(emmeans::emmeans(m3, pairwise ~ treat|phage)$contrasts)

mens2 <- mens2[,c(-4,-5)]

mens2 <- mutate(mens2,
                estimate = round(estimate, digits = 3),
                t.ratio = round(t.ratio, digits = 3),
                p.value = round(p.value, digits = 3))

mens2$p.value[mens2$p.value == 0] <- "<0.001"

library(flextable)

men_flex2 <- flextable(mens2) %>%
  set_header_labels(contrast = "Contrast", phage = "Phage", estimate = "Estimate", t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 11, part = 'all') %>%
  autofit() %>%
  align(j = c(1:4), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header') %>%
  hline(i = c(10,20), border = officer::fp_border(color="black")) 

mens3 <- data.frame(emmeans::emmeans(m3, pairwise ~ time|phage)$contrasts)

mens3 <- mens3[,c(-4,-5)]

mens3 <- mutate(mens3,
                estimate = round(estimate, digits = 3),
                t.ratio = round(t.ratio, digits = 3),
                p.value = round(p.value, digits = 3))

mens3$p.value[mens3$p.value == 0] <- "<0.001"

mens3$contrast <- gsub("t", "", mens3$contrast)

men_flex3 <- flextable(mens3) %>%
  set_header_labels(contrast = "Contrast", phage = "Phage", estimate = "Estimate", t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 11, part = 'all') %>%
  autofit() %>%
  align(j = c(1:4), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header') %>%
  hline(i = c(15,30), border = officer::fp_border(color="black")) %>%
  line_spacing(i = NULL, j = NULL, space = 0.5, part = "body")


#
res <- read.csv('data/staph_res_data.csv',header=T)

res2 <- res[,c(1,2,8)]

bact6 <- filter(bact, time == 't6' & phage == 'Y')

resb <- merge(res2, bact6, by = c('treat', 'rep'))
resb2 <- merge(bact6, res2, by = c('treat', 'rep'))

resb2$majority[resb2$majority == 'resistant / unculturable'] <- 'unculturable'

a <- ggplot(resb2, aes(x = treat2, y = lden)) +
  geom_point(aes(fill = majority), size = 3, position = position_jitter(0.1), shape = 21, col = 'black') +
  annotate('text', x = 1, y = 2, label = 'Extinct = 17') +
  annotate('text', x = 2, y = 2, label = 'Extinct = 1') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 13, colour = "black"), axis.title = element_text(size = 14), axis.text.x = element_text(size = 12, colour = "black"), title = element_text(size = 12), legend.text = element_text(size = 11), legend.title = element_text(size = 11), legend.position = 'bottom') +
  ylab('Bacteria density on day six\nlog10(CFU/mL)') +
  xlab('Treatment') +
  scale_x_discrete(labels = c('Isolate', 'Mixed-isolate')) +
  labs(fill = 'Dominant\nPhenotype', title = '(a)') +
  scale_fill_manual(values = c('#489E00', 'black', '#FF8751'), labels = c('Phage\nresistant', 'Phage\nsusceptible', 'Unculturable')) 

a

#phage

pden <- read.csv('data/staph_mixeddiversity_phageden_040624.csv',header=T)

pden <- na.omit(pden)

pden_m <- group_by(pden, treat, treat2, time) %>%
  summarise(., m_den = mean(lden),
            se = sd(lden)/sqrt(length(lden)))

pden$rep2 <- interaction(pden$treat, pden$rep)

m1 <- lmer(lden ~ treat * time + (1|rep2), data = pden)
m2 <- lmer(lden ~ treat + time + (1|rep2), data = pden)
anova(m1,m2) #sig
emmeans::emmeans(m1, pairwise ~ treat | time) #some sig comparisons between clone and mix at t4 (clone 2 and 4) and t5 (clone 4 and mix)

input_d <- log10((6.6 * 10^5 * 21.8)/6)

ggplot() +
  theme_bw() +
  geom_point(data = pden, aes(x = time, y = lden, group = interaction(time, treat), col = treat2), alpha = 0.5) +
  geom_line(data = pden, aes(x = time, y = lden, group = interaction(treat, rep), col = treat2), alpha = 0.2) +
  geom_point(data = pden_m, aes(x = time, y = m_den, group = interaction(time, treat), col = treat2), size = 2, position = position_dodge(0.2)) +
  geom_line(data = pden_m, aes(x = time, y = m_den, group = treat, col = treat2), linewidth = 1, position = position_dodge(0.2)) +
  geom_errorbar(data = pden_m, aes(x = time, ymin = m_den-se, ymax = m_den+se, group = interaction(time, treat), col = treat2), size = 1, width = 0.1, position = position_dodge(0.2)) +
  theme(axis.text.y = element_text(size = 13, colour = "black"), axis.title = element_text(size = 14), axis.text.x = element_text(size = 12, colour = "black"), title = element_text(size = 14),legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 11), strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0)) +
  ylab("Phage density log10(PFU/mL)") +
  xlab("Time (days)") +
  scale_x_discrete(labels = c("1","2","3","4","5","6")) +
  palettetown::scale_color_poke(pokemon = "charizard", spread = 5, labels = c('Isolate', 'Mixed-isolate')) +
  geom_hline(yintercept = input_d, linetype = 'dashed') +
  facet_wrap(~treat, labeller = labeller(treat = fac_labs)) +
  labs(color = 'Treatment')

pden$treat[pden$treat == 'clone1'] <- 'Isolate1'
pden$treat[pden$treat == 'clone2'] <- 'Isolate2'
pden$treat[pden$treat == 'clone3'] <- 'Isolate3'
pden$treat[pden$treat == 'clone4'] <- 'Isolate4'

pden_m$treat[pden_m$treat == 'clone1'] <- 'Isolate1'
pden_m$treat[pden_m$treat == 'clone2'] <- 'Isolate2'
pden_m$treat[pden_m$treat == 'clone3'] <- 'Isolate3'
pden_m$treat[pden_m$treat == 'clone4'] <- 'Isolate4'

d_mods <- pden %>%
  nest(-treat, -time) %>%
  mutate(., mod = map(data, ~ t.test(.x$lden, mu = input_d)))

d_mods_summary <- d_mods %>%
  unnest(d_mods, mod %>% map(broom::tidy)) %>%
  mutate(padjust = p.adjust(p.value, "fdr"))

d_mods_summary$sig <- NA

d_mods_summary$sig[d_mods_summary$padjust < 0.05] <- '*'

fac_labs <- c("(a) Isolate 1", "(b) Isolate 2", "(c) Isolate 3", "(d) Isolate 4", "(e) Mix of 24 Isolates")
names(fac_labs) <- c("Isolate1", "Isolate2", "Isolate3", "Isolate4", "mix")

ggplot() +
  theme_bw() +
  geom_point(data = pden, aes(x = time, y = lden, group = interaction(time, treat), col = treat2), alpha = 0.5) +
  geom_line(data = pden, aes(x = time, y = lden, group = interaction(treat, rep), col = treat2), alpha = 0.2) +
  geom_point(data = pden_m, aes(x = time, y = m_den, group = interaction(time, treat), col = treat2), size = 2, position = position_dodge(0.2)) +
  geom_line(data = pden_m, aes(x = time, y = m_den, group = treat, col = treat2), linewidth = 1, position = position_dodge(0.2)) +
  geom_errorbar(data = pden_m, aes(x = time, ymin = m_den-se, ymax = m_den+se, group = interaction(time, treat), col = treat2), size = 1, width = 0.1, position = position_dodge(0.2)) +
  theme(axis.text.y = element_text(size = 13, colour = "black"), axis.title = element_text(size = 14), axis.text.x = element_text(size = 12, colour = "black"), title = element_text(size = 14),legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 11)) +
  ylab("Phage density log10(PFU/mL)") +
  xlab("Time (days)") +
  scale_x_discrete(labels = c("1","2","3","4","5","6")) +
  palettetown::scale_color_poke(pokemon = "charizard", spread = 5) +
  geom_hline(yintercept = input_d, linetype = 'dashed') +
  facet_wrap(~treat,labeller = labeller(treat = fac_labs)) +
  geom_text(data = d_mods_summary, aes(x = time, y = 9.5, label = sig)) +
  labs(color = 'Treatment')

ggplot() +
  theme_bw() +
  geom_point(data = pden, aes(x = time, y = lden, group = interaction(time, treat), col = treat2), alpha = 0.5) +
  geom_line(data = pden, aes(x = time, y = lden, group = interaction(treat, rep), col = treat2), alpha = 0.2) +
  geom_point(data = pden_m, aes(x = time, y = m_den, group = interaction(time, treat), col = treat2), size = 2, position = position_dodge(0.2)) +
  geom_line(data = pden_m, aes(x = time, y = m_den, group = treat, col = treat2), linewidth = 1, position = position_dodge(0.2)) +
  geom_errorbar(data = pden_m, aes(x = time, ymin = m_den-se, ymax = m_den+se, group = interaction(time, treat), col = treat2), size = 1, width = 0.1, position = position_dodge(0.2)) +
  theme(axis.text.y = element_text(size = 13, colour = "black"), axis.title = element_text(size = 14), axis.text.x = element_text(size = 12, colour = "black"), title = element_text(size = 14),legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 11), strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0)) +
  ylab("Phage density log10(PFU/mL)") +
  xlab("Time (days)") +
  scale_x_discrete(labels = c("1","2","3","4","5","6")) +
  palettetown::scale_color_poke(pokemon = "charizard", spread = 5, labels = c('Clone', 'Mixed-clone')) +
  geom_hline(yintercept = input_d, linetype = 'dashed') +
  facet_wrap(~treat, labeller = labeller(treat = fac_labs))+
  geom_text(data = d_mods_summary, aes(x = time, y = 9.5, label = sig)) +
  labs(color = 'Treatment')

d_sum2 <- d_mods_summary
d_sum2 <- d_sum2[,c(1,2,9:11,17)]

d_sum2$time <- gsub("t", "", d_sum2$time)

d_sum2 <- mutate(d_sum2,
                estimate = round(estimate, digits = 3),
                statistic = round(statistic, digits = 3),
                p.value = round(p.value, digits = 3),
                padjust = round(padjust, digits = 3))

d_sum2$p.value[d_sum2$p.value == 0] <- "<0.001"
d_sum2$padjust[d_sum2$padjust == 0] <- "<0.001"

d_sum2_flex3 <- flextable(d_sum2) %>%
  set_header_labels(treat = 'Treatment', time = 'Time', estimate = "Estimate", statistic = "t-value", p.value = "p-value", padjust = 'p-value (fdr)') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 11, part = 'all') %>%
  autofit() %>%
  align(j = c(1:4), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header') %>%
  hline(i = c(5,10,15,20,25,30), border = officer::fp_border(color="black")) %>%
  line_spacing(i = NULL, j = NULL, space = 0.5, part = "body")


#Can 'population state' predict phage density as some variability and suggestion of phage densities increasing?

pden$majority <- NA

pden6 <- filter(pden, time == 't6')

res2$majority[res2$majority == 'resistant / unculturable'] <- 'unculturable'

res2$rep2 <- interaction(res2$treat, res2$rep)

pden6$majority[pden6$rep2 == 'Isolate1.1'] <- 'susceptible'
pden6$majority[pden6$rep2 == 'Isolate1.3'] <- 'resistant'
pden6$majority[pden6$rep2 == 'Isolate2.3'] <- 'unculturable'
pden6$majority[pden6$rep2 == 'Isolate2.5'] <- 'resistant'
pden6$majority[pden6$rep2 == 'Isolate2.6'] <- 'resistant'
pden6$majority[pden6$rep2 == 'Isolate3.3'] <- 'resistant'
pden6$majority[pden6$rep2 == 'Isolate3.6'] <- 'unculturable'
pden6$majority[pden6$rep2 == 'mix.1'] <- 'unculturable'
pden6$majority[pden6$rep2 == 'mix.2'] <- 'unculturable'
pden6$majority[pden6$rep2 == 'mix.3'] <- 'unculturable'
pden6$majority[pden6$rep2 == 'mix.4'] <- 'susceptible'
pden6$majority[pden6$rep2 == 'mix.5'] <- 'unculturable'
pden6$majority[pden6$rep2 == 'mix.6'] <- 'resistant'
pden6$majority[pden6$rep2 == 'mix.8'] <- 'resistant'
pden6$majority[pden6$rep2 == 'mix.9'] <- 'resistant'
pden6$majority[pden6$rep2 == 'mix.10'] <- 'unculturable'
pden6$majority[pden6$rep2 == 'mix.11'] <- 'resistant'
pden6$majority[pden6$rep2 == 'mix.12'] <- 'resistant'

pden6$majority[is.na(pden6$majority)] <- 'extinct'

#omit susceptible populations as only two so not meaningful analysis

pden6 <- filter(pden6, ! majority == 'susceptible')

m1 <- lm(lden ~ majority, data = pden6)
m2 <- lm(lden ~ 1, data = pden6)
anova(m1,m2)
emmeans::emmeans(m1, pairwise ~ majority)

#change majority to dominant phenotype?

b <- ggplot(pden6, aes(x = majority, y = lden)) +
  theme_bw() +
  MicrobioUoE::geom_pretty_boxplot(col = 'black', fill = 'black')+
  geom_point(shape = 21, size = 3, position = position_jitter(0.1), aes(fill = majority))  +
  scale_x_discrete(labels = c('Extinct', 'Recovered\nand resistant', 'Low-density\nand unculturable')) +
  xlab('Bacteria population phenotype')+
  theme(axis.text.y = element_text(size = 13, colour = "black"), axis.title = element_text(size = 14), axis.text.x = element_text(size = 12, colour = "black"), title = element_text(size = 12), legend.position = 'none') +
  ylab('Phage density on day six\nlog10(PFU/mL)') +
  scale_fill_manual(values = c('white', '#489E00', '#FF8751'), labels = c('Extinct', 'Phage\nresistant', 'Unculturable')) +
  labs(fill = 'Dominant\nPhenotype', title = '(b)') +
  geom_hline(yintercept = input_d, linetype = 'dashed') 

library(patchwork)

a + b + plot_layout()

#extinct vs resistant are not significantly different. Unculturable has sig higher phage densities overall. 
#is this higher than baseline inoculum?

d_mods <- pden6 %>%
  nest(-majority) %>%
  mutate(., mod = map(data, ~ t.test(.x$lden, mu = input_d)))

d_mods_summary <- d_mods %>%
  unnest(d_mods, mod %>% map(broom::tidy)) %>%
  mutate(padjust = p.adjust(p.value, "fdr"))

#unculturable is significantly higher than others suggesting that phages are replicating
