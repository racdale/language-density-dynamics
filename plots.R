library(lme4)
library(tidyverse)
library(ggplot2)
library(haven)
library(MuMIn)
library(corrr)
library(MetBrewer)
library(wesanderson)
library(cowplot)
library(broom.mixed)

########################################################################
#
# load prep data from A&E & population
#
########################################################################
 
step2 <- read_dta("data/step2.dta") %>% mutate(langCode = as_factor(language_cat),corpus = as_factor(corpus), fam = as_factor(fam))
step3 <- read_dta("data/step3.dta") %>% mutate(langCode = as_factor(language_cat), fam = as_factor(fam))

bibles <- read.csv("data/bibles_used_jun2019.txt") %>% 
  mutate(compressibility = 1-I(gzip/origSize)) %>% filter(compressionLevel==1) %>% 
  group_by(langCode,script) %>% summarize_if(is.numeric,mean,na.rm=TRUE)

eth <- read.csv("data/Table_of_Languages_recoded.txt",sep="\t") %>% 
  mutate(log_pop_eth = log10(All_Users+10),
         l1_prop = L1_Users/All_Users)

step2_eth <- left_join(step2,eth,by=c("langCode"="ISO_639"))

step2_eth_bibles <- left_join(step2_eth,bibles,by=c("langCode"="langCode"))

bible_sizes <- step2 %>% filter(corpus=="Bible_NT") %>% select(langCode,corpus_size)

step3_eth <- left_join(step3,eth,by=c("langCode"="ISO_639"))
step3_eth <- left_join(step3_eth,bible_sizes)

step3_eth_bibles <- left_join(step3_eth,bibles,by=c("langCode"="langCode"))

cor_mat <- step2_eth_bibles %>% filter(grepl("Bible_NT",corpus)) %>% group_by(langCode) %>%
	mutate(log10_corpus_size = log10(corpus_size)) %>% 
	summarize_if(is.numeric,mean,na.rm=TRUE) %>% 
	select(contains("density"),
				 log10_corpus_size,
				 log_pop_eth, # add pop size!!!!!
				 numUniqueWords
	) %>% 
	select_if(is.numeric) %>% 
	cor(m="s",use="p")

variable_order <- rownames(cor_mat)
	
########################################################################
#
# setup correlation data
#
########################################################################

cor_mat[lower.tri(cor_mat)] <- NA
cor_frame <- cor_mat %>% 
  data.frame() %>% 
  rownames_to_column(var = 'Variable_1') %>% 
  as_tibble() %>% 
  pivot_longer(-Variable_1, names_to = 'Variable_2', values_to = 'Correlation') 

cor_frame$Variable_1 <- factor(cor_frame$Variable_1, 
															 levels = variable_order,
															 labels = c('Semantic Density', 'Information Density',
															 					 'Info Density Chars', 'Corpus Size (log)',
															 					 'Population (log)', 'Num Unique Words'))
cor_frame$Variable_2 <- factor(cor_frame$Variable_2, levels = variable_order,
															 labels = c('Semantic Density', 'Information Density',
															 					 'Info Density Chars', 'Corpus Size (log)',
															 					 'Population (log)', 'Num Unique Words'))
cor_frame <- drop_na(cor_frame)

########################################################################
#
# correlation plot (lower triangle)
#
########################################################################

cor_plot <- cor_frame %>% 
  ggplot(aes(x = Variable_1, y = Variable_2, fill = Correlation)) +
  geom_tile(linewidth = 0.75, color = 'black') +
  geom_text(aes(label = round(Correlation, digits= 2)), 
  					) +
  theme_classic() +
  scale_fill_gradientn(colors = wes_palette('Zissou1', 100, type = 'continuous')) +
  scale_y_discrete(limits = rev) +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(face = 'bold',size=12),
        axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1,size=12),
  			plot.margin = unit(c(1, 1, 0.5, 1), 'cm'),
  			axis.line = element_blank(),
  			axis.ticks = element_blank()
  			) 

########################################################################
#
# build model for info density DV then plot
#
########################################################################

m3 <- step2_eth_bibles %>% filter(grepl("Bible_NT",corpus)) %>% 
	group_by(langCode,fam) %>% summarize_if(is.numeric,mean,na.rm=TRUE) %>%
  lmer(scale(information_density)~scale(semantic_density)*scale(log_pop_eth)+scale(numUniqueWords)+(1|fam),data=.)

m3_pair <- step2_eth_bibles %>% filter(grepl("Bible_NT",corpus)) %>%  # to check with 'sister' languages
	group_by(langCode,pair) %>% summarize_if(is.numeric,mean,na.rm=TRUE) %>%
	lmer(scale(information_density)~scale(semantic_density)*scale(log_pop_eth)+scale(numUniqueWords)+(1|pair),data=.)

all_stats_m3 <- broom.mixed::tidy(m3, conf.int = TRUE)
all_stats_m3$typ = 'fam'

all_stats_m3_pair <- broom.mixed::tidy(m3_pair, conf.int = TRUE) 
all_stats_m3_pair$typ = 'pair'

all_stats_m3 = rbind(all_stats_m3, all_stats_m3_pair) # get together for plot

coefs_m3 <- all_stats_m3 %>% 
  filter(effect == 'fixed') %>% 
  select(estimate, statistic, term, conf.low, conf.high, typ) %>% 
  filter(term != '(Intercept)') %>% 
  mutate(term = str_remove_all(term, pattern = 'scale\\('),
         term = str_remove_all(term, pattern = '\\)'),
  			 line_stop = ifelse(estimate < 0, conf.high, conf.low),
  			 line_stop = ifelse(0 > conf.low & 0 < conf.high, NA, line_stop),
  			 statistic_text = str_c('t = ', round(statistic, digits = 2)),
  			 position_text = ifelse(estimate < 0, -0.25, 0.25),
  			 term = case_match(term,
  			 									'log_pop_eth' ~ 'Population (log)',
  			 									'totalMatches' ~ 'Total Matches',
  			 									'semantic_density' ~ 'Semantic Density',
  			 									'information_density' ~ 'Information Density',
  			 									'numUniqueWords' ~ 'Num Unique Words',
  			 									'meanDistanceBackMatch' ~ 'Avg Dist Back Match',
  			 									'avgLengthMatch' ~ 'Avg Length Match',
  			 									'semantic_density:log_pop_eth' ~ 'Semantic Density x Population (log)',
  			 									'information_density:log_pop_eth' ~ 'Information Density x Population (log)'))


########################################################################
#
# plot info density DV
#
########################################################################

coefs_m3 <- coefs_m3 %>%
	arrange(desc(estimate)) 
coefs_m3$term <- factor(coefs_m3$term, levels = unique(coefs_m3$term))

coef_plot_info <- ggplot(coefs_m3, aes(x = estimate, y = term)) +
	geom_vline(xintercept = 0) + 
	geom_point(aes(shape=typ), size=3, color='maroon', fill='maroon', position=position_dodge(0.2)) + 
	geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width=.2, position=position_dodge(0.4)) + 
	theme_classic() +
	labs(x = 'Coefficient') + 
	theme(
		panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(), 
				axis.title.y = element_blank(),
				axis.text.y = element_text(face = 'bold', size=12),
				axis.title.x = element_text(face = 'bold', size=12),
				axis.text.x = element_text(face = 'bold', size=12),
				plot.margin = unit(c(0.25, 0.5, 1, 0.5), 'cm'), 
				legend.position = 'none'
				)

########################################################################
#
# build model for semantic density DV then plot
#
########################################################################

m3 <- step2_eth_bibles %>% filter(grepl("Bible_NT",corpus)) %>% 
	group_by(langCode,fam) %>% summarize_if(is.numeric,mean,na.rm=TRUE) %>%
	lmer(scale(semantic_density)~scale(numUniqueWords)*scale(log_pop_eth)+scale(information_density)+(1|fam),data=.)

m3_pair <- step2_eth_bibles %>% filter(grepl("Bible_NT",corpus)) %>%  # to check with 'sister' languages
	group_by(langCode,pair) %>% summarize_if(is.numeric,mean,na.rm=TRUE) %>%
	lmer(scale(semantic_density)~scale(numUniqueWords)*scale(log_pop_eth)+scale(information_density)+(1|pair),data=.)

all_stats_m3 <- broom.mixed::tidy(m3, conf.int = TRUE)
all_stats_m3$typ = 'fam'

all_stats_m3_pair <- broom.mixed::tidy(m3_pair, conf.int = TRUE) 
all_stats_m3_pair$typ = 'pair'

all_stats_m3 = rbind(all_stats_m3, all_stats_m3_pair) # get together for plot

coefs_m3 <- all_stats_m3 %>% 
	filter(effect == 'fixed') %>% 
	select(estimate, statistic, term, conf.low, conf.high, typ) %>% 
	filter(term != '(Intercept)') %>% 
	mutate(term = str_remove_all(term, pattern = 'scale\\('),
				 term = str_remove_all(term, pattern = '\\)'),
				 line_stop = ifelse(estimate < 0, conf.high, conf.low),
				 line_stop = ifelse(0 > conf.low & 0 < conf.high, NA, line_stop),
				 statistic_text = str_c('t = ', round(statistic, digits = 2)),
				 position_text = ifelse(estimate < 0, -0.25, 0.25),
				 term = case_match(term,
				 									'log_pop_eth' ~ 'Population (log)',
				 									'totalMatches' ~ 'Total Matches',
				 									'semantic_density' ~ 'Semantic Density',
				 									'information_density' ~ 'Information Density',
				 									'numUniqueWords' ~ 'Num Unique Words',
				 									'meanDistanceBackMatch' ~ 'Avg Dist Back Match',
				 									'avgLengthMatch' ~ 'Avg Length Match',
				 									'semantic_density:log_pop_eth' ~ 'Semantic Density x Population (log)',
				 									'information_density:log_pop_eth' ~ 'Information Density x Population (log)',
				 									'numUniqueWords:log_pop_eth' ~ 'N. Unique Words x Population (log)'))

########################################################################
#
# plot semantic density DV
#
########################################################################

coefs_m3 <- coefs_m3 %>%
	arrange(desc(estimate))
coefs_m3$term <- factor(coefs_m3$term, levels = unique(coefs_m3$term))

coef_plot_sem <- ggplot(coefs_m3, aes(x = estimate, y = term)) +
	geom_vline(xintercept = 0) + 
	geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.1) + 
	geom_point(aes(shape=typ), size=3, color='orange', fill='orange') + 
	#geom_segment(aes(x = 0, xend = line_stop), linetype = 'dotted') + 
	#geom_text(aes(label = statistic_text, x = position_text), vjust = 2, fontface = 'bold') + 
	theme_classic() +
	labs(x = 'Coefficient') + 
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		axis.title.y = element_blank(),
		axis.text.y = element_text(face = 'bold', size=12),
		axis.title.x = element_text(face = 'bold', size=12),
		axis.text.x = element_text(face = 'bold', size=12),
		plot.margin = unit(c(0.25, 0.5, 1, 0.5), 'cm'),
		legend.position = 'none'
	)

########################################################################
#
# assemble full panel
#
########################################################################

right_column_plot <- plot_grid(coef_plot_sem, coef_plot_info, ncol = 1)
final_plot <- plot_grid(cor_plot, right_column_plot, ncol = 2, rel_widths = c(1, .75))
ggsave(filename = 'plot_panel.png', plot = final_plot, 
			 width = 12, height = 6, units = 'in', dpi = 500, bg = 'white')

