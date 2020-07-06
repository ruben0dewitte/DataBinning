# SETUP -------------------------------------------------------------------

#Load libraries
library(tidyverse)
library(here)
library(ggrepel) # Automatically repel text in graphs
library(distributionsrd)

#Styles
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Variables
firms_nr = 1e5 #Nr. of firms
bins_nr = 15 #Nr. of bins


# DATA --------------------------------------------------------------------

set.seed(1234)

x = data.frame(firms=1:firms_nr, sales = rpareto(n=firms_nr,xmin=0.01,k=1.01), dist="Pareto (shape = 1.01)")
y = data.frame(firms=1:firms_nr, sales = rlnorm(n=firms_nr,meanlog=log(mpareto(r=1,xmin=0.01,k=1.01,lower.tail=FALSE)),sdlog=1.74), dist="Lognormal (sdlog = 1.74)")
z = data.frame(firms=1:firms_nr, sales = rlnorm(n=firms_nr,meanlog=log(mpareto(r=1,xmin=0.01,k=1.01,lower.tail=FALSE)),sdlog=1.74), dist="Left-Truncated Lognormal (sdlog=1.74, trunc at 40\\%)") 

z = z %>%
  filter(sales>quantile(sales,0.4))

x = rbind.data.frame(y,z,x)
rm(y,z)
x$dist = factor(x$dist,levels=unique(x$dist),ordered=TRUE)


# Binning -----------------------------------------------------------------

# 1. Linear Binning
x = x %>%
  group_by(dist) %>%
  mutate(bin_lin = cut(sales,breaks=seq(from=min(sales),to=max(sales),length.out=bins_nr),include.lowest=TRUE))

test = x %>%
  filter(dist=="Pareto (shape = 1.01)") 

# 2. Logarithmic Binning
x = x %>%
  group_by(dist) %>%
  mutate(log_sales = log(sales)) %>%
  mutate(bin_log = cut(sales,breaks=c(exp(seq(from=min(log_sales),to=max(log_sales),length.out=bins_nr-1)),Inf),include.lowest=TRUE))

# 3. Variable binning
x = x %>%
  group_by(dist) %>%
  mutate(bin_var = cut(sales,breaks=c(exp(seq(from=min(log_sales),to=quantile(log_sales,0.8),length.out=(bins_nr-5))),quantile(sales,c(0.825,0.875,0.95,0.99,1))),include.lowest=TRUE))

# Gather all three methods
x = x %>%
  pivot_longer(cols=c("bin_lin","bin_log","bin_var"),
               names_to = "bin_method",
               values_to = "bin")
x$bin_method = factor(x$bin_method,levels=unique(x$bin_method),labels=c("Linear Binning","Logarithmic Binning", "Variable Binning"))

#Get lower and upper bounds + binwidth
x = x %>%
  mutate(binbounds = as.character(bin)) %>%
  separate(binbounds,into=c("one","bin_lower","bin_upper","four"),sep="[\\[|\\(|,|\\]|\\)]") %>%
  select(-c("one","four")) %>%
  mutate(bin_lower = as.numeric(bin_lower), bin_upper = as.numeric(bin_upper),binwidth = bin_upper - bin_lower)

# Visualization CCDF -----------------------------------------------------------

# Empirical CCDF
x = x %>%
  group_by(dist,bin_method) %>%
  mutate(CCDF = pempirical(q=sales,data=sales,lower.tail=FALSE))

# Summarize bins
x = x %>%
  group_by(dist,bin_method,bin) %>%
  mutate(bin_CCDF = min(CCDF), bin_sales = mean(sales)) 

# Vertical line for variable bins
x = x %>%
  group_by(dist) %>%
  mutate(xintercept = ifelse(bin_method=="Variable Binning",quantile(sales,0.8),NA))

# CCDF plot
fig_out = ggplot(data=x) + 
  geom_line(aes(x=sales,y=CCDF)) + 
  geom_point(aes(x=bin_sales,y=bin_CCDF)) +
  geom_vline(aes(xintercept=xintercept),linetype = "dashed") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap(~dist + bin_method,scales="free") + 
  theme_minimal() +
  theme(panel.grid = element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),legend.position="bottom",legend.title=element_blank()) +
  ylab(TeX("ln(CCDF)")) +
  xlab(TeX("ln(sales)")) 
ggsave(filename =here("paper","fig_CCDF.png"),plot=fig_out,width=30,height=30/2,units="cm")


# Visualization PDF -------------------------------------------------------

#Frequency
x = x %>%
  group_by(dist,bin_method,bin) %>%
  mutate(freq = n()/binwidth, freq_unadjust=n()) %>%
  group_by(dist,bin_method) %>%
  mutate(freq = freq/sum(freq_unadjust),freq_unadjust=freq/sum(freq_unadjust))

# Frequency plot
fig_out = ggplot(data=x) + 
  geom_point(aes(x=bin_sales,y=freq)) +
  geom_vline(aes(xintercept=xintercept),linetype = "dashed") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap(~dist + bin_method,scales="free") + 
  theme_minimal() +
  theme(panel.grid = element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),legend.position="bottom",legend.title=element_blank()) +
  scale_color_manual(cbbPalette) + 
  ylab(TeX("ln(Frequency)")) +
  xlab(TeX("ln(sales)")) 
ggsave(filename =here("paper","fig_PDF.png"),plot=fig_out,width=30,height=30/2,units="cm")


# Estimation CDF --------------------------------------------------------------

estimation <- function(y=y,x=x,weights=NULL){

  est_out = lm(y~x,weights=weights)
  
  return_list = tibble(constant = est_out$coefficients[[1]],k = est_out$coefficients[[2]], rsquare = summary(est_out)$r.squared) 
  return(return_list)
}

x_normal = x %>%
  group_by(dist,bin_method) %>%
  mutate(log_CCDF = log(CCDF)) %>%
  filter(is.finite(log_CCDF)) %>%
  do(estimation(y=.[["log_CCDF"]],x=.[["log_sales"]])) %>%
  unnest() %>%
  mutate(weighted="unbinned data")

# Summarize bins
y = x %>%
  group_by(dist,bin_method,bin) %>%
  summarise(bin_CCDF = min(CCDF), bin_sales = mean(sales),freq_unadjust = min(freq_unadjust)) %>%
  mutate(log_bin_CCDF = log(bin_CCDF), log_bin_sales = log(bin_sales)) %>%
  filter(is.finite(log_bin_CCDF)) 

# Unweighted
y_unweighted = y %>%
  group_by(dist,bin_method) %>%
  do(estimation(y=.[["log_bin_CCDF"]],x=.[["log_bin_sales"]])) %>%
  unnest() %>%
  mutate(weighted="unweighted data")

y_weighted = y %>%
  group_by(dist,bin_method) %>%
  do(estimation(y=.[["log_bin_CCDF"]],x=.[["log_bin_sales"]],weights=.[["freq_unadjust"]])) %>%
  unnest() %>%
  mutate(weighted="weighted data")

y_weights = rbind.data.frame(x_normal,y_weighted,y_unweighted)
rm(x_normal,y_weighted,y_unweighted)
y_weights$weighted = factor(y_weights$weighted,levels=unique(y_weights$weighted),ordered=TRUE)

y = y %>%
  left_join(y=y_weights,by=c("dist","bin_method")) %>%
  mutate(fit = exp(constant + k*log_bin_sales))

fig_out = ggplot(data=y) +
  geom_point(aes(x=bin_sales,y=bin_CCDF)) +
  geom_line(aes(x=bin_sales,y=fit,color=weighted)) + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap(~dist + bin_method,scales="free") + 
  geom_text_repel(
    data=y_weights,
    mapping=aes(x=Inf,y=Inf,label=paste0("k=",formatC(-k,format="f",digits=2), " (R^2 = ",formatC(rsquare,format="f",digits=2),")"),color=weighted),segment.color = 'transparent'
  ) + 
  scale_color_manual(values=cbbPalette) +
  theme_minimal() +
  theme(panel.grid = element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),legend.position="bottom") +
  labs(color="OLS with") +
  ylab(TeX("ln(CCDF)")) +
  xlab(TeX("ln(sales)")) 
ggsave(filename =here("paper","fig_CCDF_est.pdf"),plot=fig_out,width=30,height=30/2,units="cm")

# Estimation PDF ----------------------------------------------------------

estimation <- function(y=y,x=x,weights=NULL){
  
  est_out = lm(y~x,weights=weights)
  
  return_list = tibble(constant = est_out$coefficients[[1]],k = est_out$coefficients[[2]], rsquare = summary(est_out)$r.squared) 
  return(return_list)
}

# Summarize bins
y = x %>%
  group_by(dist,bin_method,bin) %>%
  summarise(freq = min(freq), bin_sales = mean(sales),bin_CCDF = min(CCDF),freq_unadjust = min(freq_unadjust)) %>%
  mutate(log_freq = log(freq), log_bin_sales = log(bin_sales)) %>%
  filter(is.finite(log_freq)) 

# Unweighted
y_unweighted = y %>%
  group_by(dist,bin_method) %>%
  do(estimation(y=.[["log_freq"]],x=.[["log_bin_sales"]])) %>%
  unnest() %>%
  mutate(weighted="unweighted data")

y_weighted = y %>%
  group_by(dist,bin_method) %>%
  do(estimation(y=.[["log_freq"]],x=.[["log_bin_sales"]],weights=.[["freq_unadjust"]])) %>%
  unnest() %>%
  mutate(weighted="weighted data")

y_weights = rbind.data.frame(y_weighted,y_unweighted)
rm(y_weighted,y_unweighted)
y_weights$weighted = factor(y_weights$weighted,levels=unique(y_weights$weighted),ordered=TRUE)

y = y %>%
  left_join(y=y_weights,by=c("dist","bin_method")) %>%
  mutate(fit = exp(constant + k*log_bin_sales))

fig_out = ggplot(data=y) +
  geom_point(aes(x=bin_sales,y=freq)) +
  geom_line(aes(x=bin_sales,y=fit,color=weighted)) + 
  facet_wrap(~dist + bin_method,scales="free") + 
  geom_text_repel(
    data=y_weights,
    mapping=aes(x=Inf,y=Inf,label=paste0("k=",formatC(-(k+1),format="f",digits=2), " (R^2 = ",formatC(rsquare,format="f",digits=2),")"),color=weighted),segment.color = 'transparent'
  ) +
  scale_x_log10() + 
  scale_y_log10() + 
  scale_color_manual(values=cbbPalette[-1]) +
  theme_minimal() +
  theme(panel.grid = element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),legend.position="bottom") +
  labs(color="OLS with") +
  ylab(TeX("ln(Frequency)")) +
  xlab(TeX("ln(sales)")) 
ggsave(filename =here("paper","fig_PDF_est.pdf"),plot=fig_out,width=30,height=30/2,units="cm")








