iter <- 10000 #num iterations
RTO_PVA_raw <- data.frame(iteration=1:iter, matches= NA) #create unfilled data frame
n_select_RTO_raw<- nrow(RTO_raw_loci_select)
n_select_PVA_raw<- nrow(PVA_raw_loci_select)

for(i in 1:iter){
  RTO_sample_raw <- sample(RTO_raw_loci$BAR, n_select_RTO_raw, replace = FALSE)
  PVA_sample_raw <- sample(PVA_raw_loci$BAR, n_select_PVA_raw, replace = FALSE)
  match_count <- 0
  for(z in 1:length(RTO_sample_raw)) {
    if(check_match(locus = RTO_sample_raw[z], thresh = dist, table = PVA_sample_raw))
      match_count <- match_count + 1
  }
  RTO_PVA_raw$matches[i] <- match_count
} # Runs iterations of matches-- populates RTO_PVA$matches

nrow(RTO_PVA_raw)
mean_RTO_PVA_raw<- mean(RTO_PVA_raw$matches)

##Plot and determine significance
###Observed matches
obs_loci_RTO_PVA <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(obs_loci_RTO_PVA) <- colnames(RTO_raw_loci)

for(i in 1:nrow(RTO_screen_loci_select)){
  if(check_match(locus = RTO_screen_loci_select$BAR[i], thresh = dist, table = PVA_screen_loci_select$BAR))
    obs_loci_RTO_PVA <- rbind(obs_loci_RTO_PVA, RTO_screen_loci_select[i,])
}

obs<-nrow(obs_loci_RTO_PVA)
obs

###Probability distribution function
prob_df_raw<- ggplot(RTO_PVA_raw, aes(x=matches)) +geom_histogram(binwidth=1) + geom_vline(aes(xintercept= mean_RTO_PVA_raw), linetype="dashed") + geom_vline(aes(xintercept = obs), colour="green") + (labs(title= "Probability distribution between P. vacca and R. toxotes"))

ggsave(filename = "figs/RTO_PVA_prob_df_raw.png", plot = prob_df_raw)

###Cumulative distribution function
cdf_raw<- data.frame(x = 1:(max(RTO_PVA_raw$matches)), y = NA)
cdf_raw$y<- cumsum(table(RTO_PVA_raw$matches)) / nrow(RTO_PVA_raw)
qtile_raw<- cdf_raw[obs,2]
p_raw <- as.character(1-qtile_raw)


cdf_plot_raw<- ggplot(cdf_raw, aes(x,y)) + 
  geom_point(size=3) + 
  labs (x = "Number of matches", y = "Cumulative probability") +
  ggtitle("CDF of matches between R. toxotes and P. vacca \n p =", p_raw)+
  scale_x_continuous(breaks=seq(0,nrow(cdf), by = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
  geom_segment(x=obs, y=0, xend=obs, yend=qtile_raw, color = 'red', linetype = "dashed", size= 1) +
  geom_segment(x= 0, y = qtile_raw, xend = obs, yend = qtile_raw, color = 'red', linetype = "dashed", size = 1)

ggsave("figs/RTO_PVA_cum_df_raw.png", cdf_plot_raw)


