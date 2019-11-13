library(readxl)


##### dataframe setup ##### 
nycsense <- read.csv('~/OpioidModel/nycsense.csv',header = FALSE)
nyssense <- read.csv('~/OpioidModel/nyssense.csv',header = FALSE)
nycsense <- nycsense[-1,]
nyssense <- nyssense[-1,]
colnames(nycsense) <- c('POP','TYPE','pres_rate', 'no_addict_rate', 'addict_frompres_rate','addict_fromaddict_rate','addict_rate','treat_entry_rate',
                            'rehab_success_rate','relapse_rate','death_rate','addict_od_death_rate')
colnames(nyssense) <- c('POP','TYPE','pres_rate', 'no_addict_rate', 'addict_frompres_rate','addict_fromaddict_rate','addict_rate','treat_entry_rate',
                            'rehab_success_rate','relapse_rate','death_rate','addict_od_death_rate')

zcrit <- 1.645
nycsense$POP <- as.character(nycsense$POP)
nycsense$TYPE <- as.character(nycsense$TYPE)
nycconf <- nycsense[which(nycsense$TYPE == '3' | nycsense$TYPE == '4'),]
for(i in 1:nrow(nycconf)){
  if(nycconf[i,'TYPE'] == '3'){
    nycconf[i,'TYPE'] <- '1'
  }else{
    nycconf[i,'TYPE'] <- '2'
  }
}
nycsense <- nycsense[which(nycsense$TYPE == '1' | nycsense$TYPE == '2'),]

pop_vector <- vector()
order_vector <- vector()
parameter_vector <- vector()
value_vector <- vector()
for (i in 1:nrow(nycsense)){
  for (j in 3:12){
    pop_vector <- c(pop_vector,nycsense[i,1])
    order_vector <- c(order_vector,nycsense[i,2])
    parameter_vector <- c(parameter_vector,colnames(nycsense[j]))
    value_vector <- c(value_vector,nycsense[i,j])
  }
}
nycsense.df <- data.frame(pop_vector,order_vector,parameter_vector,value_vector)
colnames(nycsense.df) <- c('pop','order','parameter','value')

pop_vector <- vector()
order_vector <- vector()
parameter_vector <- vector()
conf_vector <- vector()
for (i in 1:nrow(nycconf)){
  for (j in 3:12){
    pop_vector <- c(pop_vector,nycsense[i,1])
    order_vector <- c(order_vector,nycconf[i,2])
    parameter_vector <- c(parameter_vector,colnames(nycconf[j]))
    conf_vector <- c(conf_vector,nycconf[i,j])
  }
}
nycconf.df <- data.frame(pop_vector,order_vector,parameter_vector,conf_vector)
colnames(nycconf.df) <- c('pop','order','parameter','std_error')
nyc_total.df <- merge(nycsense.df,nycconf.df, by = c('pop','order','parameter'))

nyssense$POP <- as.character(nyssense$POP)
nyssense$TYPE <- as.character(nyssense$TYPE)
nysconf <- nyssense[which(nyssense$TYPE == '3' | nyssense$TYPE == '4'),]
for(i in 1:nrow(nysconf)){
  if(nysconf[i,'TYPE'] == '3'){
    nysconf[i,'TYPE'] <- '1'
  }else{
    nysconf[i,'TYPE'] <- '2'
  }
}

nyssense <- nyssense[which(nyssense$TYPE == '1' | nyssense$TYPE == '2'),]

pop_vector <- vector()
order_vector <- vector()
parameter_vector <- vector()
value_vector <- vector()
for (i in 1:nrow(nyssense)){
  for (j in 3:12){
    pop_vector <- c(pop_vector,nyssense[i,1])
    order_vector <- c(order_vector,nyssense[i,2])
    parameter_vector <- c(parameter_vector,colnames(nyssense[j]))
    value_vector <- c(value_vector,nyssense[i,j])
  }
}
nyssense.df <- data.frame(pop_vector,order_vector,parameter_vector,value_vector)
colnames(nyssense.df) <- c('pop','order','parameter','value')

pop_vector <- vector()
order_vector <- vector()
parameter_vector <- vector()
conf_vector <- vector()
for (i in 1:nrow(nysconf)){
  for (j in 3:12){
    pop_vector <- c(pop_vector,nysconf[i,1])
    order_vector <- c(order_vector,nysconf[i,2])
    parameter_vector <- c(parameter_vector,colnames(nysconf[j]))
    conf_vector <- c(conf_vector,nysconf[i,j])
  }
}
nysconf.df <- data.frame(pop_vector,order_vector,parameter_vector,conf_vector)
colnames(nysconf.df) <- c('pop','order','parameter','std_error')
nys_total.df <- merge(nyssense.df,nysconf.df, by = c('pop','order','parameter'))

##### z-tests #####
z_stat_list <- data.frame()
for (i in 1:nrow(nyc_total.df)){
  for(j in 1:nrow(nys_total.df)){
  same_pop <- nyc_total.df[i,'pop'] == nys_total.df[j,'pop']
  same_order <- nyc_total.df[i,'order'] == nys_total.df[j,'order']
  same_parameter <- nyc_total.df[i,'parameter'] == nys_total.df[j,'parameter']
  
  if(same_pop & same_parameter & same_order){
    pop <- as.character(nys_total.df[j,'pop'])
    if(as.character(nys_total.df[j,'order']) == '2'){
      order <- 't'
    }else{
      order <- '1'
    }
    parameter <- as.character(nys_total.df[j,'parameter'])
    diff <- abs(nyc_total.df[i,'value'] - nys_total.df[j,'value'])
    pooled_p <- (nyc_total.df[i,'value'] + nys_total.df[j,'value']) / 2
    std_dev <- (pooled_p*(1-pooled_p))
    std_error <- (std_dev/500)^0.5
    z_stat <- diff / std_error
    result_row <- data.frame(pop,order,parameter,diff,std_error,z_stat)
    z_stat_list <- rbind(z_stat_list,result_row)
    }
  }
}
z_stat_list <- data.frame(z_stat_list)
z_stat_list <- z_stat_list[order(z_stat_list[,6], decreasing = TRUE),]
z_stat_s <- z_stat_list[which(z_stat_list$pop == 's'),]
z_stat_p <- z_stat_list[which(z_stat_list$pop == 'p'),]
z_stat_a <- z_stat_list[which(z_stat_list$pop == 'a'),]
z_stat_r <- z_stat_list[which(z_stat_list$pop == 'r'),]

final_list <- data.frame()

for(z_list in list(z_stat_s,z_stat_p,z_stat_a,z_stat_r)){
  z_list <<- z_list[order(z_list[,6], decreasing = TRUE),]
  #Calculate P-values at which each group test will be evaluated
  alpha_true <<- 0.05
  p_crit <- vector()
  for (i in 1:nrow(z_list)){
    p_crit[i] <- 1 - (1 - alpha_true)^(1 / (nrow(z_list) - i + 1))
  }
  z_list <- cbind(z_list,p_crit)
  #Based on values above, set p-crit and t-crit for each individual value,
  #read off those values from a critical t table
  z_list[,'z_crit'] <- qnorm(1 - z_list$p_crit)
  final_list <- rbind(final_list,z_list)
}

final_list <- data.frame(final_list)

#Evaluate null hypothesis of no difference
final_list$reject_null <- NA
for (i in 1:nrow(final_list)){
  final_list[i,'reject_null'] <- final_list[i,6] > final_list[i,8]
}

write.csv(final_list, file = "~/Desktop/sensitivity_z-test.csv", row.names = FALSE)


