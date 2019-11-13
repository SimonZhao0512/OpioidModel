args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
main <- function(test = 0){
  #Choose to run a graphics test or just default run
  if(test == 0){
modelresults <- read.csv('~/OpioidModel/modeloutput.csv',header = FALSE)
  }else{
    if(test == 1){
modelresults <- read.csv('~/OpioidModel/TestData/modeloutput_1.csv',header = FALSE)
}else if(test == 2){
  modelresults <- read.csv('~/OpioidModel/TestData/modeloutput_2.csv',header = FALSE)
}else if(test == 3){
  modelresults <- read.csv('~/OpioidModel/TestData/modeloutput_3.csv',header = FALSE)
}else (
  modelresults <- read.csv('~/OpioidModel/TestData/modeloutput_4.csv',header = FALSE)
)
}
model_type <- modelresults[1,1]
actual_addict_death <- as.vector(modelresults[1,4:10])
model_output <- modelresults[-1,]
if (model_type != 4){
  if(model_type == 2){
    colnames(model_output) <- c("Year","S->P","P->S", "A->R", "S", "P", "A", "R", "Opioid OD Deaths", "Abs Diff",'Treatment Entry')
  }else{
  colnames(model_output) <- c("Year","S->P","P->S", "A->R", "S", "P", "A", "R", "Opioid OD Deaths", "Abs Diff")
  }
}else{
  colnames(model_output) <- c('POP','TYPE','pres_rate', 'no_addict_rate', 'addict_frompres_rate','addict_fromaddict_rate','addict_rate','treat_entry_rate',
                             'rehab_success_rate','relapse_rate','death_rate','addict_od_death_rate')
}
  if (model_type == 1){
    
      #OPTION 1
      model_output$`P->S` <- as.character(model_output$`P->S`)
      YEAR <- vector()
      ACTUAL <- vector()
      for (i in 1:7){
        YEAR[i] <- 2010 + i
        ACTUAL[i] <- as.numeric(actual_addict_death[i])
      }
      
      actual_values <- as.data.frame(cbind(YEAR,ACTUAL))
      colnames(actual_values) <- c("Year","ACTUAL")
      plot <- ggplot(data = model_output, aes(x = model_output$Year, y = model_output$`Opioid OD Deaths`, colour = model_output$`P->S`)) + geom_polygon() 
      plot <- plot + xlab('Year') + ylab('Prescription Opioid Overdoses') 
      plot <- plot + ggtitle(label = 'Simulated Prescription Opioid Overdoses 2011-2017') + labs(caption = 'Varying alpha,epsilon,whatever') + theme(plot.caption = element_text(hjust = 0))  
      plot <- plot + geom_line(data = actual_values, aes(x = actual_values$Year, y = actual_values$ACTUAL), colour = 'black') + geom_point(data = actual_values, aes(x = actual_values$Year, y = actual_values$ACTUAL), colour = 'yellow')
      plot
   }else if(model_type == 2){
     
     #OPTION 2
      model_output$`S->P` <- as.character(format(model_output$`S->P`,digits = 2))
      colour_vector <- c('red','yellow','green','blue','purple')
      plot <- ggplot(data = model_output, aes(x = model_output$`P->S`, y = model_output$`A->R`)) + geom_line(aes(colour = model_output$`S->P`)) 
      plot <-  plot + xlab('Epsilon') + ylab('whatitcalled') 
      plot <- plot + ggtitle(label = 'epsiolon, alpha, and whatsitcalled Plot') + labs(caption = 'Combinations for S-ORD(2017)\nwithin desired neighborhood from A-ORD(2017)') 
      plot <- plot + theme(plot.caption = element_text(hjust = 0)) 
      plot
   }else if (model_type == 3){
     
     #OPTION 3
     parameter_table <- model_output[,c('S->P','P->S','A->R')]
     gap_table <- c('10-11:','11-12:','12-13:','13-14:','14-15:','15-16:','16-17:')
     parameter_plot <- 'Parameters (alpha,epsilon,whatitcalled):\n'
     for (i in 1:nrow(model_output)){
       if(i != 4){
     parameter_plot <- paste(parameter_plot,gap_table[i],' (',as.character(format(parameter_table[i,1], digits = 3)),',',as.character(format(parameter_table[i,2], digits = 3)),',',as.character(format(parameter_table[i,3], digits = 3)),'), ', sep = '')
       }else{
     parameter_plot <- paste(parameter_plot,gap_table[i],' (',as.character(format(parameter_table[i,1], digits = 3)),',',as.character(format(parameter_table[i,2], digits = 3)),',',as.character(format(parameter_table[i,3], digits = 3)),'),\n', sep = '')
       }
     }
     
     dataforplot <- as.data.frame(cbind(model_output$Year,model_output$`Opioid OD Deaths`))
     dataforplot$'ACTUAL' <- NA
     for (i in 1:length(actual_addict_death)){
       dataforplot$'ACTUAL'[i] = actual_addict_death[i]
     }
     dataforplot$ACTUAL <- as.numeric(dataforplot$ACTUAL)
     colnames(dataforplot) <- c("YEAR","PROJECTED","ACTUAL")
     plot <- ggplot(data = dataforplot) + geom_line(aes(x = dataforplot$YEAR, y = dataforplot$PROJECTED),colour = 'red', linetype = 'dashed')
     plot <- plot + geom_line(aes(x = dataforplot$YEAR, y = dataforplot$ACTUAL)) + geom_point(aes(x = dataforplot$YEAR, y = dataforplot$ACTUAL), colour = 'blue')
     plot <- plot + xlab('Year') + ylab('Prescription Opioid Overdoses') + ggtitle(label = 'Predicted vs. Actual Overdoses (2011-2017)') 
     plot <- plot + labs(caption = parameter_plot) + theme(plot.caption = element_text(hjust = 0))
     plot
   } else{
     #option 4
     allorder <- model_output
     allorder$POP <- as.character(allorder$POP)
     allorder$TYPE <- as.character(allorder$TYPE)
     
     pop_vector <- vector()
     order_vector <- vector()
     parameter_vector <- vector()
     value_vector <- vector()
     for (i in 1:nrow(allorder)){
     for (j in 3:12){
       pop_vector <- c(pop_vector,allorder[i,1])
       order_vector <- c(order_vector,allorder[i,2])
       parameter_vector <- c(parameter_vector,colnames(allorder[j]))
       value_vector <- c(value_vector,allorder[i,j])
      }
     }
     allorder.df <- data.frame(pop_vector,order_vector,parameter_vector,value_vector)
      colnames(allorder.df) <- c('pop','order','parameter','value')
     plot <- ggplot(data = allorder.df, aes( x = allorder.df$parameter, fill = allorder.df$pop)) + geom_col(aes(y = allorder.df$value))
     plot <- plot +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(~allorder.df$order) 
   }
return(plot)
}

main()


