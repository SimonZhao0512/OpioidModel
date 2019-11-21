args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
main <- function(){
  
modelresults <- read.csv('~/OpioidModel/modeloutput.csv',header = FALSE)

model_type <- modelresults[1,1]
if(model_type != 3){
actual_addict_death <- as.vector(modelresults[1,4:10])
}else{
  actual_addict_death <- as.vector(modelresults[1,8:14])
}
model_output <- modelresults[-1,]
if (model_type != 3){
  if(model_type == 2){
    colnames(model_output) <- c("Year","S->P","P->S", "A->R", "S", "P", "A", "R", "Opioid OD Deaths", "Abs Diff")
 }else {  
  colnames(model_output) <- c("Year","S->P","P->S", "A->R", "S", "P", "A", "R", "Opioid OD Deaths", "Abs Diff")
  }
}else{
  colnames(model_output) <- c('POP','TYPE','pres_rate', 'no_addict_rate', 'treat_entry_rate','addict_rate','illict_addict/pres','illicit_pres/pres',
                             'relapse_rate','illict_addict/illicit_pres','death_rate','addict_od_death_rate')
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
      plot <- plot + ggtitle(label = 'epsilon, alpha, and whatsitcalled Plot') + labs(caption = 'Combinations for S-ORD(2017)\nwithin desired neighborhood from A-ORD(2017)') 
      plot <- plot + theme(plot.caption = element_text(hjust = 0)) 
      plot
   }else if (model_type == 3){
     
     #OPTION 3
     
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


