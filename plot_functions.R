# This is the first plot function 
plot_function1 <- function(true_f, melt_y, melt_result){
  ggplot(data=melt_y) +
    # all the trials 
    geom_line(aes(x = time, y = value, group = factor(trial), color= "Trials"), alpha = 0.25) +
    # empirical mean
    geom_line(data = melt_y %>% dplyr::group_by(time)
              %>% dplyr::summarize(mean = mean(value)), 
              aes(x = time, y = mean, color= "Empirical Mean"), linewidth = 1, linetype =8) + 
    # true_f 
    geom_line(data =data.frame(time = 1:n_time, f = dat$f), 
              aes(x = time, y = f,color="True f"), linetype = 1, size = 1.5) + 
    # all the sample_fs (draws) to estimate f
    #geom_line(data=melt_result,
    # aes(x=time,y=value,group=factor(draws),color = "Draws of f"),alpha=0.25) + 
    # median (final estimate) of f
    geom_line(data= melt_result %>% dplyr::group_by(time)
              %>% dplyr::summarize(median=median(value)),
              aes(x=time,y=median, color ="Final Estimate of f"),linewidth= 1,linetype=1) + 
    scale_color_manual(values = c("Trials" = "darkgrey", 
                                  "Empirical Mean" = "darkgreen", 
                                  "True f" = "red",
                                  "Draws of f" = "darkgrey",
                                  "Final Estimate of f" = "blue")) +
    labs(title = "Estimate of f with Trials from Data",
         x = "Time",
         y = "Values",
         color = "Line Names")  + # Legend title 
    #guides(colour = guide_legend(override.aes = list(pch = c(16, 21, 16), fill = c("red", "white", "red")))) +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
}


# This is the second plot function 
plot_function2 <- function(true_f, melt_y, melt_result){
  ggplot(data=melt_y) +
    # all the trials 
    #geom_line(aes(x = time, y = value, group = factor(trial), color= "Trials"), alpha = 0.25) +
    # empirical mean
    geom_line(data = melt_y %>% dplyr::group_by(time)
              %>% dplyr::summarize(mean = mean(value)), 
              aes(x = time, y = mean, color= "Empirical Mean"), linewidth = 1, linetype =8) + 
    # true_f 
    geom_line(data =data.frame(time = 1:n_time, f = dat$f), 
              aes(x = time, y = f,color="True f"), linetype = 1, size = 1.5) + 
    # all the sample_fs (draws) to estimate f
    geom_line(data=melt_result, aes(x=time,y=value,group=factor(draws),color = "Draws of f"),alpha=0.25) + 
    # median (final estimate) of f
    geom_line(data= melt_result %>% dplyr::group_by(time)
              %>% dplyr::summarize(median=median(value)),
              aes(x=time,y=median, color ="Final Estimate of f"),linewidth= 1,linetype=1) + 
    scale_color_manual(values = c("Trials" = "darkgrey", 
                                  "Empirical Mean" = "darkgreen", 
                                  "True f" = "red",
                                  "Draws of f" = "darkgrey",
                                  "Final Estimate of f" = "blue")) +
    labs(title = "Estimate of f with Draws from f",
         x = "Time",
         y = "Values",
         color = "Line Names")  + # Legend title 
    #guides(colour = guide_legend(override.aes = list(pch = c(16, 21, 16), fill = c("red", "white", "red")))) +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
}


###################################################################
# This function is used to calculate the MSE along with the two plots

Calculate_MSE <- function(true_f, melt_result){
  library(dplyr)
  median_estimate <- melt_result %>% 
    group_by(time) %>% 
    summarise(median = median(value)) %>%
    dplyr::select(median)
  MSE <- mean((true_f-unlist(median_estimate))^2)
  return(MSE)
}
