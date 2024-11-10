##libraries
library(broom)
library(dplyr)
library(Biobase)
library(preprocessCore)
library(devtools)
library(caret)
library(dabestr)
library(pwr)
library(AICcmodavg)
library(BayesFactor)
library(genefilter)
#options
colramp = colorRampPalette(c(3,"white",2))(20)
dir()
##Compute sample sizes needed based on medium effect size expected  
pwr.t2n.test(n1 = , n2=65 , d =0.5 , sig.level =0.1, power=0.85)
##live imaging hca 
## analysis per time condition
## t-test 
#read data by time condition
t_data = read.csv('D:/MiguelW12/Documents/stats_script_try/gsd1a_new_stats/frames/frame_72_ee.csv', header=TRUE)
p_data = read.csv('D:/MiguelW12/Documents/stats_script_try/gsd1a_new_stats/frames/72_ee_p.csv', header=TRUE)
dim(p_data)
dim(t_data)
# sort matrixes 
rownames(t_data2) <- t_data2$group_id_pc
t_data2 = subset(t_data2, select = -c(group_id_pc))
# define group to test 
group = p_data$group_with_pc
#run test, save results
res_24 = rowttests(as.matrix(t_data),factor(group))
res_48 = rowttests(as.matrix(t_data),factor(group))
res_72 = rowttests(as.matrix(t_data),factor(group))
##export res
write.csv(as.data.frame(res_72_ee), file="res_72_ee_folowproc.csv")
# run these lines to evaluate results in excel 
# =IF(D2<0.1, TRUE)
# =IF(E2=TRUE,1,0)
# =IF(E2=FALSE,1,0)
# visalize res 
h_24 = hist(res_24$p.value)
h_48 = hist(res_48$p.value)
h_72 = hist(res_72$p.value)


####### data stats summary function

summarize_dataframe <- function(df) {
  # data frame for the summary data
  summary_df <- data.frame(
    column_name = character(),
    mean = numeric(),
    variance = numeric(),
    num_samples = integer(),
    group_means = numeric(),
    group_variances = numeric(),
    group_num_samples = integer(),
    total_mean = numeric(),
    total_variance = numeric(),
    total_num_samples = integer(),
    calculated_effect_size = numeric(),
    S = numeric(),
    stringsAsFactors = FALSE
  )
  #calculate summary statistics
  for (col_name in names(df)) {
    if (col_name != "index") {
      # mean, variance & number of samples 
      col_mean <- mean(df[[col_name]])
      col_var <- var(df[[col_name]])
      col_num_samples <- length(df[[col_name]])
      group_means <- tapply(df[[col_name]], df$index, mean)
      group_variances <- tapply(df[[col_name]], df$index, var)
      group_sd <- sqrt(group_variances)
      group_num_samples <- tapply(df[[col_name]], df$index, length)
      total_mean <- mean(df[[col_name]])
      total_var <- var(df[[col_name]])
      total_num_samples <- length(df[[col_name]])
      mean_dif <- group_means[1] - group_means[2]
      # S value 
      s_value <- sqrt(
        ((group_num_samples[1] - 1) * group_variances[1] + (group_num_samples[2] - 1) * group_variances[2]) /
          (total_num_samples - 2)
      )
      v_value <- (
        ((group_num_samples[1] - 1) * group_variances[1] + (group_num_samples[2] - 1) * group_variances[2]) /
          (total_num_samples - 2)
      )
      marg <- qt(0.975, df=total_num_samples -1) * sqrt(v_value/group_num_samples[1] + v_value/group_num_samples[2])
      low_int <- mean_dif - marg
      up_int <- mean_dif + marg
      ci_up <- up_int + group_means
      # effect size
      effect_size <- abs(group_means[1] - group_means[2]) / s_value
      # new data frame with summary statistics 
      summary_col <- data.frame(
        column_name = col_name,
        mean = col_mean,
        variance = col_var,
        num_samples = col_num_samples,
        group_means = group_means,
        group_variances = group_variances,
        group_num_samples = group_num_samples,
        total_mean = total_mean,
        total_variance = total_var,
        total_num_samples = total_num_samples,
        calculated_effect_size = effect_size,
        S = s_value,
        V = v_value,
        mean_dif = mean_dif,
        marg = marg,
        low_int = low_int,
        up_int = up_int,
        group_sd = group_sd
        
      )
      
      summary_df <- rbind(summary_df, summary_col)
    }
  }
  summary_df
}

##initiate
summary_df <- summarize_dataframe(t_data)
write.csv(as.data.frame(summary_df), file="gsd1a_frame_72-summary.csv")

##95% CI Plots based on mean-diff
two_row_ci_plot <- function(df, path) {
  
  # empty list to contain plots
  plot_list <- list()
  n <- nrow(df)
  # loop  the data frame two rows at a time
  for (i in seq(1, n, by=2)) {
    
    # get the two rows
    df_small <- df[i,]
    
    # relevant columns
    up_int <- df_small$up_int
    low_int <- df_small$low_int
    effect_size <- df_small$mean_dif
    column_name <- df_small$column_name[1]
    mean_dif <- df_small$mean_dif
    plot <- ggplot(df_small, aes(x = 1, y = effect_size)) +
      geom_point(x = 1, y = mean_dif[1], shape = 21, colour = "black", fill = "white", size = 5, stroke = 5) +
      geom_errorbar(aes(ymin = low_int, ymax = up_int), width = 0.2) +
      geom_hline(yintercept = effect_size[1], linetype = "dotted") +
      ggtitle(column_name)
    

    plot_list[[i]] <- plot
    
  }
  
  for (i in seq_along(plot_list)) {
    ggsave(paste0(path, "/plot_", i, ".png"), plot_list[[i]], width = 8, height = 6)
  }
}

##choose path and execute
two_row_ci_plot(summary_df,'D:/MiguelW12/Documents/stats_script_try/gsd1a_stats/plot_try_CI_72_ee' )

## add statistical parameters to summary data frame
#####power
add_power_to_df <- function(df) {
  
  # empty vector for power vals
  power_vals <- rep(NA, nrow(df))
  
  for (i in 1:nrow(df)) {
    d <- df$calculated_effect_size[i]
    sig.level <- 0.1
    power <- pwr.t2n.test(n1 = 65, n2 = 92, d = d, sig.level = sig.level)$power
    power_vals[i] <- power
  }
  
  df$power <- power_vals
  
  return(df)
}

##execute
p <- add_power_to_df(summary_df)
write.csv(as.data.frame(p), file="gsd1a_frame_72_ee-summary_withhpower.csv")

#####
###linear_reg with plots 
compareGroups <- function(data, plot_path, summary_path) {
  

  feature_cols <- names(data)[-1]
  

  model_list <- list()
  

  summary_df <- data.frame()
  
  for (feature in feature_cols) {
    model <- lm(data[,1] ~ data[,feature])
    model_name <- paste0("model_", feature)
    model_list[[model_name]] <- model
    plot_file <- paste0(plot_path, "/", model_name, ".png")
    png(plot_file)
    plot(model)
    dev.off()
    summary_row <- data.frame(Model = model_name,
                              R_squared = summary(model)$r.squared,
                              Adjusted_R_squared = summary(model)$adj.r.squared)
    summary_df <- rbind(summary_df, summary_row)
  }
  
  for (i in 2:length(feature_cols)) {
    combo_list <- combn(feature_cols, i, simplify = FALSE)
    for (combo in combo_list) {
      combo_string <- paste0(combo, collapse = "+")
      model <- lm(data[,1] ~ eval(parse(text = combo_string)), data = data)
      model_name <- paste0("model_", combo_string)
      model_list[[model_name]] <- model
      plot_file <- paste0(plot_path, "/", model_name, ".png")
      png(plot_file)
      plot(model)
      dev.off()
      summary_row <- data.frame(Model = model_name,
                                R_squared = summary(model)$r.squared,
                                Adjusted_R_squared = summary(model)$adj.r.squared)
      summary_df <- rbind(summary_df, summary_row)
    }
  }
  
  write.csv(summary_df, file = summary_path, row.names = FALSE)
  return(model_list)
}


##data 
## convert index to numbers
t_data[t_data=="GSD1ANONE72"] <- 2
t_data[t_data=="HCNONE72"] <- 1
## execute function & add path for plots, csv
model_list<- compareGroups(t_data,'D:/MiguelW12/Documents/stats_script/gsd1a_stats/lm_72/plots', 
                                      'D:/MiguelW12/Documents/stats_script/gsd1a_stats/lm_72/res.csv')
#### AIC
aic_model <- aictab(cand.set = model_list)
write.csv(as.data.frame(aic_model), file="AIC-gsd1a.csv")


######
### bayes upper bound- for stats frame 

add_bub_column <- function(df) {
  df$B.U.B <- (-1) / ((2.718 * df$p.value) * log(df$p.value))
  return(df)
}

####read_data 
frame = read.csv('D:/MiguelW12/Documents/stats_script/gsd1a_stats/gsd1a_stat_frames/res_stats_folowproc.csv', header=TRUE)
data_with_bub <- add_bub_column(frame)
head(data_with_bub)
write.csv(as.data.frame(data_with_bub), file="data_with_bub-gsd1a.csv")

#### Bayes factor

####split frames
split_dataframe <- function(df, save_files = FALSE) {

  index_col <- names(df)[1]
  

  df_col_names <- names(df)[-1]
  

  split_df_list <- list()
  

  for (col_name in df_col_names) {
    new_df <- df[c(index_col, col_name)]
    
    split_df_list[[col_name]] <- new_df
    
    if (save_files) {
      file_name <- paste0(col_name, ".csv")
      write.csv(new_df, file_name, row.names = FALSE)
    }
  }
  
  return(split_df_list)
}
split_dfs <- split_dataframe(t_data, save_files = TRUE)
split_dfs

# function to perform t-test using BayesFactor package
#ADJUST NAMES OF GROUPS
ttest_bf <- function(df) {
  group1 <- df[[2]][df[[1]] == "GSD1ANONE72"]
  group2 <- df[[2]][df[[1]] == "HCNONE72"]
  
  result <- ttestBF(x = group1, y = group2, rscale="medium")
  
  return(result)
}

save_ttest_results <- function(df_list, path) {
  for (df_name in names(df_list)) {
    res <- ttest_bf(df_list[[df_name]])
    res_df <- as.data.frame(res)
    file_path <- paste0(path, "/", df_name, "_ttest_results .csv")
    write.csv(res_df, file = file_path, row.names = FALSE)
  }
}

names(split_dfs)

# replace with names
names(split_dfs) <- c("calc_area1", "calc_area2", "calc_area3", "calc_intensity", "calc_text1",     
                      "lyso_area1", "lyso_area2", "lyso_area3", "lyso_intensity1", "lyso_text1", "lyso_text2",     
                      "nuc_area1", "nuc_area2", "nuc_intensity", "nuc_text1",      
                      "tmre_area1", "tmre_area2", "tmre_area3", "tmre_intensity1", "tmre_text1", "tmre_text2" ) 
my_path <- "D:/MiguelW12/Documents/stats_script/gsd1a_new_stats/BF_72/" 
save_ttest_results(split_dfs, my_path)

#####
##dabaset - For IF as well
split_df = read.csv('D:/MiguelW12/Documents/calc_area1.csv', header=TRUE)
two.group.unpaired <- 
  try_split_df %>%
  dabest(index, calc_area1,  
         idx = c("HCNONE", "GSD1ANONE"), 
         paired = FALSE)

two.group.unpaired
two.group.unpaired.meandiff <- mean_diff(two.group.unpaired)
two.group.unpaired.meandiff

## plot
plot(two.group.unpaired.meandiff, 
     color.column = group,
     rawplot.markersize = 3,
     rawplot.groupwidth = 0.6,
     rawplot.ylabel = "feature x",
     effsize.ylabel = "mean difference",
     axes.title.fontsize = 16,
     palette = c("dodgerblue3", "firebrick3"),
     theme = ggplot2::theme_gray()
     
)

### IF
pwr.t2n.test(n1 = , n2=35 , d =0.6, sig.level =0.1, power=0.8)


two.group.unpaired <- 
  frame_try %>%
  dabest(group, measuerment,  
         idx = c("HC", "GSD1A"), 
         paired = FALSE)


  
two.group.unpaired.meandiff <- mean_diff(two.group.unpaired)
two.group.unpaired.meandiff


## PLOT
plot(two.group.unpaired.meandiff, 
     color.column = group,
     rawplot.markersize = 3,
     rawplot.groupwidth = 0.6,
     rawplot.ylabel = "marker x",
     effsize.ylabel = "mean difference",
     axes.title.fontsize = 16,
     palette = c("dodgerblue3", "firebrick3"),
     theme = ggplot2::theme_gray()
     
)

#TREATMENT
two.group.unpaired <- 
  frame_try %>%
  dabest(group, measuerment,  
         idx = c("HC", "GSD1A", "GSD1A+GHF201"), 
         paired = FALSE)


two.group.unpaired

two.group.unpaired.meandiff <- mean_diff(two.group.unpaired)

two.group.unpaired.meandiff

## PLOT
plot(two.group.unpaired.meandiff, 
     color.column = group,
     rawplot.markersize = 3,
     rawplot.groupwidth = 0.6,
     rawplot.ylabel = "SIRT1 Intensity",
     effsize.ylabel = "Mean Difference",
     axes.title.fontsize = 16,
     palette = c("forestgreen", "firebrick3", "dodgerblue3"),
     theme = ggplot2::theme_classic()
     
)

#####################################################################