library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

suppvec_hist <- function(df, caller, outfile) {
  df$SUPP_VEC_STRING = str_pad(as.character(df$SUPP_VEC), 3, "left", "0")
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts$TYPE = "INS"
  ggplot(df, aes(x = SUPP_VEC_STRING, y = 1, fill = TYPE)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste("SVs by Support Vector (", caller, ")")) +
    xlab("Support Vector") +
    ylab("Count") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 16),
    ) +
    scale_fill_discrete(name = "Type") +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75)
  
    ggsave(outfile, width= 6, height = 8)
}

plot_table <- function(fn, outdir, prefix, caller, discvec){

  df <- read.table(fn, sep = "\t", header = TRUE)

  df$LenCategory = "0"
  df$LenCategory = ifelse(df$LEN > 0 & df$LEN <= 50, "INS_1-50", df$LenCategory)
  df$LenCategory = ifelse(df$LEN < 0 & df$LEN >= -50, "DEL_1-50", df$LenCategory)
  df$LenCategory = ifelse(df$LEN > 50 & df$LEN <= 200, "INS_50-200", df$LenCategory)
  df$LenCategory = ifelse(df$LEN < -50 & df$LEN >= -200, "DEL_50-200", df$LenCategory)
  df$LenCategory = ifelse(df$LEN > 200 & df$LEN <= 500, "INS_200-500", df$LenCategory)
  df$LenCategory = ifelse(df$LEN < -200 & df$LEN >= -500, "DEL_200-500", df$LenCategory)
  df$LenCategory = ifelse(df$LEN > 500 & df$LEN <= 1000, "INS_500-1000", df$LenCategory)
  df$LenCategory = ifelse(df$LEN < -500 & df$LEN >= -1000, "DEL_500-1000", df$LenCategory)
  df$LenCategory = ifelse(df$LEN > 1000 & df$LEN <= 5000, "INS_1000-5000", df$LenCategory)
  df$LenCategory = ifelse(df$LEN < -1000 & df$LEN >= -5000, "DEL_1000-5000", df$LenCategory)
  df$LenCategory = ifelse(df$LEN > 5000, "INS_5000+", df$LenCategory)
  df$LenCategory = ifelse(df$LEN < -5000, "DEL_5000+", df$LenCategory)
  df$LenCategory <- factor(df$LenCategory,levels = c("DEL_5000+", "DEL_1000-5000", "DEL_500-1000", "DEL_200-500", "DEL_50-200", "DEL_1-50",
                                                     "INS_1-50", "INS_50-200", "INS_200-500", "INS_500-1000", "INS_1000-5000", "INS_5000+" ))
  
  df %>% group_by(df$LenCategory) %>% tally()
  
  labellist = c("5k+", "1k-5k", "500-1k", "200-500", "50-200", "1-50", 
                "1-50", "50-200", "200-500", "500-1k", "1k-5k", "5k+")
  
  counts <- df %>% filter(df$TYPE == "INS" | df$TYPE == "DEL") %>% group_by(LenCategory) %>% count()
  counts$TYPE = "INS"
  
  ggplot(df %>% filter(TYPE == "INS" | TYPE == "DEL"), aes(x = LenCategory, fill = TYPE, y = 1)) +
         geom_bar(position = "stack", stat = "identity")+
         scale_x_discrete(labels=labellist) +
         xlab("Length") +
         ylab("Count") +
         labs(title = paste("Indels by Length (", caller, ")")) +
         theme(plot.title = element_text(size = 18, hjust = 0.5),
               axis.text.x = element_text(size = 12),
               axis.text.y = element_text(size = 12),
               axis.title.x = element_text(size = 16),
               axis.title.y = element_text(size = 16),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 16),
         ) +
         scale_fill_discrete(name = "Type") +
         geom_text(data = counts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75)
  outfile <- paste(outdir, prefix, "_indellenhist_", caller, ".png", sep="")
  ggsave(outfile, width= 12, height = 8)
  
  ggplot(df %>% count(TYPE), aes(x = reorder(TYPE, -n), y = n)) +
      geom_bar(stat = "identity", fill = "lightblue2") +
      labs(title = paste("SVs by Type (", caller, ")")) +
      xlab("Type") +
      ylab("Count") +
      theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 16),
      ) +
    geom_text(aes(x = TYPE, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75)
  
  outfile <- paste(outdir, prefix, "_typehist_", caller, ".png", sep="")
  ggsave(outfile, width= 6, height = 8)

  suppvec_hist(df, caller, paste(outdir, prefix, "_suppvechist_", caller, ".png", sep=""))
  suppvec_hist(df %>% filter(SPECIFIC_FLAG == 1), paste(caller, "Specific"), paste(outdir, prefix, "_specificsuppvechist_", caller, ".png", sep=""))
  suppvec_hist(df %>% filter(SPECIFIC_FLAG == 1 & PRECISE_FLAG == 1), paste(caller, "Precise+Specific"), paste(outdir, prefix, "_precisespecificsuppvechist_", caller, ".png", sep=""))
  

  df$startspan = df$MAX_START - df$MIN_START
  df$endspan = df$MAX_END - df$MIN_END
  
  multi <- df %>% filter(NUMVARS > 1)
  
  allspans = c(multi$startspan, multi$endspan)
  allspans = round(log(allspans + 1, base = 2))
  ggplot() + aes(allspans) + geom_histogram(binwidth = 1, fill = "lightblue2", color = "darkblue") + 
    xlab(expression(Breakpoint~range:~log[2](max~-~min~+~1))) + ylab("Number of merged breakpoints") +
    labs(title = paste("Breakpoint Spans (", caller, ")")) +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
    ) +     stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1) 

  outfile <- paste(outdir, prefix, "_breakpointspanhist_", caller, ".png", sep="")
  ggsave(outfile, width= 8, height = 8)
  
  filtered = df %>% filter(SPECIFIC_FLAG == 1 & PRECISE_FLAG == 1 & SUPP_VEC == discvec)
  filtered$caller <- caller
  
  totalcount <- sum(df$NUMVARS)
  mergedcount <- nrow(df)
  #return(c(filtered, totalcount, mergedcount))
  return (list(filtered, totalcount, mergedcount))
}

plot_all_callers <- function(outdir, prefix, discvec){
  jasmineresults <- plot_table(paste(outdir, prefix, ".jasmine_augmented.txt", sep = ""), outdir, prefix, "Jasmine", discvec)
  survivorresults <- plot_table(paste(outdir, prefix, ".survivor_augmented.txt", sep = ""), outdir, prefix, "Survivor", discvec)
  svtoolsresults <- plot_table(paste(outdir, prefix, ".svtools_augmented.txt", sep = ""), outdir, prefix, "svtools", discvec)
  svimmerresults <- plot_table(paste(outdir, prefix, ".svimmer_augmented.txt", sep = ""), outdir, prefix, "svimmer", discvec)
  
  ncol(jasmineresults[0])
  
  discordant <- data.frame()
  discordant <- rbind(discordant, jasmineresults[[1]])
  discordant <- rbind(discordant, survivorresults[[1]])
  discordant <- rbind(discordant, svtoolsresults[[1]])
  discordant <- rbind(discordant, svimmerresults[[1]])
  
  totals <- c(jasmineresults[[2]], survivorresults[[2]], svtoolsresults[[2]], svimmerresults[[2]],
              jasmineresults[[3]], survivorresults[[3]], svtoolsresults[[3]], svimmerresults[[3]])

  disccounts <- discordant %>% group_by(caller) %>% summarise(counts=n(), sums=sum(NUMVARS))
  #disccounts
  totals
}
  

#fn = "/home/mkirsche/eclipse-workspace/SvPopulationAnalysis/augment.txt"
#fn = "/home/mkirsche/eclipse-workspace/SvPopulationAnalysis/survaugment.txt"
#fn = "/home/mkirsche/eclipse-workspace/SvPopulationAnalysis/svtoolsaugment.txt"


outdir = "/home/mkirsche/eclipse-workspace/SvPopulationAnalysis/"

plot_all_callers(outdir, "pur_ccs", 100)
plot_all_callers(outdir, "md50", 001)



