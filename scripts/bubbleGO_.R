###############################################
#             Requires                        #
###############################################
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(Cairo)) install.packages("Cairo")
if(!require(tidyverse)) install.packages("tidyverse")

###############################################
#             Includes                        #
###############################################
library(ggplot2)
library(Cairo)
library(tidyverse)

###############################################
# Begin Bubble plots                          #
###############################################

# set base directory  ### modify accordingly ###
basedir <- "c:/leaf3_dev/comp"
setwd(basedir)

# BP Up
df.bp<-read_tsv("DE/GO/bubble_plot_selected_GOs_UP.txt")

# assign a representation factor of 0.01 to those terms not present in the data
df.bp[which(is.na(df.bp$p.value)),"rf"] <- 0.01 
# assign NA to non-significant terms
df.bp[which(df.bp$p.value > 0.05),"p.value"] <- NA

# obtain the p.values and assign a score
score <- df.bp$p.value
score[is.na(df.bp$p.value)]<-NA
score[df.bp$p.value < 0.05]<-1
score[df.bp$p.value <= 0.01]<-2
score[df.bp$p.value <= 0.001]<-3
score[df.bp$p.value <= 0.0001]<-4

go.bp.up<-ggplot(df.bp, aes(x = x_var, y = y_var, labels = Term, size =RF, col = score)) +
  geom_point(alpha = 0.7) +
  scale_colour_gradient(limits= c(1, 4), low = "#ffc000",  high = "red", na.value = "light gray") +     # pink "#f6e7e4"
  scale_size_continuous(range = c(2,20)) +
  scale_y_continuous("GO Terms (Biological Process)", breaks = df.bp$y_var, labels = df.bp$Term) +
  scale_x_discrete("D. A. S.") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey"), text=element_text(size=12))

go.bp.up

# export the plot to a png file
png(filename = "DE/GO/bubble_plot_selected_BP_UP.png",
    width = 600 * 7.5,        
    height = 600 * 6,
    units = "px",
    res = 600)  
plot(go.bp.up)
dev.off()
rm(df.bp, score, go.bp.up)

# BP Down
df.bp<-read_tsv("DE/GO/bubble_plot_selected_GOs_DOWN.txt")

# assign a representation factor of 0.01 to those terms not present in the data
df.bp[which(is.na(df.bp$p.value)),"rf"] <- 0.01 
# assign NA to non-significant terms
df.bp[which(df.bp$p.value > 0.05),"p.value"] <- NA

# obtain the p.values and assign a score
score <- df.bp$p.value
score[is.na(df.bp$p.value)]<-NA
score[df.bp$p.value < 0.05]<-1
score[df.bp$p.value <= 0.01]<-2
score[df.bp$p.value <= 0.001]<-3
score[df.bp$p.value <= 0.0001]<-4

go.bp.down<-ggplot(df.bp, aes(x = x_var, y = y_var, labels = Term, size =RF, col = score)) +
  geom_point(alpha = 0.7) +
  scale_colour_gradient(limits= c(1, 4), low = "#ffc000",  high = "red", na.value = "light gray") +     # pink "#f6e7e4"
  scale_size_continuous(range = c(2,20)) +
  scale_y_continuous("GO Terms (Biological Process)", breaks = df.bp$y_var, labels = df.bp$Term) +
  scale_x_discrete("D. A. S.") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey"), text=element_text(size=12))

go.bp.down

# export the plot to a png file
png(filename = "DE/GO/bubble_plot_selected_BP_DOWN.png",
    width = 600 * 7.5,        
    height = 600 * 6,
    units = "px",
    res = 600)  
plot(go.bp.down)
dev.off()
rm(df.bp, score, go.bp.down)

#######################
# End of bubble plots #
#######################
