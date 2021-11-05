##Install the Hudson Plot Package
library(devtools)
devtools::install_github('anastasia-lucas/hudson')
library(hudson)

disc_new_group <- disc_new_group[-c(1)]
rep_new_group <- rep_new_group[-c(1)]

vector_var <- c('alpha-tocopherol (Vitamin E) - Triglycerides', 'gamma-tocopherol (Vitamin E) - Triglycerides', 'Retinol (Vitamin A) - Triglycerides', 'Retinol (Vitamin A) - Uric acid', 'Retinol (Vitamin A) - Blood urea nitrogen')

#Generate a plot and highlight by p-value threshold
library(hudson)
data(disc_new_group)
data(rep_new_group)
emirror2(top=disc_white_plot, bottom=rep_white_plot, tline = c(0.05/nrow(disc_white_plot), 0.01052426), bline = c(0.05/nrow(rep_white_plot), 0.05013249), annotate_var = vector_var, highlight_var = vector_var, highlighter="green", color1 ="#AAAAAA", color2 = "#4D4D4D", freey = FALSE, rotatelabels = TRUE, labelangle = 90, 
        toptitle = "Discovery", bottomtitle = "Replication", hgt = 10, res = 300, file = 'white_eman_plot')


19.3/17930
[1] 0.001076408
> 0.05/nrow(rep_mexican_plot)
[1] 0.0001157407
> 43.2/17930
[1] 0.00240937
> 19.3/432
[1] 0.0446759