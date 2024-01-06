library(arrow)
library(dplyr)
library(ggplot2)


###############################################################################
## Author: Andrea Riebler [andrea *.* riebler *a*t* math *.* ntnu *. no]
##         
## Time-stamp: <2009> revised <04/2014>
##
## Description:
## Compute nonrandomized PIT values for count data according to 
## 
## Section 2.1 in: 
## Czado, C., Gneiting, T. and Held, L. (2009), 
## Predictive Model Assessment for Count Data. 
## Biometrics, 65: 1254–1261. doi: 10.1111/j.1541-0420.2009.01191.x
###############################################################################

#####################################################
## non-randomized version of the PIT histogram
##
## Params: 
## u - real number 0<=u<=1 
## Px - P(X <= x)
## Pxm1 - P(X <= x-1)
####################################################

pit.one <- function(u,x,Px,Pxm1){

   F_u <- ifelse(u <= Pxm1 , 0, pmin(1,(u-Pxm1)/(Px-Pxm1) ) )
   F_u[x==0] <- pmin(1,u/Px)[x==0] # needless!?
   if(u == 1){
    F_u <- 1
   }
   if(u == 0){
    F_u <- 0
   }
   #print(F_u)
   return(mean(F_u))
}

####################################################
## Params: 
## J - number of bins
#####################################################
pit <- function(J=10, x, Px, Pxm1){
  F_u.bar <- sapply((0:J)/J,pit.one, x=x, Px=Px, Pxm1=Pxm1)
  f_j <- J*diff(F_u.bar)
  
  erg <- list(breaks=(0:J)/J,counts=f_j, density=f_j,mids=(0:(J-1))/J+diff((0:J)/J)/2,
             xname="PIT",equidist=TRUE)
  class(erg) <- "histogram"
  return(erg)
}

#######################################################
## EXAMPLES
#######################################################

## test the function for binomial data
if(TRUE){
    n <- 1000
    prob <- runif(n)
    size <- rpois(n,100)
    x <- rbinom(n, prob=prob,size=size)
    Px <- pbinom(x, prob=prob,size=size)
    Pxm1 <- pbinom(x-1, prob=prob,size=size)
    
    plot(pit(J=20,x=x, Px=Px,Pxm1=Pxm1),ylab="Relative frequency", ylim=c(0,2),main="")
}

## From INLA output
if(TRUE){
    library(INLA)
    data(Germany)
    g = system.file("demodata/germany.graph", package="INLA")
    source(system.file("demodata/Bym-map.R", package="INLA"))
    summary(Germany)

    ## just make a duplicated column
    Germany = cbind(Germany,region.struct=Germany$region)

    # standard BYM model (without covariates)
    formula = Y ~ f(region.struct,model="besag",graph=g) +
                f(region,model="iid")


    result  =  inla(formula,family="poisson",data=Germany,E=E,
        control.compute=list(cpo=TRUE))



    ## we have three failures => recompute those manually
    result  = inla.cpo(result)
    vpit <- result$cpo$pit
    ## compute Pxm1 
    vcpo <- result$cpo$cpo
    Pxm1help <- vpit - vcpo
    ## be sure to avoid negative PITs
    Pxm1 <- ifelse(Pxm1help<0,0,Pxm1help)
    plot(pit(J=20, x=Germany$Y, Px=vpit, Pxm1=Pxm1), ylim=c(0,1.5), ylab="Relative frequency", xlab="PIT")
}

################################################################################
################################################################################



# Define the list of model identifiers
model_original_identifier <- c('baseline_f_m', 'baseline_f', 'model13_nb_feature_set1', 'model3_zinb_feature_set1', 
                               'model1_zinb_feature_set1', 'model4_zinb_feature_set1', 'model15_zinb_feature_set1', 
                               'model19_zinb_feature_set3', 'model23_zinb_feature_set4')

# Define identifiers as a list of expressions
scientific_model_identifiers <- list(
  expression(B[HMV]), 
  expression(B[HV]), 
  expression(M[1]), 
  expression(M[2]), 
  expression(M[3]), 
  expression(M[4]), 
  expression(M[5]), 
  expression(M[6]), 
  expression(M[7])
)
scientific_model_identifiers_str <- c('B[HMV]', 'B[HV]', 'M[1]', 'M[2]', 'M[3]', 'M[4]', 'M[5]', 'M[6]', 'M[7]')


# Read the parquet files into a list
pps_list <- lapply(model_original_identifier, function(model_identifier) {
  file_path <- paste0('C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/', model_identifier, '_posterior_predicitve_samples_wt.parquet')
  read_parquet(file_path)
})
# pps_list = list(pps_list[[1]])
# Loop over each data frame in the list and their corresponding identifiers
for (i in 1:length(pps_list)) {
  pps <- pps_list[[i]]
  model_identifier <- model_original_identifier[i]
  
  # Specify the file path for the PNG file
  file_path <- paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Plots/plots_for_results_section/pit_histograms/nonrandom_pit_histogram_", scientific_model_identifiers_str[i], ".png")
  
  # Set the dimensions and resolution of the output image
  png(file_path, width = 2000, height = 1200, res = 300, units = "px")
  
  # Extract the 'ged_sb' column as a numeric vector
  x <- as.numeric(pps$ged_sb)
  
  # Extract only the columns that start with "draw_"
  draw_columns <- grep("^draw_", names(pps), value = TRUE)
  
  # Initialize vectors to store the probabilities
  P_X_leq_x <- numeric(nrow(pps))
  P_X_leq_x_minus_1 <- numeric(nrow(pps))
  
  # Loop through each row to calculate probabilities
  for (j in 1:nrow(pps)) {
    row <- pps[j, ]
    x_j <- as.numeric(row[["ged_sb"]])  # Extract as a single numeric value
    draws <- as.numeric(row[draw_columns])
    draws <- na.omit(draws)  # This will remove NA values
    
    P_X_leq_x[j] <- mean(draws <= x_j)
    P_X_leq_x_minus_1[j] <- mean(draws <= (x_j - 1))
  }
  
  # Plot using the pit function
  # Assuming 'pit' is a function you have defined or loaded from a library
  hist_data <- pit(J=20, x, P_X_leq_x, P_X_leq_x_minus_1)
  
  plot(hist_data, main=scientific_model_identifiers[[i]] , xlab="PIT", ylab="Relative Frequency", col=rgb(161/255, 203/255, 226/255), border="black", cex.lab=1.25, cex.axis=1.25, cex.main=1.75, cex.sub=1.25)
  abline(h=1.0, lty=2, col="black")
  
  # Create a dataframe from the histogram object
  # ggplot_hist_data <- data.frame(mids = hist_data$mids, density = hist_data$density)
  
  # ggplot_hist_data <- as.data.frame(ggplot_hist_data)
  # ggplot(ggplot_hist_data, aes(x=mids, y=density)) +
  #  geom_bar(stat='identity', fill=rgb(161/255, 203/255, 226/255), colour="black") +
  #  theme_minimal() +
  #  theme(axis.text.x=element_text(size=9),
  #        axis.text.y=element_text(size=9),
  #        axis.title.x=element_text(size=11, face="bold"),
  #        axis.title.y=element_text(size=10, face="bold")) +
  #  ggtitle(scientific_model_identifiers[[i]]) +
  #  xlab("PIT") +
  #  ylab("Relative Frequency") +
  #  geom_hline(yintercept=1.0, linetype="dashed", color = "black")
  
  # Close the PNG device
  dev.off()
}


################################# EXPERIMENTAL ################################


pps <- arrow::read_parquet("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/model15_zinb_feature_set1_posterior_predicitve_samples_wt.parquet")
x <- as.numeric(pps$ged_sb)
# Extract only the columns that start with "draw_"
draw_columns <- grep("^draw_", names(pps), value = TRUE)

# Initialize vectors to store the probabilities
P_X_leq_x <- numeric(nrow(pps))
P_X_leq_x_minus_1 <- numeric(nrow(pps))

# Loop through each row to calculate probabilities
for (i in 1:nrow(pps)) {
  row <- pps[i, ]
  x_j <- row[["ged_sb"]]
  draws <- as.numeric(row[draw_columns])
  draws <- na.omit(draws)  # This will remove NA values
  
  P_X_leq_x[i] <- mean(draws <= x_j)
  P_X_leq_x_minus_1[i] <- mean(draws <= (x_j - 1))
}
# P_X_leq_x_df <- as.data.frame(P_X_leq_x)
# P_X_leq_x_minus_1_df <- as.data.frame(P_X_leq_x_minus_1)

plot(pit(J=20, x, P_X_leq_x, P_X_leq_x_minus_1))

# Load the ggplot2 package

# Example usage of the pit function
# Assuming x, Px, Pxm1 are defined
hist_data <- pit(J=20, x, P_X_leq_x, P_X_leq_x_minus_1)

plot(hist_data, main=expression(M[4]), xlab="PIT", ylab="Frequency", col=rgb(161/255, 203/255, 226/255), border="black", cex.lab=1.25, cex.axis=1.25, cex.main=1.75, cex.sub=1.25)
abline(h=0.5, lty=2, col="black")
# Convert the histogram object to a dataframe for ggplot
plot(hist_data, main=expression(M[4]), xlab="PIT", ylab="Frequency", col="lightblue", border="black")

hist_df <- data.frame(mids=hist_data$mids, counts=hist_data$counts)

# Create the plot using ggplot2
pit_plot <-ggplot(hist_df, aes(x=mids, y=density)) +
  geom_bar(stat='identity', fill=rgb(161/255, 203/255, 226/255), colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11, face="bold"),
        axis.title.y=element_text(size=10, face="bold")) +
  ggtitle("PIT Histogram")+
  xlab("PIT") +
  ylab("Relative Frequency") +
  geom_hline(yintercept=1.0, linetype="dashed", color = "black")

# pit_plot <-ggplot(hist_df, aes(x=mids, y=counts)) +
#   geom_bar(stat="identity", fill="lightblue", color="black") +
#   theme_minimal() +
#   labs(x="PIT", y="Frequency") +
#   ggtitle("PIT Histogram")

print(pit_plot)
##################################################################################
file_name <- paste0("pit_histogram_", model_original_identifier[3], ".png")
file_path <- file.path("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Plots/plots_for_results_section/pit_histograms", file_name)

# Set the dimensions and resolution of the output image
png(file_path, width = 3000, height = 1800, res = 300, units = "px")

# Plot the histogram using the 'hist' function or your custom 'pit' plotting
hist_data <- pit(J=20, x, Px, Pxm1) # You need to define x, Px, Pxm1 before this line
plot(hist_data, main=expression(B[HV]), xlab="PIT", ylab="Frequency", col="lightblue", border="black")

# Add more customizations to the plot here if needed

# Close the PNG device to save the file
dev.off()
###################################################################################
i = 4200
row <- pps[i, ]
x_j <- row[["ged_sb"]]
draws <- as.numeric(row[draw_columns])

P_X_leq_x_i <- mean(draws <= x_j)
P_X_leq_x_minus_1_i <- mean(draws <= (x_j - 1))
