## ---------------------------------------------------------------------------
## gibbs_hierarchical.R -- Gibbs sampler for a hierarchical multivariate normal
## model of RPPA protein expression across four FAB subtypes of AML.
## Run from the repository root:  source("R/gibbs_hierarchical.R")
## Reads: data/leukemia.csv
## ---------------------------------------------------------------------------

library("MASS")
library("stats")
library(mvtnorm)
data <- read.csv("data/leukemia.csv", sep = ";")
M0 <- data[data$FAB == "M0", -1]
M1 <- data[data$FAB == "M1", -1]
M2 <- data[data$FAB == "M2", -1]
M4 <- data[data$FAB == "M4", -1]

M0 <- cbind(Group = "M0", M0)
M1 <- cbind(Group = "M1", M1)
M2 <- cbind(Group = "M2", M2)
M4 <- cbind(Group = "M4", M4)

y <- rbind(M0, M1, M2, M4)
y = as.matrix(y[, -1])


# Define parameters
p <- ncol(y) 
n <- nrow(y)
D <- 4
omega <- solve(cov(y))
T0 <- diag(100, p)
theta0 <- rep(0, p)
U <- diag(0.01, p)
a <- p
b <- p
W <- diag(0.01, p)
m0 <- rep(0, p)
k0 <- diag(1,p) 
S <- 2000

#initializing thetaj as the mean in group j
thetaj <- list()

# Calcola la media per ogni gruppo
for (group in list(M0, M1, M2, M4)) {
  group <- group[, -which(names(group) == "Group")]
  media_group <- colMeans(group)
  thetaj[[length(thetaj) + 1]] <- media_group
}

thetaj <- lapply(thetaj, function(x) as.numeric(x))



Gibbs <- function(n, omega, T0, y_bar, theta0, thetaj, U, W, b, a, m0, k0, y, D, S) {
  
  set.seed(1234)
  # Initialize the output list, inizializza la lista di uoutput in maniera diversa
  
  omega_post = array(NA, c(p, p, S)) #S matrici 18x18 (array:matrice tridimensionale)
  thetaj_post = array(NA,c(S,p,D))  #D matrici Sx18
  theta0_post = matrix(NA,S,p)
  T0_post = array(NA, c(p, p, S))
  
  for (s in 1:S) {
    
    # Updating of omega
    an <- a + n
    Sn = (t(y)-c(theta0))%*%t(t(y)-c(theta0))
    Un <- U + Sn
    omega <- rWishart(1, an, solve(Un))[,,1] 
    
    # Updating of theta0
    thetaJ_bar = Reduce("+", thetaj) / length(thetaj)
    A0 <- solve(k0)
    b0 <- A0 %*% m0
    mp <- solve(A0 + D * T0) %*% (b0 + D * T0 %*% matrix(thetaJ_bar, ncol = 1))
    kp <- solve(A0 + D * T0)
    theta0 <- rmvnorm(1, mp, kp)
    
    # Updating of T0
    bn <- b + D
    thetaj_vec <- as.numeric(thetaj[[1]])
    theta0_vec <- as.numeric(theta0[1,])
    differenza <- thetaj_vec - theta0_vec
    differenza_matrix <- matrix(differenza, ncol = 1)
    S2 <- differenza_matrix %*% t(differenza_matrix)  
    epsilon <- 1e-6  
    S2_reg <- S2 + epsilon * diag(nrow(S2))
    Wn <- solve(S2_reg) + W
    T0 <- rWishart(1, bn, solve(Wn))[,,1]
    
    #Updating of Thetaj
    for (d in 1:D) {
      An1 <- solve(n * omega + T0)
      bn <- n * omega %*% thetaj[[d]] + T0 %*% matrix(theta0, ncol = 1)
      thetaJ_post <- rmvnorm(1, An1 %*% bn, An1)
      thetaJ_post = matrix(c(thetaJ_post), nrow=1,byrow = T)
      
    #thetaj: media delle p variabili nel gruppo d-esimo (con n=numero di osservazioni nel gruppo d-esimo)
    
    
    #Output storage
    thetaj_post[s,,d]  = thetaJ_post 
    }
    omega_post[,,s] = omega
    theta0_post[s,] = theta0
    T0_post[,,s] = T0
    
    
  }  
  
  return(posterior = list(thetaj_post = thetaj_post,
                          Omega_post = omega_post,
                          theta0_post = theta0_post,
                          T0_post = T0_post))
  
}

# Run the Gibbs sampler
results <- Gibbs(n=n, omega=omega, T0=T0, y_bar=y_bar, theta0=theta0,thetaj=thetaj, U=U, W=W, b=b, a=a, m0=m0, k0=k0, y=y, D=D, S=S)


#Inspect the results
h <- hist(results$theta0_post, plot=FALSE)

hist(results$theta0_post, main="Theta0_post", xlab="Values", ylab="Frequency", col="skyblue")

mids <- h$mids
counts <- h$counts

lines(smooth.spline(mids, counts), col="red", lwd=2)


hist(results$T0_post, xlim = c(-200,300))

hist(results$Omega_post)





#Confronto tra medie e box plot di confronto tra 4 proteine (AKT,BAD,BAX,MYC)
#
# nel nostro caso:
# theta_group1 = thetaJ_post[,,1]
# mi prendo la distribuzione della colonna 4 per 2000 iterazion (proteina BAD) #results$theta_post[,,1][,4]
#



#Qual è la probabilità che una persona scelta a random per il gruppo M0 abbia
#una concentrazione della proteina BAD nel sangue maggiore di uno del gruppo M2?

mean(results$thetaj_post[,,1][,4]>results$thetaj_post[,,2][,4])

# Ce lo aspettiamo questo valore? Si pèrchè c'è poca differenza tra i due  Valori molto simili


#Shrinkage for protein BAD, by taking the mean of the columns of thetaj_post[,,1]

exp_post=colMeans(results$thetaj_post[,,1])
plot(thetaj[[1]],exp_post, lwd = 2, xlab = "ybar for group M0", ylab = "ThetaJ Post", main = "Shrinkage for M0")



#Boxplots per il confronto tra medie
######AKT########
# Install and load the necessary packages
library(ggplot2)
library(patchwork)

# Create the data frame combining all vectors and their labels
theta_vecd1 <- c(results$thetaj_post[,,1][,1])
theta_vecd2 <- c(results$thetaj_post[,,2][,1])
theta_vecd3 <- c(results$thetaj_post[,,3][,1])
theta_vecd4 <- c(results$thetaj_post[,,4][,1])

data <- data.frame(
  value = c(theta_vecd1, theta_vecd2, theta_vecd3, theta_vecd4),
  group = factor(rep(c("AKT in M0", "AKT in M1", "AKT in M2", "AKT in M4"), 
                     times = c(length(theta_vecd1), length(theta_vecd2), length(theta_vecd3), length(theta_vecd4))))
)

# Create the boxplot and zoom it to see the mean differences
main_plot_limited <- ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("pink", "skyblue", "lightgoldenrod", "orange")) +
  coord_cartesian(ylim = c(-3, 3)) +  # Limit the y-axis to the range -5 to 5
  labs(title = "Boxplots of Theta Values (Limited Y-Axis on -3,3)", x = "Model", y = "Theta Values") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13))
main_plot_limited

#######BAD########
theta_vecd1 <- c(results$thetaj_post[,,1][,4])
theta_vecd2 <- c(results$thetaj_post[,,2][,4])
theta_vecd3 <- c(results$thetaj_post[,,3][,4])
theta_vecd4 <- c(results$thetaj_post[,,4][,4])

data <- data.frame(
  value = c(theta_vecd1, theta_vecd2, theta_vecd3, theta_vecd4),
  group = factor(rep(c("BAD in M0", "BAD in M1", "BAD in M2", "BAD in M4"), 
                     times = c(length(theta_vecd1), length(theta_vecd2), length(theta_vecd3), length(theta_vecd4))))
)

# Create the boxplot and zoom it to see the mean differences
main_plot_limited2 <- ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("coral", "limegreen", "mediumpurple", "gold")) +
  coord_cartesian(ylim = c(-5, 5)) +  # Limit the y-axis to the range -5 to 5
  labs(title = "Boxplots of Theta Values (Limited Y-Axis)", x = "Model", y = "Theta Values") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13))
main_plot_limited2
