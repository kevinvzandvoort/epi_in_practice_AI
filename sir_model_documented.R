# SIR Model for Covid-19 Transmission in R (with ggplot2 plotting)
# ----------------------------------------------------------------
# This script implements a simple SIR (Susceptible-Infected-Recovered) epidemic model
# using ordinary differential equations (ODEs), parameterized by R0.
# Epidemic curves are plotted using ggplot2 for high-quality graphics.
#
# Author: kevinvzandvoort
# Date: 2025-10-15
#
# Requirements:
#   - R (https://cran.r-project.org/)
#   - deSolve package for ODE integration
#   - ggplot2 package for plotting
#
# Usage:
#   - Adjust R0_value and gamma as needed to explore different scenarios.
#   - Run the script to simulate and plot the epidemic curves (S, I, R).
#
# References:
#   - Kermack, W.O. & McKendrick, A.G. (1927). A contribution to the mathematical theory of epidemics.
#   - https://github.com/kevinvzandvoort/epi_in_practice_AI

# Install and load required packages
if (!requireNamespace("deSolve", quietly = TRUE)) {
  install.packages("deSolve")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(deSolve)
library(ggplot2)

# -----------------------------------------------------------------------------
# SIR model function
# -----------------------------------------------------------------------------
# Arguments:
#   time      : Current time (numeric, from ODE solver)
#   state     : Named vector of current states: S (Susceptible), I (Infected), R (Recovered)
#   parameters: Named vector of model parameters: R0, gamma, N (population size)
# Returns:
#   List of derivatives (dS, dI, dR)
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Transmission rate (beta) is calculated from R0 and gamma
    beta <- R0 * gamma
    
    # SIR differential equations
    dS <- -beta * S * I / N         # Susceptible individuals becoming infected
    dI <- beta * S * I / N - gamma * I  # Infected individuals becoming recovered
    dR <- gamma * I                 # Recovered individuals
    
    # Return list of derivatives for ODE solver
    return(list(c(dS, dI, dR)))
  })
}

# -----------------------------------------------------------------------------
# Initial conditions and parameters
# -----------------------------------------------------------------------------
N     <- 1e6      # Total population size (e.g., 1 million)
I0    <- 1        # Initial number of infected individuals
R0    <- 0        # Initial number of recovered individuals (usually 0)
S0    <- N - I0   # Initial number of susceptible individuals
init  <- c(S = S0, I = I0, R = R0) # Named vector of initial state

# -----------------------------------------------------------------------------
# Model parameters
# -----------------------------------------------------------------------------
R0_value  <- 2.5      # Basic reproduction number (typical for Covid-19: 2-3)
gamma     <- 1/5      # Recovery rate (1/mean infectious period, e.g., 5 days)
parameters <- c(R0 = R0_value, gamma = gamma, N = N) # Named vector of parameters

# -----------------------------------------------------------------------------
# Time sequence for simulation
# -----------------------------------------------------------------------------
times <- seq(0, 160, by = 1) # Simulate for 160 days in daily steps

# -----------------------------------------------------------------------------
# Run SIR model simulation
# -----------------------------------------------------------------------------
out <- ode(y = init, times = times, func = sir_model, parms = parameters)

# Convert output to a data frame for plotting and analysis
out <- as.data.frame(out)

# -----------------------------------------------------------------------------
# Reshape data for ggplot2
# -----------------------------------------------------------------------------
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}
library(reshape2)
out_long <- melt(out, id.vars = "time", measure.vars = c("S", "I", "R"),
                 variable.name = "Compartment", value.name = "Count")

# -----------------------------------------------------------------------------
# Plot epidemic curves: Susceptible, Infected, Recovered using ggplot2
# -----------------------------------------------------------------------------
ggplot(out_long, aes(x = time, y = Count, color = Compartment)) +
  geom_line(size = 1.2) +
  labs(title = "SIR Model: Covid-19 Transmission",
       x = "Time (days)",
       y = "Number of people",
       color = "Compartment") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c(S = "blue", I = "red", R = "green"))

# -----------------------------------------------------------------------------
# Output key results (peak infected, total infected at end)
# -----------------------------------------------------------------------------
cat("Peak number of infected individuals:", max(out$I), "\n")
cat("Total number of recovered individuals at epidemic end:", tail(out$R, 1), "\n")

# -----------------------------------------------------------------------------
# End of script
# -----------------------------------------------------------------------------
