#### Global variables and configurations ####

NSeries <- 50000 # Sets the number of series used in the MonteCarlo Simulation

numberBTCMiners <- 1750 # Used to set BTC production cap and calculate investment

BTCStartType <- "triangular" # Used to determine the starting value of BTC diffusion process
# Possible values: "deterministic", "triangular", "regression"
# Fine tuning can be done in the BTCStart section; overriden if excelVersion = TRUE

extraProduction <- TRUE # Set TRUE to sell excedent of power production at PLD prices
# Always FALSE if excelVersion is TRUE

testBTCReturns <- FALSE # Measure mu and sigma of BTC diffusions statistics
# Can slow down calculations

savePDF <- FALSE # Set TRUE to save histograms in PDF format
# See other variables at the chart generation command line

# Global Chart Configurations
require(repr) # Used to resize plots
options(repr.plot.width=8, repr.plot.height=3)

#### Wind Power Data and Variables ####

wind.velocity = c(2.3,2.6,2.4,2.1,2.1,2.95,3.5,3.6,3.8,3.1,3,2.5) # Based on (Lira and Moita Neto, 2017)

windToPower <- 1871.208247 # Given by the technology used

power.production <- wind.velocity * windToPower

## Chart ##

barplot(power.production, col="blue", xlab="Average Monthly Output", ylab="MWh", names.arg=month.abb, ylim=c(0,8000))

#### Stochastic variables estimation ####

## PLD Data ##

PLDMax <- 150 # As defined by ONS + expected future increase adjustment

PLDStart <- 75 # Defined in paper

PLDEta <- 0.08038 # Defined in paper

PLDRevertingMean <- 86.30 # Defined in paper

PLDVolatility <- 0.557 # Defined in paper

## BTC Data ##

BTCMu <- -0.0366 # Defined in paper

BTCSigma <- 0.2228 # Defined in paper

## Diffusion Functions ##

PLD_Diffusion <- function(start, len, n, RevertingMean, Volatility, Eta, Max) {
    # Diffusion1 <- RevertingMean*exp(-Volatility^2/(4*Eta))
    Diffusion2 <- exp(-Eta)
    Diffusion3 <- (log(RevertingMean)-Volatility^2/(2*Eta))*(1-Diffusion2)
    Diffusion4 <- Volatility*sqrt((1-Diffusion2^2)/(2*Eta))
    
    x = matrix(NA, nrow=(len+1), ncol=n)
    x[1, ] = start
    
    for(i in 2:(len+1)){
        x[i, ] = exp(log(x[i-1, ])*Diffusion2+Diffusion3+Diffusion4*rnorm(n,0,1))
        x[i, x[i, ] > Max] = Max
    }
    return(x)
}

PLD_Deterministic <- function(start, len, RevertingMean, Volatility, Eta, Max) {
    Diffusion1 <- log(RevertingMean)-Volatility^2/(2*Eta)
    Diffusion2 <- Volatility^2/(4*Eta)
    
    x = rep(NA, len)
    
    for(i in 1:len){
        x[i] = min(Max, exp(log(start)*exp(-Eta*i)+Diffusion1*(1-exp(-Eta*i))+Diffusion2*(1-exp(-2*Eta*i))))
    }
    return(c(start, x))
}

BTC_Diffusion <- function(start, len, n, mu, sigma) {
    x = matrix(NA, nrow=(len+1), ncol=n)
    x[1, ] = start
    
    for(i in 2:(len+1)){
        x[i, ] = x[i-1, ]*exp(rnorm(n, mu-sigma^2/2, sigma))
    }
    return(x)
}

BTC_Deterministic <- function(start, len, mu, sigma) {
    x = rep(NA, (len+1))
    
    if (length(start) > 1) start <- mean(start)
    x[1] = start
    
    for(i in 2:(len+1)){
        x[i] = x[i-1]*exp(mu)
    }
    return(x)
}

## Printing AUxiliar Function ##

format_text <- function(text, decPlaces) { 
    format(round(as.numeric(text), decPlaces), nsmall=decPlaces, big.mark=",") # Display results with comma separator and without decimals
}

## BTC Mining Data ## 

BTCInitialPrice <- 8105 # as of 2019-11-19 
# Used if BTCStart is set to "deterministic"

BTCInitialDifficulty <- 12720005267390.5 # as of 2019-11-19
# Used for Deterministic and Triangular starts

BTCMinerHash <- 56*10^12 # Antminer S17 Pro - Used to calculate consumption costs

BTCMinerPower <- 2212 # Watts # Antminer S17 Pro - Used to calculate consumption costs

BTCMinerCost <- 1900 # Antminer S17 Pro - Used to calculate investment

BTCRegression <- 0.089497807 # Used if BTCStart is set to "regression"

BTCTriangle <- c(5000, 7000, 9000) # Used if BTCStart is set to "triangular"

## BTC Starting point ##

BTCConsumptionFactor <- (2^32*BTCMinerPower)/(BTCMinerHash*3600*12.5*1000) # kWh/BTC

BTCStartDeterministic <- BTCInitialPrice / (BTCConsumptionFactor * BTCInitialDifficulty)

BTCStartRegression <- BTCRegression / (BTCConsumptionFactor*10^8)

require(extraDistr)

BTCStartTriangle <- rtriang(NSeries, BTCTriangle[1], BTCTriangle[3], BTCTriangle[2]) / (BTCConsumptionFactor * BTCInitialDifficulty)

BTCStart <- switch(BTCStartType, "deterministic" = BTCStartDeterministic, "triangular" = BTCStartTriangle, "regression" = BTCStartRegression)

if (BTCStartType == "triangular") hist(BTCStart, border="darkred", col="red", xlab="BTC Price/Difficulty Starting Point", main="", freq=FALSE)

## PLD Series ##

PLD.matrix <- PLD_Diffusion(PLDStart, 72, NSeries, PLDRevertingMean, PLDVolatility, PLDEta, PLDMax)

PLD.period.means <- apply(PLD.matrix, 1, mean) # Means for each period from all simulations

PLD.deterministic.series <- PLD_Deterministic(PLDStart, 72, PLDRevertingMean, PLDVolatility, PLDEta, PLDMax)

## Chart ##

matplot(1:73, cbind(PLD.matrix[, 1000], PLD.period.means, PLD.deterministic.series), type='l', xlab='Periods', ylab='series')
legend("top", inset=.02, legend=c("Random Series","MonteCarlo Mean","Deterministic Series"),col=c("black", "red", "green"), lty=1:3, cex=0.6, horiz=TRUE, bty="n")

## BTC Series and Charts ##

BTC.1.deterministic.series <- c(rep(0, 22), BTC_Deterministic(BTCStart, 26, BTCMu, BTCSigma), rep(0, 24))

BTC.2.deterministic.series <- c(rep(0, 22), BTC_Deterministic(BTCStart, 26, BTCMu, BTCSigma), BTC_Deterministic(BTCStart, 26, BTCMu, BTCSigma)[-(1:3)])

plot(BTC.1.deterministic.series, type="l", xlab="Period", ylab="Price/Difficulty")

plot(BTC.2.deterministic.series, type="l", xlab="Period", ylab="Price/Difficulty")

BTC.1.matrix <- BTC_Diffusion(BTCStart, 26, NSeries, BTCMu, BTCSigma)
BTC.1.matrix <- rbind(matrix(0, nrow=22, ncol=NSeries), BTC.1.matrix, matrix(0, nrow=24, ncol=NSeries))

BTC.1.period.means <- apply(BTC.1.matrix, 1, mean) # Means for each period from all simulations

matplot(1:73, cbind(BTC.1.matrix[, 2], BTC.1.period.means, BTC.1.deterministic.series), type='l', xlab='Periods', ylab='series')
legend("top", inset=.02, legend=c("Random Series","MonteCarlo Mean","Deterministic Series"),col=c("black", "red", "green"), lty=1:3, cex=0.6, horiz=TRUE, bty="n")

BTC.2.matrix <- BTC_Diffusion(BTCStart, 26, NSeries, BTCMu, BTCSigma)

BTC.2.matrix <- BTC.2.matrix[-c(1:3), ]
BTC.2.matrix <- rbind(BTC.1.matrix[1:49, ], BTC.2.matrix)

BTC.2.period.means <- apply(BTC.2.matrix, 1, mean) # Means for each period from all simulations

matplot(1:73, cbind(BTC.2.matrix[, 2], BTC.2.period.means, BTC.2.deterministic.series), type='l', xlab='Periods', ylab='series')
legend("top", inset=.02, legend=c("Random Series","MonteCarlo Mean","Deterministic Series"),col=c("black", "red", "green"), lty=1:3, cex=0.6, horiz=TRUE, bty="n")

## Tests for BTC Diffusion Process ##

create_tests_results <- function(diffusion.table, rowname, mu, sigma) {
    
    diffusion.diff <- apply(log(diffusion.table), 2, diff)
    diffusion.SDs <- apply(diffusion.diff, 2, sd)
    diffusion.SDs.mean <- mean(diffusion.SDs)
    diffusion.SDs.sd <- sd(diffusion.SDs)
    
    diffusion.means <- apply(diffusion.diff, 2, mean)
    diffusion.means.mean <- mean(diffusion.means)
    diffusion.means.sd <- sd(diffusion.means)
    
    delta.sd <- mean(diffusion.SDs) - sigma
    delta.mu <- mean(diffusion.means) - (mu - sigma^2/2)
    
    x.vec = format_text(c(diffusion.SDs.mean, diffusion.SDs.sd, diffusion.means.mean, diffusion.means.sd, delta.sd, delta.mu), 4)
    
    x.matrix <- matrix(c(x.vec), nrow=1)
    colnames(x.matrix) <- c("Mean of SDs", "SD of SDs", "Mean of means", "SD of means", "SD test", "Mu test")
    rownames(x.matrix) <- rowname
    
    return(list(Table = x.matrix, MuVec = diffusion.means))
}

if (testBTCReturns) {
    test_results_1 <- create_tests_results(BTC.1.matrix[23:49, ], "BTC 1", BTCMu, BTCSigma)
    test_results_2 <- create_tests_results(BTC.2.matrix[50:73, ], "BTC 2", BTCMu, BTCSigma)
    
    rbind(test_results_1$Table, test_results_2$Table)
}

if (testBTCReturns) hist(test_results_1$MuVec, border="darkred", col="red", xlab="Mu BTC 1", main="", freq=FALSE)

if (testBTCReturns) hist(test_results_2$MuVec, border="darkred", col="red", xlab="Mu BTC 2", main="", freq=FALSE)

#### Obtaining Cash Flow and NPV ####

## Project Cash Flow Data ##

BRLUSDEx <- 4 # R$/US$

InitialInvestment <- 9379943/BRLUSDEx # Anticipation of the investment based on (Fontanet, 2012) # Defined in paper

FixedCosts <- 52757/BRLUSDEx # Based on (Fontanet, 2012) # Defined in paper

VariableCosts <- 0.14 # % # Defined in paper

WACC <- 0.08 # Annual # Defined in paper

RF <- 0.05 # Annual # Defined in paper

Taxes <- 0.25*0.08+0.09*0.12 # Simplified tax regime # Defined in paper

Refrig <- 0.85 # Defined in paper

## Project Cash Flow Functions ##

P_and_L <- function(revenue, variableCosts, fixedCost, taxesPerc, depreciation, taxOnRevenue) {
    if (taxOnRevenue) {
        PandL <- revenue*(1-taxesPerc)-variableCosts-fixedCost
    } else {
        PandL <- (revenue-variableCosts-fixedCost-depreciation)*(1-taxesPerc)+depreciation
    }
    return(PandL)
}

select_greater <- function(input1, input2) {
    # The best results from each scenario are chosen a posteriori 
    # This is a simplified version of an instant and more granular choice
    
    results <- input1
    index.vec <- which(input2 > input1)
    results[index.vec] <- input2[index.vec]
    return(results)
}

npv <- function(cf, r) {
    len <- length(cf)
    r.vec <- rep(NA, len)
    for (i in 1:len){
        r.vec[i] <- (1+r)^(i-1)
    }
    return(sum(cf / r.vec))
}

## Auxiliary variables and calculations ##

BTCMax <- (numberBTCMiners*BTCMinerPower/1000)/Refrig # Set BTC production cap

RFMonthly <- (1+RF)^(0.0833333333333333)-1

WACCMonthly <- (1+WACC)^(0.0833333333333333)-1

BTCInvestment.1 <- numberBTCMiners*BTCMinerCost*((1+RFMonthly)/(1+WACCMonthly))^22

BTCInvestment.2 <- numberBTCMiners*BTCMinerCost*((1+RFMonthly)/(1+WACCMonthly))^46

production.series <- c(NA, rep(power.production, 6))

capped.production.series <- production.series
capped.production.series[production.series > BTCMax] <- BTCMax

PLD.output <- PLD.matrix[26:73, ] * production.series[26:73]

BTC.1.output <- BTC.1.matrix[26:49, ] * capped.production.series[26:49] * Refrig * 1000
BTC.1.output <- rbind(BTC.1.output, PLD.output[25:48, ])

BTC.2.output <- BTC.2.matrix[26:73, ] * capped.production.series[26:73] * Refrig * 1000

if (extraProduction) { # Adds extra output if TRUE
  extra.production.series <- production.series - capped.production.series
  extra.output.2 <- PLD.matrix[26:73, ] * extra.production.series[26:73]
  extra.output.1 <- rbind(extra.output.2[1:24, ], matrix(0, nrow=24, ncol=NSeries))
  BTC.1.output <- BTC.1.output + extra.output.1
  BTC.2.output <- BTC.2.output + extra.output.2
}

## Cash Flow Results ##

PLD.cashflow <- P_and_L(PLD.output, PLD.output*VariableCosts, FixedCosts, Taxes, 0 , TRUE)
PLD.cashflow <- rbind(matrix(0, nrow=25, ncol=ncol(PLD.cashflow)), PLD.cashflow)

BTC.1.cashflow <- P_and_L(BTC.1.output, PLD.output*VariableCosts, FixedCosts, Taxes, 0, TRUE)
BTC.1.cashflow <- rbind(matrix(0, nrow=25, ncol=ncol(BTC.1.cashflow)), BTC.1.cashflow)
BTC.1.cashflow <- select_greater(PLD.cashflow, BTC.1.cashflow)

BTC.2.cashflow <- P_and_L(BTC.2.output, PLD.output*VariableCosts, FixedCosts, Taxes, 0, TRUE)
BTC.2.cashflow <- rbind(matrix(0, nrow=25, ncol=ncol(BTC.2.cashflow)), BTC.2.cashflow)
BTC.2.cashflow <- select_greater(PLD.cashflow, BTC.2.cashflow)

PLD.mean.cashflow <- apply(PLD.cashflow, 1, mean)
BTC.1.mean.cashflow <- apply(BTC.1.cashflow, 1, mean)
BTC.2.mean.cashflow <- apply(BTC.2.cashflow, 1, mean)

## Charts ##

i = 4

matplot(1:73, cbind(PLD.cashflow[, i], BTC.1.cashflow[, i], BTC.2.cashflow[, i]), type='l', xlab='Periods', ylab='series')
legend("topleft", inset=.1, legend=c("PLD","BTC2","BTC2+2"),col=c("black", "red", "green"), lty=1:3, cex=0.6, horiz=TRUE, bty="n", title=paste("MonteCarlo Series",i))

matplot(1:73, cbind(PLD.mean.cashflow, BTC.1.mean.cashflow, BTC.2.mean.cashflow), type='l', xlab='Periods', ylab='series')
legend("topleft", inset=.08, legend=c("PLD","BTC2","BTC2+2"),col=c("black", "red", "green"), lty=1:3, cex=0.6, horiz=TRUE, bty="n", title="MonteCarlo Means")

## Obtaining Deterministic Results for Debugging ##

PLD.output.series <- PLD.period.means * production.series
PLD.output.series[1:25] <- 0

BTC.1.output.series <- c(BTC.1.deterministic.series[1:49]*capped.production.series[1:49]*Refrig*1000, PLD.output.series[50:73])

BTC.2.output.series <- BTC.2.deterministic.series*capped.production.series*Refrig*1000

if(extraProduction) { # Adds extra output if TRUE
    extra.output.2 <- extra.production.series[26:73]*PLD.deterministic.series[26:73]
    extra.output.1 <- c(extra.output.2[1:24], rep(0, 24))
    BTC.1.output.series[26:73] <- BTC.1.output.series[26:73] + extra.output.1
    BTC.2.output.series[26:73] <- BTC.2.output.series[26:73] + extra.output.2
}

PLD.deterministic.cashflow <- PLD.output.series*(1-VariableCosts-Taxes)-FixedCosts

BTC.1.deterministic.cashflow <- BTC.1.output.series*(1-Taxes)-PLD.output.series*VariableCosts-FixedCosts

BTC.2.deterministic.cashflow <- BTC.2.output.series*(1-Taxes)-PLD.output.series*VariableCosts-FixedCosts

## Chart ##

matplot(1:73, cbind(PLD.deterministic.cashflow, BTC.1.deterministic.cashflow, BTC.2.deterministic.cashflow), type='l', xlab='Periods', ylab='series')
legend("topleft", inset=.08, legend=c("PLD","BTC2","BTC2+2"),col=c("black", "red", "green"), lty=1:3, cex=0.6, horiz=TRUE, bty="n", title="Deterministic Series")

## Calculating Project NPV ##

## NPV for mean results - used for debugging only ##

NPV.PLD.mean.cashflow <- PLD.mean.cashflow
NPV.PLD.mean.cashflow[1] <- -InitialInvestment

NPV.PLD.mean <- npv(NPV.PLD.mean.cashflow, r=RFMonthly)

NPV.BTC.1.mean.cashflow <- BTC.1.mean.cashflow
NPV.BTC.1.mean.cashflow[1] <- -InitialInvestment
NPV.BTC.1.mean.cashflow[23] <- BTC.1.mean.cashflow[23]-BTCInvestment.1

NPV.BTC.1.mean <- npv(NPV.BTC.1.mean.cashflow, r=RFMonthly)

NPV.BTC.2.mean.cashflow <- BTC.2.mean.cashflow
NPV.BTC.2.mean.cashflow[1] <- -InitialInvestment
NPV.BTC.2.mean.cashflow[23] <- BTC.2.mean.cashflow[23]-BTCInvestment.1
NPV.BTC.2.mean.cashflow[47] <- BTC.2.mean.cashflow[47]-BTCInvestment.2

NPV.BTC.2.mean <- npv(NPV.BTC.2.mean.cashflow, r=RFMonthly)

## Project NPV for each scenario ##

NPV.PLD.cashflow <- PLD.cashflow
NPV.PLD.cashflow[1, ] <- NPV.PLD.cashflow[1, ]-InitialInvestment

NPV.PLD <- apply(NPV.PLD.cashflow, 2, npv, r=RFMonthly)

NPV.BTC.1.cashflow <- BTC.1.cashflow
NPV.BTC.1.cashflow[1, ] <- NPV.BTC.1.cashflow[1, ]-InitialInvestment
NPV.BTC.1.cashflow[23, ] <- NPV.BTC.1.cashflow[23, ]-BTCInvestment.1

NPV.BTC.1 <- apply(NPV.BTC.1.cashflow, 2, npv, r=RFMonthly)

NPV.BTC.2.cashflow <- BTC.2.cashflow
NPV.BTC.2.cashflow[1, ] <- NPV.BTC.2.cashflow[1, ]-InitialInvestment
NPV.BTC.2.cashflow[23, ] <- NPV.BTC.2.cashflow[23, ]-BTCInvestment.1
NPV.BTC.2.cashflow[47, ] <- NPV.BTC.2.cashflow[47, ]-BTCInvestment.2

NPV.BTC.2 <- apply(NPV.BTC.2.cashflow, 2, npv, r=RFMonthly)

## Third scenario not included in paper ##

choice.vector <- (NPV.PLD > NPV.BTC.1)

BTC.3.cashflow <- BTC.2.cashflow
BTC.3.cashflow[49:73, choice.vector] <- PLD.cashflow[49:73, choice.vector]

NPV.BTC.3.cashflow <- BTC.3.cashflow
NPV.BTC.3.cashflow[1, ] <- NPV.BTC.3.cashflow[1, ]-InitialInvestment
NPV.BTC.3.cashflow[23, ] <- NPV.BTC.3.cashflow[23, ]-BTCInvestment.1
NPV.BTC.3.cashflow[47, !choice.vector] <- NPV.BTC.3.cashflow[47, !choice.vector]-BTCInvestment.2

NPV.BTC.3 <- apply(NPV.BTC.3.cashflow, 2, npv, r=RFMonthly)

#### Displaying Results ####

## Auxiliary Functions ##

create_more_results_table <- function(NPV, rowname) {
    
    meanNPV = mean(NPV)
    meanNegNPV = mean(NPV[NPV < 0])
    meanPosNPV = mean(NPV[NPV > 0])
    
    nElements = round(length(NPV)*0.05,0)
    CVaRNPV = mean(sort(NPV)[1:nElements])
    VaRNPV =  sort(NPV)[nElements]
    
    lambda <- 0.5
    ECPNPV = (1-lambda)*meanNPV + lambda*CVaRNPV
    EqECPNPV = meanNPV + lambda/(1-lambda) * (CVaRNPV - VaRNPV)
    
    upsideIndex = meanPosNPV/abs(meanNegNPV)
    
    x.vec = format_text(c(meanNPV, CVaRNPV, ECPNPV, EqECPNPV), 0)
    x.vec = c(x.vec, format_text(upsideIndex, 2))
    x.vec = c(x.vec, format_text(c(meanNegNPV, meanPosNPV), 0))
        
    x.matrix <- matrix(x.vec, nrow=1)
    colnames(x.matrix) <- c("Mean of NPV", "CVaR at 95%", "ECP at 0.5", "Certainty Eq.", "Upside Potential", "Mean of -NPV", "Mean of +NPV")
    rownames(x.matrix) <- rowname
    
    return(x.matrix)
}

create_results_table <- function(NPV, rowname, BaseNPV = NULL) {
    
    meanNPV = mean(NPV)
    
    nNegNPV = sum(NPV < 0)
    nElements = length(NPV)
    PercNeg = nNegNPV/nElements
    
    if (!is.null(BaseNPV)) {
        compareMeans = meanNPV/mean(BaseNPV) - 1
        comparePerc = PercNeg - sum(BaseNPV < 0)/nElements
        compare.vec = format_text(c(compareMeans, comparePerc), 2) 
    } else {
        compare.vec = rep("-", 2)
    }
    
    x.vec = c(format_text(meanNPV, 0), format_text(PercNeg, 2))
    x.vec = c(x.vec, compare.vec)
        
    x.matrix <- matrix(x.vec, nrow=1)
    colnames(x.matrix) <- c("Mean of NPV", "Perc. Neg.", "NPV/Base", "Perc. - Base")
    rownames(x.matrix) <- rowname
    
    return(x.matrix)
}

create_results_hist <- function(NPV, breaksLen, plotXLim, textXAdj, textYAdj, subTextAdj, lineAdj, plotWidth, plotHeight) {
    options(repr.plot.width=plotWidth, repr.plot.height=plotHeight)
    
    NPV = NPV/1000000
    meanNPV = mean(NPV)
    minNPV = min(NPV)
    maxNPV = max(NPV)
    percNeg <- round(as.numeric(sum(NPV < 0)/length(NPV)*100), 0)
    
    breaks.vec <- seq(minNPV, maxNPV+breaksLen, by=breaksLen)
    
    quant95 <- quantile(NPV, 0.95)
    text.pos.x <- c(-textXAdj, quant95/2, quant95 + textXAdj)
    
    text.labels <- c(paste(percNeg, "%", sep=""), paste(95 - percNeg, "%", sep=""), "5%")
    text.labels <- paste("<", text.labels, ">")
    text.colors <- c("black", "red", "black")
    
    hist.obj <- hist(NPV, breaks=breaks.vec, border="darkred", col="red", xlab="NPV in Millions", main="", xlim=plotXLim, freq=FALSE)
    abline(v=c(0, quant95), lty=2)

    textYPos <- max(hist.obj$density) + textYAdj
    text(meanNPV, subTextAdj, labels=paste("| Mean:", round(meanNPV, 2)), adj=0, cex=0.8)
    text(quant95, subTextAdj, labels=paste(" Q95:", round(quant95, 2)), adj=0, cex=0.8)

    text(text.pos.x, textYPos, labels=text.labels, col=text.colors, cex=0.8)
    text(-textXAdj, textYPos - lineAdj, labels="Min:", cex=0.8)
    text(quant95+textXAdj, textYPos - lineAdj, labels="Max:", cex=0.8)
    text(-textXAdj, textYPos - 2 * lineAdj, labels=round(minNPV, 2), cex=0.8)
    text(quant95+textXAdj, textYPos - 2 * lineAdj, labels=round(maxNPV, 2), cex=0.8)
    
    return(hist.obj)
}

## Mean results - for debugging purposes only ##

paste("NPV of PLD Mean Series: ", format_text(mean(NPV.PLD.mean), 0)) 

paste("NPV of BTC_1 Mean Series: ", format_text(mean(NPV.BTC.1.mean), 0)) 

paste("NPV of BTC_2 Mean Series: ", format_text(mean(NPV.BTC.2.mean), 0)) 


#### Base Scenario Results ####

results.matrix.PLD = create_results_table(NPV.PLD, "Base")
results.matrix.PLD

## Chart ##

if (savePDF) pdf(file="./Plots/Base.pdf", width=8, height=5)

PLD.hist <- create_results_hist(NPV.PLD, 0.625, c(-10,30), 3.5, 0.002, -0.0036, 0.008, 8, 6)

if (savePDF) dev.off()

#### First Scenario Results ####

results.matrix.BTC.1 = create_results_table(NPV.BTC.1, "2 Years", NPV.PLD)
results.matrix.BTC.1

## Chart ##

if (savePDF) pdf(file="./Plots/2+0_Scenario.pdf", width=8, height=5)

BTC.1.hist <- create_results_hist(NPV.BTC.1, 0.625, c(-10,30), 3.5, 0.002, -0.0022, 0.004, 8, 6)

if (savePDF) dev.off()

#### Second Scenario Results ####

results.matrix.BTC.2 = create_results_table(NPV.BTC.2, "2+2 Years", NPV.PLD)
results.matrix.BTC.2

## Chart ##

if (savePDF) pdf(file="./Plots/2+2_Scenario.pdf", width=8, height=5)

BTC.2.hist <- create_results_hist(NPV.BTC.2, 0.625, c(-10,30), 3.5, 0.002, -0.0017, 0.004, 8, 6)

if (savePDF) dev.off()

#### Third Scenario Results - not included in paper ####

ChoiceScenarios <- NSeries - sum(choice.vector)
PercChoice <- ChoiceScenarios/NSeries*100
paste("Scenarios with Choice:", ChoiceScenarios, "out of", NSeries, "(", PercChoice, "%)") 

results.matrix.BTC.3 = create_results_table(NPV.BTC.3, "2 or 4 Years", NPV.PLD)
results.matrix.BTC.3

## Chart ##

if (savePDF) pdf(file="./Plots/2+2_Scenario+Choice.pdf", width=8, height=5)

BTC.3.hist <- create_results_hist(NPV.BTC.3, 0.625, c(-10,30), 3.5, 0.002, -0.0017, 0.004, 8, 6)

if (savePDF) dev.off()

#### Aggregated Results ####

rbind(results.matrix.PLD, results.matrix.BTC.1, results.matrix.BTC.2, results.matrix.BTC.3)
