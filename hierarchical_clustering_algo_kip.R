covMat <- read.csv('C:/Users/anany/Downloads/Cov_Mat - Cov_mat.csv', header = TRUE)
corMat <- read.csv('C:/Users/anany/Downloads/Cov_Mat - Cor_Mat.csv', header = TRUE)
clustOrder <- hclust(dist(corMat), method = 'single')$order
getIVP <- function(covMat) {
  invDiag <- 1/diag(as.matrix(covMat))
  weights <- invDiag/sum(invDiag)
  return(weights)
}
getClusterVar <- function(covMat, cItems) {
  covMatSlice <- covMat[cItems, cItems]
  weights <- getIVP(covMatSlice)
  cVar <- t(weights) %*% as.matrix(covMatSlice) %*% weights
  return(cVar)
}
getRecBipart <- function(covMat, sortIx) {
  w <- rep(1,ncol(covMat))
  w <- recurFun(w, covMat, sortIx)
  return(w)
}
recurFun <- function(w, covMat, sortIx) {
  subIdx <- 1:trunc(length(sortIx)/2)
  cItems0 <- sortIx[subIdx]
  cItems1 <- sortIx[-subIdx]
  cVar0 <- getClusterVar(covMat, cItems0)
  cVar1 <- getClusterVar(covMat, cItems1)
  alpha <- 1 - cVar0/(cVar0 + cVar1)
  
  # scoping mechanics using w as a free parameter
  w[cItems0] <- w[cItems0] * alpha
  w[cItems1] <- w[cItems1] * (1-alpha)
  
  if(length(cItems0) > 1) {
    w <- recurFun(w, covMat, cItems0)
  }
  if(length(cItems1) > 1) {
    w <- recurFun(w, covMat, cItems1)
  }
  return(w)
}
out <- getRecBipart(covMat, clustOrder)
out
require(tseries)
require(PerformanceAnalytics)
require(quantmod)
require(Quandl)
require(RTL)
Quandl.api_key("Shvn3L7jAJjZexzssGLR") # not displaying my own api key, sorry 


p<-q$Adj.Close
print(p)
# function to append missing (I.E. assets not selected) asset names and sort into original order
appendMissingAssets <- function(wts, allAssetNames, wtsDate) {
  absentAssets <- allAssetNames[!allAssetNames %in% names(wts)]
  absentWts <- rep(0, length(absentAssets))
  names(absentWts) <- absentAssets
  wts <- c(wts, absentWts)
  wts <- xts(t(wts), order.by=wtsDate)
  wts <- wts[,allAssetNames]
  return(wts)
}
symbols <- c("SPY", "VGK",    "EWJ",    "EEM",    "VNQ",    "RWX",    "IEF",    "TLT",    "DBC",    "GLD")    
rets <- list()

q<-read.csv("C:/Users/anany/Downloads/ETF Quotes/SPY.csv")
price<-xts(q$Adj.Close,order.by=as.Date(q$Date))
pr<-xts::lag.xts(price,k=1,na.pad=T)
spyret<- Return.calculate(pr,method="discrete")
quotes<-data.frame(Dates=q$Date)
t<-quotes[,-1]
for(i in 1:length(symbols)) {
  
  # quandl command to download from EOD database. Free users should use write.zoo in this loop.
  p<-paste("C:/Users/anany/Downloads/Quotes_ETF/",symbols[i],".csv",sep="")
  q<-read.csv(p)
  x<-symbols[i]

  
  price<-xts(q$Adj.Close,order.by=as.Date(q$Date))
  pr<-xts::lag.xts(price,k=1,na.pad=T)
  returns <- Return.calculate(pr,method="discrete")
  colnames(returns) <- symbols[i]
  rets[[i]] <- returns
}

as.Date(q$Date[0])
rets <- na.omit(do.call(cbind, rets))
invVolWts <- list()
minVolWts <- list()
hrpWts <- list()
ep <- endpoints(rets, on =  "months")
nMonths = 6 # month lookback (6 as per parameters from allocateSmartly)
nVol = 20 # day lookback for volatility (20 ibid)
for(i in 1:(length(ep)-nMonths)) {
  
  # get returns subset and compute absolute momentum
  retSubset <- rets[c(ep[i]:ep[(i+nMonths)]),]
  retSubset <- retSubset[-1,]
  moms <- Return.cumulative(retSubset)
  
  # select top performing assets and subset returns for them
  highRankAssets <- rank(moms) >= 6 # top 5 assets
  posReturnAssets <- moms > 0 # positive momentum assets
  selectedAssets <- highRankAssets & posReturnAssets # intersection of the above
  selectedSubset <- retSubset[,selectedAssets] # subset returns slice
  
  if(sum(selectedAssets)==0) { # if no qualifying assets, zero weight for period
    
    wts <- xts(t(rep(0, ncol(retSubset))), order.by=last(index(retSubset)))
    colnames(wts) <- colnames(retSubset)
    invVolWts[[i]] <- minVolWts[[i]] <- hrpWts[[i]] <- wts
    
  } else if (sum(selectedAssets)==1) { # if one qualifying asset, invest fully into it
    
    wts <- xts(t(rep(0, ncol(retSubset))), order.by=last(index(retSubset)))
    colnames(wts) <- colnames(retSubset)
    wts[, which(selectedAssets==1)] <- 1
    invVolWts[[i]] <- minVolWts[[i]] <- hrpWts[[i]] <- wts
    
  } else { # otherwise, use weighting algorithms
    
    cors <- cor(selectedSubset) # correlation
    volSubset <- tail(selectedSubset, nVol) # 20 day volatility
    vols <- StdDev(volSubset)
    covs <- t(vols) %*% vols * cors
    
    # minimum volatility using portfolio.optim from tseries
    minVolRets <- t(matrix(rep(1, sum(selectedAssets))))
    minVolWt <- portfolio.optim(x=minVolRets, covmat = covs)$pw
    names(minVolWt) <- colnames(covs)
    minVolWt <- appendMissingAssets(minVolWt, colnames(retSubset), last(index(retSubset)))
    minVolWts[[i]] <- minVolWt
    
    # inverse volatility weights
    invVols <- 1/vols 
    invVolWt <- invVols/sum(invVols) 
    invNames <- colnames(invVolWt)
    invVolWt <- as.numeric(invVolWt) 
    names(invVolWt) <- invNames
    invVolWt <- appendMissingAssets(invVolWt, colnames(retSubset), last(index(retSubset)))
    invVolWts[[i]] <- invVolWt
    
    # hrp weights
    clustOrder <- hclust(dist(cors), method = 'single')$order
    hrpWt <- getRecBipart(covs, clustOrder)
    names(hrpWt) <- colnames(covs)
    hrpWt <- appendMissingAssets(hrpWt, colnames(retSubset), last(index(retSubset)))
    hrpWts[[i]] <- hrpWt
  }
}
invVolWts <- round(do.call(rbind, invVolWts), 3) # round for readability
minVolWts <- round(do.call(rbind, minVolWts), 3)
hrpWts <- round(do.call(rbind, hrpWts), 3)

# allocate to cash if no allocation made due to all negative momentum assets
invVolWts$cash <- 0; invVolWts$cash <- 1-rowSums(invVolWts)
hrpWts$cash <- 0; hrpWts$cash <- 1-rowSums(hrpWts)
minVolWts$cash <- 0; minVolWts$cash <- 1-rowSums(minVolWts)

# cash value will be zero
rets$cash <- 0

# compute backtest returns
invVolRets <- Return.portfolio(R = rets, weights = invVolWts)
minVolRets <- Return.portfolio(R = rets, weights = minVolWts)
hrpRets <- Return.portfolio(R = rets, weights = hrpWts)
compare <- cbind(invVolRets, minVolRets, hrpRets)
colnames(compare) <- c("invVol", "minVol", "HRP")
charts.PerformanceSummary(compare)
rbind(table.AnnualizedReturns(compare), maxDrawdown(compare), CalmarRatio(compare))  
compare2<-cbind(invVolRets,minVolRets,hrpRets,spyret)
charts.PerformanceSummary(compare2)
rbind(table.AnnualizedReturns(compare2), maxDrawdown(compare2), CalmarRatio(compare2))
q<-read.csv("C:/Users/anany/Downloads/Quotes_ETF/OBXD.OL.csv")
price3<-xts(q$Adj.Close,order.by=as.Date(q$Date))
pr3<-xts::lag.xts(price3,k=1,na.pad=T)
obxret<- Return.calculate(as.numeric(pr3),method="discrete")
q2<-read.csv("C:/Users/anany/Downloads/Quotes_ETF/ACWI.csv")
price2<-xts(q2$Adj.Close,order.by=as.Date(q2$Date))
pr2<-xts::lag.xts(price2,k=1,na.pad=T)
acwiret<- Return.calculate(pr2,method="discrete")
compare3<-cbind(invVolRets,minVolRets,hrpRets,spyret,acwiret)
charts.PerformanceSummary(compare3)
rbind(table.AnnualizedReturns(compare3), maxDrawdown(compare3), CalmarRatio(compare3))
colnames(compare3)=c("invVol","minVol","HRP","SPY","ACWI")
chart.CumReturns(compare3,legend.loc = "topleft")
compare4<-cbind(compare3$HRP,compare3$minVol,compare3$ACWI)
compare5<-compare4*10000
chart.CumReturns(compare5,legend.loc = "topleft",main="Performance comparison of strategy against Global Stock Index ( MSCI ACWI)")
chart.CumReturns(compare5,legend.loc = "topleft", main="Performance comparison of strategy against Global Stock Index ( MSCI ACWI)")
colnames(compare4)<-c("Strategy_1","Strategy_2","Benchmark")
chart.CumReturns(compare4,legend.loc = "topleft")
plot(compare4)
cumret<-cumprod(1+compare4)
cumret<-cumret*100
compare4[is.na
         (compare4)]<-0
?plot
plot(cumret, main="Performance comparison of strategy against Global Stock Index (ACWI)", ylab="Returns (in %)", fill="red","black","green")
legend("bottomright",legend=c("Strategy_1","Strategy_2","Benchmark"))
 ?chart.CumReturns
results1<-rbind(table.AnnualizedReturns(compare4),maxDrawdown(compare4),CalmarRatio(compare4))
write.csv(results1,"D:/PB COLLAB/HRP.csv",row.names = TRUE)
write.csv(compare4,"D:/PB COLLAB/Returns_Dataset.csv",row.names = TRUE)
dr_percent<-Drawdowns(compare4)*100
plot(dr_percent, main="Drawdown comparison of strategy against Global Stock Index (ACWI) ",ylab="Drawdown in %age",legend("bottomright",legend=c("Strategy_1","Strategy_2","Benchmark")))
legend("topleft",legend=c("Strategy_1","Strategy_2","Benchmark"),col=c("balck","red","green"), bg="lightblue")

chart.CumReturns(compare4,legend)
fig<-plotly::plot_ly
library("ggplot2")
plt1<-ggplot2::ggplot()+
  geom_line(data=cumret,aes(x=rownames(cumret),y=Strategy_1),color="Black")+
  geom_line(data=cumret,aes(x=rownames(cumret),y=Strategy_2),color="Red")+
  geom_line(data=cumret,aes(x=rownames(cumret),y=Benchmark),color="Green")+
  xlab("Dates")+
  ylab("Returns")
print(plt1)
plt1<-autoplot(cumret,facets=FALSE,ylab="Returns in %",main = "Performance comparison of strategy against Global Stock Index (ACWI)",le)
plt1+labs(fill="Indices")
Plt2<-autoplot(dr_percent, facets=FALSE,main="Drawdown comparison of strategy against Global Stock Index (ACWI) ",ylab="Drawdown in %age")
print(Plt2)