
arima <- function(phi, theta, len, c, media = 0, sigma = 1)
{
    nn<-length(phi)
    mm<-length(theta)
    
    ruid <- rnorm(n = len, mean = media, sd = sigma )
    ruido<- c(numeric(mm),ruid) #inserting zeros for initial condition
    x <- numeric(len+nn) #inserting zeros for initial condition
    
    for( i in  nn:(len-1) ) {
        
        print(i)
        x[i+1]<- c + x[i+1] + ruido[i+1]
        
        for (j in 1:nn){
            #print(j)
            x[i+1]<-  x[i+1] + phi[j] * x[i-j+1]
            
        }
        
        for(j in 1:mm){
            #print(j)
            
            x[i+1]<- x[i+1] + theta[j]*ruido[i+1-j]
            
        }
        
    }
    
    return(tail(x,len - nn))
}




rednoise <- function(phi, len, media = 0, sigma = 1)
{
    ruido <- rnorm(n = len,mean = media, sd = sigma )
    x <- numeric(len)
    for( i in  1:(len-1) ) {x[i+1]<- phi*x[i]+ruido[i]}
    return(x)
}


ColouredGaussian<-function( len, media = 0, sigma = 1, alpha)
{
    # alpha = 1/tau
    
    ruido<-rnorm(n = len,mean = media, sd = sigma )
    
    x<-numeric(len)
    
    f <- exp(- alpha )
    
    x[1] <- ruido[1]
    
    
    for(i in 2:len)
    {
        x[i]<-  f*x[i-1] + sqrt(1-f**2)*ruido[i]
    }
    
    return(x)
    
}

ModelHeitorNorm<-function( len, media = 0, sigma = 1, alpha)
{
    # alpha = 1/tau
    
    ruido<-rnorm(n = len,mean = media, sd = sigma )
    
    x<-numeric(len)
    
    f <- exp(- alpha )
    
    x[1] <- ruido[1]
    
    
    for(i in 2:len)
    {
        x[i]<-  f*x[i-1] + (1-f)*ruido[i]
    }
    
    return(x)
    
}

ModelHeitorUnif<-function( len, minimo = -1, maximo = 1, alpha)
{
    # alpha = 1/tau
    
    ruido<-runif(n = len, min = minimo, max = maximo )
    
    x<-numeric(len)
    
    f <- exp(- alpha )
    
    x[1] <- ruido[1]
    
    
    for(i in 2:len)
    {
        x[i]<-  f*x[i-1] + (1-f)*ruido[i]
    }
    
    return(x)
    
}

ModelHeitorDiscrete<-function( len,  alpha)
{
    # alpha = 1/tau
    #ruido<- sample(0:1,size = 100 ,replace = TRUE)
    ruido<- sample(0:1,size = len ,replace = TRUE)
    ruido[which(ruido==0)]<- -1
    
    x<-numeric(len)
    
    f <- exp(- alpha )
    
    x[1] <- ruido[1]
    
    
    for(i in 2:len)
    {
        x[i]<-  f*x[i-1] + (1-f)*ruido[i]
    }
    
    return(x)
}

ColouredGaussianns<-function(ns, len, media = 0, sigma = 1, alpha)
{
    # alpha = 1/tau
    x<-matrix(ncol = ns , nrow = len)
    
    f <- exp(- alpha )
    
    for(j in 1:ns)
    {
        ruido<-rnorm(n = len,mean = media, sd = sigma )
        
        x[1,j] <- ruido[1]
        
        
        for(i in 2:len)
        {
            x[i,j]<-  f*x[i-1,j] + sqrt(1-f**2)*ruido[i]
        }
        
    }
    
    return(x)
    
}



RuidoColorido<- function(alfa,len, media = 0, sigma = 1, dts=0.0001)
{
    
    xi<-0
    xii<-numeric(len)
    ruido<-rnorm(mean = media, sd = sigma,n = len)*dts
    
    for( i in 1:999)
    {
        for(j in 1:999)
        {
            xi<-xi + exp(-alfa*i)*0.5*(ruido[j]+ruido[j+1])
        }
        xii[i]<-xi
    }
    return(xii)
}

DiscreteDerivative<-function(a,b){
    
    
    nn<-length(a)-1
    
    fim<-numeric(nn)
    
    for(i in 1: nn )
    {
        fim[i]<-(a[i+1]-a[i])/(b[i+1]-b[i])
        #  derivhard[i]<-(fim$S[i+1]-fim$S[i])/(fim$E[i+1]-fim$E[i])
    }
    
    return(fim)
    
}
