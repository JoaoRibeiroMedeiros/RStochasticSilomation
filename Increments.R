


incrementsD<-function(x,D) {
    
    for(i in 1:D)
    {
        x0<-head(x,length(x)-1)
        x1<-tail(x,length(x)-1)
        x<-(x1-x0)
    }

    return(x)
    
}


sidemean<-function(matrix, axis = 1) {
    
    sidesize <- dim(matrix)[axis]
    
    finalmean <- numeric(sidesize)
    
    for(i in 1:sidesize)
    {
        if(axis==1)
        {
            finalmean[i]<- mean(matrix[i,])
        }
        else
        {
            finalmean[i]<- mean(matrix[,i])
        }
    }
    
    return(finalmean)
}

lseq <- function(from = 1, to = 100000, length.out = 6) {

    exp(seq(log(from), log(to), length.out = length.out))
}

incrementos<-function(x) {
    
    x0<-head(x,length(x)-1)
    x1<-tail(x,length(x)-1)
    
    v<-(x1-x0)
    v0<-head(v,length(x)-2)
    v1<-tail(v,length(x)-2)
    
    a<-(v1-v0)
    a0<-head(a,length(x)-3)
    a1<-tail(a,length(x)-3)
    
    da<-(a1-a0)
    da1<-tail(da,length(x)-4)
    da0<-head(da,length(x)-4)
    
    dda<-(da1-da0)
    dda1<-tail(dda,length(x)-5)
    dda0<-head(dda,length(x)-5)
    
    ddda<-(dda1-dda0)
    
    
    
    final<-list("v"= v, "a"= a, "da" = da, "dda" = dda, "ddda" = ddda )
    
    return(final)
}
