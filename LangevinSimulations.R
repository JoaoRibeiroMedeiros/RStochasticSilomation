SimulacaoLangevin <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    
    it<- 0
    
    
    
    for(i in 1:ns)
    {
        
        print(i)
        
        x<-0
        v<-0
        u<-0
        y<-0
        
        ji1<-0
        jd1<-0
        
        
        it<-0
        ic<-0
        ir<-0
        ig<-0
        
        forceest0<-0
        forcedet0<-0
        
        while(it < ntc)
        {
            
            it <- it + 1
            
            forceest <- sqrt(2*gamma*temp*dt)*rnorm(1)
            
            #forceeststrato<- (forceest + forceest0)/2
            #forceest0<-forceest
            
            forcedet <-  (-gm*u -c1*y)*dt
            
            #forcedetstrato <- (forcedet + forcedet0)/2
            #forcedet0 <- forcedet
            
            v <- u + forcedet + forceest
            x <- y + v*dt
            
            ji1 <- ji1 + (v+u)*forceest/2
            jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            # caso bimodal
            # ji1<-ji1+ (v+u)*eta*dt/2
            # jd1<-jd1+ gamma*(v^2+u^2)*dt/2
            
            ig <- ig + 1
            
            u <- v
            y <- x
            
            if(ig==nc)
            {
                
                xx[it/nc,i]<- x
                vv[it/nc,i]<- v
                ee[it/nc,i]<- (ji1-jd1)
                ji[it/nc,i]<- ji1
                jd[it/nc,i]<- jd1
                
                ig<- 0
                
            }
        }
        
    }
    
    
    final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp)
    
    return(final)
}


SimulacaoLangevinColorida <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp, dtemp, alfa, dalfa){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    ruido<-matrix( nrow = nt, ncol = ns)
    
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    
    it<- 0
    
    alfazero<- alfa
    
    tempzero <- temp
    
    alfa<- alfa- dalfa
    
    temp<- temp - dtemp
    
    for(i in 1:ns)
    {
        alfa<- alfa + dalfa
        temp<- temp + dtemp
        
        f<-exp(-alfa*dt)
        
        #force2<- force2 * exp(- alfa * dt  )  - gm * u * exp(- alfa * dt  ) * dt
        #sera q precisa declarar ?
        #forceestpre <- numeric(ntc)
        #forceestpos <- numeric(ntc)
        #forceest <- sqrt(2*gamma*temp*dt)*rnorm(1)
        
        forcedet0<- 0
        forceestpos<- 0
        
        print(i)
        
        x<-0
        v<-0
        u<-0
        y<-0
        
        ji1<-0
        jd1<-0
        
        
        it<-0
        ic<-0
        ir<-0
        ig<-0
        
        
        
        
        
        while(it < ntc)
        {
            it <- it + 1
            
            forceestpre <- sqrt(2*gamma*temp*dt)*rnorm(1)
            
            forceestpos <-  forceestpos*f + sqrt(1-f**2)*forceestpre
            
            forcedet <-  (-gm*u -c1*y)*dt
            
            v <- u + forcedet + forceestpos
            x <- y + v*dt
            
            ji1 <- ji1 + (v+u)*forceestpos/2
            jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            # caso bimodal
            # ji1<-ji1+ (v+u)*eta*dt/2
            # jd1<-jd1+ gamma*(v^2+u^2)*dt/2
            
            ig <- ig + 1
            
            u <- v
            y <- x
            
            if(ig==nc)
            {
                
                xx[it/nc,i]<- x
                vv[it/nc,i]<- v
                ee[it/nc,i]<- (ji1-jd1)
                ji[it/nc,i]<- ji1
                jd[it/nc,i]<- jd1
                
                ruido[it/nc,i]<-forceestpos
                
                ig<- 0
                
            }
        }
        
    }
    
    
    final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee, "ruido" = ruido, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp, "dalfa" = dalfa, "alfa"= alfazero, "dtemp" = dtemp)
    
    return(final)
}

SimulacaoLangevinColoridaInternal <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    
    it<- 0
    
    
    
    for(i in 1:ns)
    {
        
        print(i)
        
        x<- 0
        v<- 0
        u<- 0
        y<- 0
        
        ji1<- 0
        jd1<- 0
        
        
        it<- 0
        ic<- 0
        ir<- 0
        ig<- 0
        
        forceest0<- 0
        forcedet0<- 0
        force2<- 0
        
        while(it < ntc)
        {
            
            it <- it + 1
            
            forceest <- sqrt(2*gamma*temp*dt)*rnorm(1)
            
            #forceeststrato<- (forceest + forceest0)/2
            #forceest0<-forceest
            
            force2<- force2 * exp(- alfa * dt  )  - gm * (v+u) * exp(- alfa * dt  ) * dt/(2* alfa)
            
            
            forcedet <- force2 + ( -c1*y)*dt
            
            #forcedetstrato <- (forcedet + forcedet0)/2
            #forcedet0 <- forcedet
            
            v <- u + forcedet + forceest
            x <- y + v*dt
            
            ji1 <- ji1 + (v+u)*forceest/2
            jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            # caso bimodal
            # ji1<-ji1+ (v+u)*eta*dt/2
            # jd1<-jd1+ gamma*(v^2+u^2)*dt/2
            
            ig <- ig + 1
            
            u <- v
            y <- x
            
            if(ig==nc)
            {
                
                xx[it/nc,i]<- x
                vv[it/nc,i]<- v
                ee[it/nc,i]<- (ji1-jd1)
                ji[it/nc,i]<- ji1
                jd[it/nc,i]<- jd1
                
                ig<- 0
                
            }
        }
        
    }
    
    
    final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp)
    
    return(final)
}

SimulacaoLangevinKickF <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp, deltaf, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    
    it<- 0
    
    
    
    for(i in 1:ns)
    {
        
        print(i)
        
        x<-0
        v<-0
        u<-0
        y<-0
        
        ji1<-0
        jd1<-0
        
        
        it<-0
        ic<-0
        ir<-0
        ig<-0
        
        forceest0<-0
        forcedet0<-0
        
        indicadorkick<- 0
        
        while(it < ntc)
        {
            
            it <- it + 1
            
            forceest <- sqrt(2*gamma*temp*dt)*rnorm(1)
            
            
            #forceeststrato<- (forceest + forceest0)/2
            #forceest0<-forceest
            
            forcedet <-  (-gm*u -c1*y)*dt
            
            #forcedetstrato <- (forcedet + forcedet0)/2
            #forcedet0 <- forcedet
            
            
            
            
            force <- forcedet + forceest
            
            
            if( it*dt >= tt & indicadorkick==0)
            {
                force <- force + deltaf
                print(i)
                indicadorkick<- 1
            }
            
            v <- u + force
            x <- y + v*dt
            
            ji1 <- ji1 + (v+u)*forceest/2
            jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            # caso bimodal
            # ji1<-ji1+ (v+u)*eta*dt/2
            # jd1<-jd1+ gamma*(v^2+u^2)*dt/2
            
            ig <- ig + 1
            
            u <- v
            y <- x
            
            if(ig==nc)
            {
                
                xx[it/nc,i]<- x
                vv[it/nc,i]<- v
                ee[it/nc,i]<- (ji1-jd1)
                ji[it/nc,i]<- ji1
                jd[it/nc,i]<- jd1
                
                ig<- 0
                
            }
        }
        
    }
    
    
    final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp, "kick" = deltaf, "tkick" = tt)
    
    return(final)
}



analiseentropicaLangevin <- function(listasim ){
    
    #xai<-listasim$xa
    
    mediaE<-numeric(listasim$ns)
    
    
    for(i in 1:listasim$ns){
        mediaE[i]<-mean(listasim$energia[,i])
    }
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        hhx<-hist(listasim$posicao[,i], breaks=100, plot= FALSE)
        hhv<-hist(listasim$velocidade[,i], breaks=100, plot= FALSE)
        
        ssx<- - hhx$density*log(hhx$density)*(hhx$breaks[2]-hhx$breaks[1])
        ssx[which(is.na(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(is.na(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        
        stoth[i]<- sx + sv
    }
    
    
    
    fim<-list("S"=stoth, "E"=mediaE,  "ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt)
    
    
    return(fim)
    
    
}


