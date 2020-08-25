


SimulacaoBimodal <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
    
    tempo<-seq(dt*nc, nt*dt*nc, dt*nc)
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    it<- 0
    
    for(i in 1:ns)
    {
        
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                v<- u + force*dt
                x<- y+ v*dt
                
                
                #  u<-v
                #y<-x
                #force<-  -gm*u -c1*y + eta/massa
                #v<- u + force*dt
                #x<- y+ v*dt
                
                
                
                
                ji1 <- ji1 + (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt   #!!!!!
                
                
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ji[it/nc,i]<- ji1
                    jd[it/nc,i]<- jd1
                    ruido[it/nc,i]<- eta
                    
                    
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "ruido" = ruido, "t" = tempo)
    
    return(final)
}

SimulacaoBimodalInternal <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
    
    tempo<-seq(dt*nc, nt*dt*nc, dt*nc)
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    #velocoarse<-numeric(nc*nt)
    
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    it<- 0
    
    for(i in 1:ns)
    {
        #
        print(i)
        
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        
        force2<-0
        
        while(it < ntc)
        {
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            
            ie<- 0
            
            #force2<-0
            
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                #velocoarse[it]<-v
                
                u<-v
                y<-x
                
                # verificar FDT
                #force1<-  -c1*y + (eta*sqrt(gamma))/massa
                # for(ij in 1:it)
                #{
                #force2<- force2 * exp(- alfa * dt  )  - gm * u * exp(- alfa * dt  ) * dt
                #force02<- force2
                
                force1<-  -c1*y + eta/massa
                
                force<- force1 + force2
                
                v<- u + force*dt
                
                x<- y+ v*dt
                
                force2<- force2 * exp(- alfa * dt  )  - gm * (v+u) * exp(- alfa * dt  ) * dt/(2* alfa)
                
                # force2<- force2 * exp(- alfa * dt  )  - ( xa**2) *(v+u) * exp(- alfa * dt  ) * dt/(2* alfa)
                
                
                ji1 <- ji1 + (v+u)*eta*dt/2
                
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                
                #jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt   #!!!!!
                #jd1 <- jd1 + (force2*(v+u)/2)*dt
                #incluí o termo de massa que ta dividido lai em cima!!!
                
                jd1 <- jd1 + (massa*force2*(v+u)/2)*dt   #!!!!!
                
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1 +jd1)
                    ji[it/nc,i]<- ji1
                    jd[it/nc,i]<- jd1
                    ruido[it/nc,i]<- eta
                    
                    
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "ruido" = ruido, "t" = tempo)
    
    return(final)
}




SimulacaoBimodalTcte <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, Temp){
    
    tempo<-seq(dt*nc, nt*dt*nc, dt*nc)
    
    ntc<-nt*nc
    
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    alfai<-alfa
    
    ie<-0
    
    it<-0
    
    xai <-numeric(ns)
    
    
    for(i in 1:ns)
    {
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfai<- alfai+ dalfa
        #AQUI A LINHA QUE GARANTE TEMPERATURA CONSTANTE = 1, para g=1 cm=1 we=1 alfa=1 xa=1
        
        xa<- sqrt( Temp*(gamma*(ctelastica+alfai*(gamma+massa*alfai)))/(massa*alfai)   )
        
        xai[i]<-xa
        
        xb<- -xa
        ua<-alfai/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                
                it<-it+1
                ie<-ie+1
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                
                v<- u + force*dt
                x<- y+ v*dt
                #print(x)
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    #print(x)
                    #print(v)
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ji[it/nc,i]<- (ji1)
                    jd[it/nc,i]<- (jd1)
                    
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    
    
    
    
    #function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, Temp)
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "ji" = ji, "jd" = jd , "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa, "dalfa"= dalfa, "xa"=xai, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "t"=tempo)
    return(final)
}


SimulacaoBimodalXActe <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa){
    
    
    ntc<-nt*nc
    
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    alfai<-alfa
    
    ie<-0
    
    it<-0
    
    
    #Temp<- numeric(ns)
    
    for(i in 1:ns)
    {
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfai<- alfai+ dalfa
        
        xb<- -xa
        ua<-alfai/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                
                it<-it+1
                ie<-ie+1
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                
                v<- u + force*dt
                x<- y+ v*dt
                #print(x)
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    #print(x)
                    #print(v)
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    
    
    
    
    #function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, Temp)
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp)
    return(final)
}

SimulacaoBimodalODlite <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
    
    tempo<-seq(dt*nc, nt*dt*nc, dt*nc)
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    # ee<-matrix( nrow = nt, ncol = ns)
    
    # ji<-matrix( nrow = nt, ncol = ns)
    #jd<-matrix( nrow = nt, ncol = ns)
    #ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica
    
    gm<- gamma
    
    ie<- 0
    
    it<- 0
    
    for(i in 1:ns)
    {
        print(i)
        
        #  Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<- 0
        v<- 0
        v12<- 0
        ji1<- 0
        jd1<- 0
        alfa<- alfa + dalfa
        xa<- xa + dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        
        it<- 0
        ic<- 0
        ir<- 0
        ig<- 0
        ie<- 0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                
                dx<- - ( c1 * y  + eta)*(dt / gm)
                
                x <-  y + dx
                
                v <- dx/dt
                
                
                #ji1 <- ji1 - (v+u)*eta*dt/2 -------> original era com sinal + mas troquei pra acertar (????)
                
                #ji1 <- ji1 - (v+u)*eta*dt/2
                
                #jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig <- ig+1
                
                if(ig==nc)
                {
                    
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    #  ee[it/nc,i]<- (ji1-jd1)
                    #ji[it/nc,i]<- ji1
                    #jd[it/nc,i]<- jd1
                    #ruido[it/nc,i]<- eta
                    
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    #final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "ruido" = ruido, "t" = tempo)
    
    final<-list( "posicao"=xx, "velocidade"=vv,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt,  "t" = tempo)
    
    return(final)
}
