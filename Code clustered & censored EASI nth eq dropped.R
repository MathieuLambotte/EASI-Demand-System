#### Clustered & Censored EASI #### dispo=select(df,ends_with("Dispo"))
easi <-  function(shares=shares,dispo=dispo,log.price=log.price,var.soc=var.soc,log.exp=log.exp,y.power=FALSE,labels.share=FALSE,labels.soc=FALSE,py.inter=FALSE,zy.inter=FALSE,pz.inter=FALSE,interpz=FALSE){
    #### Setting (data transformation) ####    
    ## y.power = hightest power of y
    ## nsoc = mumber of demographics variables
    ## neq =  number of equations (without the last item)
    if (!y.power){ny <- 3} else ny <- y.power
    nsoc <- ncol(var.soc)
    neq <- ncol(shares)-1
    
    # Number of observations
    n <- length(log.exp)
    
    ## Matrix of socio-demographic variables
    z <- matrix(0,n,nsoc)
    for (i in 1:nsoc)
      z[,i] <- var.soc[,i]
    LABELS.Z <- c()
    for (i in 1:nsoc)
      LABELS.Z <-c(LABELS.Z,paste("z",i,sep=""))
    colnames(z) <- LABELS.Z
    
    # Labels or names of the budget shares: s1 - sneq
    LABELS.W <- c()
    for (i in 1:(neq+1))
      LABELS.W <-c(LABELS.W,paste("s",i,sep=""))
    
    LABELS.P <- c()
    for (i in 1:(neq+1))
      LABELS.P <-c(LABELS.P,paste("p",i,sep=""))
    
    ## Matrix of budget shares, "s", and the matrix of log.prices, "p"
    s <- matrix(0,n,neq+1)
    p <- matrix(0,n,neq+1)
    for (i in 1:(neq+1)){
      s[,i] <- shares[,i]
      p[,i] <- log.price[,i]
    }
    colnames(s) <- LABELS.W
    colnames(p) <- LABELS.P
    
    ## labels.share
    if (length(labels.share)==1){
      labels.share <- LABELS.W}
    
    ## labels.price
    labels.price <- rep(0,neq)
    for (i in 1:(neq)){
      labels.price[i] <- paste("p",labels.share[i],sep="")}
    
    ## labels.soc
    if (length(labels.soc)==1){
      labels.soc <- LABELS.Z}
    
    ## interpz
    interpz <- ifelse((length(interpz)>1),interpz,c(1:nsoc))
    
    
    ## Convergence criteria for iterated method. If model includes     
    ## interactive terms,  we also define the the convergence variable 
    ## (implicit utility or parameter)                              
    
    interact <- ((zy.inter)||(py.inter)||(pz.inter))
    if (interact){
      conv_param <- 1
      conv_y     <- 0
      conv_crit  <- 0.000001 
    } else conv_crit  <- 0.000001 
    
    #### Price normalization & interactions prices & z/y if wanted ####       
    ## Price normalisation: "np[,i]-np[,neq]" not needed, just fasten computations
    np <- matrix(0,n,neq)
    for (i in 1:neq)
      np[,i] <- p[,i]-p[,neq+1]
    
    LABELS.np <- c()
    for (i in 1:neq){
      LABELS.np <- c(LABELS.np,paste("np",i,sep=""))
    }
    colnames(np) <- LABELS.np
    
    ## Backup of price matrix "np_backup" for future reference
    np_backup <- matrix(0,n,neq)
    for (i in 1:neq)
      np_backup[,i] <- np[,i]
    tempo <- c()
    for (i in 1:neq)
      tempo <- c(tempo,paste("np",i,"_backup",sep=""))
    colnames(np_backup) <- tempo
    
    ## Initialization of matrix Ap to 0
    Ap <- matrix(0,1,neq)
    for (i in 1:neq)
      Ap[i] <- 0
    
    ## Initialization of matrix  Bp and vector pBp to 0
    if (interact){
      Bp <- matrix(0,1,neq)
      for (i in 1:(neq))
        Bp[i] <- 0
      pBp <- 0
    }
    
    ## Initialization of matrix pAp to 0
    pAp <- 0
    
    
    
    ## Creation of interaction variables "np[,i]*z[,j]"                   ##
    ## Not all variables need to be interacted.  Those we wish to         ##
    ## interact with prices will be specified by interpz=c().             ##
    ## For example, writing interpz=c(1,2) means that we wish to interact ##
    ## the first two demographic variables with the prices.               ##
    
    
    if (pz.inter){
      npz <- matrix(0,n,neq*length(interpz))
      for (i in 1:neq){
        for (j in interpz){
          npz[,match(j,interpz)+(i-1)*length(interpz)]  <- np[,i]*z[,j]}}
      
      ## Creation of the names of the interaction variables "np[,i]*z[,j]"
      tempo <- c()
      for (i in 1:neq){
        for (j in interpz){
          tempo <- c(tempo,paste("np",i,"z",j,sep=""))}}
      
      LABELS.npz <- tempo
      colnames(npz) <- tempo
      
      tempo <- c()
      for (j in interpz){
        for (i in 1:neq){
          tempo <- c(tempo,paste("np",i,"z",j,sep=""))}}
      LABELS.npz2 <- tempo
    }  else npz <- LABELS.npz <- LABELS.npz2 <- c()
    
    #### Y Stone & Y titda & interaction Y & z/p if wanted ####    
    ## Computation of  y_stone=x-p'w,
    ## and of the instrument, y_tilda=x-p'w^bar
    y_stone <- log.exp; y_tilda <- log.exp;
    
    mean_s <- matrix(0,n,neq)
    for (i in 1:(neq))
      mean_s[,i] <- mean(s[,i])
    
    for (i in 1:(neq)){
      y_tilda <- y_tilda - mean_s[,i]*p[,i]
      y_stone <- y_stone - s[,i]*p[,i]
    }
    
    ## Creation of y^i et y_inst^i
    y      <- y_stone
    y_inst <- y_tilda
    
    YY <- matrix(0,n,ny)
    Yinst <- matrix(0,n,ny)
    for (i in 1:ny){
      YY[,i] <- y^i
      Yinst[,i] <- y_inst^i
    }
    
    LABELS.YY <- c()
    LABELS.Yinst <- c()
    for (i in 1:ny){
      LABELS.YY <- c(LABELS.YY,paste("y",i,sep=""))
      LABELS.Yinst <- c(LABELS.Yinst,paste("y_inst",i,sep=""))
    }
    
    colnames(YY) <- LABELS.YY
    colnames(Yinst) <- LABELS.Yinst
    
    ## Creation of y*z and z*y_inst only if required.
    if (zy.inter){
      yz <- matrix(0,n,nsoc)
      yzinst <- matrix(0,n,nsoc)
      for (i in 1:nsoc){
        yz[,i]  <- y*z[,i]
        yzinst[,i]  <- y_inst*z[,i]}
      
      tempo <- c()
      tempo2 <- c()
      for (i in 1:nsoc){
        tempo <- c(tempo,paste("yz",i,sep=""))
        tempo2 <- c(tempo2,paste("yzinst",i,sep=""))
      }
      
      colnames(yz) <- tempo
      colnames(yzinst) <- tempo2
    } else yz <- yzinst <- c()
    
    ## Creation of  y*p and y_inst*p only if required
    if (py.inter){
      ynp <- matrix(0,n,neq)
      for (i in 1:neq)
        ynp[,i]  <- y*np[,i]
      
      tempo <- c()
      for (i in 1:neq)
        tempo <- c(tempo,paste("ynp",i,sep=""))
      
      LABELS.ynp <- tempo
      colnames(ynp) <- tempo
      
      ynpinst <- matrix(0,n,neq)
      for (i in 1:neq)
        ynpinst[,i]  <- y_inst*np[,i]
      
      tempo <- c()
      for (i in 1:neq)
        tempo <- c(tempo,paste("ynpinst",i,sep=""))
      
      colnames(ynpinst) <- tempo
    } else ynp <- ynpinst <- c()
    #### Equation formula & restrictions (symmetry ect) ####    
    # List of covariates for each equation (labels.share)
    form5 <- c("cste")
    for (i in 1:ny)
      form5 <- c(form5,colnames(YY)[i])
    form6 <- form5
    for (i in 1:nsoc)
      form6 <- c(form6,colnames(z)[i])
    form7 <- form6
    if (zy.inter){
      for (i in 1:nsoc)
        form7 <- c(form7,colnames(yz)[i])}
    form8 <- form7
    for (i in 1:neq)
      form8 <- c(form8,colnames(np)[i])
    form9 <- form8
    if (py.inter){
      for (i in 1:neq)
        form9 <- c(form9,colnames(ynp)[i])}
    form10 <- form9
    if (pz.inter){
      for (i in 1:(neq*length(interpz)))
        form10 <- c(form10,LABELS.npz2[i])
    }
    
    varlist <- form10
    
    ## Number of variables per equation       
    dim_varlist <- length(varlist)
    
    ## Constraints of symmetry on price coefficients
    TT <- matrix(0,2,neq*(neq-1)/2)
    k <- 0
    for (i in 1:(neq-1)){
      for (j in ((i+1):(neq))){
        k <- k+1
        aa <- paste("eq",i,"_np",j,"-","eq",j,"_np",i,"=0",sep="")
        TT[1,k] <- aa
      }}
    TT <- t(TT)
    TT <- TT[,1]
    
    ## Constraints of symmetry on y*p coefficients only if required
    if (py.inter){
      TT2 <- matrix(0,2,neq*(neq-1)/2)
      k <- 0
      for (i in 1:(neq-1)){
        for (j in ((i+1):(neq))){
          k <- k+1
          aa <- paste("eq",i,"_ynp",j,"-","eq",j,"_ynp",i,"=0",sep="")
          TT2[1,k] <- aa
        }}
      TT2 <- t(TT2)
      TT2 <- TT2[,1]
    } else TT2 <- c()
    
    ## Constraints of symmetry on p*z coefficients only if required
    if (pz.inter){
      TT3 <- matrix(0,2,neq*(length(interpz)*(neq-1))/2)
      k <- 0
      
      for (t in interpz){
        for (i in 1:(neq-1)){
          for (j in ((i+1):(neq))){
            k <- k+1
            aa <- paste("eq",i,"_np",j,"z",t,"-","eq",j,"_np",i,"z",t,"=0",sep="")
            TT3[1,k] <- aa
          }}}
      TT3 <- t(TT3)
      TT3 <- TT3[,1]} else TT3 <- c()
    
    ## full list of constraints
    Rmat <- c(TT,TT2,TT3)
    
    ## Creation of system of equations
    form1 <- c(0)
    for (i in 1:ny)
      form1 <- paste(form1,"+",colnames(YY)[i])
    form2 <- form1
    for (i in 1:nsoc)
      form2 <- paste(form2,"+",colnames(z)[i])
    form3 <- form2
    if (zy.inter){
      for (i in 1:nsoc)
        form3 <- paste(form3,"+",colnames(yz)[i])}
    form4 <- form3
    for (i in 1:neq)
      form4 <- paste(form4,"+",colnames(np)[i])
    form5 <- form4
    if (py.inter){
      for (i in 1:neq)
        form5 <- paste(form5,"+",colnames(ynp)[i])}
    form6 <- form5
    if (pz.inter){
      for (i in 1:(neq*length(interpz)))
        form6 <- paste(form6,"+",LABELS.npz2[i])
    }
    
    
    system <- list()
    for (i in 1:neq){
      system <- c(system, list(formula(paste(paste("eqS",i,sep=""), "<-", paste("s",i,sep=""), "~", form6))))
    }
    
    ## Creation of the list of instruments for the 3SLS estimation
    tempo <- c()
    for (i in 1:nsoc)
      tempo <- c(tempo,paste("z",i,sep=""))
    colnames(z) <- tempo
    form11 <- c(0)
    for (i in 1:nsoc)
      form11 <- paste(form11,"+",colnames(z)[i])
    form22 <- form11
    for (i in 1:neq)
      form22 <- paste(form22,"+",colnames(np)[i])
    form33 <- form22
    if (pz.inter){
      for (i in 1:(neq*length(interpz)))
        form33 <- paste(form33,"+",LABELS.npz2[i])}
    form44 <- form33
    for (i in 1:ny)
      form44 <- paste(form44,"+",colnames(Yinst)[i])
    form55 <- form44
    if (zy.inter){
      for (i in 1:nsoc)
        form55 <- paste(form55,"+",colnames(yzinst)[i])}
    form66 <- form55
    if (py.inter){
      for (i in 1:neq)
        form66 <- paste(form66,"+",colnames(ynpinst)[i])}
    
    zlist <- form66
    zlist <- paste("~",zlist)
    
    
    ## Creation of the internal database to estimate the complete system
    new.data <- new.data_dep <- data.frame(cbind(s,YY,z,yz,np,ynp,npz,Yinst,yzinst,ynpinst,np_backup))
    
    ## Initialization of vectors and criteria before creating the instruments
    y <- y_stone
    y_backup <- y_stone
    y_old <- y_stone
    y_change <- 0
    crit_test <- 1
    iter <- 0
    conv_y <- 1
    
    #### Creation of instruments (iterative linear 3sls) ####
    cat("\n","*** Please wait during the creation of final instruments... *** ","\n")
    while (crit_test>conv_crit) {
      
      iter <- iter+1
      fit3sls <- systemfit( system, "3SLS", inst = formula(zlist), data = new.data,restrict.matrix=Rmat)
      
      if (interact){  
        if (iter>1) params_old <- params}
      
      params=fit3sls$coefficients
      
      pAp <- 0
      pBp <- 0
      y_old <- y
      
      
      if (!interact){
        #********* Predicted Values ********
        pred <- predict(fit3sls,new.data)
        shat <- matrix(0,n,neq)
        for (i in 1:neq)
          shat[,i] <- pred[,i]
        
        for (i in 1:neq){
          new.data[,paste("np",i,sep="")] <- 0
        }
        #********* Predicted Values with p=0 and no interactions ********
        pred <- predict(fit3sls,new.data)
        shat_p0 <- matrix(0,n,neq)
        for (i in 1:neq)
          shat_p0[,i] <- pred[,i]
      }
      
      
      if (interact){
        #******** y^i = 1 *********
        for (j in 1:ny)
          new.data[,paste("y",j,sep="")] <- 1}
      
      if (py.inter){  
        #********* y*p = p ******** 
        for (j in 1:neq)
          new.data[,paste("ynp",j,sep="")] <- new.data[,paste("np",j,sep="")]}
      
      if (zy.inter){
        #********* y*z = z ********
        for (j in 1:nsoc)
          new.data[,paste("yz",j,sep="")] <- new.data[,paste("z",j,sep="")]}
      
      #********* Predicted Values with y=1 and interactions ********
      pred <- predict(fit3sls,new.data)
      shat_y1 <- matrix(0,n,neq)
      for (i in 1:neq)
        shat_y1[,i] <- pred[,i]
      
      #********* all p variables (p, p*z, p*y) are set to 0 ********
      for (i in 1:neq){
        new.data[,paste("np",i,sep="")] <- 0
        if (py.inter) new.data[,paste("ynp",i,sep="")] <- 0
      }
      
      if (pz.inter){
        for (i in 1:neq){
          for (j in interpz){
            new.data[,paste("np",i,"z",j,sep="")] <- 0
          }}}
      
      #********* Predicted Values with y=1 and p variables are set to 0 and interactions********
      pred <- predict(fit3sls,new.data)
      shat_y1_p0 <- matrix(0,n,neq)
      for (i in 1:neq)
        shat_y1_p0[,i] <- pred[,i]
      
      #******** p and p*z are restored  ******** p=p_backup & p*z=p_backup*z *****
      for (i in 1:neq)
        new.data[,paste("np",i,sep="")] <-new.data[,paste("np",i,"_backup",sep="")]
      
      if (pz.inter){
        for (i in 1:neq){
          for (j in interpz){
            new.data[,paste("np",i,"z",j,sep="")] <- new.data[,paste("np",i,"_backup",sep="")]*new.data[,paste("z",j,sep="")]
          }}}
      
      #******** y variables and its interactions are set to 0 ********* y^i=y*p=y*z=0 ****
      if (interact){
        for (i in 1:ny)
          new.data[,paste("y",i,sep="")] <- 0
        
        if (py.inter){
          for (i in 1:neq)
            new.data[,paste("ynp",i,sep="")] <- 0}
        
        if (zy.inter){
          for (i in 1:nsoc)
            new.data[,paste("yz",i,sep="")] <- 0}
        
        
        #******** Predicted Values with y=0 and its interactions ********
        pred <- predict(fit3sls,new.data)
        shat_y0 <- matrix(0,n,neq)
        for (i in 1:neq)
          shat_y0[,i] <- pred[,i]
        
        #*********** p variables are set to 0  *************
        for (i in 1:neq){
          new.data[,paste("np",i,sep="")] <- 0
          new.data[,paste("ynp",i,sep="")] <- 0
        }
        
        for (i in 1:neq){
          for (j in interpz){
            new.data[,paste("np",i,"z",j,sep="")] <- 0
          }}
        
        #******** Predicted Values with p=0 and its interactions ****
        pred <- predict(fit3sls,new.data)
        shat_y0_p0 <- matrix(0,n,neq)
        for (i in 1:neq)
          shat_y0_p0[,i] <- pred[,i]
      }
      
      #*********** p variables are restored : p=p_backup
      for (i in 1:neq)
        new.data[,paste("np",i,sep="")] <-new.data[,paste("np",i,"_backup",sep="")]
      
      #********** Ap, pAp, Bp, pBp ******
      if (!interact){shat_y0 <- shat; shat_y0_p0 <- shat_p0}
      Ap <- matrix(0,n,neq)
      for (i in 1:neq)
        Ap[,i] <- shat_y0[,i]-shat_y0_p0[,i]
      
      pAp <- 0
      #********** pAp=pAp+p*Ap ******
      for (i in 1:neq)
        pAp <- pAp+new.data[,paste("np",i,sep="")]*Ap[,i]
      
      if (interact){
        Bp <- matrix(0,n,neq)
        for (i in 1:neq)
          Bp[,i] <- shat_y1[,i]-shat_y1_p0[,i]-(shat_y0[,i]-shat_y0_p0[,i])
        
        pBp <- 0
        #********** pBp=pBp+p*Bp ******
        for (i in 1:neq)
          pBp <- pBp+new.data[,paste("np",i,sep="")]*Bp[,i]
      } else Bp <- pBp<- 0
      
      
      pAp <- round(1000000*pAp+0.5)/1000000
      pBp <- round(1000000*pBp+0.5)/1000000
      
      #********* Update of y *****
      y <- (y_stone+0.5*pAp)/(1-0.5*pBp)
      
      #********* Update of y^i *****
      for (i in 1:ny)
        new.data[,paste("y",i,sep="")] <- y^i
      
      if (zy.inter){
        #********* Update of y*z *****
        for (i in 1:nsoc)
          new.data[,paste("yz",i,sep="")] <- y*new.data[,paste("z",i,sep="")]}
      
      if (py.inter){
        #********* y*p=y*p_backup ***********
        for (i in 1:neq)
          new.data[,paste("ynp",i,sep="")] <-y*new.data[,paste("np",i,"_backup",sep="")]}
      
      if (pz.inter){
        #********* z*p=z*p_backup ***********
        for (i in 1:neq){
          for (j in interpz){
            new.data[,paste("np",i,"z",j,sep="")] <- new.data[,paste("np",i,"_backup",sep="")]*new.data[,paste("z",j,sep="")]
          }}}
      
      if (interact){
        #********* Update of crit_test if conv_param=1 *******
        if ((iter>1)&(conv_param==1)){
          params_change <- (params-params_old)
          crit_test <- sum(params_change^2)
        }}
      
      #********* Update of y_change ********
      y_change <- abs(y-y_old)
      
      #********* Update of crit_test if conv_y=1 *******
      if (conv_y==1) crit_test <- max(y_change)
      
      #********* Update of crit_test if conv_y=1 (no interactions) *******
      if (!interact) crit_test <- max(y_change)
      
      #********* After the first iteration, conv_param replace conv_y *******
      conv_y <- 0
      
      #********* Interface => convergence *******
      cat("iteration = ", iter, "\n")
      cat("crit_test = ", crit_test, "\n")
      
      #********* Creation of instruments  yinst yinst*p and yinst*z *****
      y_inst <- (y_tilda+0.5*pAp)/(1-0.5*pBp)
      Yinst <- matrix(0,n,ny)
      for (i in 1:ny){
        Yinst[,i] <- y_inst^i
      }
      colnames(Yinst) <- LABELS.Yinst
      if (py.inter){
        #********* ypinst = y_inst*p **********
        for (i in 1:neq)
          new.data[,paste("ynpinst",i,sep="")] <- y_inst*new.data[,paste("np",i,sep="")]}
      
      if (zy.inter){  
        #********* yzinst = y_inst*z **********
        for (i in 1:nsoc)
          new.data[,paste("yzinst",i,sep="")] <- y_inst*new.data[,paste("z",i,sep="")]
      }
      
    }
    
    
    
    #### Final IV estimation with previous instruments ####
    new.data <- new.data_dep
    
    #*** Initialization of vectors, parameters and criteria *****
    if (interact){y_old <- y} else y_old <- y_stone
    y_change <- 0
    iter <- 0
    crit_test <- 1
    
    cat("\n","*** Creation of final instruments successfully completed... ***","\n")
    cat("\n","*** Please wait during the estimation... ***","\n")
    while (crit_test>conv_crit) {
      iter <- iter+1
      fit3sls <- systemfit( system, "3SLS", inst = formula(zlist), data = new.data,restrict.matrix=Rmat)
      
      if (interact){  
        if (iter>1) params_old <- params}
      
      params=fit3sls$coefficients
      
      pAp <- 0
      pBp <- 0
      y_old <- y
      
      if (!interact){
        #********* Predicted Values ********
        pred <- predict(fit3sls,new.data)
        shat <- matrix(0,n,neq)
        for (i in 1:neq)
          shat[,i] <- pred[,i]
        
        for (i in 1:neq){
          new.data[,paste("np",i,sep="")] <- 0
        }
        
        #********* Predicted Values with p=0 and no interactions ********
        pred <- predict(fit3sls,new.data)
        shat_p0 <- matrix(0,n,neq)
        for (i in 1:neq)
          shat_p0[,i] <- pred[,i]
      }
      
      if (interact){
        #******** y^i = 1 *********
        for (j in 1:ny)
          new.data[,paste("y",j,sep="")] <- 1}
      
      if (py.inter){  
        #********* y*p = p ******** 
        for (j in 1:neq)
          new.data[,paste("ynp",j,sep="")] <- new.data[,paste("np",j,sep="")]}
      
      if (zy.inter){
        #********* y*z = z ********
        for (j in 1:nsoc)
          new.data[,paste("yz",j,sep="")] <- new.data[,paste("z",j,sep="")]}
      
      #********* Predicted Values with y=1 and interactions ********
      pred <- predict(fit3sls,new.data)
      shat_y1 <- matrix(0,n,neq)
      for (i in 1:neq)
        shat_y1[,i] <- pred[,i]
      
      #********* all p variables (p, p*z, p*y) are set to 0 ********
      for (i in 1:neq){
        new.data[,paste("np",i,sep="")] <- 0
        if (py.inter) new.data[,paste("ynp",i,sep="")] <- 0
      }
      
      if (pz.inter){
        for (i in 1:neq){
          for (j in interpz){
            new.data[,paste("np",i,"z",j,sep="")] <- 0
          }}}
      
      #********* Predicted Values with y=1 and p variables are set to 0 and interactions********
      pred <- predict(fit3sls,new.data)
      shat_y1_p0 <- matrix(0,n,neq)
      for (i in 1:neq)
        shat_y1_p0[,i] <- pred[,i]
      
      #******** p and p*z are restored  ******** p=p_backup & p*z=p_backup*z *****
      for (i in 1:neq)
        new.data[,paste("np",i,sep="")] <-new.data[,paste("np",i,"_backup",sep="")]
      
      if (pz.inter){
        for (i in 1:neq){
          for (j in interpz){
            new.data[,paste("np",i,"z",j,sep="")] <- new.data[,paste("np",i,"_backup",sep="")]*new.data[,paste("z",j,sep="")]
          }}}
      
      #******** y variables and its interactions are set to 0 ********* y^i=y*p=y*z=0 ****
      if (interact){
        for (i in 1:ny)
          new.data[,paste("y",i,sep="")] <- 0
        
        if (py.inter){
          for (i in 1:neq)
            new.data[,paste("ynp",i,sep="")] <- 0}
        
        if (zy.inter){
          for (i in 1:nsoc)
            new.data[,paste("yz",i,sep="")] <- 0}
        
        #******** Predicted Values with y=0 and its interactions ********
        pred <- predict(fit3sls,new.data)
        shat_y0 <- matrix(0,n,neq)
        for (i in 1:neq)
          shat_y0[,i] <- pred[,i]
        
        #*********** p variables are set to 0  *************
        for (i in 1:neq){
          new.data[,paste("np",i,sep="")] <- 0
          new.data[,paste("ynp",i,sep="")] <- 0
        }
        
        for (i in 1:neq){
          for (j in interpz){
            new.data[,paste("np",i,"z",j,sep="")] <- 0
          }}
        
        #******** Predicted Values with p=0 and its interactions ****
        pred <- predict(fit3sls,new.data)
        shat_y0_p0 <- matrix(0,n,neq)
        for (i in 1:neq)
          shat_y0_p0[,i] <- pred[,i]
      }
      
      #*********** p variables are restored : p=p_backup
      for (i in 1:neq)
        new.data[,paste("np",i,sep="")] <-new.data[,paste("np",i,"_backup",sep="")]
      
      #********** Ap, pAp, Bp, pBp ******
      if (!interact){shat_y0 <- shat; shat_y0_p0 <- shat_p0}
      Ap <- matrix(0,n,neq)
      for (i in 1:neq)
        Ap[,i] <- shat_y0[,i]-shat_y0_p0[,i]
      
      pAp <- 0
      #********** pAp=pAp+p*Ap ******
      for (i in 1:neq)
        pAp <- pAp+new.data[,paste("np",i,sep="")]*Ap[,i]
      
      if (interact){
        Bp <- matrix(0,n,neq)
        for (i in 1:neq)
          Bp[,i] <- shat_y1[,i]-shat_y1_p0[,i]-(shat_y0[,i]-shat_y0_p0[,i])
        
        pBp <- 0
        #********** pBp=pBp+p*Bp ******
        for (i in 1:neq)
          pBp <- pBp+new.data[,paste("np",i,sep="")]*Bp[,i]
      } else Bp <- pBp<- 0
      
      
      pAp <- round(1000000*pAp+0.5)/1000000
      pBp <- round(1000000*pBp+0.5)/1000000
      
      #********* Update of y *****
      y <- (y_stone+0.5*pAp)/(1-0.5*pBp)
      
      #********* Update of y^i *****
      for (i in 1:ny)
        new.data[,paste("y",i,sep="")] <- y^i
      
      if (zy.inter){
        #********* Update of y*z *****
        for (i in 1:nsoc)
          new.data[,paste("yz",i,sep="")] <- y*new.data[,paste("z",i,sep="")]}
      
      if (py.inter){
        #********* y*p=y*p_backup ***********
        for (i in 1:neq)
          new.data[,paste("ynp",i,sep="")] <-y*new.data[,paste("np",i,"_backup",sep="")]}
      
      if (pz.inter){
        #********* z*p=z*p_backup ***********
        for (i in 1:neq){
          for (j in interpz){
            new.data[,paste("np",i,"z",j,sep="")] <- new.data[,paste("np",i,"_backup",sep="")]*new.data[,paste("z",j,sep="")]
          }}}
      
      if (interact){
        #********* Update of crit_test if conv_param=1 *******
        if ((iter>1)&(conv_param==1)){
          params_change <- (params-params_old)
          crit_test <- sum(params_change^2)
        }}
      
      #********* Update of y_change ********
      y_change <- abs(y-y_old)
      
      #********* Update of crit_test if conv_y=1 *******
      if (conv_y==1) crit_test <- max(y_change)
      
      #********* Update of crit_test if conv_y=1 (no interactions) *******
      if (!interact) crit_test <- max(y_change)
      
      #********* After the first iteration, conv_param replace conv_y *******
      conv_y <- 0
      
      #********* Interface => convergence *******
      cat("iteration = ", iter, "\n")
      cat("crit_test = ", crit_test, "\n")
      
    }
    cat("\n","*** Estimation successfully completed ***","\n")
    #### Sample Selection ####
    
    ## Matrix of dispo
    LABELS.D <- c()
    for (i in 1:(neq+1))
      LABELS.D <-c(LABELS.D,paste("d",i,sep=""))
    
    
    d<- matrix(0,n,neq+1)
    for (i in 1:(neq+1)){
      d[,i] <- dispo[,i]
    }
    colnames(d) <- LABELS.D
    
    ## Formula probits
    form1 <- c("0")
    for (i in 1:nsoc)
      form1 <- paste(form1,"+",colnames(z)[i])
    form_probit <- list()
    for (i in 1:(neq+1)){
      form_probit <- c(form_probit, list(formula(paste(paste("eqS",i,sep=""), "<-", paste("s",i,sep=""), "~", paste(form1,paste0("d",i),sep="+")))))
    }
    
    #Probits
    selection.data<-data.frame(cbind(s,z,d))
    for (i in 1:(neq+1)){
      for (k in 1:n){
        selection.data[k,i]<-ifelse(selection.data[k,i]>0,1,0)
      }
    }
    #export the share of censored budjet shares per products' family
    #amount_zero<-as.data.frame(matrix(nrow=neq,ncol=2))
    ######for (i in 1:neq){
    #####  amount_zero[i,]<-table(as.factor(selection.data[,i]))/n*100
   #### }
   ### rownames(amount_zero)<-labels.share
   ## colnames(amount_zero)<-levels(as.factor(selection.data[,i]))
    #write.xlsx(amount_zero,file=paste(Sys.Date(),"amout_zero.xlsx"))
    
    
    #heteroske probits
    probit<-list()
    for (i in 1:(neq+1)){
      probit[[i]]<-glm(form_probit[[i]],selection.data,family=binomial(link="probit"))}
    for( i in 1:(neq+1)){
      selection.data[,i]<-as.factor(as.character(selection.data[,i]))
    }
  
    
    #export the probit results in case
    #stargazer(probit,column.labels=str_replace(labels.share,"_"," "),column.separate=c(rep(1,14)),
          #    covariate.labels=str_replace_all(c(labels.soc,LABELS.D),"_"," "),title="Probits_Results",type="html",out=paste(Sys.Date(),"Het_Probits_Results.htm"),align=TRUE)
    #PDF & CDF from probits
    pdf<-as.data.frame(matrix(ncol=(neq+1),nrow=n))
    cdf<-as.data.frame(matrix(ncol=(neq+1),nrow=n))
    for (i in 1:(neq+1)){
      pdf[,i]<-dnorm(predict(probit[[i]],type="link"))
      cdf[,i]<-pnorm(predict(probit[[i]],type="link"))
    }
    
    #Label pdf & cdf
    LABELS.PDF <- c()
    for (i in 1:(neq+1))
      LABELS.PDF <-c(LABELS.PDF,paste("pdf",i,sep=""))
    
    LABELS.CDF <- c()
    for (i in 1:(neq+1))
      LABELS.CDF <-c(LABELS.CDF,paste("cdf",i,sep=""))
    
    colnames(pdf)<-LABELS.PDF    
    colnames(cdf)<-LABELS.CDF   
    new.data<-cbind(new.data,pdf,cdf)
    
    #### final estimate  ####
    ## Constraints of symmetry on price coefficients
    TT <- matrix(0,2,(neq)*(neq-1)/2)
    k <- 0
    for (i in 1:(neq-1)){
      for (j in ((i+1):(neq))){
        k <- k+1
        aa <- paste("eq",i,"_cdf",i,":np",j,"-","eq",j,"_cdf",j,":np",i,"=0",sep="")
        TT[1,k] <- aa
      }}
    TT <- t(TT)
    TT <- TT[,1]
    
    
    
    
    ## Constraints of symmetry on y*p coefficients only if required
    if (py.inter){
      TT2 <- matrix(0,2,neq*(neq-1)/2)
      k <- 0
      for (i in 1:(neq-1)){
        for (j in ((i+1):neq)){
          k <- k+1
          aa <- paste("eq",i,"_ynp",j,"-","eq",j,"_ynp",i,"=0",sep="")
          TT2[1,k] <- aa
        }}
      TT2 <- t(TT2)
      TT2 <- TT2[,1]
    } else TT2 <- c()
    
    ## Constraints of symmetry on p*z coefficients only if required
    if (pz.inter){
      TT3 <- matrix(0,2,(neq-1)*(length(interpz)*(neq))/2)
      k <- 0
      
      for (t in interpz){
        for (i in 1:(neq-1)){
          for (j in ((i+1):(neq))){
            k <- k+1
            aa <- paste("eq",i,"_np",j,"z",t,"-","eq",j,"_np",i,"z",t,"=0",sep="")
            TT3[1,k] <- aa
          }}}
      TT3 <- t(TT3)
      TT3 <- TT3[,1]} else TT3 <- c()
    
    
    ## full list of constraints
    Rmat <- c(TT,TT2,TT3)
    
    form1 <- c("y1")
    for (i in 2:ny)
      form1 <- paste(form1,"+",colnames(YY)[i])
    form2 <- form1
    for (i in 1:nsoc)
      form2 <- paste(form2,"+",colnames(z)[i])
    form3 <- form2
    if (zy.inter){
      for (i in 1:nsoc)
        form3 <- paste(form3,"+",colnames(yz)[i])}
    form4 <- form3
    for (i in 1:neq)
      form4 <- paste(form4,"+",colnames(np)[i])
    form5 <- form4
    if (py.inter){
      for (i in 1:neq)
        form5 <- paste(form5,"+",colnames(ynp)[i])}
    form6 <- form5
    if (pz.inter){
      for (i in 1:(neq*length(interpz)))
        form6 <- paste(form6,"+",LABELS.npz2[i])
    }
    
    
    censored_system <- list()
    for (i in (1:(neq))){
      censored_system <- c(censored_system, list(formula(paste(paste("eqS",i,sep=""), "<-", paste("s",i,sep=""), "~",paste0("0 +cdf",i,":(",form6,")"),paste0("+ pdf",i)))))
    }
    
    
    form11 <- c("z1")
    for (i in 2:nsoc)
      form11 <- paste(form11,"+",colnames(z)[i])
    form22 <- form11
    for (i in 1:neq)
      form22 <- paste(form22,"+",colnames(np)[i])
    form33 <- form22
    if (pz.inter){
      for (i in 1:(neq*length(interpz)))
        form33 <- paste(form33,"+",LABELS.npz2[i])}
    form44 <- form33
    for (i in 1:ny)
      form44 <- paste(form44,"+",colnames(Yinst)[i])
    form55 <- form44
    if (zy.inter){
      for (i in 1:nsoc)
        form55 <- paste(form55,"+",colnames(yzinst)[i])}
    form66 <- form55
    if (py.inter){
      for (i in 1:(neq+1))
        form66 <- paste(form66,"+",colnames(ynpinst)[i])}
    
    zlist_censored<-list()
    for (i in 1:(neq)){
      zlist_censored<-c(zlist_censored,list(formula(paste("~",paste0("0 + cdf",i,":(",form66,") +pdf",i)))))
    }
    
    fit3sls <- systemfit(censored_system, "3SLS", inst = (zlist_censored), data = new.data,restrict.matrix=Rmat,useMatrix=T)
    
  
 
  
  
  ## Calculation of w_j
    W=matrix(0,n,neq)
    for (i in 1:neq) W[,i] <- predict(fit3sls)[,i]
    
    colnames(W) <- labels.share[1:neq]
    
    
    ## Preparation of the display of results
    VARS=c("pdf1")
    for (i in 1:ny)
      VARS <- c(VARS,paste("y^",i,sep=""))
    for (i in 1:nsoc)
      VARS <- c(VARS,labels.soc[i])
    if (zy.inter)
      for (i in 1:nsoc)
        VARS <- c(VARS,paste("y*",labels.soc[i],sep=""))
    for (i in 1:neq)
      VARS <- c(VARS,labels.price[i])
    if (py.inter){
      for (i in 1:neq)
        VARS <- c(VARS,paste("y*",labels.price[i],sep=""))}
    if (pz.inter){
      for (i in interpz){
        for (j in 1:neq){
          VARS <- c(VARS,paste(labels.soc[i],"*",labels.price[j],sep=""))
        }}
    }
    
    VARS2 <- c("~")
    for (i in 1:length(VARS))
      VARS2 <- paste(VARS2,"+",VARS[i])
    
    tempo <- c()
    for (i in 1:nsoc)
      tempo <- c(tempo,paste("z",i,sep=""))
    colnames(z) <- tempo
    form11 <- c()
    for (i in 1:nsoc)
      form11 <- c(form11,labels.soc[i])
    form22 <- form11
    for (i in 1:neq)
      form22 <- c(form22,labels.price[i])
    form33 <- form22
    if (pz.inter){
      for (i in interpz){
        for (j in 1:neq){
          form33 <- c(form33,paste(labels.soc[i],"*",labels.price[j],sep=""))
        }}}
    form44 <- form33
    for (i in 1:ny)
      form44 <- c(form44,paste("inst_y^",i,sep=""))
    form55 <- form44
    if (zy.inter){
      for (i in 1:nsoc)
        form55 <- c(form55,paste("inst_y*",labels.soc[i],sep=""))}
    form66 <- form55
    if (py.inter){
      for (i in 1:neq)
        form66 <- c(form66,paste("y*",labels.price[i],sep=""))}
    
    VARINST <- form66
    
    VARINST2 <- c("~")
    for (i in 1:length(VARINST))
      VARINST2 <- paste(VARINST2,"+",VARINST[i])
    
    a=summary(fit3sls)
    for (i in 1:neq){
      a$eq[[i]][1] <- labels.share[i]
      rownames(a$eq[[i]][8][[1]]) <- VARS
      a$eq[[i]][3] <- VARS2
      a$eq[[i]][4] <- VARINST2
    }
    
    Residuals=summary(fit3sls)$residuals
    colnames(Residuals) <- labels.share[1:neq]
    
    
    
    
    # The EASI function returns "Result", an object of class "easi", namely  #
    # a list of results that includes:                                       # 
    # - Residuals                                                            #
    # - The covariance matrix of estimated parameters                        #
    # - The matrix of the fitted budget shares                               #
    # - a summary of the result (a)				                                   #
    # - the implicit utility (y)					                                   #
    # - the list of all variables                                            #
    # - others objects used in the estimation : var.soc (demographic varia-  #
    # bles), shares (matrix of budget shares), log.price (matrix of prices), #
    # neq (number of equations without the last equation), y.power (hightest #
    # power of y), nsoc (number of demographic variables), interact (a dum-  #
    # my equal to one if there are interaction), py.inter (a dummy equal to  #
    # one if there are interactions between Prices and y, pz.inter (a dummy  #
    # equal to one if there are interactions between Prices and demographic  #
    # variables, zy.inter (a dummy equal to one if there are interactions    #
    # between demographics and y, interpz (that allows to choose the varia-  #
    # bles Z to crossed with the price), fit3sls (object of class systemfit) #
    # log.exp (the logaritm of total expenditure), labels.price (the names of#
    # prices), labels.soc (the names of demographic variables), labels (the  #
    # names of budget shares), dim_varlist (nombre de variables)             #
    
    
    Result <-  list(
      Residuals=Residuals,
      CoefCov=summary(fit3sls)$coefCov,
      fitted.w=W,
      summary=a,
      y=y,
      varlist=varlist,
      var.soc=var.soc,
      shares=shares,
      log.price=log.price,
      neq=neq,
      y.power=ny,
      nsoc=nsoc,
      interact=interact,
      py.inter=py.inter,
      zy.inter=zy.inter,
      pz.inter=pz.inter,
      interpz=interpz,
      fit3sls=fit3sls,
      log.exp=log.exp,
      labels.price=labels.price,
      labels.soc=labels.soc,
      labels.share=labels.share,
      dim_varlist=dim_varlist,
      VARS=VARS,
      new.data=new.data
    
    )
    class(Result) <- c("easi","systemfit")
    Result
}
#### SY elasticities (Hicksian/compensated) ####   Castellon 2012  

  SY_hicksian<-function(object=object){
    fit3sls <- object$fit3sls
    varlist <- object$varlist
    var.soc <- object$var.soc
    shares <- object$shares
    log.price <- object$log.price
    neq <- object$neq
    y.power <- object$y.power
    nsoc <- object$nsoc
    interact <- object$interact
    py.inter <- object$py.inter
    zy.inter <- object$zy.inter
    pz.inter <- object$pz.inter
    interpz <- object$interpz
    log.exp <- object$log.exp
    labels.price <- object$labels.price
    labels.soc <- object$labels.soc
    labels.share <- object$labels.share
 coef<-fit3sls$coefficients
 cdf<-Result$new.data[,grepl("cdf",colnames(Result$new.data))]
 #get the neq-1 x neq-1 directly estimated elasticites
     Comp_Elas<-matrix(0,neq,neq)
    for(i in 1:neq){
      for(j in 1:neq){
        Comp_Elas[i,j]<-mean(cdf[,i])*coef[paste0("eq",i,"_cdf",i,":np",j)]/mean(shares[,i]) +mean(shares[,j])-as.numeric((i==j))
      }
    }
     
     #use budgjet constraint to get the last equation's elastictiies (sum wi*eij=0)
    dropped<-matrix(0,1,neq)
    for(i in 1:neq){
      for(j in 1:neq){
        dropped[,i]<-0-sum(mean(shares[,j])*Comp_Elas[j,i])
      }
    }
    
    temp<-rbind(Comp_Elas,dropped)
    
  
    sum_col<-matrix(0,neq,1)
    for(i in 1:neq){
      for(j in 1:(neq+1)){
        sum_col[i,]<-sum(mean(shares[,j])*temp[j,i])
      }
    }
    
    dropped2<-matrix(0,neq+1,1)
    for(i in 1:(neq+1)){
      for(j in 1:neq){
        dropped2[i,]<-sum(mean(shares[,j])*temp[i,j])
      }
    }
    temp2<-cbind(temp,dropped2)
    sum_row<-matrix(0,neq+1,1)
    for(i in 1:(neq+1)){
      for(j in 1:(neq+1)){
        sum_row[i,]<-sum(mean(shares[,j])*temp2[i,j])
      }
    }
      rownames(temp2)<-labels.share
      colnames(temp2)<-labels.share
    res<-list(temp2,sum_row,sum_col)
    return(res)
}
#### SY Income elasticties #### Castellon 2015
  SY_income<-function(object=object){
    fit3sls <- object$fit3sls
    varlist <- object$varlist
    var.soc <- object$var.soc
    shares <- object$shares
    log.price <- object$log.price
    neq <- object$neq
    y.power <- object$y.power
    nsoc <- object$nsoc
    interact <- object$interact
    py.inter <- object$py.inter
    zy.inter <- object$zy.inter
    pz.inter <- object$pz.inter
    interpz <- object$interpz
    log.exp <- object$log.exp
    labels.price <- object$labels.price
    labels.soc <- object$labels.soc
    labels.share <- object$labels.share
    coef<-fit3sls$coefficients 
    np<-object$log.price
    cdf<-Result$new.data[,grepl("cdf",colnames(Result$new.data))]
    
    
       
    #get the neq-1  directly estimated elasticites (castellon 2015 p11)
  bjr=matrix(0,y.power,neq)
  for (i in 1:neq){
    for (j in 1:y.power){
      bjr[j,i] <- coef[paste("eq",i,"_cdf",i,":y",j,sep="")]
    }}
  sum_bjr<-matrix(0,1,neq)
  for (i in 1:(neq)){
    for (j in 1:y.power){
      tempo <- j*bjr[j,i]*mean(y)^{j-1}
      sum_bjr[,i] <- sum_bjr[,i]+tempo}}   

  
Elas_exp<-matrix(0,neq,1)
for(i in 1:neq){
  Elas_exp[i,]<-(1/mean(shares[,i]))*((mean(cdf[,i])*sum_bjr[,i]))+1
}    

# get last one income elasticity (yen 2003 footnote 9)
last_one<-matrix(0,neq,1)
for (i in 1:neq){
  last_one[i,]<-mean(shares[,i])*Elas_exp[i]
}
final_one<-c(last_one,(1-sum(last_one)))
final_one<-final_one/colMeans(shares)
sum<-sum(final_one*colMeans(shares))
final<-list(final_one,sum)
return(final)
}
  #### SY elasticities (Marshallian) ####   slutsky equation 
  SY_marshallian<-function(object=object,elas_hick=elas_hick,income_elas=income_elas){
    
    elas_hick<-elas_hick[[1]]
    shares <- object$shares
    labels.share<-object$labels.share
    neq<-object$neq+1
    income_elas<-income_elas[[1]]
    
    elas_marshall<-matrix(0,neq,neq) 
    for(i in 1:neq){
      for(j in 1:neq){
        elas_marshall[i,j]<-elas_hick[i,j]-mean(shares[,i])*income_elas[i] 
      }
    }
    
   
    rownames(elas_marshall)<-labels.share
    colnames(elas_marshall)<-labels.share
    res<-list(elas_marshall,sum_row,sum_col)
    return(res)
  }
#### fonctions annexes ####
print.easi<-function(object){
  cat("\n","*** systemfit results ***","\n")
  cat("\n",object$summary$method,"\n")
  cat("\n","*** Coefficient ***","\n")
  print(coef(object))
}
coef.easi <-
  function(object=object,...){
    neq <- object$neq
    VARS <- object$VARS
    VARS2 <- c()
    for (i in 1:neq)
      VARS2 <- c(VARS2,paste(paste("eq",i,sep=""),VARS,sep="_"))
    tp <- object$summary$coefficients
    rownames(tp) <- VARS2
    tp
  }
vcov.easi <-
  function(object=object,...){
    neq <- object$neq
    VARS <- object$VARS
    VARS2 <- c()
    for (i in 1:neq)
      VARS2 <- c(VARS2,paste(paste("eq",i,sep=""),VARS,sep="_"))
    tp <- object$CoefCov
    colnames(tp) <- rownames(tp) <- VARS2
    tp
  }
residuals.easi <-
  function(object=object,...){
    tp <- object$summary$residuals
    colnames(tp) <- object$labels.share[1:object$neq]
    tp
  }
summary.easi <-
  function(object=object,...){
    tp <- object$summary
    tp
  }
predict.easi <-
  function(object=object,...){
    tp <- object$fitted.w
    tp
  }
#### intermediate.blocs ####
intermediate.blocs <-
  function(object=object,log.price=log.price,var.soc=var.soc,log.exp=log.exp){
    
    shares <- object$shares
    fit3sls <- object$fit3sls
    varlist <- object$varlist
    neq <- object$neq
    y.power <- object$y.power
    nsoc <- object$nsoc
    interact <- object$interact
    py.inter <- object$py.inter
    zy.inter <- object$zy.inter
    pz.inter <- object$pz.inter
    interpz <- object$interpz
    log.exp <- object$log.exp
    labels.price <- object$labels.price
    labels.soc <- object$labels.soc
    labels.share <- object$labels.share
    dim_varlist <- object$dim_varlist
    y <- object$y
    
    n <- length(log.exp)
    nsoc <- ncol(var.soc)
    neq <- ncol(shares)-1
    
    # Number of observations
    n <- length(log.exp)
    
    ## Matrix of socio-demographic variables
    z <- matrix(0,n,nsoc)
    for (i in 1:nsoc)
      z[,i] <- var.soc[,i]
    LABELS.Z <- c()
    for (i in 1:nsoc)
      LABELS.Z <-c(LABELS.Z,paste("z",i,sep=""))
    colnames(z) <- LABELS.Z
    var.soc=z
    #** Recovery of the estimated coefficients***
    Estimates <- summary(fit3sls)$coefficients
    
    coef=Estimates[,1]
    n <- length(log.exp)
    
    
    ## Labels or names of the equations:
    noms <- c()
    for (i in 1:neq)
      noms <- c(noms, paste("eq",i,sep=""))
    
    #***** new price Matrix *******
    P=log.price
    
    
    #****  new sociodemographic matrix *****
    Z=cbind(1,z)
    
    #**** logarithm of total expenditure *******
    lnx <- log.exp
    
    #**** Budget shares matrix *******
    w=matrix(0,n,neq+1)
    for (i in 1:(neq+1))
      w[,i] <- shares[,i]
    
    
    ## Recovery of coefficients for the variables p and p * z Note: the first element of Z is a constant to capture the direct price effects
    
    my.array <-array (0, dim = c (nsoc+1,(neq+1),(neq+1)))
    
    ## Step 1: Recovery of coefficients of p variables
    a0 <- matrix(0,neq,neq)
    for (i in 1:neq){
      for (j in 1:neq){
        a0[i,j] <- coef[paste("eq",i,"_cdf",i,":np",j,sep="")]
        my.array[,,i][1,j] <- a0[i,j]
      }}
    my.array[,,neq+1][1,1:neq] <- 0-apply(a0,2,sum)
    
    for (i in 1:neq){
      my.array[,,i][1,neq+1] <- my.array[,,neq+1][1,i] 
    }
    my.array[,,neq+1][1,neq+1] <- 0-sum(my.array[,,neq+1][1,1:neq])
    
    
    
    ## Step 2: Recovery of coefficients of p*z variables only if required
    if (pz.inter){
      for (i in 1:neq){
        for (j in interpz){
          for (k in 1:neq){
            my.array[,,i][j+1,k] <- coef[paste("eq",i,"_np",k,"z",j,sep="")]
          }}}
    }
    
    for (t in interpz){
      for (i in 1:neq){
        for (j in 1:neq){
          my.array[,,neq+1][t+1,i] <- my.array[,,neq+1][t+1,i]-my.array[,,j][t+1,i]
        }
      }
    }
    
    for (t in interpz){
      for (i in 1:neq){
        my.array[,,i][t+1,neq+1] <- my.array[,,neq+1][t+1,i]
      }
    }
    
    for (t in interpz){
      for (i in 1:neq){
        my.array[,,neq+1][t+1,neq+1] <- my.array[,,neq+1][t+1,neq+1]-my.array[,,i][t+1,neq+1]
      }
    }
    
    ## construction of the sum "sum_j sum_k sum_t a_jkt z_t p_j p_k" (calculation of y)
    ## "EASI made EASIER" (PENDAKUR 2008 -  page 11 formula 22)
    a <- my.array
    tot = 0
    for (j in 1:neq){
      for (k in 1:neq){
        for (t in c(1:(nsoc))){
          tempo <- a[t,k,j]*P[,k]*P[,j]*Z[,t]
          tot <- tot+tempo
        }}}
    
    ## Recovery of coefficients of p*y variables only if required
    bjk=matrix(0,neq+1,neq+1)
    tot2 = 0
    if (py.inter){
      for (i in 1:neq){
        for (j in 1:neq){
          bjk[j,i] <- coef[paste("eq",i,"_ynp",j,sep="")]
        }}
      
      for (j in 1:neq){
        bjk[j,neq+1] <- 0-sum(bjk[j,1:neq])
      }
      
      for (j in 1:neq){
        bjk[neq+1,j] <- bjk[j,neq+1]
      }
      
      bjk[neq+1,neq+1] <- 0-sum(bjk[neq+1,1:neq])
      
      colnames(bjk)<- c(noms,"Others")
      
      ## construction of the sum "sum_j sum_k  b_jk p_j p_k" (calculation of y)
      ## "EASI made EASIER" (PENDAKUR 2008 -  page 11 formula 22)
      for (j in 1:neq){
        for (k in 1:neq){
          tempo <- bjk[j,k]*P[,j]*P[,k]
          tot2 <- tot2+tempo
        }}
    }
    
    ## construction of the sum "sum_j w_j p_j" (calculation of y)
    ## "EASI made EASIER" (PENDAKUR 2008 -  page 11 formula 22)
    tot0 = 0
    for (j in 1:neq){
      tempo <- w[,j]*P[,j]
      tot0 <- tot0+tempo
    }
    
    ## Calculation of y
    ## "EASI made EASIER" (PENDAKUR 2008 -  page 11 formula 22)
    if (interact){
      y <- (lnx - tot0+1/2*tot)/(1-1/2*tot2)} else y <- (lnx-tot0+1/2*tot)
    
    ## Recovery of coefficients of y^r variables (calculation of w_j)
    bjr=matrix(0,y.power,neq+1)
    for (i in 1:neq){
      for (j in 1:y.power){
        bjr[j,i] <- coef[paste("eq",i,"_cdf",i,":y",j,sep="")]
      }}
    
    for (j in 1:y.power){
      bjr[j,neq+1] <- 0-sum(bjr[j,1:neq])
    }
    
    
    colnames(bjr)<- c(noms,"others")
    
    ## Recovery of coefficients of z variables (calculation of w_j)
    gjt=matrix(0,nsoc,neq+1)
    for (i in 1:neq){
      for (j in 1:nsoc){
        gjt[j,i] <- coef[paste("eq",i,"_cdf",i,":z",j,sep="")]
      }}
    
    for (j in 1:nsoc){
      gjt[j,neq+1] <- 0-sum(gjt[j,1:neq])
    }
    
    colnames(gjt)<- c(noms,"others")
    
    
    ## Recovery of coefficients of y*z variables (calculation of w_j) only if required
    hjt=matrix(0,nsoc,neq+1)
    if (zy.inter){
      for (i in 1:neq){
        for (j in 1:nsoc){
          hjt[j,i] <- coef[paste("eq",i,"_yz",j,sep="")]
        }}
      
      for (j in 1:nsoc){
        hjt[j,neq+1] <- 0-sum(hjt[j,1:neq])
      }
      
    }
    colnames(hjt) <- c(noms,"others")
    
    ## Recovery of the constants (calculation of w_j)
    cc=c()
    for (i in 1:neq)
      cc=cbind(cc,coef[paste("eq",i,"_cdf",i,":z1",sep="")])
    
    cc=cbind(cc,1-sum(cc))
    
    
    
    Result <-  list(
      CoefCov=fit3sls$coefCov,
      a=a,
      y=y,
      varlist=varlist,
      var.soc=var.soc,
      shares=shares,
      log.price=log.price,
      neq=neq,
      y.power=y.power,
      nsoc=nsoc,
      interact=interact,
      py.inter=py.inter,
      zy.inter=zy.inter,
      pz.inter=pz.inter,
      interpz=interpz,
      fit3sls=fit3sls,
      log.exp=log.exp,
      labels.price=labels.price,
      labels.soc=labels.soc,
      labels.share=labels.share,
      dim_varlist=dim_varlist,
      n=n,
      coef=coef,
      my.array=my.array,
      tot=tot,
      tot2=tot2,
      tot0=tot0,
      bjk=bjk,
      P=P,
      w=w,
      Z=Z,
      bjr=bjr,
      gjt=gjt,
      hjt=hjt,
      cc=cc,
      lnx=lnx
      
    )
    
    Result
    
  }

#### elastic ####
elastic <-
  function(object=object,type = c("price", "income", "demographics"),sd=FALSE){
    
    
    type <- match.arg(type)
    
    EPDELTA <- ifelse(((type=="price")&(sd==TRUE)),TRUE,FALSE)
    ERDELTA <- ifelse(((type=="income")&(sd==TRUE)),TRUE,FALSE)
    EZDELTA <- ifelse(((type=="demographics")&(sd==TRUE)),TRUE,FALSE)
    
    
    fit3sls <- object$fit3sls
    varlist <- object$varlist
    var.soc <- object$var.soc
    shares <- object$shares
    log.price <- object$log.price
    neq <- object$neq
    y.power <- object$y.power
    nsoc <- object$nsoc
    interact <- object$interact
    py.inter <- object$py.inter
    zy.inter <- object$zy.inter
    pz.inter <- object$pz.inter
    interpz <- object$interpz
    log.exp <- object$log.exp
    labels.price <- object$labels.price
    labels.soc <- object$labels.soc
    labels.share <- object$labels.share
    
    
    n <- length(log.exp)
    
    temp <- intermediate.blocs(object,log.price=log.price,var.soc=var.soc,log.exp=log.exp)
    my.array <- temp$my.array
    tot <- temp$tot
    tot2 <- temp$tot2
    tot0 <- temp$tot0
    bjk <- temp$bjk
    P <- temp$P
    w <- temp$w
    Z <- temp$Z
    bjr <- temp$bjr
    gjt <- temp$gjt
    hjt <- temp$hjt
    cc <- temp$cc
    noms <- object$noms
    lnx <- object$log.exp
    y <- object$y
    
    
    if (type == "price"){
      ### Calculation of log.price elasticities
      
      #*** semi-elasticities with respect to log.prices
      ## page 13 formula 23 EASI made EASIER (Pendakur 2008)
      EP=matrix(0,neq+1,neq+1)
      a <- my.array 
      for(i in 1:(neq+1)){
        for (k in 1:(neq+1)){
          tot10 = 0
          for (t in (1:(nsoc+1))){
            tempo <- a[t,k,i]*Z[,t]
            tot10 <- tot10+tempo
          }
          tot10 <- tot10+bjk[k,i]*y
          EP[k,i]=mean(tot10)
        }}
      
      
      colnames(EP) <- rownames(EP) <- c(labels.price[1:neq],"pothers")
      
      
      ## Calculation of standard deviations of log.price elasticities
      ## (Delta method) if EPDELTA=TRUE 
      if (EPDELTA){
        ttt <- colnames(summary(fit3sls)$coefCov)
        EP_SE=matrix(0,neq+1,neq+1)
        ELASTPRICE_SE=matrix(0,neq+1,neq+1)
        for (i in 1:neq){
          for (j in 1:neq){
            tt <- paste("eq",i,"_cdf",i,":np",j,sep="")
            
            if (pz.inter){
              for (t in interpz){
                tt <- c(tt,paste("eq",i,"_cdf",i,":np",j,"z",t,sep=""))
              }}
            
            if (py.inter){
              tt <- c(tt,paste("eq",i,"_cdf",i,":ynp",j,sep=""))
            }
            
            
            tnum <- match(tt,ttt)
            
            DD <- summary(fit3sls)$coefCov[tnum,tnum]
            
            MAT <- Z[,1]
            if (pz.inter)
              MAT <- cbind(MAT,Z[,interpz+1])
            if (py.inter)
              MAT <- cbind(MAT,y)
            
            
            EP_SE[i,j] <- median(sqrt(diag(as.matrix(MAT)%*%as.matrix(DD)%*%t(as.matrix(MAT)))))
            ELASTPRICE_SE[i,j] <- EP_SE[i,j]/mean(shares[,j])
            rm(DD)}}
        
        
        EP_SE[neq+1,1:neq] <- sqrt(apply(EP_SE[1:neq,1:neq]^2,1,sum))
        EP_SE[,neq+1] <- EP_SE[neq+1,]
        EP_SE[neq+1,neq+1] <- sqrt(sum(EP_SE[neq+1,1:neq]))
        
        for (i in (1:(neq+1))){
          ELASTPRICE_SE[i,neq+1] <- EP_SE[i,neq+1]/mean(shares[,neq+1])}
        
        ELASTPRICE_SE[neq+1,1:neq] <- ELASTPRICE_SE[1:neq,neq+1]
        ELASTPRICE_SE[neq+1,neq+1] <- sqrt(sum(ELASTPRICE_SE[1:neq,neq+1]))
        
        
        colnames(ELASTPRICE_SE) <- rownames(ELASTPRICE_SE) <- colnames(EP_SE) <- rownames(EP_SE) <- c(labels.price[1:neq],"pothers")
        
      }
      
      
      
      #*** Normalised Slutsky matrix
      #** matrix of compensated quantity derivatives with respect to unlogged log.prices)
      #** own-log.price Slutsky terms in Pendakur
      ## page 849  "Tricks with Hicks : The EASI demand system" (Lewbel & Pendakur 2008)
      EPS=EP+apply(shares[,1:(neq+1)],2,mean)%*%t(apply(shares[,1:(neq+1)],2,mean))-matrix(diag(apply(shares[,1:(neq+1)],2,mean)),neq+1,neq+1)
      colnames(EPS) <- rownames(EPS) <- c(labels.price[1:neq],"pothers")
      
      
      #*** Compensated (good-specific) expenditures elasticities with respect to log.prices
      #** own-log.price Quant elast in Pendakur
      ## page 849 "Tricks with Hicks : The EASI demand system" (Lewbel & Pendakur 2008)
      EPQ=solve(diag(apply(shares[,1:(neq+1)],2,mean)))%*%(EP+apply(shares[,1:(neq+1)],2,mean)%*%t(apply(shares[,1:(neq+1)],2,mean)))
      colnames(EPQ) <- rownames(EPQ) <- c(labels.price[1:neq],"pothers")
      
      
      #*** calculation of elasticity of good j with respect to the log.price of good i
      #*** calculation of elasticity of good j with respect to income
      ajk <- my.array
      ELASTPRICE <- matrix(0,neq+1,neq+1)
      for (i in 1:(neq+1)){
        for (q in 1:(neq+1)){
          
          C <- 0
          for (j in 1:y.power){
            tempo <- j*bjr[j,i]*y^{j-1}
            C <- C+tempo}
          
          D <- 0
          if (zy.inter){
            for (j in 1:nsoc){
              tempo <- hjt[j,i]*Z[,j+1]
              D <- D+tempo}}
          
          E <- 0
          for (t in (1:(nsoc))){
            tempo <- ajk[t,q,i]*Z[,t]
            E <- E+tempo
          } 
          
          G <- 0
          if (py.inter){
            for (k in 1:(neq+1)){
              tempo <- bjk[k,i]*P[,k]
              G <- G+tempo }
            
            F <- bjk[q,i]*y
          }
          
          U <- 0
          for (d in 1:(neq+1)){
            for (t in (1:(nsoc))){
              tempo <- ajk[t,d,i]*Z[,t]*P[,d]
              U <- U+tempo
            } }
          
          B <- -mean(shares[,q]+U)/mean(1-1/2*tot2)-mean(y)/mean(1-1/2*tot2)*mean(G)
          
          H <- B*(C+D+G)+E+F
          
          ELASTPRICE[q,i] <- mean(H) / mean(shares[,i]) - as.numeric((i==q))
          
        }
      }
      
      colnames(ELASTPRICE) <- rownames(ELASTPRICE) <- c(labels.price[1:neq],"pothers")
      
    }
    
    
    
    
    if (type=="income"){
      
      #*** calculation of elasticity of good j with respect to the log.price of good i
      #*** calculation of elasticity of good j with respect to income
      ajk <- my.array
      ELASTINCOME <- matrix(0,1,neq+1)
      for (i in 1:(neq+1)){
        for (q in 1:(neq+1)){
          
          C <- 0
          for (j in 1:y.power){
            tempo <- j*bjr[j,i]*y^{j-1}
            C <- C+tempo}
          
          D <- 0
          if (zy.inter){
            for (j in 1:nsoc){
              tempo <- hjt[j,i]*Z[,j+1]
              D <- D+tempo}}
          
          E <- 0
          for (t in (1:(nsoc+1))){
            tempo <- ajk[t,q,i]*Z[,t]
            E <- E+tempo
          } 
          
          G <- 0
          if (py.inter){
            for (k in 1:(neq+1)){
              tempo <- bjk[k,i]*P[,k]
              G <- G+tempo }
            
            F <- bjk[q,i]*y
          }
          
          U <- 0
          for (d in 1:(neq+1)){
            for (t in (1:(nsoc+1))){
              tempo <- ajk[t,d,i]*Z[,t]*P[,d]
              U <- U+tempo
            } }
          
          B <- -mean(shares[,q]+U)/mean(1-1/2*tot2)-mean(y)/mean(1-1/2*tot2)*mean(G)
          
          H <- B*(C+D+G)+E+F
          
          ELASTINCOME[1,i] <- 1 + (mean(C+D+G))/mean(shares[,i])
          
        }
      }
      
      colnames(ELASTINCOME) <- c(labels.share[1:neq],"others")
 #### castellon 2012 doest work     
      a<-matrix(0,neq,neq)
      
      for (i in 1:neq){
        a[i,1]<-mean(shares[,i])*(1/(mean(cdf[,i])*sum_bjr[,i])+mean(shares[,i])*mean(p[,i]))
      }
      
      for (i in 1:neq){
        k<-1
        for (j in setdiff(1:neq,i)){
          k<-k+1
          a[i,k]<-mean(shares[,j])*mean((p[,j]))
        }
      }
      
      b<-matrix(0,nrow=neq,ncol=neq+1)
      for (i in 1:neq){
        b[i,1]<-mean(shares[,i])/(mean(cdf[,i])*sum_bjr[,i])+1
      }
      
      for (i in 1:neq){
        for (k in 2:(neq+1)){
          b[i,k]<-mean(shares[,k])*mean((p[,k]))
        }
      }
      b[,2]<-rowSums(b[,c(2:neq+1)])
      b<-b[,1]+b[,2]
      missing_one<-solve(as.matrix(a),as.matrix(b))
      last_one<-matrix(0,neq,1)
      for (i in 1:neq){
        last_one[i,]<-mean(shares[,i])*missing_one[i]
      }
      last_one<-c(last_one,(1-sum(last_one))/mean(shares[,neq+1]))
      last_one<-last_one/colMeans(shares)
    }   
      
      ## Calculation of income elasticities of budget shares
      ## page 13 formula 23 "EASI made EASIER" (Pendakur 2008)
      ER=matrix(0,1,(neq+1))
      for(i in 1:(neq+1)){
        tot11 = 0
        tempo1=tempo2=tempo3=0
        
        for (t in 1:y.power){
          tempo0 <- bjr[t,i]*t*y^(t-1)
          tempo2 <- tempo2+tempo0
        }
        
        if (zy.inter){   for (t in 1:nsoc){
          tempo0 <- hjt[t,i]*Z[,t+1]
          tempo3 <- tempo3+tempo0
        }}
        
        if (py.inter){  for (k in 1:(neq+1)){
          tempo0 <- bjk[k,i]*P[,k]
          tempo1 <- tempo1+tempo0
        }}
        
        tot11 <- tot11+tempo1+tempo2+tempo3
        ER[i]=mean(tot11)
      }
      
      colnames(ER) <- c(labels.share[1:neq],"others")
      
      ## Calculation of standard deviations of income elasticities
      ## (delta method) if "ERDELTA=TRUE"
      if (ERDELTA){
        ttt <- colnames(summary(fit3sls)$coefCov)
        ER_SE=matrix(0,1,neq+1)
        ELASTINCOME_SE=matrix(0,1,neq+1)
        for (i in (1:neq)){
          tt <- c()
          for (j in (1:y.power)){
            tt <- c(tt,paste("eq",i,"_y",j,sep=""))}
          
          
          if (zy.inter){
            for (t in (1:nsoc)){
              tt <- c(tt,paste("eq",i,"_yz",t,sep=""))
            }}
          if (py.inter){
            for (j in 1:neq){
              tt <- c(tt,paste("eq",i,"_ynp",j,sep=""))
            }
          }
          
          tnum <- match(tt,ttt)
          
          DD <- summary(fit3sls)$coefCov[c(tnum),c(tnum)]
          
          MAT <- c()
          for(r in 1:y.power){
            MAT <- cbind(MAT,r*y^{r-1})}
          if (zy.inter) MAT <- cbind(MAT,Z[,-1])
          if (py.inter) MAT <- cbind(MAT,P[,(1:neq)])
          
          ER_SE[1,i] <- median(sqrt(diag(as.matrix(MAT)%*%as.matrix(DD)%*%t(as.matrix(MAT)))))
          ELASTINCOME_SE [1,i] <- ER_SE[1,i]/mean(shares[,i])
          
          
          
        }
        ER_SE[1,neq+1] <- -sqrt(sum(ER_SE[1,1:neq]^2))
        ELASTINCOME_SE [1,neq+1] <- ER_SE[1,neq+1]/mean(shares[,neq+1])
        
        ER_SE <- as.matrix(ER_SE)
        ELASTINCOME_SE <- as.matrix(ELASTINCOME_SE)
        
        colnames(ER_SE) <-colnames(ELASTINCOME_SE) <- c(labels.share[1:neq],"others")
      }
      
      
      
    }
    
    
    if (type=="demographics"){
      ## Calculation of sociodemographic elasticities of budget shares 
      ## page 13 formula 23 "EASI made EASIER" (Pendakur 2008)
      EZ=matrix(0,nsoc,(neq+1))
      a <- my.array
      for(i in 1:(neq+1)){
        tempo4=tempoo=0
        for (t in 1:nsoc){
          
          if (interact){
            for (k in 1:(neq+1)){
              if (t %in% interpz){
                tempoo <- a[t,k,i]*P[,k]
                tempo4 <- tempo4+tempoo   
              }
            }
          }   
          tot12 <- gjt[t,i]+hjt[t,i]*y
          
          tot12 <- tot12+tempo4
          EZ[t,i]=mean(tot12)
        }
      }
      
      colnames(EZ) <- c(labels.share[1:neq],"others")
      rownames(EZ) <- labels.soc
      
      ## Calculation of standard deviations of sociodemographic elasticities 
      ## (delta method) if "EZDELTA=TRUE"
      if (EZDELTA){
        ttt <- colnames(summary(fit3sls)$coefCov)
        EZ_SE=matrix(0,nsoc,neq+1)
        for (i in 1:neq){
          for (j in 1:nsoc){
            tt <- c()
            tt <- c(tt,paste("eq",i,"_z",j,sep=""))
            if (zy.inter){
              tt <- c(tt,paste("eq",i,"_yz",j,sep=""))
            }
            if (pz.inter){
              if (j  %in% interpz){
                for (t in 1:neq){
                  tt <- c(tt,paste("eq",i,"_np",t,"z",j,sep=""))
                }}}
            
            
            tnum <- match(tt,ttt)
            
            DD <- summary(fit3sls)$coefCov[tnum,tnum]
            
            MAT <- Z[,1]
            
            if (zy.inter)
              MAT <- cbind(MAT,y)
            if (pz.inter){
              if (j  %in% interpz){
                MAT <- cbind(MAT,P[,(1:neq)])
              }}
            
            
            EZ_SE[j,i] <- median(sqrt(diag(as.matrix(MAT)%*%as.matrix(DD)%*%t(as.matrix(MAT)))))
            
          }}
        for (j in 1:nsoc) EZ_SE[j,neq+1] <- sqrt(sum(EZ_SE[j,1:neq]^2))
        
        colnames(EZ_SE) <- c(labels.share[1:neq],"others")
        rownames(EZ_SE) <- labels.soc
        
        
      }
      
    }
    
    if (type=="demographics") ER=ELASTINCOME=EP=EPS=EPQ=ELASTPRICE="Not calculated"
    if (type=="income") EZ=EP=EPS=EPQ=ELASTPRICE="Not calculated"
    if (type=="price") ER=ELASTINCOME=EZ="Not calculated"
    if (!EPDELTA)EP_SE=ELASTPRICE_SE="Not calculated"
    if (!ERDELTA) ER_SE=ELASTINCOME_SE="Not calculated"
    if (!EZDELTA) EZ_SE="Not calculated"
    
    
    Result <-  list(
      #EL1=`*** Semi-elasticities of budget shares in respect to log.prices ***`
      EP=EP,
      #EL2="*** Standard deviations of Semi-elasticities of budget shares in respect to log.prices ***",
      EP_SE=EP_SE,
      #EL3="*** Matrix of compensated quantity derivatives with respect to unlogged log.prices ***",
      EPS=EPS,
      #EL4="*** Compensated (good-specific) expenditures with respect to log.prices ***"
      EPQ=EPQ,
      #EL5="*** Elasticities of quantities in respect to log.prices ***",
      ELASTPRICE=ELASTPRICE,
      #EL6="*** Elasticities of quantities in respect to income ***",
      ELASTINCOME=ELASTINCOME,
      #EL7="*** Standard deviations of Elasticities of quantities in respect to log.prices ***",
      ELASTPRICE_SE=ELASTPRICE_SE,
      #EL8="*** Standard deviations of Elasticities of quantities in respect to income ***",
      ELASTINCOME_SE=ELASTINCOME_SE,
      #EL9="*** Semi-elasticities of budget shares in respect to real expenditures ***",
      ER=ER,
      #EL10="*** Standard deviations of Semi-elasticities of budget shares in respect to real expenditures ***",
      ER_SE=ER_SE,
      #EL11="*** Semi-elasticities of budget shares in respect to demographics ***",
      EZ=EZ,
      #EL12="*** Standard deviations of Semi-elasticities of budget shares in respect to demographics ***",
      EZ_SE=EZ_SE
    )
    
    Result
    
#### concavity ####
concavity <-
  function(object=object){
    
    cat("\n","*** Please wait during the test of local concavity *** ","\n")
    
    fit3sls <- object$fit3sls
    varlist <- object$varlist
    var.soc <- object$var.soc
    shares <- object$shares
    log.price <- object$log.price
    neq <- object$neq
    y.power <- object$y.power
    nsoc <- object$nsoc
    interact <- object$interact
    py.inter <- object$py.inter
    zy.inter <- object$zy.inter
    pz.inter <- object$pz.inter
    interpz <- object$interpz
    log.exp <- object$log.exp
    labels.price <- object$labels.price
    labels.soc <- object$labels.soc
    labels.share <- object$labels.share
    
    n <- length(log.exp)
    
    temp <- intermediate.blocs(object,log.price=log.price,var.soc=var.soc,log.exp=log.exp)
    my.array <- temp$my.array
    tot <- temp$tot
    tot2 <- temp$tot2
    tot0 <- temp$tot0
    bjk <- temp$bjk
    P <- temp$P
    w <- temp$w
    Z <- temp$Z
    bjr <- temp$bjr
    gjt <- temp$gjt
    hjt <- temp$hjt
    cc <- temp$cc
    noms <- object$noms
    lnx <- object$log.exp
    y <- object$y
    
    #** Recovery of the estimated coefficients***
    Estimates <- summary(fit3sls)$coefficients
    
    coef=Estimates[,1]
    n <- length(log.exp)
    
    ## Labels or names of the equations:
    noms <- c()
    for (i in 1:neq)
      noms <- c(noms, paste("eq",i,sep=""))
    
    #***** new price Matrix *******
    P=log.price
    
    #****  new sociodemographic matrix *****
    Z=cbind(1,var.soc)
    
    #**** logarithm of total expenditure *******
    lnx <- log.exp
    
    #**** Budget shares matrix *******
    w=matrix(0,n,neq+1)
    for (i in 1:(neq+1))
      w[,i] <- shares[,i]
    
    ## Calculation of first and second derivative of y with respect to p_j and p_k
    ## uprime
    h1 <- my.array
    toth1 <- matrix(0,n,neq)
    u1p <- matrix(0,n,neq)
    for (j in 1:neq){
      toth1[,j] = 0
      for (k in 1:neq){
        for (t in c(1:(nsoc+1))){
          tempo <- h1[t,k,j]*P[,k]/exp(P[,j])*Z[,t]
          toth1[,j] <- toth1[,j]+tempo
        }
      }
      u1p[,j]=-w[,j]/exp(P[,j])+toth1[,j]
    }
    
    ## usecond
    toth2 <- array(0,c(neq,neq,n))
    u2p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        toth2[j,k,] = 0
        for (t in c(1:(nsoc+1))){
          tempo <- h1[t,k,j]/(exp(P[,k])*exp(P[,j]))*Z[,t]
          toth2[j,k,] <- toth2[j,k,]+tempo
        }
      }
    }
    u2p=toth2
    
    ## vprime
    toth3 <- matrix(0,n,neq)
    v1p <- matrix(0,n,neq)
    for (k in 1:neq){
      tempo <- bjk[j,k]/exp(P[,j])*P[,k]
      toth3[,k] <- toth3[,k]+tempo
    }
    for (j in 1:neq){
      v1p[,j]=-toth3[,j]}
    
    
    ## vsecond
    toth4 <- array(0,c(neq,neq,n))
    v2p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        toth4[j,k,] = bjk[j,k]/(exp(P[,j])*exp(P[,k]))
      }
    }
    v2p=toth4
    
    S=y+tot0-1/2*tot-1/2*tot2*y
    
    u=lnx - tot0+1/2*tot
    v=1-1/2*tot2
    
    w <- array(0,c(neq,neq,n))
    wp <- array(0,c(neq,neq,n))
    lp <- array(0,c(neq,neq,n))
    l <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        w[j,k,] <- u1p[,j]*v-v1p[,j]*u
        l[j,k,] <- v*v
        wp[j,k,] <- u2p[j,k,]*v + u1p[,j]*v1p[,j]-v1p[,j]*u1p[,j]-v2p[j,k,]*u  
        lp[j,k,] <- 2*v1p[,j]*v
      }
    }
    
    Hy1 <- wp/l
    Hy2 <- (wp*l-lp*w)/(l)^2
    
    ## Calculation of of first and second derivative of S0 with respect to p_j and p_k :
    S0 <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S0[j,k,] <- tot0
      }
    }
    
    S01p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S01p[j,k,] <- shares[,j]/exp(P[,j])
      }
    }
    
    S02p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S02p[j,k,] <- -(j==k)*shares[,j]/exp(P[,j])^2
      }
    }
    
    ## Calculation of of first and second derivative of S1 with respect to p_j and p_k :
    S1 <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S1[j,k,] <- -1/2*tot
      }
    }
    
    S11p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S11p[j,k,] <- -toth1[,j]
      }
    }
    
    S12p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S12p[j,k,] <- -toth2[j,k,]
      }
    }
    
    ## Calculation of of first and second derivative of S2 with respect to p_j and p_k :
    S2 <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S2[j,k,] <- -1/2*tot2
      }
    }
    
    S21p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S21p[j,k,] <- -toth3[,j]
      }
    }
    
    S22p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S22p[j,k,] <- -toth4[j,k,]
      }
    }
    
    ## First derivate of S in respect with p_j
    SSp <- Hy1 + S01p + S11p + S21p*y + S2*Hy1
    
    ## Second derivate of S
    SS2p <- Hy2 + S02p + S12p + S22p*y + S21p*Hy1 + S21p*Hy1+S2*Hy2
    
    w <- shares
    
    ## Calculation of first and second derivative of y with respect to p_j and p_k
    ## uprime
    h1 <- my.array
    toth1 <- matrix(0,n,neq)
    u1p <- matrix(0,n,neq)
    for (j in 1:neq){
      toth1[,k] = 0
      for (k in 1:neq){
        for (t in c(1:(nsoc+1))){
          tempo <- h1[t,k,j]*P[,k]/exp(P[,k])*Z[,t]
          toth1[,k] <- toth1[,k]+tempo
        }
      }
      u1p[,k]=-w[,k]/exp(P[,k])+toth1[,k]
    }
    
    ## usecond
    toth2 <- array(0,c(neq,neq,n))
    u2p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        toth2[,k,j] = 0
        for (t in c(1:(nsoc+1))){
          tempo <- h1[t,k,j]/(exp(P[,j])*exp(P[,k]))*Z[,t]
          toth2[k,j,] <- toth2[k,j,]+tempo
        }
      }
    }
    u2p=toth2
    
    ## vprime
    toth3 <- matrix(0,n,neq)
    v1p <- matrix(0,n,neq)
    for (j in 1:neq){
      tempo <- bjk[j,k]/exp(P[,k])*P[,j]
      toth3[,j] <- toth3[,j]+tempo
    }
    for (k in 1:neq){
      v1p[,k]=-toth3[,k]}
    
    ## vsecond
    toth4 <- array(0,c(neq,neq,n))
    v2p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        toth4[k,j,] = bjk[j,k]/(exp(P[,k])*exp(P[,j]))
      }
    }
    v2p=toth4
    
    S=y+tot0-1/2*tot-1/2*tot2*y
    
    u=lnx - tot0+1/2*tot
    v=1-1/2*tot2
    
    wp <- array(0,c(neq,neq,n))
    l <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        l[k,j,] <- v*v
        wp[k,j,] <- u2p[k,j,]*v + u1p[,k]*v1p[,k]-v1p[,k]*u1p[,k]-v2p[k,j,]*u  
        
      }
    }
    
    Hy1 <- wp/l
    
    ## Calculation of of first and second derivative of S0 with respect to p_j and p_k :
    S0 <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S0[k,j,] <- tot0
      }
    }
    
    S01p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S01p[k,j,] <- shares[,k]/exp(P[,k])
      }
    }
    
    ## Calculation of of first and second derivative of S1 with respect to p_j and p_k :
    S1 <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S1[k,j,] <- -1/2*tot
      }
    }
    
    S11p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S11p[k,j,] <- -toth1[,k]
      }
    }
    
    ## Calculation of of first and second derivative of S2 with respect to p_j and p_k :
    S2 <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S2[k,j,] <- -1/2*tot2
      }
    }
    
    S21p <- array(0,c(neq,neq,n))
    for (j in 1:neq){
      for (k in 1:neq){
        S21p[k,j,] <- -toth3[,k]
      }
    }
    
    ## First derivate of S in respect with p_k
    SSpp <- Hy1 + S01p + S11p + S21p*y + S2*Hy1
    
    
    
    hess <- exp(S)*(SS2p+SSp*SSpp)
    
    testconcav <- matrix(0,n,1)
    for (i in 1:n){
      testconcav[i] <- semidefiniteness(hess[,,i],positive=TRUE)
      if (!testconcav[i]) testconcav[i] <- 1 else testconcav[i] <- 0
    }
    
    
    for (i in (10:1)/10){ 
      if (mean(testconcav)<=i)
        a <- paste(" Result: The cost function is concave on more than ",(i-0.1)*100,"% of the sample. ",sep="")
    }
    
    cat("\n")
    a
  }

#### engel curves ####
engel <-
  function(object=object,file=file,sd=FALSE,limY=FALSE){
    
    WDELTA <- ifelse(sd,TRUE,FALSE) 
    fit3sls <- object$fit3sls
    varlist <- object$varlist
    var.soc <- object$var.soc
    shares <- object$shares
    log.price <- object$log.price
    neq <- object$neq
    y.power <- object$y.power
    nsoc <- object$nsoc
    interact <- object$interact
    py.inter <- object$py.inter
    zy.inter <- object$zy.inter
    pz.inter <- object$pz.inter
    interpz <- object$interpz
    log.exp <- object$log.exp
    labels.price <- object$labels.price
    labels.soc <- object$labels.soc
    labels.share <- object$labels.share
    dim_varlist <- object$dim_varlist
    y <- object$y
    
    n <- length(log.exp)
    
    
    temp <- intermediate.blocs(object,log.price=log.price,var.soc=var.soc,log.exp=log.exp)
    my.array <- temp$my.array
    tot <- temp$tot
    tot2 <- temp$tot2
    tot0 <- temp$tot0
    bjk <- temp$bjk
    P <- temp$P
    w <- temp$w
    Z <- temp$Z
    bjr <- temp$bjr
    gjt <- temp$gjt
    hjt <- temp$hjt
    cc <- temp$cc
    noms <- temp$noms
    lnx <- temp$log.exp
    y <- temp$y
    
    ## Calculation of w_j
    W=matrix(0,n,neq)
    
    ajk <- my.array
    for (i in 1:neq){
      
      tot3 <- tot4 <- tot5 <- tot6 <- tot7 <- 0
      
      for (j in 1:y.power){
        tempo <- bjr[j,i]*y^j
        tot3 <- tot3+tempo}
      
      for (j in 1:nsoc){
        tempo <- gjt[j,i]*Z[,j+1]
        tot4 <- tot4+tempo}
      
      if (zy.inter){
        for (j in 1:nsoc){
          tempo <- hjt[j,i]*Z[,j+1]*y
          tot5 <- tot5+tempo}}
      
      if (pz.inter){
        for (k in 1:neq){
          for (t in (1:(nsoc+1))){
            tempo <- ajk[t,k,i]*Z[,t]*P[,k]
            tot6 <- tot6+tempo
          }}} 
      
      if (py.inter){
        for (k in 1:neq){
          tempo <- bjk[k,i]*P[,k]*y
          tot7 <- tot7+tempo }}
      
      
      W[,i] <- cc[i]+tot3+tot4+tot5+tot6+tot7
      
    }
    colnames(W) <- labels.share[1:neq]
    
    
    ## Calculation of standard deviations of the fitted budget shares (if WDELTA=TRUE) - (delta method)
    if (WDELTA){
      nb <- dim_varlist
      MAT <- rep(1,n)
      for (i in 1:y.power)
        MAT <- cbind(MAT,y^i)
      for (i in 1:nsoc)
        MAT <- cbind(MAT,Z[,i+1])
      if (zy.inter){
        for (i in 1:nsoc)
          MAT <- cbind(MAT,y*Z[,i+1])
      }
      for (i in 1:neq)
        MAT <- cbind(MAT,P[,i])
      if (py.inter){
        for (i in 1:neq)
          MAT <- cbind(MAT,y*P[,i])
      }
      if (pz.inter){
        for (i in interpz){
          for (j in 1:neq){
            MAT <- cbind(MAT,Z[,i+1]*P[,j])}}
      }
      
      nn <- 1
      W_ecart <- matrix(0,n,neq)
      
      for (i in 1:neq){
        DD <- summary(fit3sls)$coefCov[nn:(nn+nb-1),nn:(nn+nb-1)]
        
        W_e <- MAT%*%DD%*%t(MAT)
        
        W_ecart[,i] <- sqrt(diag(W_e))
        rm(W_e)
        
        nn <- nn+nb
      }
    }
    
    
    ## Labels of W matrix
    colnames(W)=noms
    
    ## Engel Curves
    quant=quantile(lnx,seq(0,1,0.01))
    
    ee=rep(1,n)
    for (j in 1:n) {
      for (i in 1:100)
        if (lnx[j]>quant[i]) 
          ee[j]=i
    }
    
    Wm=matrix(0,100,neq)
    for (i in 1:100){
      for (j in 1:neq)
        Wm[i,j]=median(W[ee==i,j])
    }
    
    Wm_autres=1-apply(Wm,1,sum)
    Wm=cbind(Wm,Wm_autres)
    
    
    ## Calculation of confidence intervals for fitted budget shares (if WDELTA=TRUE)
    if (WDELTA){
      Wme=matrix(0,100,neq+1)
      for (i in 1:100){
        for (j in 1:neq)
          Wme[i,j]=median(W_ecart[ee==i,j])
      }
      
      for (i in 1:100){
        Wme[i,neq+1]=sqrt(sum(Wme[i,1:neq]^2))
      }
      
      Wmep <- Wm+1.96*Wme
      Wmem <- Wm-1.96*Wme
      
      Wmep <- Wmep[(1:20)*5,]
      Wmem <- Wmem[(1:20)*5,]
    }
    
    ### management of labels.share
    if (length(labels.share)<2) labels.share <- noms
    
    limYY <- c()
    if (length(limY)<2){
      for (i in 1:neq) 
        limYY <- c(limYY,c(0,summary(w[,i])[5]))} else {
          limYY <- limY}
    
    ss <- seq(1,neq*2,by=2)
    
    ## Export of Engel curves in the parent folder under the name "file". pdf
    ## File name is entered on the command line
    
    pdf(paste("./",file,".pdf"))
    for(i in 1:neq)
    {
      # smoothing cubic
      sp <- smooth.spline(c(1:100),Wm[,i], spar = 0.9)
      y.loess <- loess(Wm[,i] ~ c(1:100), span=0.75, data.frame(xxx=c(1:100), yyy=Wm[,i]))
      y.predict <- predict(y.loess, data.frame(xxx=c(1:100)))
      plot(c(1:100),Wm[,i],xlab="Percentiles of total expenditure",ylab="Budget shares",col="green",ylim=c(limYY[ss[i]],limYY[ss[i]+1]))
      title(main=labels.share[i])
      xx <- seq(1,100, len=20)
      # plot of the adjustment curve 
      lines(predict(sp,xx), col = "red")
      lines(c(1:100),y.predict,col="blue")
      lines(ksmooth(c(1:100),Wm[,i], "normal", bandwidth=10), col="black")
      if (WDELTA){  points(c((1:20)*5),Wmep[,i],pch="+",cex=1,col="violet")
        points(c((1:20)*5),Wmem[,i],pch="+",cex=1,col="violet")}
    }
    sp <- smooth.spline(c(1:100),Wm[,i+1], spar = 0.9)
    plot(c(1:100),Wm[,i+1],xlab="Percentiles of total expenditure",ylab="Budget shares",col="green")
    title(main ="Others Goods")
    xx <- seq(1,100, len=20)
    y.loess <- loess(Wm[,i+1] ~ c(1:100), span=0.75, data.frame(xxx=c(1:100), yyy=Wm[,i]))
    y.predict <- predict(y.loess, data.frame(xxx=c(1:100)))
    # plot of the adjustment curve
    lines(predict(sp,xx), col = "red")
    lines(c(1:100),y.predict,col="blue")
    lines(ksmooth(c(1:100),Wm[,i+1], "normal", bandwidth=10), col="black")
    if (WDELTA){  points(c((1:20)*5),Wmep[,i+1],pch="+",cex=1,col="violet")
      points(c((1:20)*5),Wmem[,i+1],pch="+",cex=1,col="violet")}
    dev.off()
    
    
    ## refreshment
    if (WDELTA){rm(MAT); rm(W_ecart); rm(DD)}
    
    result <- W
    result
  }




    
#### simulation ####
simulations <-
      function(object=object,log.price_new=log.price_new,var.soc_new=var.soc_new,log.exp_new=log.exp_new){
        
        fit3sls <- object$fit3sls
        varlist <- object$varlist
        shares <- object$shares
        neq <- object$neq
        y.power <- object$y.power
        nsoc <- object$nsoc
        interact <- object$interact
        py.inter <- object$py.inter
        zy.inter <- object$zy.inter
        pz.inter <- object$pz.inter
        interpz <- object$interpz
        labels.price <- object$labels.price
        labels.soc <- object$labels.soc
        labels.share <- object$labels.share
        dim_varlist <- object$dim_varlist
        log.exp <- object$log.exp
        
        n <- length(log.exp)
        
        log.price=log.price_new 
        var.soc=var.soc_new
        log.exp <- log.exp_new 
        
        
        
        temp <- intermediate.blocs(object,log.price=log.price,var.soc=var.soc,log.exp=log.exp)
        my.array <- temp$my.array
        tot <- temp$tot
        tot2 <- temp$tot2
        tot0 <- temp$tot0
        bjk <- temp$bjk
        P <- temp$P
        w <- temp$w
        Z <- temp$Z
        bjr <- temp$bjr
        gjt <- temp$gjt
        hjt <- temp$hjt
        cc <- temp$cc
        noms <- object$noms
        lnx <- object$log.exp
        y <- object$log.exp
        
        
        
        ## Calculation of w_j
        W=matrix(0,n,neq)
        
        ajk <- my.array
        for (i in 1:neq){
          
          tot3 <- tot4 <- tot5 <- tot6 <- tot7 <- 0
          
          for (j in 1:y.power){
            tempo <- bjr[j,i]*y^j
            tot3 <- tot3+tempo}
          
          for (j in 1:nsoc){
            tempo <- gjt[j,i]*Z[,j]
            tot4 <- tot4+tempo}
          
          if (zy.inter){
            for (j in 1:nsoc){
              tempo <- hjt[j,i]*Z[,j]*y
              tot5 <- tot5+tempo}}
          
          if (pz.inter){
            for (k in 1:neq){
              for (t in (1:(nsoc))){
                tempo <- ajk[t,k,i]*Z[,t]*P[,k]
                tot6 <- tot6+tempo
              }}} 
          
          if (py.inter){
            for (k in 1:neq){
              tempo <- bjk[k,i]*P[,k]*y
              tot7 <- tot7+tempo }}
          
          
          W[,i] <- cc[i]+tot3+tot4+tot5+tot6+tot7
          
        }
        
        
        colnames(W) <- labels.share[1:neq]
        
        Result <-  list(
          CoefCov=fit3sls$coefCov,
          fit3sls=fit3sls,
          y=y,
          varlist=varlist,
          var.soc=var.soc,
          shares=shares,
          log.price=log.price,
          neq=neq,
          y.power=y.power,
          nsoc=nsoc,
          interact=interact,
          py.inter=py.inter,
          zy.inter=zy.inter,
          pz.inter=pz.inter,
          interpz=interpz,
          log.exp=log.exp,
          labels.price=labels.price,
          labels.soc=labels.soc,
          labels.share=labels.share,
          dim_varlist=dim_varlist,
          fitted.w=W
        )
        
        Result
        
      }
    