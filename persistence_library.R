
require(ncdf4); require(lubridate); require(MASS)

### EOF ANALYSIS ###

comp.EOF <- function(lon,lat,data,th,anomaly=T,weight=T,maps=T){
  
  # lon, lat: to get coordinates
  # data: lon x lat x time
  # th: threshold for cumulative explained variance (e.g., 90%)
  # anomaly: calculate EOFs from anomalies [True/False]
  # weight: apply latitudinal weighting [True/False]
  # maps: return EOF patterns [True/False]
  
  grid <- expand.grid(x=lon,y=lat)
  steps <- dim(data)[3]
  numCells <- length(lon)*length(lat)
  
  data.df <- data.frame(x=rep(grid$x,steps),y=rep(grid$y,steps),
                        timestep=rep(c(1:steps),each=numCells),
                        d=as.vector(data))
  
  # Weigh data
  if (weight==T){
    data.df$d <- cos(pi*data.df$y/180)*data.df$d 
  }
  
  data.EOF <- acast(data.df,x+y~timestep,value.var='d')
  
  if (anomaly==T){
    data.EOF <- sweep(data.EOF,1,rowMeans(data.EOF),"-")
  }
  
  EOF <- prcomp(t(data.EOF),scale.=F,center=T)
  n <- which((cumsum(EOF$sdev^2)/sum(EOF$sdev^2))>=th)[1]
  
  if (maps){
    # Back to lon/lat format
    eof.df <- as.data.frame(EOF$rotation[,1:n])
    eof.df <- melt(eof.df)
    # Get lon, lat coordinates
    coordLon <- as.numeric(sapply(strsplit(rownames(EOF$rotation),"\\_"),`[[`,1))
    coordLat <- as.numeric(sapply(strsplit(rownames(EOF$rotation),"\\_"),`[[`,2))
    eof.df$lon <- rep(coordLon,n)
    eof.df$lat <- rep(coordLat,n)
    # Reorder columns
    eof.df <- eof.df[,c(3,4,1,2)]
    # Reapply weighting
    if (weight==T){
      eof.df$value <- eof.df$value/cos(pi*eof.df$lat/180)
    }
    PC <- acast(eof.df,lon~lat~variable) 
    return(list(EOF$x[,1:n],(cumsum(EOF$sdev^2)/sum(EOF$sdev^2))[1:n],PC))
  }else{
    return(list(EOF$x[,1:n],res[1:n]))
  }
}

### AUTOCORRELATION ###

# Decorrelation time
decorrelation.time <- function(x,lag.max){
  return(sum(as.numeric(acf(x,lag.max=lag.max,na.action=na.pass,plot=F)[[1]])))
}

# Characteristic time
characteristic.time <- function(x,lag.max){
  y <- as.numeric(acf(x,lag.max=lag.max,na.action=na.pass,plot=F)[[1]])
  return(sum(y*(1-(0:lag.max)/lag.max)))
}

### QUASISTATIONARY STATES ###

# Weighting function
phi <-function(d,epsilon){return(0.5*(1+cos(pi*(d/epsilon)**2)))}

# Calculate average tendencies
tendency <- function(x,data,epsilon){
  
  # x: data point at which to calculate tendency
  # data: data matrix (time x EOF/space dimension)
  # epsilon: distance threshold to select neighbours
  
  distVec <- as.numeric(sqrt(rowSums(sweep(data,2,x,'-')**2)))
  inds <- which(distVec<=epsilon)
  inds <- inds[inds>1 & inds<nrow(data)] # avoid extremities
  w <- phi(distVec[inds],epsilon)/sum(phi(distVec[inds],epsilon))
  return(as.numeric(colSums(sweep((data[inds+1,]-data[inds-1,])/2,1,w,'*'))))
}

# Solve {|T(x)|**2=0} and remove non-significant/divergent solutions
optimise <- function(data,frac.n=0.1){
  
  # data: data matrix (time x EOF/space dimension)
  # frac.n: percentile of pairwise distances to select neighbours in tendency calculation
  #   [Vautard 1990 uses frac.n=0.1]
  
  # Distance to select neighbours
  epsilon <- as.numeric(quantile(dist(data,method='euclidean'),frac.n))
  
  # Optimise starting from all data points
  optim_params <- array(NA,c(nrow(data),ncol(data)))
  optim_tendency <- NA*1:nrow(data)
  for (t in 1:nrow(data)){
    tryCatch({
      res <- optim(par=as.numeric(data[t,]),
                   function(x){return(sum(tendency(x,data=data,epsilon)**2))},
                   control=list(ndeps=0.02*apply(data,2,sd)),
                   method="BFGS")
      optim_params[t,] <- res$par
      optim_tendency[t] <- res$value
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  # Tendencies on observations
  obs_tendency <- apply(data,1,function(x){return(sum(tendency(x,data=data,epsilon)**2))},data=data,epsilon)
  
  # Enforce sufficient number of neighbours among observations
  # + low enough tendency
  num_n <- apply(optim_params,1,function(x){return(sum(sqrt(rowSums(sweep(data,2,x,'-')**2))<=epsilon))})
  inds <- which((!is.na(num_n)) & (num_n>=frac.n*nrow(data)) & (optim_tendency<=quantile(obs_tendency,0.1)))
  
  # Statistical tests
  T_all <- t(apply(data,1,tendency_vec,data=data,epsilon))
  J_all <- quantile(apply(T_all,1,function(x){return(sum(x**2))}),0.1)
  
  # Loop
  select <- 0*inds
  for (i in 1:length(inds)){
    
    x_s <- optim_params[inds[i],] # x*
    distVec <- as.numeric(sqrt(rowSums(sweep(data,2,x_s,'-')**2)))
    jnds <- which(distVec<=epsilon) # phase-space neighbours
    
    # Tendency
    T_s <- tendency_vec(x_s,data,epsilon) # T(x*)
    
    # Covariance matrix
    w <- phi(distVec[jnds],epsilon)/sum(phi(distVec[jnds],epsilon))
    N_s <- 1/sum(w**2)
    C_s <- matrix(0,ncol(data),ncol(data))
    for (t in 1:length(jnds)){
      C_s <- C_s + w[t]*(T_s-T_all[jnds[t],])%*%t(T_s-T_all[jnds[t],])
    }
    C_s <- C_s/N_s
    
    # Generate Gaussian distributed tendencies using C_s
    T_mc <- MASS::mvrnorm(1000, T_s, C_s)
    if (quantile(apply(T_mc,1,function(x){return(sum(x**2))}),0.9)<J_all){
      T2_mc <- MASS::mvrnorm(1000, 0*T_s, C_s)
      if (quantile(apply(T2_mc,1,function(x){return(sum(x**2))}),0.95)>sum(T_s**2)){
        select[i] <- 1
      }
    }
  }
  inds <- inds[select==1]
  
  return(list(optim_tendency[inds],optim_params[inds,]))
}

### TRANSITION PROBABILITIES ###

tr.prob <- function(y){
  
  # y: integer time series (categorical data)
  
  states <- sort(unique(y))
  probs <- 0*states
  for (n in 1:length(states)){
    probs[n] <- mean(y[which(y==states[n])+1]==states[n])
  }
  return(list(states,probs))
}

### RESIDENCE TIMES ###

# Calculate the ensemble of residence times for a given category/state
res.time <- function(y,y0){
  
  # y: integer time series (categorical data)
  # y0: selected system state
  
  w <- rle(1*(y==y0))
  return(w$lenghts[w$values==1])
}

### DISPERSION COEFFICIENT ###

# Runs declustering
decluster <- function(v,th,r){
  
  # v: data series
  # th: threshold to define (extreme) event occurrences
  # r: run length for declustering
  
  x <- 1*(v>=th)
  y <- x
  w <- rle(x)
  
  # Set runs of zeros with length<r to 1
  final_pos <- cumsum(w$lengths)
  for (uu in 1:(length(w$lengths)-1)){
    if ((w$lengths[uu]<r) & (w$values[uu]==0)){
      y[(final_pos[uu]-w$lengths[uu]+1):final_pos[uu]] <- 1
    }
  }
  if (x[length(x)]==1){
    uu <- w$lengths[length(w$lengths)]
    if ((w$lengths[uu]<r) & (w$values[uu]==0)){
      y[(final_pos[uu]-w$lengths[uu]+1):final_pos[uu]] <- 1
    }
  }
  w <- rle(y)
  
  # Maximum element of each run of '1's
  final_pos <- cumsum(w$lengths)
  out_max_idx <- out_max <- out_clust_len <- out_clust_evts <- out_clust_tot <- 0*1:sum(w$values==1)
  index <- 1
  for (i in 1:length(final_pos)){
    if (w$values[i]==1){
      vec <- v[(final_pos[i]-w$lengths[i]+1):final_pos[i]]
      vec2 <- x[(final_pos[i]-w$lengths[i]+1):final_pos[i]]
      out_max[index] <- max(vec)
      out_max_idx[index] <- which.max(vec)+final_pos[i]-w$lengths[i]
      out_clust_len[index] <- length(vec)
      out_clust_evts[index] <- sum(vec2)
      out_clust_tot[index] <- sum(vec)
      index <- index+1
    }
  }
  res <- 0*v
  res[out_max_idx] <- 1
  return(res)
}

dispersion <- function(x,window.size){
  
  # x: binary time series of event occurrences
  # window.size: window width
  
  n <- 0*1:floor(length(x)/window.size)
  t <- 1
  for (i in seq(1,window.size*floor(length(x)/window.size),window.size)){
    n[t] <- sum(x[i:(i+window.size-1)])
    t <- t+1
  }
  return(var(n)/mean(n)-1)
}

### RIPLEY'S K FUNCTION ###

ripleyK <- function(x,window.size){
  
  # x: binary time series of event occurrences
  # window.size: width of Ripley's K window
  
  S <- 0
  inds <- which(x==1)
  for (i in inds){
    S <- S+sum(x[max((i-window.size),0):(i-1)])+sum(x[(i+1):min((i+window.size),length(x))])
  }
  return(S/sum(x[inds]))
}

### RECURRENCE PLOTS ###

# Number of diagonal lines
diag.line.num <- function(rp,l){
  
  # rp: binary recurrence plot [n x n]
  # l: length of diagonal line
  
  n <- nrow(rp)
  res <- 0
  # i=1
  for (j in 3:(n-l)){
    if (l>1){
      res <- res+as.numeric((1-rp[1+l,j+l])*prod(diag(rp[1:(1+l-1),j:(j+l-1)]))) 
    }else{
      res <- res+as.numeric((1-rp[1+l,j+l])*rp[1,j])
    }
  }
  # i>1
  for (i in 2:(n-l-2)){
    for (j in (i+2):(n-l)){
      if (l>1){
        res <- res+as.numeric((1-rp[i-1,j-1])*(1-rp[i+l,j+l])*prod(diag(rp[i:(i+l-1),j:(j+l-1)])))
      }else{
        res <- res+as.numeric((1-rp[i-1,j-1])*(1-rp[i+1,j+1])*rp[i,j])
      }
    }
  }
  return(res)
}

# Number of vertical lines
vert.line.num <- function(rp,v){
  
  # rp: binary recurrence plot [n x n]
  # v: length of vertical line
  
  n <- nrow(rp)
  res <- 0
  # for i=1
  if (v>1){
    res <- res+as.numeric((1-rp[1,1+v])*prod(rp[1,1:v])) 
  }else{
    res <- res+as.numeric((1-rp[1,1+v])*rp[1,1])
  }
  # for i>1
  for (i in 2:(n-v)){
    if (v>1){
      res <- res+as.numeric((1-rp[i,i+v])*prod(rp[i,i:(i+v-1)])) 
    }else{
      res <- res+as.numeric((1-rp[i,i+v])*rp[i,i])
    }
  }
  return(res)
}

# Canonical indices: recurrence rate
calc.RR <- function(rp){
  
  # rp: binary recurrence plot [n x n]
  
  return(mean(rp))
}

calc.DET <- function(rp,lmin){
  
  # rp: binary recurrence plot [n x n]
  # lmin: minimum diagonal line length
  
  DET <- 0
  for (l in lmin:(nrow(rp)-1)){ # maximum diagonal length is nrow(rp)-1
    DET <- DET + l * diag.line.num(rp,l)
  }
  N <- 0
  for (l in 1:(lmin-1)){
    N <- N + l * diag.line.num(rp,l)
  }
  return(DET/(DET+N))
}

calc.LAM <- function(rp,vmin){
  
  # rp: binary recurrence plot [n x n]
  # vmin: minimum vertical line length
  
  # Laminarity
  LAM <- 0
  for (v in vmin:(nrow(rp)-2)){ # maximum vertical length is nrow(rp)-1
    LAM <- LAM + v * vert.line.num(rp,v)
  }
  N <- 0
  for (l in 1:(vmin-1)){
    N <- N + l * vert.line.num(rp,l)
  }
  return(LAM/(LAM+N))
}

# Number of recurrences
calc.num.rec <- function(rp,p,q,r.decl){
  
  # rp: binary recurrence plot [n x n]
  # p: length of recurrence line
  # q: recurrence horizon
  # r.decl: runs declustering window
  
  res <- 0*1:nrow(rp)
  
  for (i in 1:(nrow(rp)-q+1)){
    
    # Submatrix
    X <- RP.bin[i:(i+p-1),i:(i+q-1)]
    
    if (sum(X)>p){ # at least 1 recurrence outside the diagonal
      
      # Loop through to detect diagonals
      select_diag <- sum_diag <- c()
      for (j in (r.decl+1):(N-p+1)){
        
        w <- rle(c(diag(X[,j:(j+p-1)]),0))
        
        # If zeros are separated by less than r.decl steps
        if (max(w$lengths[w$values==0])<r.decl){
          
          select_diag <- c(select_diag,j)
          sum_diag <- c(sum_diag,sum(diag(X[,j:(j+p-1)])))
        }
      }
      
      # Remove neighbouring diagonals
      if (length(select_diag)>1){
        x <- 0*(1:N)
        x[select_diag] <- 1
        y <- 0*(1:N)
        y[select_diag] <- sum_diag
        res <- myDecluster(y,x,r=r.decl)
        # Make sure diagonals are separated by >=r.decl day
        select_diag <- res[[1]]
        num <- 1
        if (length(select_diag)>1){
          for (u in 1:(length(select_diag)-1)){
            ok <- 0
            for (ii in 1:p){
              x <- as.numeric(X[ii,(select_diag[u]+ii-1):(select_diag[u+1]+ii-1)])
              x[1] <- 1
              x[length(x)] <- 1
              # At least r.decl zeros?
              if (sum(x==0)>(r.decl-1)){
                ok <- ok+1
              }
            }
            if (ok==p){
              num <- num+1
            }
          }
        }
        res[i] <- num
      }
    }
  }
}








