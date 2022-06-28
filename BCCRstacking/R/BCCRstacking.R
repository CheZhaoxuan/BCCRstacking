#' @title An Efficient Genome Selection Program for Big Data Based on Machine Learning
#'
#' @description
#' \code{BCCRstacking}calculation of GEBV based on gradient boosting compressed component regression algorithm and stacking ensemble learning algorithm
#'
#' @param y.train a phenotype vector of type numeric
#' @param X.train a genotype matrix of type numeric
#' @param y.test a phenotype vector of type numeric, the phenotype of verification group
#' @param X.test a genotype matrix of type numeric, the verification group(or no phenotype group) required to calculate GEBV
#' @param CV a logical variable, decide whether to include covariates when compressing components
#' @param cc_num the number of compressed components
#'
#' @return a list including compressed components and GEBV(need provide X.test)
#' @export
#'
#' @examples
#' data(example_dataset)
#' library(Rfast)
#' y.train <- phenotype[c(1:200),1]
#' X.train <- genotype[c(1:200),]
#' X.test <- genotype[-c(1:200),]
#' BCCRS <- BCCRstacking(y.train = y.train, X.train = X.train, X.test = X.test, CV=FALSE, cc_num=NULL)
#'

BCCRstacking <- function(
  y.train = y.train,
  X.train = X.train,
  y.test = NULL,
  X.test = NULL,
  CV = FALSE,
  cc_num = NULL)
{
  BCCRstacking.Version <- function(){
    BCCRstacking.Version="BCCRstacking v1.00, June 27, 2022"
    return(BCCRstacking.Version)
  }

  print("------------------------------------- Welcome to BCCRstacking --------------------------------------")
  print("--------------------------- Checking input data format for  BCCRstacking ---------------------------")
  if(is.null(y.train)){stop("y.train must be required")}
  if(!is.numeric(y.train)){stop("y.train must be numeric")}
  if(!is.vector(y.train)){stop("y.train must be a vector")}

  if(is.null(X.train)){stop("X.train must be required")}
  if(!is.numeric(X.train)){stop("X.train must be numeric")}
  if(!is.matrix(X.train)){stop("X.train must be a matrix")}
  if(sum(is.na(X.train)) != 0){stop("NAs are not allowed in X.train")}

  if(length(y.train) != nrow(X.train)){stop("The individual dimensions of phenotype and genotype should be consistent")}

  if(!is.null(y.test) & length(y.test)!=nrow(X.test)){stop("The individual dimensions of y.test and X.test should be consistent")}

  if(!is.null(cc_num) & !is.integer(cc_num)){stop("cc_num should be NULL or integer")}
  if(is.integer(cc_num) & !length(cc_num)==1){stop("cc_num should be NULL or one number")}
  if(!is.null(cc_num)){
    print("------------------------------------- cc_num has been provided -------------------------------------")
    print("----------------------- The number of CC will be less than or equal to cc_num ----------------------")
  }

  if(!is.logical(CV)){stop("CV should be logical")}

  if(!is.null(X.test) & !is.numeric(X.test)){stop("X.test must be numeric")}
  if(!is.null(X.test) & !is.matrix(X.test)){stop("X.test must be a matrix")}
  if(!is.null(X.test) & !is.matrix(X.test)){stop("X.test must be a matrix")}
  if(sum(is.na(X.test)) != 0){stop("NAs are not allowed in X.test")}
  if(is.null(X.test)){
    print("-------------------------------------- X.test is not provided --------------------------------------")
    print("----------------------------------------- Just caculate CC -----------------------------------------")
  }
  if(!is.null(X.test)){
    print("------------------------------------- X.test has been provided -------------------------------------")
    print("--------------------------------------- Caculate CC and GEBV ---------------------------------------")
  }

  normalize <- function(x){
    return((x-min(x))/(max(x)-min(x)))
  }


  covar<-function(y,x){
    denom<-dim(x)[1]-1;
    (Rfast::eachcol.apply(x,y) - Rfast::colmeans(x)*sum(y))/denom
  }


  Boosted_CC <- function(
    y.train = y.train,
    X.train = X.train,
    X.test = NULL,
    CV1 = NULL,
    cc_num = NULL)
  {
    YY <- normalize(y.train)
    n <- length(YY)
    m <- Rfast::colmeans(X.train)
    sx <- Rfast::colVars(X.train, suma = n * m)
    N <- length(y.train)
    CVm <- matrix(1, nrow=N, ncol=1)
    CV1 <- CV1
    CV <- as.matrix(cbind(CVm,CV1))
    U1 <- crossprod(CV, YY)
    U2 <- solve(crossprod(CV), U1)
    CVt <- CV %*% U2
    if(!is.null(cc_num)){
      S1M = matrix(0, nrow=N, ncol=cc_num)
      if(!is.null(X.test)){N2 <- dim(X.test)[1]}
      if(!is.null(X.test)){S2M = matrix(0, nrow=N2, ncol=cc_num)}
      for (i in 1:cc_num){
        yt <- YY - CVt
        yt <- as.numeric(yt)
        r <- covar(yt, X.train)
        beta <- r/sx
        beta[is.na(beta)] <- Inf
        delSNP <- which(abs(beta) == Inf)
        if(length(delSNP)>0){
          beta <- beta[-c(delSNP)]
          S1M[, i] <- (X.train[,-c(delSNP)]%*%beta)/length(beta)
          if(!is.null(X.test)){S2M[, i] <- (X.test[,-c(delSNP)]%*%beta)/length(beta)}
          YY <- normalize(YY - (CVt + X.train[,-c(delSNP)]%*%beta))
        }
        else{
          S1M[, i] <- (X.train%*%beta)/length(beta)
          if(!is.null(X.test)){S2M[, i] <- (X.test%*%beta)/length(beta)}
          YY <- normalize(YY - (CVt + X.train%*%beta))
        }
        if(i > 1){
          condition = abs(cor(S1M[,i-1],S1M[,i]))
          if(i > 5 && condition >= 0.9999){break}
        }
      }
    }
    else{
      S1M = matrix(0, nrow=N, ncol=1)
      if(!is.null(X.test)){N2 <- dim(X.test)[1]}
      if(!is.null(X.test)){S2M = matrix(0, nrow=N2, ncol=1)}
      i = 1
      condition = 0
      while(!(i > 5 && condition >= 0.9999)){
        yt <- YY - CVt
        yt <- as.numeric(yt)
        r <- covar(yt, X.train)
        beta <- r/sx
        beta[is.na(beta)] <- Inf
        delSNP <- which(abs(beta) == Inf)
        if(length(delSNP)>0){
          beta <- beta[-c(delSNP)]
          S1M[, i] <- (X.train[,-c(delSNP)]%*%beta)/length(beta)
          if(!is.null(X.test)){S2M[, i] <- (X.test[,-c(delSNP)]%*%beta)/length(beta)}
          YY <- normalize(YY - (CVt + X.train[,-c(delSNP)]%*%beta))
        }
        else{
          S1M[, i] <- (X.train%*%beta)/length(beta)
          if(!is.null(X.test)){S2M[, i] <- (X.test%*%beta)/length(beta)}
          YY <- normalize(YY - (CVt + X.train%*%beta))
        }
        if(i > 1){
          condition = abs(cor(S1M[,i-1],S1M[,i]))
        }
        i = i+1
        S1M = cbind(S1M, rep( 0, N))
        if(!is.null(X.test)){S2M = cbind(S2M, rep(0, N2))}
      }
      i = i-1
    }

    S1M <- S1M[,1:i]
    if(!is.null(X.test)){S2M <- S2M[,1:i]}
    name.col <- unlist(lapply(1:i, function(x){paste0("CC", paste0(x, collapse=""), sep = "")}))
    colnames(S1M) <- name.col
    if(is.null(X.test)){S2M = NULL}else{colnames(S2M) <- name.col}
    return(list(CCtrain = S1M, CCtest = S2M))

  }

  Blink.LDRemove <- function(
    GDneo=NULL,
    LD=0.7,
    Porder=NULL,
    block=1000,
    LD.num =50)
  {
    Blink.LDRemoveBlock<-function(GDneo=NULL,LD=NULL,Porder=NULL){
      GDneo=as.matrix(GDneo)
      n=nrow(GDneo)
      corr=cor(GDneo)
      corr[is.na(corr)]=1
      corr[abs(corr)<=LD]=0
      corr[abs(corr)>LD]=1
      Psort=as.numeric(matrix(1,1,ncol(corr)))
      # print(ncol(corr))
      for(i in 2:ncol(corr)){
        p.a=Psort[1:(i-1)]
        p.b=as.numeric(corr[1:(i-1),i])
        index=(p.a==p.b)
        index[(p.a==0)&(p.b==0)]=FALSE
        if(sum(index)!=0) Psort[i]=0
      }
      seqQTN=Porder[Psort==1]
      #seqQTN=Porder[Psort==1,c(1,2),drop=F]
      return(seqQTN)
    }

    ##############
    GDneo = as.matrix(GDneo)
    SNP.index = apply(GDneo, 2, sd) != 0
    GDneo = GDneo[, SNP.index]
    Porder = Porder[SNP.index]
    l = block
    seqQTN=NULL
    lp=length(Porder)
    k=ceiling(lp/l)
    GDneo=as.matrix(GDneo)
    n=nrow(GDneo)

    for(i in 1:k){
      bottom=(i-1)*l+1
      up=l*i
      if(up>lp) up = lp
      Porderb=Porder[bottom:up]
      index = seq(bottom:up)
      GDneob = GDneo[,index]
      # cat("i is ",i,"\n")
      # print(length(index))
      seqQTNs = Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb)
      # print(seqQTN)
      seqQTN = append(seqQTN,seqQTNs)
      if(k >1){
        index1 = which(Porder %in% seqQTN)
        Porderb = Porder[index1]
        GDneob = GDneo[,index1]
        if(length(index1)>1){
          seqQTN = Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb)
        }else{
          seqQTN = Porderb
        }
      }
      if(LD.num < length(seqQTN)) break
    }
    rm(GDneob,Porderb)
    return(seqQTN)
  }

  ##remove NA
  missing = is.na(y.train) #index for missing phenotype
  if(length(which(missing==TRUE))!=0){
    y.train <- as.matrix(y.train[!missing,])
    X.train <- X.train[!missing,]
  }
  rm(missing)


  #BoostedCCR
  if(is.null(X.test)){
    print("--------------------------------------------------")
    if(CV == FALSE){
      if(!is.null(cc_num)){
        gg <- Boosted_CC(y.train = y.train, X.train = X.train, X.test = NULL, cc_num = cc_num)
        CCtrain0 <- gg$CCtrain[,]
        rownames(CCtrain0) <- rownames(X.train)
      }else{
        gg <- Boosted_CC(y.train = y.train, X.train = X.train, X.test = NULL)
        CCtrain0 <- gg$CCtrain[,]
        rownames(CCtrain0) <- rownames(X.train)
      }
    }
    if(CV == TRUE){
      CV1 <- myPC0[rownames(X.train),1]
      if(!is.null(cc_num)){
        gg <- Boosted_CC(y.train = y.train, X.train = X.train, X.test = NULL, CV1 = CV1, cc_num = cc_num)
        CCtrain0 <- gg$CCtrain[,]
        rownames(CCtrain0) <- rownames(X.train)
      }else{
        gg <- Boosted_CC(y.train = y.train, X.train = X.train, X.test = NULL, CV1 = CV1)
        CCtrain0 <- gg$CCtrain[,]
        rownames(CCtrain0) <- rownames(X.train)
      }
    }
    results <- list(CCtrain=CCtrain0)
    print("----------------------------------------------------------------------------------------------------")
  }


  if(!is.null(X.test)){
    X <- rbind(X.train, X.test)

    ##GLM___train___first SNP-selected
    CVglm <- cbind(rep(1,length(y.train)))
    colnames(CVglm) <- c("intercept")
    P <- apply(X.train, 2, function(x){
      Model <- RcppEigen::fastLmPure(y = y.train, X = as.matrix(cbind(x, CVglm), nrow = length(y.train)))
      return(2*pnorm(abs(Model$coefficients[1]/Model$se[1]), lower.tail = FALSE))})
    PP <- sort(P[which(-log10(P) >- log10(0.01/ncol(X.train)))])
    xzp <- length(PP)
    print("----------")

    ##PCA___geno___xzSNP
    if(xzp < 100 & xzp > 0){
      PCAall <- prcomp(X[,names(PP[1:length(PP)])])
      myPC0 <- PCAall$x
    }

    ##bin___train___second SNP-selected
    if(xzp >= 100){
      Porder <- as.matrix(PP)
      Psort <- Blink.LDRemove(GDneo=as.matrix(X.train[,rownames(Porder)]), Porder = Porder, LD = 0.7)
      indexQTN <- Porder[,1]%in%Psort
      PPnew <- Porder[indexQTN,]
    }
    print("--------------------")

    ##PCA___geno___bin
    if(xzp >= 100){
      PCAall <- prcomp(X[,names(PPnew[1:length(PPnew)])])
      myPC0 <- PCAall$x
    }
    print("------------------------------")


    if(xzp == 0){
      PCAnum=0
    }else{
      PCAnum <- ncol(myPC0)
    }

    #stepwise___PCA-selected
    if(PCAnum == 0){
      index <- NULL
    }
    if(PCAnum == 1){
      index <- colnames(myPC0)
    }
    if(PCAnum > 1 & PCAnum < 100){
      res <- StepReg::stepwise(y.train~. ,data = as.data.frame(cbind(y.train,myPC0[rownames(X.train),])), selection = "bidirection", sle = 0.01, sls = 0.01)
      index <- res$Varaibles[-1]
    }
    if(PCAnum >= 100){
      res<-StepReg::stepwise(y.train~. ,data=as.data.frame(cbind(y.train,myPC0[rownames(X.train), 1:100])), selection = "bidirection", sle = 0.01, sls = 0.01)
      index <- res$Varaibles[-1]
    }
    print("----------------------------------------")

    #BoostedCCR
    if(CV == FALSE){
      if(!is.null(cc_num)){
        gg <- Boosted_CC(y.train = y.train, X.train = X.train, X.test = X.test, cc_num = cc_num)
        CCtrain0 <- gg$CCtrain[,]
        CCtest0 <- gg$CCtest[,]
        rownames(CCtrain0) <- rownames(X.train)
        rownames(CCtest0) <- rownames(X.test)
      }else{
        gg <- Boosted_CC(y.train = y.train, X.train = X.train, X.test = X.test)
        CCtrain0 <- gg$CCtrain[,]
        CCtest0 <- gg$CCtest[,]
        rownames(CCtrain0) <- rownames(X.train)
        rownames(CCtest0) <- rownames(X.test)
      }
    }
    if(CV == TRUE){
      CV1 <- myPC0[rownames(X.train),1]
      if(!is.null(cc_num)){
        gg <- Boosted_CC(y.train = y.train, X.train = X.train, CV1 = CV1, X.test = X.test, cc_num = cc_num)
        CCtrain0 <- gg$CCtrain[,]
        CCtest0 <- gg$CCtest[,]
        rownames(CCtrain0) <- rownames(X.train)
        rownames(CCtest0) <- rownames(X.test)
      }else{
        gg <- Boosted_CC(y.train = y.train, X.train = X.train, X.test = X.test, CV1 = CV1)
        CCtrain0 <- gg$CCtrain[,]
        CCtest0 <- gg$CCtest[,]
        rownames(CCtrain0) <- rownames(X.train)
        rownames(CCtest0) <- rownames(X.test)
      }
    }
    print("--------------------------------------------------")

    #Merge PCA and CC
    if(xzp>0){
      trainPC <- as.matrix(myPC0[rownames(CCtrain0),index])
      testPC <- as.matrix(myPC0[rownames(CCtest0),index])
      colnames(trainPC) = index
      colnames(testPC) = index
    }else{
      trainPC <- NULL
      testPC <- NULL
    }
    CCtrain <- cbind(trainPC,CCtrain0)
    CCtest <- cbind(testPC,CCtest0)

    #Prepare stacking strategy
    X1 <- cbind(y.train, CCtrain)
    fold1 <- cut(seq(1,nrow(X1)), breaks = 5, labels = FALSE)
    Mt1 <- vector()
    Mt2 <- vector()
    Mt3 <- vector()
    Mt4 <- vector()
    Mtt1 <- matrix(0, 5, nrow(X.test))
    Mtt2 <- matrix(0, 5, nrow(X.test))
    Mtt3 <- matrix(0, 5, nrow(X.test))
    Mtt4 <- matrix(0, 5, nrow(X.test))
    print("------------------------------------------------------------")

    ##stacking base learners
    for (j in 1:5) {
      testIndex <- which(fold1 == j, arr.ind = TRUE)
      y.test1 <- X1[testIndex, 1]
      X.test1 <- X1[testIndex, -1]
      y.train1 <- X1[-testIndex, 1]
      X.train1 <- X1[-testIndex, -1]
      train1 <- cbind(y.train1, X.train1)
      colnames(train1) <- c("y.train1", colnames(X.train1))

      #Bridge
      Bridge <- grpreg::gBridge(as.matrix(X.train1), train1[,1], group=1:ncol(X.train1), family="gaussian")
      T1 <- predict(Bridge, as.matrix(X.test1), type="response", lambda= grpreg::select(Bridge, "BIC", df.method = "active")$lambda)
      Mt1 <- c(Mt1, T1)
      Mtt1[j,] <- predict(Bridge, as.matrix(CCtest), type="response",lambda= grpreg::select(Bridge, "BIC", df.method = "active")$lambda)

      #EN
      EN <- glmnet::glmnet(as.matrix(X.train1), train1[, 1], family = "gaussian", alpha = 0.5)
      T2 <- predict(EN, newx = as.matrix(X.test1), type = "response")
      n <- which(cor(y.test1, T2[,-1]) == max(cor(y.test1, T2[,-1])))
      T2 <- T2[, n+1]
      Mt2 <- c(Mt2, T2)
      TT2 <- predict(EN, newx = as.matrix(CCtest), type = "response")
      if(is.null(y.test)){
        Mtt2[j,] <- rowMeans(TT2)
      }
      if(!is.null(y.test)){
        n2 <- which(cor(y.test,TT3[,-1])==max(cor(y.test, TT3[,-1])))
        n2 <- N1[1]
        Mtt2[j,] <- TT2[, n2+1]
      }

      #grlasso
      grlasso <- grpreg::grpreg(as.matrix(X.train1), y.train1, group = 1:ncol(X.train1), penalty = "grLasso", family = "gaussian")
      T3 <- predict(grlasso, as.matrix(X.test1), type = "response", lambda = grpreg::select(grlasso)$lambda)
      Mt3 <- c(Mt3, T3)
      Mtt3[j,] <- predict(grlasso, as.matrix(CCtest), type = "response", lambda = grpreg::select(grlasso)$lambda)

      #gbm
      gbm=gbm::gbm(y.train1 ~ ., data = as.data.frame(train1), distribution = 'gaussian', interaction.depth = 5, n.trees = 250, shrinkage = 0.1, n.minobsinnode = 5)
      T4=predict(gbm, as.data.frame(X.test1), n.trees = 250)
      Mt4 <- c(Mt4, T4)
      Mtt4[j,] <- predict(gbm, as.data.frame(CCtest), n.trees = 250)
    }
    print("--------------------------------------------------------------------------------")

    ##stacking meta learner
    #上面循环结束，产生新meta变量
    M1 <- as.matrix(Matrix::colMeans(Mtt1))
    M2 <- as.matrix(Matrix::colMeans(Mtt2))
    M3 <- as.matrix(Matrix::colMeans(Mtt3))
    M4 <- as.matrix(Matrix::colMeans(Mtt4))

    #合并最后一轮进行预测的meta训练和测试数据集
    MX.train <- cbind(Mt1,Mt2,Mt3,Mt4)
    MX.test <- cbind(M1,M2,M3,M4)
    colnames(MX.test) <- colnames(MX.train)

    MX.data = as.data.frame(cbind(y.train, MX.train))
    colnames(MX.data) <- c("y.train",colnames(MX.train))

    #meta learner
    lasso.model = elasticnet::enet(as.matrix(MX.data[,-1]), MX.data[,1], lambda = 0)
    MX.T <- predict(lasso.model,MX.test,type='fit')$fit
    if(is.null(y.test)){
      MX.T <- rowMeans(MX.T)
    }
    if(!is.null(y.test)){
      N <- which(cor(y.test,MX.T[,-1]) == max(cor(y.test, MX.T[,-1])))
      MX.T <- MX.T[, N+1]
    }
    GEBV <- MX.T
    print("----------------------------------------------------------------------------------------------------")

    results <- list(CCtrain=CCtrain0, CCtest=CCtest0, GEBV=GEBV)
  }

  print("-------------------------------------- BCCRstacking runs ends --------------------------------------")
  return(results)

}
