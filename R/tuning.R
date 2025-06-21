# one way to compute gamma: quartiles of between class euclidean distance
get.gamma0=function(formula, dat.train){
    tmp=model.frame(formula, dat.train)
    x1=model.matrix(formula, tmp[tmp[,1]==1,])[,-1,drop=FALSE]
    x2=model.matrix(formula, tmp[tmp[,1]==0,])[,-1,drop=FALSE]
    dist. = getK(X=x1,"euclidean",X2=x2)
    gamma0=.5/quantile(dist., c(.9,.75,.5))
    gamma0
}


# return optimal tuning param
tune.it=function(formula, dat.train, dat.tune, method, kernel, verbose=TRUE, step.size=2){
    
    kernel=substr(kernel,1,1)
    
    gammas=get.gamma0(formula, dat.train)
    if(verbose) print(gammas)
    
    f1=function(params) {
        sapply(1:nrow(params), simplify="array", function(row.){
            param=params[row.,]
            if (verbose) print(param)
            if (method=="rauc") {
                fit = rauc (formula, dat=dat.train, s=1, lambda=param$lamcost, kernel=kernel, para=param$gamma, verbose=verbose)
                auc.tune = fast.auc(predict(fit, dat.tune), dat.tune$y)
# the following is blocked out b/c svmw is not on cran, and calling svmw::svml is not allowed when rcmdcheck aucm
#            } else if (method=="svm"){
#                if (kernel=="q") {
#                    fit = svmw::svml (formula, dat.train, cost=param$lamcost, kernel="p", degree=2, fitted=FALSE)
#                } else {
#                    fit = svmw::svml (formula, dat.train, cost=param$lamcost, kernel=kernel, gamma=param$gamma, fitted=FALSE)
#                }                
#                auc.tune = ifelse (is.null(dat.tune), -svmw::gacv(fit), fast.auc(predict(fit, dat.tune), dat.tune$y))
            } else if (method=="dcsauc"){
                fit = dcsauc (formula, dat.train, lambda=param$lamcost, kernel=kernel, para=param$gamma, verbose=verbose)
                auc.tune = attr(predict(fit, dat.tune), "auc")
            } else if (method=="srauc"){
                fit = srauc (formula, dat.train, lambda=param$lamcost, kernel=kernel, para=param$gamma, verbose=verbose)
                auc.tune = attr(predict(fit, dat.tune), "auc")
            } else stop("tune.it: something is wrong")
            res = c(lamcost=param$lamcost, gamma=unname(param$gamma), tune=auc.tune)
            if (verbose) print(res)
            res
                  
        })    
    }
    
            
    # fit model on an initial grid of 7 x 1 tuning parameters
    lamcosts=step.size^seq(-3,3,by=1)     
    params=expand.grid(lamcost=lamcosts, gamma=gammas[2])    
    res=f1(params)
    best.lamcost = lamcosts[which.max(res["tune",])]
    
    # search for best lamcost if not in the initial grid
    if (verbose) print("extend lambda")
    if (best.lamcost==min(lamcosts) | best.lamcost==max(lamcosts)) {
        incr = ifelse (best.lamcost==min(lamcosts), 1/step.size, step.size)
        while (TRUE){
            params=expand.grid(lamcost=c(best.lamcost*incr), gamma=gammas[2])
            res2=f1(params)
            if (res2["tune",1]<=max(res["tune",])) {
                break
            } else {
                best.lamcost=best.lamcost*incr
                res=cbind(res, res2)
            }
        } 
    }
    
    # construct a 3 x 3 grid of tuning parameters around best.lamcost and gammas and fit model
    if (verbose) print("3x3 grid")
    lamcosts=c(best.lamcost/sqrt(step.size), best.lamcost, best.lamcost*sqrt(step.size))
    if (kernel=="r") {
        params=expand.grid(lamcost=lamcosts, gamma=gammas)
    } else if (kernel=="l" | kernel=="q") {
        params=expand.grid(lamcost=lamcosts, gamma=gammas[2])
    } else stop("kernel not supported: "%.%kernel)
    res=f1(params)
    
    tmp=res[,which.max(res["tune",])]
    if(kernel=="r") {
        res = tmp[c("lamcost","gamma")]
    } else if (kernel=="l" | kernel=="q") {
        res = tmp["lamcost"] # add drop=FALSE does not change things
    } else stop("kernel not supported: "%.%kernel)
    
    # doing the following violates some programming principle, dat.test is no longer an argument
#    if (!is.null(dat.test)) {
#        if (method=="svm"){
#            fit = svmw::svml (formula, dat.train, cost=res["lamcost"], kernel=kernel, gamma=res["gamma"], fitted=FALSE)
#            auc.train = fast.auc(predict(fit, dat.train), dat.train$y)
#            auc.test = fast.auc(predict(fit, dat.test), dat.test$y)
#        } else stop("method not supported: "%.%method)       
#        res=c(res, "auc.train"=auc.train, "auc.test"=auc.test)
#    }
    
    
    return (res)
}


# dat is a data frame, one of the columns has to be "y": 0/1, 1 for case
sample.for.cv=function(dat, v, seed){
    set.seed(seed)
    n1=sum(dat$y==1)
    n2=sum(!dat$y==1)
    test.set = c(sample(which (dat$y==1), round(n1/v)), sample(which (dat$y==0), round(n2/v)))
    train.set = setdiff(1:nrow(dat), test.set)
    list("train"=train.set, "test"=test.set)
}
