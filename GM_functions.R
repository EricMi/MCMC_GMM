#######################'
#' 1.  Common functions 
#######################' 

gmllk <- function(x , Mu , Sigma , p) {
    #' Log-likelihood in the Gaussian mixture model. 
    #' x: dataset: a n*d matrix for n points with d features each.
    #' Mu: a k*d matrix with k the number of components: the centers
    #' Sigma: a d*d*k array:: the convariance matrices.
    #' p: a vector of length k: the mixture weights
    #' returns:  the log-likelihood (single number)
    k <- length(p)
    if(is.vector(x)) {
        x <- matrix(x, nrow=1)}
    n <- nrow(x)
    mat_dens <- vapply(1:k, function(j) {
        dmnorm(x, mean = Mu[j,], varcov = Sigma[,,j], log=FALSE)
    }, FUN.VALUE = numeric(n)) ##  n rows, k columns.
    if(is.vector(mat_dens)) {
        mat_dens <- matrix(mat_dens, nrow = 1)
    }
    vect_dens <-   mat_dens%*%matrix(p,ncol=1) ## vector of size n
    return(sum(log(vect_dens)))
}

gmcdf <- function(x , Mu , Sigma , p) {
    #' multivariate cumulative distribution function in a GMM. 
    #' x: a single point (vector of size d)
    #' Mu, Sigma, p: see gmllk.
    #' returns: the cdf at point x. 
    k <- length(p)
    vect_cdf <- vapply(1:k, function(j){
        pmnorm(x, mean = Mu[j,], varcov = Sigma[,,j])
    }, FUN.VALUE = numeric(1))
    return(sum(p*vect_cdf))
}


initem <- function(x , k) {
    #' Initialisation for EM and VB based on kmeans. 
    #'x: dataset: a n*d matrix for n points with d features each.
    #'k: number of components for the inferred mixture
    #' returns: a list with entries p, Mu, Sigma: respectively a vector of size k (weights), a k*d matrix (centers) and a d*d*k array (empirical covariance matrix)
    init <- kmeans(x = x, centers = k, iter.max = 100, nstart = 1,
                   algorithm = c("Hartigan-Wong"), trace=FALSE)
    Mu <- init$centers
    d <- ncol(x)
    Sigma <- array(dim=c(d,d,k))
    p <- rep(0,k)
    for( i in (1:k)) {
        inds = which(init$cluster==i)
        n = length(inds)
        tildeX = t(t(x[inds,]) -Mu[i,])  
        sig = 1/n * t(tildeX) %*% tildeX
        Sigma[,,i] <- sig
        p[i] <-  n/nrow(x)
    }
    return(list(p = p, Mu = Mu, Sigma = Sigma ))
}


draw_sd <- function(mu , sigma) {
    #' draws  an ellipsoid  around the mean of a gaussian distribution
    #' which corresponds to the density level set of the univariate
    #' 0.95 quantile.
    #' mu: vector of size d the dimension
    #' sigma: a d*d covariance matrix.
    #' returns: a 2*100 matrix containing abscissas and ordinates of
    #' the ellipsoid to be drawn. 
    L <-  chol(sigma)
    angles <- seq(0, 2*pi, length.out=100)
    U <- 1.64* rbind(cos(angles), sin(angles))
    X <- mu + t(L) %*% U
    return(X)
}

nanDetector <- function(X) {
    #' returns TRUE if X contains NaNs
   # examine data frames
   if(is.data.frame(X)) { 
       return(any(unlist(sapply(X, is.nan))))
   }
   #  examine vectors, matrices, or arrays
   if(is.numeric(X)) {
       return(any(is.nan(X)))
   }
   #  examine lists, including nested lists
   if(is.list(X)) {
       return(any(rapply(X, is.nan)))
   }
   return(FALSE)
}

wrapper <- function(x , y , FUN, ...) {
    #' applies a function on a grid with abscissas x, y.
    #' x, y: vectors of same length.
    sapply(seq_along(x), FUN = function(i){FUN(x[i], y[i],...)})
}

########################'
##' 2.  EM functions 
########################'

estep <- function(x , Mu , Sigma , p) {
    #' E-step in the EM algorithm.
    #' x, Mu, Sigma, p: see gmllk. 
    #' returns:  a n*k matrix (responsibilities)
    k <- length(p)
    N <- nrow(x)
    Respons <- matrix(0, ncol=k, nrow=N)
    for(i in 1:k) {
        Respons[,i] <- p[i] * dmnorm(x, mean=Mu[i,], varcov=Sigma[,,i])
    }
    respons <- Respons / apply(Respons, 1, sum)
    return(respons)
}

mstep <- function(x , respons) {
    #' M-step in the EM algorithm.
    #' x: the dataset, see gmllk.
    #' respons: a n*k matrix as returned by estep.
    #' returns the optimized parametres (Mu, Sigma, p) as a list 
    N <- nrow(respons)
    d <- ncol(x)
    k <- ncol(respons)
    NK <- apply(respons, 2, sum)
    p <- NK / sum(NK)
    mustar <- t(respons) %*% x / NK
    Sigstar=array(dim=c(d,d,k))
    for( i in (1:k)){
        x_centered = t(t(x) - mustar[i,])
        weighted_cov_k <- vapply(1:N, function(j) {
            respons[j,i] * x_centered[j,] %*% t(x_centered[j,])
        }, FUN.VALUE = numeric(d^2)) ##  array of dimension (N,d,d)
        Sigstar[,,i] <- apply(weighted_cov_k, 1, sum) / NK[i]
    }
    return(list(Mu = mustar, Sigma = Sigstar, p = p))    
}

emalgo <- function(x, k, tol=1e-6) {
    #' EM algorithm.
    #' x: the dataset, see gmllk.
    #' The tolerance parameter used for the stopping rule:
    #' the algorithm stops when the log-likelihood increase is less than tol.
    #' returns: if  T=number of iteration,  a list composed of
    #' - Muarray: a k*d*T array (centers)
    #' - Sigmaarray: a d*d*k*T array (covariance matrices)
    #' - pmat: a k*d matrix (weights)
    #' - objective : a vector of size T: log-likelihood values.
    #' - last: a list composed of the  parameters  returned by the latest M-step. 
    init <- initem(x = x, k = k)
    d <- ncol(x)
    res <- list(Muarray = array(dim=c(k, d, 0) ),
                Sigmaarray = array(dim=c(d, d,k, 0) ),
                pmat = matrix(nrow = k, ncol=0))
    obj <- gmllk( x = x, Mu = init$Mu, Sigma = init$Sigma, p = init$p)
    ## obj: objective function (to be maximized)  used as a stopping criterion. 

    
    ## initialisation:
    current <- init
    ## storing initial values:
    res$Muarray <- abind(res$Muarray, current$Mu, along =3)
    res$Sigmaarray <- abind(res$Sigmaarray, current$Sigma, along=4)
    res$pmat <- cbind(res$pmat, current$p)
    
    niter <- 0
    continue <- TRUE 
    while(continue) {
        niter <- niter+1
        respons <- estep(x=x, Mu=current$Mu, Sigma=current$Sigma, p=current$p)
        current <- mstep(x=x, respons=respons)
                
        res$Muarray <- abind(res$Muarray, current$Mu, along = 3)
        res$Sigmaarray <- abind(res$Sigmaarray, current$Sigma, along = 4)
        res$parray <- cbind(res$pmat, current$p)
 
        newobj <- gmllk( x=x, Mu = current$Mu, Sigma = current$Sigma, p=current$p)
        if(newobj - obj[niter] < tol){continue =FALSE}
        obj <- c(obj,newobj)
    }    
    res$objective <- obj
    res$last <- current
    return(res)
}

################################'
##' 3. Variational Bayes functions 
################################' 

vbMstep <- function(x , respons , alpha0 ,  W0inv , nu0 , m0 , beta0)
    #' x: the data. A n*d matrix 
    #' respons: current q(z): a n*k matrix (responsibilities)
    #' alpha0: vector of size k: dirichlet prior parameter on p
    #' W0inv, nu0: parameters for the Wishart prior on Lambda.
    #' W0inv: d*d matrix, inverse of the Wishart parameter.
    #' nu0 > d-1:  is a real.
    #' m0 : mean parameter (d vector) for the Gaussian-Wishart prior on  mu
    #' beta0: scale parameter for the gaussian-wishart  prior on mu (>0)
    ##' 
    #' returns: a list made of ( alpha , Winv, Nu , M , Beta):  optimal parameters for q(p),
    #' q(mu_j, Lambda_j), j=1, ...,k: 
    #' alpha: k-vector ; Winv: d*d*k array ; Nu: a k-vector ; M: k*d matrix ;
    #' Beta: k-vector                   
{
    K <-  ncol(respons)
    NK <- apply(respons, 2, sum) # a vector of size k
    NK <- sapply(NK, function(x){max(x, 1e-300)}) ## avoids divisions by zero
    alpha <- alpha0 - 1 + NK
    Nu <- nu0 + NK
    Beta <- beta0 + NK
    barx <- (t(respons) %*% x )/NK  # k*d matrix: weighted mean of each cluster 
    M <-  (matrix(beta0 * m0, nrow=K, ncol=d, byrow=TRUE) + barx * NK) / Beta
    d <- ncol(x)
    Winv <- array(dim=c(d,d,K))
    for( j in (1:K)){
        tildeX = t(t(x) - barx[j,])   
        Sj = 1/NK[j] * t(respons[,j] * tildeX) %*% tildeX
        Winv[,,j] <- W0inv + NK[j] * Sj + (barx[j,] - m0) %*% t((barx[j,] - m0)) * (beta0 * NK[j]) / (beta0 + NK[j])
    }

    return(list(alpha = alpha, Winv = Winv, Nu = Nu, M= M, Beta =Beta)) 
}


vbEstep <- function(x, alpha, Winv, Nu, M, Beta)
    #' computation of the variational responsibilities. 
    #' x: the data. A n*d matrix
    #' alpha: a k vector: current dirichlet parameter for q(p)
    #' Winv : a d*d*k array: current inverses of the W parameter for the Wishart q(Lambda)
    #' Nu: a k vector: current degrees of freedom parameter for the Wishart q(Lambda)
    #' M: a k*d matrix: current mean parameters for the Gaussian q(Mu | Lambda)
    #' Beta: a k vector: current scale parameters for the Gaussian q(Mu | Lambda)
    #' returns: a n*k matrix: the responsibilities for each data point.  
{
    d <- ncol(M)
    k <- length(alpha)
    N <- nrow(x)
    Eloglambda <-  # k vector
        sapply(1:k, function(j){
            sum(digamma( (Nu[j] + 1 - (1:d) )/2) )+ d * log(2) - log(det(Winv[,,j]))
        })
    Elogrho <- # k vector
        digamma(alpha) - digamma(sum(alpha))    
    Equadratic <- # k*N  matrix
        d / Beta  + Nu * t( sapply(1:k, function(j){ ## a N * k matrix
            Wj <- solve(Winv[,,j])
            sapply(1:N, function(n){# a N vector
                t(M[j,] -x[n, ]) %*% Wj %*% (M[j,] -x[n, ])})
        }))
    logResponsT <- Elogrho + Eloglambda / 2 - log(2*pi) * d / 2 - Equadratic / 2
    logRespons <- t(logResponsT) ## N * k
    logRespons <- logRespons - apply(logRespons, 1, max) #' avoids numerical precision loss. 
    respons <- exp(logRespons) ##  N * k matrix
    Z <- apply(respons, 1 , sum ) # N vector
    respons <- respons / Z ##N * k matrix
    return(respons)
}

lowerBound <- function(x, respons, alpha, Winv, Nu, M, Beta,
                       alpha0, W0inv, nu0, m0, beta0)
    #' Lower bound on the log likelihood from the classical decomposition.
    #' Used as a stopping criterion in the VB-algorithm.
    #' x, respons, alpha0, W0inv, nu0, m0, beta0: see vbMstep
    #' alpha, Winv, Nu, M, Beta: see vbEstep. 
    #' returns: a single number: the value of the lower bound. 
{
    d <- ncol(x)
    ## compute tildeLambda = Eloglambda
    k <- length(alpha)
    Eloglambda <-  # k vector
        sapply(1:k, function(j){
            sum(digamma( (Nu[j] + 1 - (1:d) )/2) )+ d * log(2) - log (det(Winv[,,j]))
        })
    Elogrho <- # k vector
        digamma(alpha) - digamma(sum(alpha))    

    NK <- apply(respons, 2, sum)  # k vector
    NK <- sapply(NK, function(x){max(x, 1e-300)})
    barx <- (t(respons) %*% x )/NK  # k*d matrix
    W <- Winv
    for(j in 1:k)
    {
        W[,,j] <- solve(Winv[,,j])
    }
    logCalpha0 <- lgamma(d*alpha0) -  d * lgamma(alpha0)
    logCalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha)) # k vector
    lBwish0 <-  (nu0/2)* log(det(W0inv)) -( (nu0 * d/2)* log(2)  +
                                           (d*(d-1)/4) * log( pi)  + 
                                         sum(lgamma((nu0 + 1 - (1:d))/2)) ) 
    lBwish <- sapply(1:k, function(j){
         (Nu[j]/2)* log(det(Winv[,,j])) -( (Nu[j] * d/2)* log(2)  +
                                           (d*(d-1)/4) * log( pi)  + 
                                         sum(lgamma((Nu[j] + 1 - (1:d))/2)) )
    }) # k vector
    
    Hqlambda <-  - lBwish - (Nu - d - 1)/2 * Eloglambda + Nu * d/2 # k vector
    llk <- 0
    for (j in 1:k)
    {
        tildeX = t(t(x) - barx[j,])   
        Sj = 1/NK[j] * t(respons[,j] * tildeX) %*% tildeX
        Wj <- W[,,j]
                    
        llk = llk + NK[j]*
            (Eloglambda[j] - d/Beta[j] - Nu[j] * sum(diag(Sj %*%Wj )) -
             Nu[j ] * t(barx[j,] - M[j,]) %*% Wj %*% (barx[j,] - M[j,]) - 
             d * log(2*pi)) 
    }
    llk <- llk/2
    llkZ <- sum(Elogrho * t(respons))
 
    llkrho <- logCalpha0 + (alpha0-1)*sum(Elogrho)

    llkmulambda <- 1/2 * sum (d * log(beta0/(2*pi)) + Eloglambda - d * beta0/Beta -
                              beta0 * sapply(1:k, function(j){
                                  Nu[j] * t(M[j,] - m0) %*% W[,,j]%*%(M[j,] - m0)})) +
        k * lBwish0 +
        (nu0 -d-1)/2 *sum(Eloglambda ) -
        1/2 * sum(sapply(1:k, function(j){Nu[j] * sum(diag(W0inv %*% W[,,j])) }))
    
    inds <- which(is.infinite(log(respons)))
    respons[inds] <-  1e-300
    llqZ <- sum(respons * log(respons))
    llqrho <- sum((alpha-1) * Elogrho) + logCalpha
    llqmulambda <-  sum( 1/2 * Eloglambda + d/2 * log(Beta/(2*pi) ) - d/2 - Hqlambda)
    res <- llk + llkZ + llkrho+ llkmulambda - llqZ - llqrho - llqmulambda
    if(nanDetector(res)) {stop("NaNs detected!\n")}
    return(res)

}


vbalgo <- function(x, k, alpha0,  W0inv, nu0, m0, beta0, tol=1e-5)
    #' x: the data. n*d matrix
    #' k: the number of mixture components. 
    #' alpha0, W0inv, nu0, m0, beta0: prior hyper-parameters, see vbMstep.
    #' returns: a list composed of (alphamat,  Winvarray, Numat, Marray, Betamat, responsarray):
    #'   optimal parameters for q(p), q(mu_j, Lambda_j), j=1, ...,k
    #'   alphamat: K* Tmatrix,  Winvarray: d*d*T array,  Numat: a k*T matrix-vector, 
    #'   Marray: k*d*T array,  Betamat: k*T matrix, responsarray: n*k*T matrix, 
    #'   where T is the number of steps.
{
    N <- nrow(x)
    init <-  initem(x=x,k=k)
    res <- list(alphamat=matrix(nrow=k, ncol=0),
                Winvarray = array(dim=c(d,d,k,0)),
                Numat = matrix(nrow=k, ncol=0),
                Marray= array(dim=c(k,d,0) ),
                Betamat = matrix(nrow=k, ncol=0),
                responsarray = array(dim=c(N,k, 0)),
                lowerbound = c()
                )
    d <- ncol(x)
    Winvstart <- array(dim=c(d,d,k))
    for(j in 1:k){
        Winvstart[,,j] <- init$p[j] * N *  init$Sigma[,,j]
        }
    current <- list( alpha = N * init$p, 
                    Winv = Winvstart, 
                    Nu = N* init$p, 
                    M = init$Mu,
                    Beta = N * init$p)
    ## current: current lit of hyper parrameters fir the variational distribution
    
    continue <- TRUE
    niter <- 0
    while(continue){        
        niter <- niter+1
        respons <- vbEstep(x=x, alpha=current$alpha, Winv=current$Winv,
                           Nu=current$Nu, M=current$M, Beta=current$Beta)
        if(nanDetector(respons)) {stop("NaNs detected!\n")}
        vbOpt <- vbMstep(x=x, respons=respons, alpha0=alpha0, W0inv=W0inv,
                         nu0=nu0, m0=m0, beta0=beta0)
        if(nanDetector(vbOpt)) {stop("NaNs detected!\n")}
        current <- vbOpt
        res$alphamat <- cbind(res$alphamat, current$alpha)
        res$Winvarray <- abind(res$Winvarray, current$Winv,along=4)
        res$Numat <- cbind(res$Numat, current$Nu)
        res$Marray <- abind(res$Marray, current$M,along=3)
        res$Betamat <- cbind(res$Betamat, current$Beta)
        res$responsarray <- abind(res$responsarray, respons,along=3)
        lastvalue <- lowerBound(x=x, respons=respons,
                                alpha=current$alpha,
                                Winv = current$Winv, Nu = current$Nu,
                                M = current$M, Beta=current$Beta,
                                alpha0 = alpha0, W0inv = W0inv,
                                nu0=nu0, m0= m0 , beta0 = beta0)
        res$lowerbound <- c(res$lowerbound, lastvalue)
        
        if(niter>=2){
            if( niter == 200  ||  lastvalue - res$lowerbound[niter-1] < tol)
            {continue <- FALSE}
        }
    }
        return(res)
        
}

vbPredictive <- function(x, alpha, Beta, M, Winv, Nu)
    #' predictive density  based on the VB approximation
    #' x: a single point (vector): where to evaluate the density. 
    #' alpha, Winv, Nu, M, Beta: see vbEstep
    #' returns the value of the posterior predictive (= posterior mean of the mixture density)
    #'    at point x. 

{
    k <- length(alpha)
    d <-  length(x)
    vect_dens <- vapply(X= 1:k, FUN= function(j){
        W <- Winv[,,j]
        L <-  (1 + Beta[j])/   ((Nu[j] + 1 - d) * Beta[j])  *   W
        L <- 1/2 * (L + t(L))
        return(dmt(x=x, mean = M[j,], S = L, df = Nu[j] + 1 - d, log=FALSE))} ,
        FUN.VALUE= numeric(1) )
    return(sum(alpha * vect_dens) / sum(alpha))
}


vbPredictiveCdf <- function(x, alpha, Beta, M, Winv, Nu)
    #' predictive cumulative distribution function based on the VB approximation.
    #' x: a single point (vector): where to evaluate the cdf. 
    #' alpha, Winv, Nu, M, Beta: the VB posterior parameters,  see vbEstep
    #' returns the value of the variational posterior predictive cdf
    #'  (= mean of the mixture cdf under the variational posterior predictive)
    #'   at point x. 
{
    k <- length(alpha)
    d <-  length(x)
    vectcdf <- vapply(X= 1:k, FUN= function(j){
        W <- Winv[,,j]
        L <-  (1 + Beta[j])/ ((Nu[j] + 1 - d) * Beta[j])  *   W
        L <- 1/2 * (L + t(L)) ## ensures symmetry despite numerical errors
        return(pmt(x=x, mean = M[j,], S = L, df = Nu[j] + 1 - d, log=FALSE))} ,
        FUN.VALUE= numeric(1) )

    return(sum(alpha * vectcdf) / sum(alpha))
}

##############################################'
###' 4. MCMC functions 
##############################################'

dprior <- function( Mu, Sigma, p,
                   hpar = list( alpha0=rep(1, length(p)),
                               m0 = rep(0, ncol(Mu)), beta0 = 1, 
                               W0 = diag(ncol(Mu)), nu0 = ncol(Mu)))
    #'log-prior density on (Mu, Sigma, p)
    #' Mu, Sigma, p: see gmllk
    #' hpar: a list of hyper-parameters composed of
    #' - alpha0: a k vector: dirichlet prior on p
    #' - m0: a d vector: mean parameter for the Gaussian-Wishart prior on Mu
    #' - beta0: a single number >0: scale parameter for the Gaussian-Wishart prior on Mu
    #' - W0: covariance parameter for the inverse-wishart distribution on Sigma
    #' - nu0: degrees of freedom >d-1 for the wishart distribution on Sigma.  
    
{
    d <- ncol(Mu)
    k <- length(p)
    prior_p <- ddirichlet(p, alpha=hpar$alpha0, log = TRUE)
    prior_MuSigma <- sum(sapply(1:k, function(j){
        dnorminvwishart(mu = Mu[j,], mu0 = hpar$m0, lambda = hpar$beta0,
                        Sigma = Sigma[,,j], S = hpar$W0, nu = hpar$nu0,
                        log = TRUE)}))
    return(prior_p + prior_MuSigma)
}

rproposal <- function( Mu, Sigma, p, ppar=list(var_Mu = 0.1,
                                               nu_Sigma = 10,
                                               alpha_p = 10))
    #' random generator according to a proposal kernel centered at the current value.
    #' Mu, Sigma, p: current mixture parameters, see gmllk.
    #' ppar: a list made of :
    #' - var_Mu: variance parameter for the gaussian kernel for Mu.
    #' - nu_Sigma: degrees of freedom for the Wihart kernel for Sigma
    #' - alpha_p: concentration aprameter for the Dirichlet kernel for p
    #' returns: a list fo proposal parameters: (Mu, Sigma, p), where 
    #'   p ~ dirichlet(alpha) with mean = alpha/sum(alpha) = current p and
    #'   concentration parameter sum(alpha) = alpha_p.
    #' Mu : a k*d matrix and Sigma: a d*d*k array: 
    #'   for j in 1:k,  Mu[j,]~ Normal(mean= current Mu[j,], covariance = var_Mu*Identity)
    #'   Sigma[,,j]~ Wishart(W = 1/nu_Sigma * current Sigma[,,j] ; nu = nu_Sigma)
    
{
    d <- ncol(Mu)
    k <- length(p)
    alphaProp <- sapply(ppar$alpha_p * p, function(x){max(x,1e-30)})
    ## this avoids numerical errors

    p <- rdirichlet(n=1, alpha = alphaProp)
    p <- sapply(p, function(x){max(x,1e-30)})
    p <- p/sum(p)
    for(j in (1:k)) {
        Mu[j,] <- rmvn(n=1, mu=Mu[j,], Sigma=ppar$var_Mu * diag(d))
        Sigma[,,j] <- rwishart(nu=ppar$nu_Sigma, S=Sigma[,,j] / ppar$nu_Sigma)
    }
    return(list(Mu = Mu, Sigma = Sigma, p = p))
}

MHsample <- function(x, k, nsample,
                     init=list(Mu = matrix(0,ncol=ncol(x), nrow=k ),
                               Sigma = array(rep(diag(ncol(x)), k),
                                             dim=c(ncol(x), ncol(x), k)),
                               p = rep(1/k, k)),
                     hpar= list( alpha0=rep(1, length(p)),
                                m0 = rep(0, ncol(Mu)), beta0 = 1, 
                                W0 = diag(ncol(Mu)), nu0 = ncol(Mu)),
                     ppar = list(var_Mu = 0.1,
                                 nu_Sigma = 10,
                                 alpha_p = 10) ) {
    #' x: the data. A n*d matrix.
    #' k: the number of mixture components.
    #' nsample: number of MCMC iterations
    #' init: starting value for the the MCMC. Format: list(Mu, Sigma, p), see gmllk for details
    #' hpar: a list of hyper-parameter for the prior: see dprior.
    #' ppar: a list of parameter for the proposal: see rproposal.
    #' returns: a sample produced by the Metropolis-Hastings algorithm, together with
    #' the log-posterior density (unnormalized) across iterations, and number of acepted proposals.  as a list composed of
    #' - p: a k*nsample matrix
    #' - Mu: a k*d*nsample array
    #' - Sigma: a d*d*k*nsample array
    #' - lpostdens: the log posterior density (vector of size nsample)
    #' - naccept! number of accepted proposals. 

    d <- ncol(x)
    output <- list(p = matrix(nrow=k, ncol=nsample),
                   Mu = array(dim = c(k, d, nsample)),
                   Sigma = array(dim = c(d, d ,k, nsample)),
                   lpostdens = rep(0, nsample),
                   naccept = 0
                   )
    current <- init
    current$lpost <- gmllk(x=x, Mu=current$Mu,
                            Sigma = current$Sigma, p=current$p) +
        dprior(Mu = current$Mu, Sigma = current$Sigma, p = current$p,
               hpar = hpar)
    ## lpost: logarithm of the unnormalized psoterior density.
    
    for (niter in 1:nsample){
        proposal <- rproposal(Mu = current$Mu, Sigma = current$Sigma, p=current$p,
                              ppar = ppar)

        proposal$lpost <- gmllk(x=x, Mu=proposal$Mu,
                                Sigma = proposal$Sigma, p=proposal$p) +
            dprior(Mu = proposal$Mu, Sigma = proposal$Sigma, p = proposal$p,
                   hpar = hpar)
        
        llkmoveSigma <- sum(vapply(1:k, FUN = function(j){
            dwishart(Omega =proposal$Sigma[,,j], nu=ppar$nu_Sigma,
                     S = 1/ppar$nu_Sigma * current$Sigma[,,j] , log=TRUE)},
            FUN.VALUE = numeric(1)))

        llkbackSigma <- sum(vapply(1:k, FUN = function(j){
            dwishart(Omega =current$Sigma[,,j], nu=ppar$nu_Sigma,
                     S = 1/ppar$nu_Sigma * proposal$Sigma[,,j] , log=TRUE)},
            FUN.VALUE = numeric(1)))
        
        alphaPropmove <- sapply(ppar$alpha_p * current$p, function(x){max(x,1e-30)})
        alphaPropback <- sapply(ppar$alpha_p * proposal$p, function(x){max(x,1e-30)})
        lacceptratio <-  ## logarithm of the acceptance ratio.
            min(0, (proposal$lpost + ddirichlet(current$p, alphaPropback, log=TRUE) + llkbackSigma) - 
                   (current$lpost + ddirichlet(proposal$p, alphaPropmove, log=TRUE) + llkmoveSigma))

        U <- runif(1)
        if(U < exp(lacceptratio)){
            current <- proposal
            output$naccept <- output$naccept + 1
        }
        output$p[,niter] <- current$p
        output$Mu[,,niter] <- current$Mu
        output$Sigma[,,,niter] <- current$Sigma
        output$lpostdens[niter] <- current$lpost            
    }
    
    return(output)
} 

cdfTrace <- function(x , sample , burnin = 0 , thin = 1) {
    #' Traces the evolution of the gmcdf at point x through the MCMC iterations.
    #'  Can be used for convergence monitoring. 
    #' x: a single point (vector of size d)
    #' burnin, thin: see MHpredictive
    #' returns: a vector of length [ (nsample - burnin )/thin ]

    k <- nrow(sample$p)
    nsample <- ncol(sample$p)
    inds <- (burnin+1):nsample
    inds <- inds[inds%%thin==0]
    output <- vapply(inds , function(niter){
        gmcdf(x = x, Mu = sample$Mu[,,niter], Sigma = sample$Sigma[,,,niter], p = sample$p[,niter])
    },FUN.VALUE = numeric(1))
    
    return(output)
}


MHpredictive <- function(x , sample , burnin=0, thin=1) {
    #' posterior predictive density computed from MH output. 
    #' x: vector size d (single point)
    #' sample: output from the MCMC algorithm should contain
    #'    entries Mu, Sigma, p as in MHsample's output
    #' burnin: length of the burn-in period
    #'   (number of sample being discarded at the beginning of the chain).
    #' thin: thinning parameter: only 1 sample out of 'thin' will be kept
    #' returns: a single numeric value

    k <- nrow(sample$p)
    nsample <- ncol(sample$p)
    inds <- (burnin+1):nsample
    inds <- inds[inds%%thin==0]
    vectllk <- vapply(inds, function(niter){
      sum(vapply(1:k, function(j) {
        sample$p[j,niter] * dmnorm(x, sample$Mu[j,,niter], sample$Sigma[,,j,niter])
      }, FUN.VALUE = numeric(1)))
    }, FUN.VALUE = numeric(1))
    
    return(mean(vectllk))
}

MHpredictiveCdf <- function(x , sample , burnin = 0, thin = 1) {
    #' posterior predictive cdf computed from MH output.
    #' arguments: see MHpredictive.
    #' returns: a single numeric value. 

    k <- nrow(sample$p)
    nsample <- ncol(sample$p)
    inds <- (burnin+1):nsample
    inds <- inds[inds%%thin==0]

    vectcdf <- vapply(inds, function(niter){
      gmcdf(x = x, Mu = sample$Mu[,,niter], Sigma = sample$Sigma[,,,niter], p = sample$p[,niter])
    }, FUN.VALUE = numeric(1))
    
    return(mean(vectcdf))
}