## demix.R - solves for sample heterogeneity
##   Function comments are in roxygen2 format to facilitate easy conversion to .Rd documentation
##
## by Artem Sokolov

#' Constrained non-negative matrix factorization
#'
#' Learns a heterogeneous sample mixture model of the form x = Bw + b,
#'  where x is a p-by-1 column vector representing the sample,
#'  B is a p-by-m matrix of m basis vectors, w is a m-by-1 vector of
#'  mixture coefficients, and b is a bias term. The model is subject to
##  several constraints, which are outlined below.
#'
#' The model is subject to non-negativity contraints and
#'  the \sum_j w_j = 1 convexity contraint.
#'
#' Furthermore, we assume that only the first m-1 bases are specified and
#'  the model is to learn the m^th basis. However, it is also assumed that
#'  w_m is known, and we must learn w_i for i = 1, ..., (m-1). The bias term
#'  b is also learned.
#'
#' @param x A p-by-1 vector of the sample in p dimensions
#' @param Bm1 A p-by-(m-1) matrix of the first (m-1) basis vectors
#' @param wm The m^th entry of the mixture coeffcient vector
#' @param nIter Number of training iterations
#' @param fix.bias If TRUE, the bias term is fixed at 0
#' @return A list with two elements:
#' \describe{
#'   \item{B}{the fully specified p-by-m matrix of basis vectors}
#'   \item{w}{the fully specified m-by-1 vector of mixture coefficients}
#'   \item{b}{the bias term in the model}
#' }
consNMF <- function( x, Bm1, wm, nIter=1000, fix.bias=FALSE )
{
    ## Verify dimensionality
    p <- length(x)
    m <- ncol(Bm1) + 1
    stopifnot( nrow(Bm1) == p )
    stopifnot( (wm > 0) & (wm < 1) )

    ## Scale the data for numerical stability
    x <- x / p
    Bm1 <- Bm1 / p
    
    ## Initial estimates
    B <- cbind( Bm1, 0 )
    f <- (1-wm) / (m-1)		## Start with equal contribution fractions
    w <- c( rep( f, m-1 ), wm )
    b <- 0

    ## Define the objective function
    f.obj <- function()
    { r <- x - B %*% w - b; sqrt(mean( r * r )) }
    
    cat( "Initial RMSE :", f.obj(), "\n" )

    for( iter in 1:nIter )
    {
        ## Estimate the missing basis
        u <- B[,-m] %*% w[-m]
        B[,m] <- (x-b-u) / w[m]
        j <- which( B[,m] < 0 )
        B[j,m] <- 0

#        cat( "RMSE after estimating the missing basis:", f.obj(), "\n" )

        ## Re-estimate the mixture coefficients
        y <- x - B[,m] * w[m] - b
        w[-m] <- convreg( B[,-m], y, 1 - w[m] )

#        cat( "RMSE after estimating the mixture coefficients:", f.obj(), "\n" )

        ## Estimate the bias term
        if( fix.bias == FALSE )
            b <- mean( x - B %*% w )
#        cat( "RMSE after estimating the bias term:", f.obj(), "\n" )

        if( iter %% as.integer(nIter / 5) == 0 )
            cat( "RMSE after iteration", iter, ":", f.obj(), "\n" )
    }
    
    list( B=B*p, w=w, b=b*p )
}

#' Linear regression with non-negative and convex constraints.
#'
#' Trains a linear model subject to the following constraints:
#'   min (y - Xw)^T (y - Xw) s.t. w^T 1 = d, and w >= 0
#'
#' @param X n-by-p matrix of n samples in p-dimensional space
#' @param y n-by-1 vector of labels
#' @param d scalar that the sum of model weights has to be equal to
#' @return A p-by-1 vector of the learned linear model weights
convreg <- function( X, y, d=1 )
{
    library( quadprog )

    ## Pass the call to a quadratic programming solver
    p <- ncol(X)
    Amat <- cbind( rep(1,p), diag(p) )
    bvec <- c( d, rep(0,p) )
    res <- solve.QP( Dmat=t(X) %*% X, dvec=t(X)%*%y, Amat=Amat, bvec=bvec, meq=1 )
    res$solution
}
