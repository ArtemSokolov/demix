## demix.R - solves for sample heterogeneity
##   Function comments are in roxygen2 format to facilitate easy conversion to .Rd documentation
##
## by Artem Sokolov

library( quadprog )

#' Constrained non-negative matrix factorization
#'
#' Learns a heterogeneous sample mixture model of the form x = Bw + b,
#'  where x is a p-by-1 column vector representing the sample,
#'  B is a p-by-m matrix of m basis vectors, w is a m-by-1 vector of
#'  mixture coefficients, and b is a bias term. The model is subject to
#'  several constraints, which are outlined below.
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
    f <- (1-wm) / (m-1)		## Start with equal contribution fractions
    w <- c( rep( f, m-1 ), wm )
    b <- 0

    ## Precompute constant factors
    BB <- t( Bm1 ) %*% Bm1
    Amat <- cbind( rep(1,m-1), diag(m-1) )
    bvec <- c( 1 - w[m], rep(0,m-1) )
    
    ## Initial variable setup
    u <- Bm1 %*% w[-m]
    xb <- x - b

    for( iter in 1:nIter )
    {
        ## Estimate the missing basis
        Bm <- xb - u
        j <- Bm < 0
        Bm[j] <- 0

        ## Re-estimate the mixture coefficients
        y <- xb - Bm
        res <- solve.QP( Dmat=BB, dvec=t(y) %*% Bm1, Amat=Amat, bvec=bvec, meq=1 )
        w[-m] <- res$solution
        u <- Bm1 %*% w[-m]

        ## Estimate the bias term
        if( fix.bias == FALSE )
            b <- mean( x - Bm - u )
        xb <- x - b

        ## Display progress
        if( iter %% as.integer(nIter / 5) == 0 )
        { r <- x - Bm - u - b
            cat( "RMSE after iteration", iter, ":", sqrt(mean( r * r )), "\n" ) }
    }

    ## Compose the final basis matrix
    B <- cbind( Bm1, Bm / w[m] )
    
    list( B=B*p, w=w, b=b*p )
}

