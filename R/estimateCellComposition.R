#' Estimate cell compositions in test samples based on reference profile and projection.
#'
#' @param test.beta Beta matrix of samples to be tested. Probes are in rows. Samples are in columns. Row names have to be probe IDs.
#' @param ref.beta.mat Matrix of mean beta values for the reference cells. Row names have to be probe IDs. For EPIC arrays,
#' ref.projection.EPIC$ref.beta.mat, pre-calculated from GSE167998, can be used.
#' @param projection Projection matrix to convert selected probes to PCA space. Row names have to be probe IDs. Columns are PCs.For EPIC arrays,
#' ref.projection.EPIC$projection, pre-calculated from GSE167998, can be used.
#' @param n.PC number of PCs to include for deconvolution.
#' @param extended Boolean variable to indicate whether to include Bcell, CD4T, CD8T and Gran in output.
#'
#' @return A matrix of cell compositions in the test samples. Samples are in rows, cell types in columns.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' test.res=estimateCellComposition(test.beta = test.beta, ref.beta.mat = ref.projection.EPIC$ref.beta.mat,projection = ref.projection.EPIC$projection, n.PC = 20,extended = F)
#' plot_celltype(test.res*100,test.pd*100,celltype = 'Bmem')
#' }
estimateCellComposition=function(test.beta, ref.beta.mat, projection, n.PC=20,extended=TRUE){
  # find overlap between test and ref
  common.pid = intersect(rownames(ref.beta.mat),rownames(test.beta))
  print(paste0(length(common.pid),' probes found in test data.'))
  # create overlapping matrices and vectors
  projection0 = projection[match(common.pid,rownames(projection)),1:n.PC]
  ref.beta.mat = data.matrix(ref.beta.mat[match(common.pid,rownames(ref.beta.mat)),])
  test.beta = t(test.beta[match(common.pid,rownames(test.beta)),])
  projection = projection[match(common.pid,rownames(projection)),]

  # estimate cell composition for each sample
  test.res=t(apply(test.beta,1,function(test0){
    y=crossprod(projection[,1:n.PC], test0)
    X=crossprod(projection[,1:n.PC], ref.beta.mat)

    cp_X_pd = Matrix::nearPD(crossprod(X), base.matrix = TRUE)

    if (!isTRUE(cp_X_pd$converged))
      stop("Could not find nearest positive definite matrix.")
    res = quadprog::solve.QP(Dmat = cp_X_pd$mat, dvec = crossprod(y,X),
                             Amat = cbind(1, diag(ncol(X))), bvec = c(1, rep(0,ncol(X))), meq = 1)
    stats::setNames(round(res$solution, 7), colnames(ref.beta.mat))
  }))

  test.res=as.data.frame(test.res)
  if(extended){
    test.res$Bcell=test.res$Bmem+test.res$Bnv
    test.res$CD4T = test.res$CD4mem + test.res$CD4nv + test.res$Treg
    test.res$CD8T = test.res$CD8mem + test.res$CD8nv
    test.res$Gran = test.res$Bas + test.res$Eos + test.res$Neu + test.res$Mono
  }
  test.res
}
