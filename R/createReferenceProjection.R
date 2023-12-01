#' Create projection matrix and reference profiles from a beta matrix of reference samples.
#'
#' @param beta.ref A matrix of beta values. Probes are in rows, samples in columns. Row names have to be probe IDs.
#' @param cell.type A vector of cell names. Should have the same length as columns of beta.ref.
#' @param n.top.probe Number of probes for each cell type to be included for the projection.
#' @param top.probes A vector of probe IDs. If provided, no probe selection will be performed and the projection will be created based on these probes.
#'
#' @return A list with 3 elements:
#' * "projection" is the projection matrix to convert selected probes to PCA space.Row names have to be probe IDs. Columns are PCs;
#' * "ref.beta.mat" is the matrix of mean beta values for the reference cells. Rows are probes, columns samples. Row names have to be probe IDs;
#' * "pca_res" is the PCA result from prcomp.
#' @export
#'
#' @examples
#' \dontrun{
#'ref.projection.EPIC = createReferenceProjection(beta.ref=ref.beta, cell.type = ref.cell,n.top.probe = 200)
#' }
createReferenceProjection=function(beta.ref, cell.type, n.top.probe=150, top.probes=NULL){
  stopifnot(is.matrix(beta.ref),ncol(beta.ref)==length(cell.type))

  # perform probe selection
  if(is.null(top.probes)){
    cells = unique(cell.type)
    # perform row t-test
    probe.stat.list = lapply(cells,function(x){
      print(paste0('Probe selection for ',x))
      tstat=genefilter::rowttests(beta.ref,as.factor(cell.type==x), tstatOnly = TRUE)
      y=tstat$statistic
      names(y)=rownames(tstat)
      y
    })

    # select top probes for PCA based on t-stat
    top.probe.list=lapply(probe.stat.list,function(x){
      order(x)[c(1:n.top.probe/2, (length(x)-n.top.probe/2+1):length(x))]
    })

    # create ref eset
    all.top.probes = unique(unlist(top.probe.list))
  }

  else
    all.top.probes = top.probes
  beta.ref = beta.ref[all.top.probes,]

  # PCA of ref eset
  pca_res = stats::prcomp(t(beta.ref), scale. = F, center = F)

  projection.id=rownames(pca_res$rotation)
  projection = pca_res$rotation
  # create ref.betas as rowMeans
  ref.betas = tapply(1:ncol(beta.ref),cell.type,function(i){
    rowMeans(beta.ref[,i])
  })
  ref.beta.mat=as.data.frame(do.call(cbind,ref.betas))
  colnames(ref.beta.mat)=names(ref.betas)

  list(projection=projection, ref.beta.mat=ref.beta.mat,pca.res=pca_res)
}
