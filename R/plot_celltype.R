#' Plot predicted against actual cell compositions for specific cell type.
#'
#' @param predicted A data.frame of predicted cell compositions. Samples are in rows. Cell types are in columns.
#' @param actual A data.frame of actual cell compositions. Samples are in rows. Cell types are in columns.
#' @param celltype A string variable containing a single cell type to plot. Must be present as a column name in both predicted and actual data.frames.
#' @param plot A boolean variable to indicate whether to plot.
#'
#' @return A list with 2 elements:
#' * A data.frame of 3 columns: Predicted, Actual, and celltype.
#' * A ggplot of the scatter plot of predicted against actual values for the cell type.
#' @export
#'
#' @examples
#' \dontrun{
#' test.res=estimateCellComposition(test.beta = test.beta, ref.beta.mat = ref.projection.EPIC$ref.beta.mat,projection = ref.projection.EPIC$projection, n.PC = 20,extended = F)
#' plot_celltype(test.res*100,test.pd*100,celltype = 'Bmem')
#' }
plot_celltype <- function(predicted, actual, celltype, plot=T) {

  # Extract the specific column for the celltype
  predicted_values <- predicted[[celltype]]
  actual_values <- actual[[celltype]]

  # Calculate R-square
  cor_val <- stats::cor(predicted_values, actual_values)
  r_square <- cor_val^2

  # Calculate RMSE
  rmse <- sqrt(mean((predicted_values - actual_values)^2))

  # Create a data frame for plotting
  plot_data <- data.frame(Predicted = predicted_values, Actual = actual_values, celltype=celltype)

  # Plot using ggplot2
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Actual, y = Predicted)) +
    ggplot2::expand_limits(x=0,y=0) +
    ggplot2::geom_point(ggplot2::aes(color = Actual)) +
    ggplot2::geom_smooth(method = 'lm', col = 'red',fullrange=TRUE, se=FALSE) +
    ggplot2::geom_abline(slope=1,intercept = 0,col='red',linetype='dashed') +
    ggplot2::ggtitle(bquote(.(celltype) ~ " " ~ R^2 ~ ":" ~ .(round(r_square, 2)) ~ " " ~ "RMSE:" ~ .(round(rmse, 2)))) +
    ggplot2::theme(legend.position = "none")

  if(plot)
    print(p)
  list(plot_data,p)
}
