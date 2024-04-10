#
# plotting_helper_functions.r
# General helper functions for plotting, primarily with ggplot2.  
#
# Andrew Roberts
# 

library(reshape2)

plot_Gaussian_pred_1d <- function(X_new, pred_mean, pred_var=NULL, include_design=!is.null(X_design), 
                                  include_CI=!is.null(pred_var), CI_prob=0.9, y_new=NULL,
                                  X_design=NULL, y_design=NULL, transformation=NULL, plot_title=NULL,
                                  xlab="x", ylab="y") {
  # Produces a Gaussian process prediction plot for one-dimensional input space. 
  #
  # Args:
  #    X_new: numeric, or one-column matrix, the input test locations. 
  #    pred_mean: numeric, or one-column matrix, the predictive mean at the test locations. 
  #    pred_var: numeric, or one-column matrix, the predictive variance at the test locations. 
  #    include_design: logical(1), whether or not to plot the design (i.e. training) points. 
  #    include_CI: logical(1), whether or not to plot confidence intervals. 
  #    CI_prob: numeric value in (0,1); e.g. `0.9` corresponds to 90% confidence interval. 
  #    y_new: numeric, or one-column matrix, the true response values at the prediction locations. 
  #    X_design: numeric, or one-column matrix, the design locations. 
  #    y_design: numeric, or one-column matrix, the response values at the design locations. 
  #    transformation: character, not yet implemented but intended to allow the options 
  #                    "LN", "truncated", or "rectified" as transformations of the Gaussian predictions. 
  # 
  # Returns: 
  #    ggplot2 object. 
  
  assert_that(is.numeric(X_new) || (ncol(X_new)==1), msg="plot_Gaussian_pred_1d() requires 1d input space.")
  if(!is.null(transformation)) .NotYetImplemented() 
  
  # Set default title, if not provided. 
  if(is.null(plot_title)) plot_title <- paste0("GP Predictions")
  if(include_CI) plot_title <- paste0(plot_title, ", ", 100*CI_prob, "% CI")
  if(!is.null(transformation)) plot_title <- paste0(plot_title, " ", transformation, "transform")
  
  # Confidence intervals. 
  if(include_CI) {
    CI_tail_prob <- 0.5 * (1-CI_prob)
    CI_upper <- qnorm(CI_tail_prob, pred_mean, sqrt(pred_var))
    CI_lower <- qnorm(CI_tail_prob, pred_mean, sqrt(pred_var), lower.tail=FALSE)
  } else {
    CI_upper <- NULL
    CI_lower <- NULL
  }
  
  plt <- plot_pred_1d_helper(X_new, pred_mean, include_design=include_design, include_CI=include_CI,
                             CI_lower=CI_lower, CI_upper=CI_upper, y_new=y_new, X_design=X_design,
                             y_design=y_design, plot_title=plot_title, xlab=xlab, ylab=ylab) 
  return(plt)
}


plot_pred_1d_helper <- function(X_new, pred_mean, include_design=!is.null(X_design), 
                                include_CI=!is.null(CI_lower), CI_lower=NULL, CI_upper=NULL, 
                                y_new=NULL, X_design=NULL, y_design=NULL, plot_title=NULL,
                                xlab="x", ylab="y") {
  # This is used by the `gpWrapper` and `llikEmulator` classes to produce plots 
  # summarizing the predictive distribution in the case of a 1d input space. 
  
  # Set default title, if not provided. 
  if(is.null(plot_title)) plot_title <- "Model Predictions"
  
  # Base plot: mean at prediction locations. 
  df_pred <- data.frame(x=drop(X_new), y_mean=drop(pred_mean))
  plt <- ggplot() + geom_line(aes(x, y_mean), df_pred, color="blue") + 
    ggtitle(plot_title) + xlab(xlab) + ylab(ylab)
  
  # Confidence intervals. 
  if(include_CI) {
    df_pred$CI_upper <- CI_upper
    df_pred$CI_lower <- CI_lower
    plt <- plt + geom_line(aes(x, CI_upper), df_pred, color="gray") + 
      geom_line(aes(x, CI_lower), df_pred, color="gray")
  }
  
  # True values at prediction locations. 
  if(!is.null(y_new)) {
    df_pred$y_true <- y_new
    plt <- plt + geom_line(aes(x, y_new), df_pred, color="red")
  }
  
  # Design points. 
  if(include_design) {
    df_design <- data.frame(x=drop(X_design), y=drop(y_design))
    plt <- plt + geom_point(aes(x,y), df_design, color="black")
  }
  
  return(plt)
  
}


plot_heatmap <- function(X, y, samples_kde=NULL, points_mat=NULL,  
                         raster=FALSE, point_coords=NULL, main_title="Heatmap", 
                         invert_colors=TRUE, legend_label="y",
                         log_transform=FALSE, log_func_str="log",
                         point_coords_shape=8, point_coords_col="black", 
                         points_mat_size=1, point_coords_size=3, 
                         samples_kde_lab="KDE", points_mat_lab="points_mat", 
                         KDE_opacity=1.0, xlab="x1", ylab="x2") {
  # Plots a 2d heatmap or contours of a scalar quantity `y`. Optionally overlays contours 
  # of a 2d kernel density estimate from `samples_kde`. The input locations are given by the 
  # M x 2 matrix `X`. If these input locations correspond to an evenly-spaced grid, 
  # then `raster = TRUE` may be set to create a classic heatmap produced over a dense 
  # grid. Alternatively, if `X` consists of more sparsely sampled or non-evenly-spaced 
  # locations, `raster = FALSE` should be set and the resulting plot will plot the 
  # individual points, which will still be colored in heatmap fashion. 
  #
  # Args:
  #    X: matrix, of dimension M x 2; the input locations used in the heatmap.
  #    y: numeric(M), the scalar output value used to determine the colors in the heatmap. 
  #    samples_kde: matrix, with 2 columns. These input points will not be part of the heatmap.
  #                 Instead, they will be used to construct a 2D KDE estimate and the
  #                 contours of this KDE will be overlaid on the heatmap. 
  #    points_mat: matrix, with 2 columns. These input points 
  #                    will be directly plotted as points on the plot and colored red. This 
  #                    argument is typically used to plot design points. 
  #    raster: logical(1), see above description. Set to TRUE when `X` is a dense grid of evenly-spaced  
  #            points FALSE when the points in `X` are not evenly-spaced or are sparse. 
  #    point_coords: numeric(2), coordinates to plot a single point. This typically 
  #                  corresponds to the location of some "true" value. The marker used for the 
  #                  point is given by `point_coords_shape` and `point_coords_col` which defaults 
  #                  to a black star. 
  #    main_title: character(1), the title of the plot. 
  #    invert_colors: logical(1), whether to invert the color scheme of the plot. 
  #    legend_label: character(1), the title for the legend which indicates the color scale. 
  #    log_transform: logical(1), if TRUE log transforms the scalar output data. 
  #    log_func_str: character(1), the name of the log function to use as the log transformation
  #                  if `log_transform` is TRUE. Default is "log" (natural log). Another option 
  #                  is "log10" for base 10. 
  #                  
  #
  # Returns:
  #    ggplot2 object. 
  
  assert_that(ncol(X)==2)
  
  # General settings. 
  color_direction <- ifelse(invert_colors, 1, -1)
  color_breaks <- c()
  color_values <- c()
  if(log_transform)  legend_label <- paste0(log_func_str, "(", legend_label, ")")
   
  # Store data as data.frame for ggplot2. 
  df <- as.data.frame(cbind(X, y))
  colnames(df) <- c("x1", "x2", "y")
  if(log_transform) df$y <- get(log_func_str)(df$y)
  
  #
  # Heatmap. 
  #
  # NOTE: it is essential that in neither case a "color scale" is added, since this causes difficulties when trying to 
  # manually set a color scale later. To make this work in the non-raster case, I set "`shape = 21`, which is a point 
  # that has both fill and color attributes. We can then use the fill attribute as the mapping, while removing 
  # the border of these points with `stroke = NA`. 
  if(raster) {
    plt <- ggplot(data=df) + 
            geom_tile(mapping=aes(x=x1, y=x2, fill=y)) + 
            scale_fill_viridis(discrete=FALSE, direction=color_direction) + 
            labs(fill=legend_label)
  } else {
    plt <- ggplot(data=df) + 
            geom_point(mapping=aes(x=x1, y=x2, fill=y), shape=21, stroke=NA) + # Color set using fill attribute. 
            scale_fill_viridis(discrete=FALSE, direction=color_direction) +
            labs(fill=legend_label)
  }
  
  # Title and axis labels. 
  plt <- plt + ggtitle(main_title) + xlab(xlab) + ylab(ylab)
  
  # Density contours from samples. 
  if(!is.null(samples_kde)) {
    assert_that(ncol(samples_kde)==2)
    samples_kde <- as.data.frame(samples_kde)
    colnames(samples_kde) <- c("x1", "x2")
    
    # Setting color to the character string `samples_kde_lab` allows us to manually add the 
    # color such that it will show up in the legend. 
    plt <- plt + geom_density_2d(data=samples_kde, mapping=aes(x=x1, y=x2, color=samples_kde_lab), alpha=KDE_opacity)
    color_breaks <- c(color_breaks, samples_kde_lab)
    color_values <- c(color_values, setNames("blue", samples_kde_lab))
  }
  
  # Plot points. 
  if(!is.null(points_mat)) {
    assert_that(ncol(points_mat)==2)
    points_mat <- as.data.frame(points_mat)
    colnames(points_mat) <- c("x1", "x2")
    
    plt <- plt + geom_point(data=points_mat, mapping=aes(x=x1, y=x2, color=points_mat_lab), size=points_mat_size)
    color_breaks <- c(color_breaks, points_mat_lab)
    color_values <- c(color_values, setNames("red", points_mat_lab))
  }
  
  # Mark specific point in plot. 
  if(!is.null(point_coords)) {
    plt <- plt + geom_point(data=data.frame(x1=point_coords[1], x2=point_coords[2]), 
                            aes(x=x1, y=x2), color=point_coords_col, 
                            shape=point_coords_shape, size=point_coords_size)
  }
  
  # Legend. 
  if(length(color_breaks) > 0) {
    plt <- plt + scale_colour_manual(aesthetics="colour", name="", breaks=color_breaks, values=color_values)
  }
  
  return(plt)
}

