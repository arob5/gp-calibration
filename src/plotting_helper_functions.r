#
# plotting_helper_functions.r
# General helper functions for plotting, primarily with ggplot2.  
#
# Andrew Roberts
# 

library(reshape2)

plot_Gaussian_pred_1d <- function(X_new, pred_mean, pred_var=NULL, include_design=!is.null(X_design), 
                                  include_interval=TRUE, interval_method="pm_std_dev",
                                  N_std_dev=1, CI_prob=0.9, y_new=NULL, X_design=NULL, y_design=NULL,
                                  transformation=NULL, plot_title=NULL, xlab="x", ylab="y",
                                  ground_truth_col="black", design_color="black", ...) {
  # Produces a Gaussian process prediction plot for one-dimensional input space. 
  #
  # Args:
  #    X_new: numeric, or one-column matrix, the input test locations. 
  #    pred_mean: numeric, or one-column matrix, the predictive mean at the test locations. 
  #    pred_var: numeric, or one-column matrix, the predictive variance at the test locations. 
  #    include_design: logical(1), whether or not to plot the design (i.e. training) points. 
  #    include_interval: logical(1), whether or not to plot confidence intervals. 
  #    interval_method: character(1), either "CI" or "pm_std_dev". If "CI", computes 
  #                     100*`CI_prob`% confidence interval. If "pm_std_dev" (pm="plus-minus"), 
  #                     the interval is defined by adding/subtracting `N_std_dev` standard 
  #                     deviations. 
  #    CI_prob: numeric value in (0,1); e.g. `0.9` corresponds to 90% confidence interval. Only used 
  #             if `include_interval` is TRUE and `interval_method` is "CI". 
  #    N_std_dev: integer, the number of standard deviations to add/subtract from the mean to 
  #               define the interval. Only used if `include_interval` is TRUE and `interval_method` 
  #               is "pm_std_dev". 
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
  if(is.null(plot_title)) {
    plot_title <- paste0("GP Predictions")
    if(include_interval && (interval_method == "CI")) plot_title <- paste0(plot_title, ", ", 100*CI_prob, "% CI")
    if(include_interval && (interval_method == "pm_std_dev")) plot_title <- paste0(plot_title, ", +/- ", N_std_dev, " std dev")
    if(!is.null(transformation)) plot_title <- paste0(plot_title, " ", transformation, "transform")
  }

  # Confidence intervals. 
  if(include_interval && (interval_method == "CI")) {
    CI_tail_prob <- 0.5 * (1-CI_prob)
    CI_upper <- qnorm(CI_tail_prob, pred_mean, sqrt(pred_var))
    CI_lower <- qnorm(CI_tail_prob, pred_mean, sqrt(pred_var), lower.tail=FALSE)
  } else if(include_interval && (interval_method == "pm_std_dev")) {
    CI_upper <- pred_mean + N_std_dev * sqrt(pred_var)
    CI_lower <- pred_mean - N_std_dev * sqrt(pred_var)
  } else {
    CI_upper <- NULL
    CI_lower <- NULL
  }
  
  plt <- plot_pred_1d_helper(X_new, pred_mean, include_design=include_design, include_CI=include_interval,
                             CI_lower=CI_lower, CI_upper=CI_upper, y_new=y_new, X_design=X_design,
                             y_design=y_design, plot_title=plot_title, xlab=xlab, ylab=ylab, 
                             ground_truth_col=ground_truth_col, design_color=design_color, ...) 
  return(plt)
}


plot_pred_1d_helper <- function(X_new, pred_mean, include_design=!is.null(X_design), 
                                include_CI=!is.null(CI_lower), CI_lower=NULL, CI_upper=NULL, 
                                y_new=NULL, X_design=NULL, y_design=NULL, plot_title=NULL,
                                xlab="x", ylab="y", ground_truth_col="black", 
                                design_color="black", design_pt_size=1.5, line_thickness=1.0) {
  # This is used by the `gpWrapper` and `llikEmulator` classes to produce plots 
  # summarizing the predictive distribution in the case of a 1d input space. It provides 
  # a generic interface to plot the fit to a curve over a one-dimensional input space 
  # that includes: 1.) predictive mean (point estimate); 2.) error bars; 
  # 3.) design points; 4.) ground truth curve. 
  
  # Set default title, if not provided. 
  if(is.null(plot_title)) plot_title <- "Model Predictions"
  
  plt <- ggplot()
  
  # Predictive mean is always plotted. 
  df_pred <- data.frame(x=drop(X_new), y_mean=drop(pred_mean))
  
  # Plot confidence intervals first, so they don't cover up other layers. 
  if(include_CI) {
    df_pred$CI_upper <- CI_upper
    df_pred$CI_lower <- CI_lower
    plt <- plt + geom_ribbon(aes(x=x, ymin=CI_lower, ymax=CI_upper), df_pred, fill="gray")
  }
  
  # True values at prediction locations. 
  if(!is.null(y_new)) {
    df_pred$y_true <- y_new
    plt <- plt + geom_line(aes(x=x, y=y_true), df_pred, color=ground_truth_col, 
                           linetype="dashed", linewidth=line_thickness)
  }
  
  # Base plot: mean at prediction locations.
  plt <- plt + geom_line(aes(x, y_mean), df_pred, color="blue", linewidth=line_thickness) + 
               ggtitle(plot_title) + xlab(xlab) + ylab(ylab)
  
  # Design points. 
  if(include_design) {
    df_design <- data.frame(x=drop(X_design), y=drop(y_design))
    plt <- plt + geom_point(aes(x,y), df_design, color=design_color, size=design_pt_size)
  }
  
  return(plt)
  
}


plot_curves_1d_helper <- function(X_new, pred, include_design=!is.null(X_design), 
                                  y_new=NULL, X_design=NULL, y_design=NULL, plot_title=NULL,
                                  xlab="x", ylab="y", ground_truth_col="black", design_color="black", 
                                  design_pt_size=1.5, line_thickness=1.0) {
  # Similar to `plot_pred_1d_helper` but does not plot prediction intervals; instead plots 
  # multiple curves. This is useful for comparing multiple approximations to a function. 
  # Can still plot design points, just like `plot_pred_1d_helper`. 
  # `pred` is a matrix of dimension `(nrow(X_new), K)`, where `K` is the number of curves
  # to plot. A baseline curve will also be plotted if `y_new` is provided. 
  
  # Set default title, if not provided. 
  if(is.null(plot_title)) plot_title <- "Approximations"
  
  # Column names of `pred` are used for plot labels. 
  if(is.null(dim(pred))) pred <- matrix(pred, ncol=1)
  if(is.null(colnames(pred))) colnames(pred) <- paste0("y", 1:ncol(pred))
  
  df_pred <- cbind(x=drop(X_new), as.data.frame(pred))
  df_pred <- melt(df_pred, id.vars="x", variable.name="approx", value.name="y")
  
  plt <- ggplot()
  
  # True values at prediction locations. 
  if(!is.null(y_new)) {
    df_true <- data.frame(x=drop(X_new), y=drop(y_new))
    plt <- plt + geom_line(aes(x=x, y=y), df_true, color=ground_truth_col, 
                           linetype="dashed", linewidth=line_thickness)  
  }
  
  plt <- plt + 
          geom_line(aes(x=x, y=y, col=approx), df_pred, linewidth=line_thickness) + 
          ggtitle(plot_title) + xlab(xlab) + ylab(ylab)
  
  # Design points. 
  if(include_design) {
    df_design <- data.frame(x=drop(X_design), y=drop(y_design))
    plt <- plt + geom_point(aes(x,y), df_design, color=design_color, size=design_pt_size)
  }
  
  return(plt)
  
}


plot_heatmap <- function(X, y, samples_kde=NULL, points_mat=NULL,  
                         raster=FALSE, point_coords=NULL, main_title="Heatmap", 
                         invert_colors=TRUE, legend_label="y",
                         log_transform=FALSE, log_func_str="log",
                         point_coords_shape=8, point_coords_col="black", 
                         points_mat_size=1, point_coords_size=3, 
                         samples_kde_lab="KDE", log_transform_kde=log_transform, 
                         points_mat_lab="points_mat", KDE_opacity=1.0, xlab="x1", ylab="x2") {
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
  #    log_transform_kde: logical(1), if TRUE computes the 2d kernel density estimate of 
  #                       `samples_kde` then takes the log of the esitmated density and plots
  #                       the contours of the result. If FALSE, then no log is taken. Defaults 
  #                       to the value `log_transform`, but it may be necessary to have these 
  #                       values differ in some cases; e.g. if `y` is already log-transformed 
  #                       before calling this function. 
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
  
  # Density contours from samples. I'd like to use stat_density_2d or geom_density_2d here 
  # but I cannot figure out how to plot the contours of the log of the estimated density 
  # in the case that `log_transform` is TRUE. Therefore, I compute the density 
  # separately using MASS::kde2d (which is what ggplot2 uses under the hood) and then 
  # use this to manually plot the contour function. It may also be nice to incorporate 
  # the package "ggdensity" for more interpretable plots, but again I don't know 
  # if this would be able to log-transform the plots as needed here. 
  if(!is.null(samples_kde)) {
    assert_that(ncol(samples_kde)==2)
    samples_kde <- as.data.frame(samples_kde)
    colnames(samples_kde) <- c("x1", "x2")
    
    # Compute 2d kernel density estimate. 
    kde <- MASS::kde2d(samples_kde$x1, samples_kde$x2, n=100)
    df_kde <- expand.grid(kde$x1, kde$x2)
    colnames(df_kde) <- c("x1", "x2")
    df_kde$z <- reshape2::melt(kde$z)$value
    if(log_transform_kde) {
      df_kde$z <- get(log_func_str)(df_kde$z)
      samples_kde_lab <- paste0(log_func_str, "(", samples_kde_lab, ")")
    }

    # Plot contours of the KDE. Manually setting the color value so that it shows up 
    # in the legend. 
    plt <- plt + geom_contour(aes(x=x1, y=x2, z=z), df_kde, alpha=KDE_opacity)
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


# -----------------------------------------------------------------------------
# ggplot themes and formatting functions. 
# -----------------------------------------------------------------------------

ggtheme_journal <- function(legend_position="none", title_size=35) {
  
  theme_journal <- theme(legend.position=legend_position, 
                         panel.grid.minor=element_blank(),  panel.grid.major=element_blank(),
                         panel.background=element_blank(), panel.border=element_blank(), 
                         axis.line.x = element_line(size=0.5, linetype="solid", colour="black"),
                         axis.line.y = element_line(size=0.5, linetype="solid", colour="black"),
                         plot.background=element_blank(), 
                         axis.title=element_text(size=22), 
                         plot.title=element_text(size=title_size))
  
  return(theme_journal)
}

ggformat_journal <- function(plt, remove_title=TRUE, xlim=NULL, ylim=NULL, ...) {
  
  if(remove_title) plt <- plt + ggtitle(NULL)
  if(!is.null(xlim)) plt <- plt + xlim(xlim)
  if(!is.null(ylim)) plt <- plt + ylim(ylim)
  
  plt + ggtheme_journal(...)
}










