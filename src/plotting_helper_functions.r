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


plot_curves_1d_helper <- function(X_new, pred, df_by=NULL, 
                                  include_design=!is.null(X_design), 
                                  y_new=NULL, X_design=NULL, y_design=NULL, 
                                  plot_title=NULL, xlab="x", ylab="y", 
                                  ground_truth_col="black", design_color="black", 
                                  design_pt_size=1.5, line_thickness=1.0, 
                                  legend=FALSE) {
  # Similar to `plot_pred_1d_helper` but does not plot prediction intervals; instead plots 
  # multiple curves. This is useful for comparing multiple approximations to a function. 
  # Can still plot design points, just like `plot_pred_1d_helper`. 
  # `pred` is a matrix of dimension `(nrow(X_new), K)`, where `K` is the number of curves
  # to plot. A baseline curve will also be plotted if `y_new` is provided. 
  # By default, if `df_by` is NULL then each line (column of `pred`) will be plotted a 
  # different color. The names of the lines in the legend will be taken from `colnames(pred)`
  # if is it non-NULL; otherwise default names will be assigned. Passing `df_by` allows 
  # specifying a finer mapping based on both color and linetype. `df_by` should have number
  # of rows equal to `ncol(pred)` and ordered in the same way as the columns in `pred`. 
  # Only columns in `df_by` called "color" or "linetype" are extracted, and these are mapped
  # to the ggplot attributes of the same names. 
 
  # Set default title, if not provided. 
  if(is.null(plot_title)) plot_title <- "Approximations"
  
  # Column names of `pred` are used for plot labels. 
  if(is.null(dim(pred))) pred <- matrix(pred, ncol=1)
  if(is.null(colnames(pred))) {
    ids <- paste0("y", 1:ncol(pred))
    colnames(pred) <- ids
  } else {
    ids <- colnames(pred)
  }

  df_pred <- cbind(x=drop(X_new), as.data.frame(pred))
  df_pred <- melt(df_pred, id.vars="x", variable.name="id", value.name="y")
  
  # Optionally add classes for each line. 
  if(!is.null(df_by)) {
    valid_by_cols <- is.element(colnames(df_by), c("color","linetype"))
    assert_that(any(valid_by_cols))
    df_by <- df_by[, colnames(df_by)[valid_by_cols]]
    df_by$id <- ids 
    df_pred <- merge(df_pred, df_by, by="id")
  }
  
  plt <- ggplot()
  
  # Plot lines. 
  by_color <- TRUE
  by_linetype <- FALSE
  if(!is.null(df_by)) {
    by_color <- ("color" %in% colnames(df_by))
    by_linetype <- ("linetype" %in% colnames(df_by))
  } else {
    df_pred$color <- df_pred$id
  }
  
  if(by_color && !by_linetype) {
    plt <- plt + 
            geom_line(aes(x=x, y=y, col=color), df_pred, linewidth=line_thickness) + 
            ggtitle(plot_title) + xlab(xlab) + ylab(ylab)
  } else if(by_color && by_linetype) {
    plt <- plt + 
            geom_line(aes(x=x, y=y, col=color, linetype=linetype), df_pred, linewidth=line_thickness) + 
            ggtitle(plot_title) + xlab(xlab) + ylab(ylab)
  } else {
    stop("No by cols identified.")
  }
  
  # True values at prediction locations. 
  if(!is.null(y_new)) {
    df_true <- data.frame(x=drop(X_new), y=drop(y_new))
    plt <- plt + geom_line(aes(x=x, y=y), df_true, color=ground_truth_col, 
                           linetype="dashed", linewidth=line_thickness)  
  }
  
  # Design points. 
  if(include_design) {
    df_design <- data.frame(x=drop(X_design), y=drop(y_design))
    plt <- plt + geom_point(aes(x,y), df_design, color=design_color, size=design_pt_size)
  }
  
  # Remove legend, if requested.
  if(!legend) plt <- plt + theme(legend.position="none")
  
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
  #                will be directly plotted as points on the plot and colored red. This 
  #                argument is typically used to plot design points. 
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


get_input_grid_1d_projection = function(x_names, x_vary=x_names, X_list=NULL, 
                                        X_fixed=NULL, X_bounds=NULL, 
                                        n_points_default=100L) {
  # A helper function for plotting 1d projections of a scalar-valued function
  # with multiple inputs. In particular, constructs sets of input points to the 
  # function where only one variable is varied at a time, with the remaining 
  # variables fixed at some values. This function provides various ways to 
  # specify the points over which each variable is varied, as well as the 
  # values of the points at which the fixed (non-varying) variables are fixed 
  # at. In particular, these sets of points may be passed explicitly via 
  # `X_list` and `X_fixed`, respectively. If not passed, then default values 
  # will be chosen based on `X_bounds`.
  #
  # Args:
  #    x_names: the full set of variable names, including all fixed and varied
  #             parameters. Required to be able to specify the full parameter 
  #             dimension `d`.
  #    x_vary: character, the vector of variable names that will be varied. 
  #             i.e., the axes onto which 1d projections will be considered.
  #             Defaults to all variables.
  #    X_list: list, providing an explicit way to specify the values at which 
  #            the non-fixed variables are varied. The names of this list must 
  #            be elements of `x_vary`. Each element is a numeric vector or 
  #            one-column matrix containing the input grid for the respective 
  #            variable. If a variable in `x_vary` is not present in 
  #            `names(X_list)`, then a default input grid for this variable will
  #            be constructed using `X_bounds`.
  #    X_fixed: matrix, of shape (m,d). Each row is a value at which to fix the 
  #             non-varying parameters. Column names must correspond to 
  #             `x_names`.Technically, if a variable is never varied then it 
  #             need not be in the matrix, but it is typically easiest to 
  #             include all parameters as columns.
  #    X_bounds: matrix, of shape (2,d) with the two rows specifying lower and 
  #              upper bounds for each parameter, respectively. If a varied 
  #              variable is not provided in `X_list` then it must have bounds
  #              specified, in which case its input grid is linearly spaced 
  #              between these bounds with `n_points_default` points. If 
  #              `X_fixed` is NULL, then the bounds must be provided, in which 
  #              case the fixed parameters are set to their midpoints between 
  #              the bounds.
  #    n_points_default: integer, the number of points to use in default grid 
  #                      construction, when the grid is not explicitly passed.
  # 
  # Returns:
  # list, with names set to `x_vary`. Each element contains a sublist of
  # length `nrow(X_fixed)` (length 1 if `X_fixed` is NULL). Specifically, 
  # for the returned list `l`, the element `l[["x"]][[j]]` is a matrix 
  # containing the input grid in which the variable "x" is varied, and all
  # other variables are fixed at the values specified by `X_fixed[,j]` (or 
  # their midpoints determined by `X_bounds` is `X_fixed` is NULL). This 
  # matrix will have shape (n,d), where `d = length(x_names)` is the full 
  # parameter dimension and `n` is the number of points given in 
  # `X_list[["x"]]`, if provided, and otherwise is equal to `n_points_default`.
  
  assert_that(all(x_vary %in% x_names))
  
  # Total parameter dimension, and number of parameters that will be varied.
  d <- length(x_names)
  n_vary <- length(x_vary)
  n_fixed <- d - 1L

  # Each variable being varied must either have a grid explicitly defined in 
  # `X_list` or have bounds provided in `X_bounds`.
  for(x_name in x_names) {
    assert_that(isTRUE(x_name %in% names(X_list)) || 
                isTRUE(x_name %in% colnames(X_bounds)))
  }
  
  # If `X_fixed` is not provided, then `X_bounds` must be given.
  assert_that(!(is.null(X_fixed) && is.null(X_bounds)))
  if(!is.null(X_bounds)) assert_that(setequal(x_names, colnames(X_bounds)))
  
  # Construct default input grids for varied parameters, if not provided.
  if(is.null(X_list)) X_list <- list()
  for(x_name in x_vary) {
    if(is.null(X_list[[x_name]])) {
      x_bounds <- X_bounds[,x_name]
      X_list[[x_name]] <- seq(x_bounds[1], x_bounds[2], 
                              length.out=n_points_default)
    }
  }
  
  # If not provided, set values of fixed parameters to midpoints with respect to 
  # `X_bounds`. At present, all variables need to be fixed in the same way, 
  # either using `X_fixed` or `X_bounds`.
  if(is.null(X_fixed)) {
    X_fixed <- matrix(apply(X_bounds, 2, mean), nrow=1)
    colnames(X_fixed) <- colnames(X_bounds)
  }
  
  # Create input grids.
  l <- list()
  for(x_name in x_vary) {
    l[[x_name]] <- list()
    x_names_fixed <- setdiff(x_names, x_name)
    x_grid <- matrix(X_list[[x_name]], ncol=1, dimnames=list(NULL,x_name))
    n_grid <- nrow(x_grid)
    
    for(j in 1:nrow(X_fixed)) {
      X_grid <- cbind(x_grid, matrix(X_fixed[j,x_names_fixed,drop=FALSE],  
                                     nrow=n_grid, ncol=n_fixed, byrow=TRUE, 
                                     dimnames=list(NULL,x_names_fixed)))
      l[[x_name]][[j]] <- X_grid[,x_names]
    }
  }
 
  return(l)
}


plot_true_pred_scatter <- function(y_pred, y_true, include_CI=FALSE, CI_lower=NULL,
                                   CI_upper=NULL, y_design=NULL, plot_title=NULL,
                                   xlab="observed", ylab="predicted", ...) {
  # Produces a scatter plot of predictions against ground truth values. Note that 
  # the plot does not display any information about the inputs at which the 
  # responses are being considered. 
  
  # Set default title/labels, if not provided. 
  if(is.null(plot_title)) plot_title <- "Predictions vs. Truth"
  if(is.null(xlab)) xlab <- "observed"
  if(is.null(ylab)) ylab <- "predicted"
  
  plt <- ggplot()
  df_pred <- data.frame(y_pred=y_pred, y_true=y_true)
  
  # Plot confidence intervals first, so they don't cover up other layers. 
  if(include_CI) {
    df_pred$CI_upper <- CI_upper
    df_pred$CI_lower <- CI_lower
    df_pred$captures_truth <- (df_pred$y_true >= df_pred$CI_lower) & 
                              (df_pred$y_true <= df_pred$CI_upper)
    plt <- plt + 
            geom_segment(aes(x=y_true, y=CI_lower, xend=y_true, yend=CI_upper, 
                             color=captures_truth), df_pred) + 
            scale_color_manual(values = c("FALSE"="orange", "TRUE"="gray"))
  }
  
  # Plot point predictions, and line y=x. 
  plt <- plt + 
          geom_point(aes(x=y_true, y=y_pred), df_pred, shape=1) + 
          geom_abline(slope=1, intercept=0, color="red") + 
          ggtitle(plot_title) + xlab(xlab) + ylab(ylab)
          
  return(plt)
  
}


# -----------------------------------------------------------------------------
# ggplot themes and formatting functions. 
# -----------------------------------------------------------------------------

ggtheme_journal <- function(legend_position="none", legend_title=element_blank(), title_size=35, 
                            legend_size=20, ...) {
  
  theme_journal <- theme(legend.position=legend_position, 
                         legend.title=legend_title,
                         legend.text=element_text(size=legend_size),
                         panel.grid.minor=element_blank(),  panel.grid.major=element_blank(),
                         panel.background=element_blank(), panel.border=element_blank(), 
                         axis.line.x = element_line(size=0.5, linetype="solid", colour="black"),
                         axis.line.y = element_line(size=0.5, linetype="solid", colour="black"),
                         plot.background=element_blank(), 
                         axis.title=element_text(size=22), 
                         plot.title=element_text(size=title_size), ...)
  
  return(theme_journal)
}

ggformat_journal <- function(plt, remove_title=TRUE, xlim=NULL, ylim=NULL, ...) {
  
  if(remove_title) plt <- plt + ggtitle(NULL)
  if(!is.null(xlim)) plt <- plt + xlim(xlim)
  if(!is.null(ylim)) plt <- plt + ylim(ylim)
  
  plt + ggtheme_journal(...)
}


get_common_lims <- function(...) {
  # Given ggplot objects as arguments, returns the smallest xlim and ylim 
  # that includes all of the plots xlims/ylims without cutting anything out. 
  # See https://stackoverflow.com/questions/7705345/how-can-i-extract-plot-axes-ranges-for-a-ggplot2-object
  # for code on extracting xlim/ylim from plot. 
  
  plt_list <- list(...)
  
  xlims <- lapply(plt_list, function(plt) layer_scales(plt)$x$range$range)
  ylims <- lapply(plt_list, function(plt) layer_scales(plt)$y$range$range)
  
  xmin <- min(sapply(xlims, function(x) x[1]))
  xmax <- max(sapply(xlims, function(x) x[2]))
  ymin <- min(sapply(ylims, function(y) y[1]))
  ymax <- max(sapply(ylims, function(y) y[2]))
  
  return(list(xlim=c(xmin, xmax), ylim=c(ymin, ymax)))
  
}


align_plots <- function(..., theme=NULL, theme_args=NULL) {
  # Aligns the axes of a set of ggplot plots using `get_common_lims()` 
  # and then  optionally adds a plot theme to each plot. 
  #
  # Args: 
  #    ...: ggplot plot objects. 
  #    theme: a function to add themes to a ggplot plot. Must be able to 
  #           be called as `theme(plt, ...)`, where the named argument 
  #           list `theme_args` will take the place of `...`. See 
  #           `ggformat_journal` for such a function. 
  #    theme_args: named list of arguments to pass to `theme`. This 
  #                should not include the `plt` argument to `theme()` 
  #                since this argument will be filled in by `...`. 
  #
  # Returns: 
  #    A list of the formatted plots.
  
  plot_list <- list(...)
  
  # Align axis limits. 
  lims_post_trim <- get_common_lims(...)
  
  # Add theme. 
  if(!is.null(theme)) {
    theme_call <- function(plt) do.call(theme, c(plt=plt, theme_args))
    plot_list <- lapply(plot_list, theme_call)
  }
  
  return(plot_list)
}




