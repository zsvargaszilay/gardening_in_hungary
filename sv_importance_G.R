### edited function from https://cran.r-project.org/web/packages/shapviz/shapviz.pdf ###

sv_importance_G <- function(object, kind = c("bar", "beeswarm", "both", "no"),
                                  max_display = 15L, fill = "#fca50a", bar_width = 2/3,
                                  bee_width = 0.4, bee_adjust = 0.5,
                                  viridis_args = getOption("shapviz.viridis_args"),
                                  color_bar_title = "Feature value",
                                  show_numbers = FALSE, format_fun = format_max,
                                  number_size = 3.2, ...) {
  stopifnot("format_fun must be a function" = is.function(format_fun))
  kind <- match.arg(kind)
  imp <- .get_imp(get_shap_values(object))
  
  if (kind == "no") {
    return(imp)
  }
  
  # Deal with too many features
  if (ncol(object) > max_display) {
    imp <- imp[seq_len(max_display)]
  }
  ord <- names(imp)
  object <- object[, ord]  # not required for kind = "bar"
  
  # ggplot will need to work with data.frame
  imp_df <- data.frame(feature = factor(ord, rev(ord)), value = imp)
  is_bar <- kind == "bar"
  if (is_bar) {
    p <- ggplot2::ggplot(imp_df, ggplot2::aes(x = value, y = feature)) +
      ggplot2::geom_bar(fill = fill, width = bar_width, stat = "identity", ...) +
      ggplot2::labs(x = "mean(|SHAP value|)", y = ggplot2::element_blank())
  } else {
    # Prepare data.frame for beeswarm plot
    S <- get_shap_values(object)
    X <- .scale_X(get_feature_values(object))
    df <- transform(
      as.data.frame.table(S, responseName = "value"),
      feature = factor(Var2, levels = rev(ord)),
      color = as.data.frame.table(X)$Freq
    )
    df$color<-as.factor(df$color)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = feature))
    if (kind == "both") {
      p <- p +
        ggplot2::geom_bar(
          data = imp_df, fill = fill, width = bar_width, stat = "identity"
        )
    }
    p <- p +
      ggplot2::geom_vline(xintercept = 0, color = "darkgray") +
      ggplot2::geom_point(
        ggplot2::aes(color = color),
        position = position_bee(width = bee_width, adjust = bee_adjust),
        ...
      ) +scale_color_brewer(palette = "Set2") +
      ggplot2::labs(
        x = "SHAP value", y = ggplot2::element_blank(), color = color_bar_title
      ) +
      ggplot2::theme(legend.box.spacing = grid::unit(0, "pt"))
  }
  if (show_numbers) {
    p <- p +
      ggplot2::geom_text(
        data = imp_df,
        ggplot2::aes(
          x = if (is_bar) value + max(value) / 60 else
            min(df$value) - diff(range(df$value)) / 20,
          label = format_fun(value)
        ),
        hjust = !is_bar,
        size = number_size
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = 0.05 + c(0.12 *!is_bar, 0.09 * is_bar))
      )
  }
  p
}


.get_imp <- function(z) {
  if (is.matrix(z)) {
    return(sort(colMeans(abs(z)), decreasing = TRUE))
  }
  # list/mshapviz
  imp <- sapply(z, function(x) colMeans(abs(x)))
  imp[order(-rowSums(imp)), ]
}

.scale_X <- function(X) {
  X_scaled <- apply(data.matrix(X), 2L, FUN = .min_max_scale)
  if (nrow(X) == 1L) t(X_scaled) else X_scaled
}

.min_max_scale <- function(z, na.rm = TRUE) {
  r <- range(z, na.rm = na.rm)
  d <- diff(r)
  if (is.na(d) || d == 0) {
    z[!is.na(z)] <- 0.5
    return(z)
  }
  (z - r[1L]) /(r[2L] - r[1L])
}

# Beeswarm position
position_bee <- function(width = NULL, adjust = NULL) {
  ggplot2::ggproto(NULL, PositionBee, width = width, adjust = adjust)
}

PositionBee <- ggplot2::ggproto(
  "PositionBee",
  ggplot2::Position,
  required_aes = c("x", "y"),
  
  setup_params = function(self, data) {
    list(
      width = if (!is.null(self$width)) self$width else
        ggplot2::resolution(data$y, zero = FALSE) * 0.4,
      adjust = if (!is.null(self$adjust)) self$adjust else 0.5
    )
  },
  
  compute_panel = function(self, data, params, scales) {
    data <- ggplot2::flip_data(data, params$flipped_aes)
    y_jit <- ave2(data$x, g = data$y, FUN = shifter, adjust = params$adjust)
    data <- ggplot2::transform_position(
      data, trans_y = function(y) y + y_jit * params$width
    )
    ggplot2::flip_data(data, params$flipped_aes)
  }
)

#===========================================================================
# Helper functions to produce the beeswarm plots
#===========================================================================

# Example
# ggplot(iris, aes(Species, Sepal.Width)) +
#   geom_point(position = position_bee(), aes(color = Species))

# Beeswarm position
position_bee <- function(width = NULL, adjust = NULL) {
  ggplot2::ggproto(NULL, PositionBee, width = width, adjust = adjust)
}

PositionBee <- ggplot2::ggproto(
  "PositionBee",
  ggplot2::Position,
  required_aes = c("x", "y"),
  
  setup_params = function(self, data) {
    list(
      width = if (!is.null(self$width)) self$width else
        ggplot2::resolution(data$y, zero = FALSE) * 0.4,
      adjust = if (!is.null(self$adjust)) self$adjust else 0.5
    )
  },
  
  compute_panel = function(self, data, params, scales) {
    data <- ggplot2::flip_data(data, params$flipped_aes)
    y_jit <- ave2(data$x, g = data$y, FUN = shifter, adjust = params$adjust)
    data <- ggplot2::transform_position(
      data, trans_y = function(y) y + y_jit * params$width
    )
    ggplot2::flip_data(data, params$flipped_aes)
  }
)

# Shift values according to their density in the unit interval by quasi-random numbers
shifter <- function(y, ...) {
  if (length(y) == 1L) {
    return(0)
  }
  dens <- stats::density(y, ...)
  dens_y <- dens[["y"]] / max(dens[["y"]])
  shift <- halton_sequence(length(y))[rank(y, ties.method = "first")] - 0.5
  2 * shift * stats::approx(dens[["x"]], dens_y, xout = y)[["y"]]
}

# "stats::ave" for grouping variable "g" and additional arguments ...
ave2 <- function(x, g = NULL, FUN = mean, ...) {
  if (is.null(g)) {
    x[] <- FUN(x, ...)
  } else {
    split(x, g) <- lapply(split(x, g), FUN, ...)
  }
  x
}

# First n values of the 1-dimensional Halton sequence (= van der Corput sequence)
# https://en.wikipedia.org/wiki/Halton_sequence
halton_sequence <- function(n, b = 2) {
  vapply(seq_len(n), halton, FUN.VALUE = 0.0)
}

# i-th element of above sequence
halton <- function(i, b = 2) {
  f <- 1
  r <- 0
  while (i > 0) {
    f <- f / b
    r <- r + f * (i %% b)
    i <- trunc(i / b)
  }
  r
}