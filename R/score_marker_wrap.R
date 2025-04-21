reals_rev1 <- function(x, panel = c(100:400), shi = 1, ploidy = 2,
                       left.cond = c(0.4, 3), right.cond = 0.2, windowL = 0.5,
                       windowR = 0.5) {

  neg <- which(x$wei <= 0)
  if (length(neg) > 0) {
    x <- list(pos = x$pos[-neg], hei = x$hei[-neg], wei = x$wei[-neg])
  }
  picos <- numeric()
  for (d in 1:length(x$wei)) {
    respi <- which(
      ((panel - x$wei[d]) <= windowL) &
        (x$wei[d] - panel) <= windowR
    )

    if (length(respi) > 0) {
      picos[d] <- 1
    } else {
      picos[d] <- 0
    }
  }
  z1 <- which(picos == 1)
  if (length(z1) > 0) {
    x2 <- list(pos = x$pos[z1], hei = x$hei[z1], wei = x$wei[z1])
    x3 <- Fragman::separate(x2, shi, type = "bp")
    highest <- which(x3$hei == max(x3$hei))
    che <- x3$hei[1:highest[1]]
    ha <- which(che >= (x3$hei[highest] * left.cond[1]))
    cha <- x3$wei[1:highest[1]]
    ha2 <- unique(c(which(abs(cha - x3$wei[highest]) >=
      left.cond[2]), highest))
    chu <- x3$hei[highest[1]:length(x3$hei)]
    hu <- (highest[1]:length(x3$hei))[which(chu > (max(x3$hei) *
      right.cond))]
    ss1 <- intersect(ha, ha2)
    ha3 <- unique(c(ss1, hu))
    if (length(ha3) > 0) {
      x3 <- list(
        pos = x3$pos[ha3], hei = x3$hei[ha3],
        wei = x3$wei[ha3]
      )
    } else {
      x3 <- x3
    }
    z2 <- length(x3$pos)
    if (z2 == 1) {
      x4 <- list(pos = rep(x3$pos, ploidy), hei = rep(
        x3$hei,
        ploidy
      ), wei = rep(x3$wei, ploidy))
    }
    if (z2 > 1) {
      toget <- sort(x3$hei, decreasing = TRUE)[1:ploidy]
      z3 <- which(x3$hei %in% toget)
      x4 <- list(
        pos = x3$pos[z3], hei = round(x3$hei[z3]),
        wei = x3$wei[z3]
      )
    }
  } else {
    x4 <- list(
      pos = rep(0, ploidy), hei = rep(0, ploidy),
      wei = rep(0, ploidy)
    )
  }
  return(x4)
}

homo_panel_rev1 <- function(x, panel, windowL = 0.49, windowR = 0.49) {
  newxxx <- numeric()
  for (i in 1:length(x$wei)) {
    v <- which(
      ((panel - x$wei[i]) <= windowL) &
        (x$wei[i] - panel) <= windowR
    )

    if (length(v) > 0) {
      y <- panel[v]
    } else {
      y <- 0
    }
    newxxx[i] <- y
  }
  x$wei <- newxxx
  return(x)
}


#' Score Markers Wrapper
#'
#' This is a revision of the Fragman script score.markers, for the original
#' instructions and parameters, run '?score.markers'. This revision designates
#' separate parameters for Left and Right search windows.
#'
#' @param my.inds The list output from the fsa_batch_imp or storing.inds
#' function that contains the channel information from the individuals that you
#' want to score.
#' @param channel The number of the channel you wish to analyze. Typically 1 is
#' blue, 2 is green, 3 yellow, and 4 red.
#' @param n.inds (optional) A vector specifying which fsa files to score.
#' @param panel A vector containing the expected allele sizes for this marker.
#' @param shift All peaks at that distance from the tallest peak will be
#' ignored and be considered noise.
#' @param ladder A vector containing the expected peaks for your ladder.
#' @param channel.ladder The channel number where your ladder can be found.
#' @param ploidy The name is a relic of the fact that [Fragman::score.markers]
#' was originally written for plants. In the context of pooled egg samples it
#' is used to specify the number of possible alleles in the marker.
#' @param left.cond The first part is a percentile (0-1) that corresponds to
#' the height that a peak to the left of the tallest peak must be in order to
#' be considered real. The second argument is a number of base pairs that a
#' peak to the left of the tallest peak must be away to be considered as real.
#' @param right.cond A percentile (0-1) that corresponds to the height that a
#' peak to the right of the tallest peak must be in order to be real.
#' @param warn TRUE/FAlSE Do you want to receive warnings when detecting the
#' ladder?
#' @param windowL the window means that all peaks closer by that distance to
#' the left of the panel peaks will be accounted as peaks.
#' @param windowR the window means that all peaks closer by that distance to
#' the right of the panel peaks will be accounted as peaks.
#' @param init.thresh A value that sets a minimum intensity in order for a
#' peak to be called.
#' @param ladd.init.thresh We don't recommend messing with this parameter
#' unless your ladder has special circumstances. See [Fragman::score.markers]
#' @param method In cases where samples weren't sized using the
#' info.ladder.attach function, this technique steps in to identify ladder
#' peaks. You have three method options using an argument: "cor" explores all
#' potential peak combinations and thoroughly searches for correlations to
#' identify the correct peaks corresponding to expected DNA weights; "ci"
#' constructs confidence intervals to identify peaks meeting specified
#' conditions from earlier arguments; "iter2" applies an iterative strategy to
#' identify the most likely peaks aligning with your ladder expectations. The
#' default method is "iter2."
#' @param env Please do not change this parameter, it is used to detect the
#' users environment.
#' @param my.palette (optional) A character vector specifying which colors
#' to use for the output RFU plots.
#' @param plotting TRUE/FALSE Do you want to create pdf output plots?
#' @param plotdir The name of the directory where output pdf plots should
#' be stored.
#' @param pref The number of plots to be drawn in the output plot.
#'
#' @importFrom stats lm
#' @importFrom Fragman big.peaks.col
#' @importFrom Fragman separate
#' @importFrom grDevices pdf
#' @importFrom Fragman transp
#' @importFrom graphics abline
#' @importFrom grDevices dev.off
#' @importFrom utils flush.console
#' @importFrom qpdf pdf_combine
#' @importFrom Fragman find.ladder
#' @importFrom graphics axis
#' @importFrom graphics legend
#' @importFrom stats predict
#' @return The score_markers_rev3 function will return a list containing three
#'  variables: $pos, $hei, and $wei. These correspond to the index position for
#'  the intensities, the intensity of each peak, and the weight in base pairs
#'  based on the ladder respectively. If plotting = TRUE, a pdf file will
#'  also have been created in the specified directory. This pdf file allows
#'  you to visually inspect how all of the peaks were scored.
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", package = "pooledpeaks")
#' mock_fsa_batch_imp_output<- fsa_batch_imp(file_path, channels = 5,
#' fourier = FALSE, saturated = FALSE, lets.pullup = FALSE,
#' plotting = FALSE, rawPlot = FALSE)
#' panel <- c(176,179,182,185,188,191,194,197,200,203,206)
#' ladder <- c( 140, 160, 180, 200, 214, 220,240, 250, 260, 280, 300, 314)
#' mock_fsa_batch_imp_output <- associate_dyes(mock_fsa_batch_imp_output,
#'                              file_path)
#' \donttest{score_markers_rev3(my.inds = mock_fsa_batch_imp_output,
#'                             channel = 1,
#'                             channel.ladder = 5,
#'                             panel = "panel",
#'                             ladder = ladder,
#'                             init.thresh = 200,
#'                             ploidy = length(panel),
#'                             plotting = FALSE)
#'                             }

score_markers_rev3 <- function(my.inds, channel= 1, n.inds = NULL, panel= NULL,
                               shift = 0.8, ladder, channel.ladder = NULL,
                               ploidy = 2, left.cond = c(0.6, 3),
                               right.cond = 0.35, warn = FALSE, windowL = 0.5,
                               windowR = 0.5, init.thresh = 200,
                               ladd.init.thresh = 200, method = "iter2",
                               env = parent.frame(), my.palette = NULL,
                               plotting = FALSE, plotdir = "plots_scoring/",
                               pref = 3) {
  if (plotting == TRUE) {
    message(paste0("You are writing plots for all samples to the directory '",
    plotdir, "'.
  For faster calculation but no plots, set 'plotting = FALSE'\n"))
  } else {
    message(paste0(" You have set 'plotting = FALSE', meaning this command will
    return
                   only a list peak scores\n"))
  }

  if (is.character(panel)) {
    mic <- panel
    if (exists(panel, envir = env, inherits = FALSE)) {
      panel <- get(panel, envir = env)
    } else {
      stop(call. = FALSE, 'Check that your panel object is specified in quotes
      e.g. panel = "my_panel" and exists in the environment')
    }
  } else if (!is.null(panel) && !is.numeric(panel)) {
    stop(call. = FALSE, 'Panel should be either a character string representing
    the variable name or a numeric vector')
  }

  dev <- 50
  thresh <- NULL

  if (length(n.inds) > length(my.inds)) {
    message(paste(
      "You are trying to examine more individuals than you actually read?
      You selected in 'n.inds' argument",
      length(n.inds), "individuals but you only provided",
      length(my.inds), " individuals. Please select a number of individuals
      smaller or same size than the ones contained in 'my.inds' argument"
    ))
    stop
  }

  message(paste("Scoring ", mic, "peaks\nPlease be patient!\n"))

  if (method == "ci") {
    message(paste("Please make sure you have used the same 'dev' value
                  you found convenient for your ladder detection or probably
                  your call will not match"))
  }

  if (is.null(channel.ladder)) {
    channel.ladder <- dim(my.inds[[1]])[2]
  } else {
    channel.ladder <- channel.ladder
  }

  if (dim(my.inds[[1]])[2] < channel.ladder) {
    warning(paste("ERROR, you have indicated an argument channel.ladder=5,
                  but your data contains less channels/colors"))
    stop
  }

  if (is.null(n.inds)) {
    n.inds <- 1:length(my.inds)
  } else {
    n.inds <- n.inds
  }

  if (is.null(thresh)) {
    thresh <- rep(list(c(1, 1, 1, 1, 1)), length(my.inds))
  }

  tot <- length(n.inds)
  my.inds2 <- list(NA)
  thresh2 <- list(NA)
  for (i in 1:length(n.inds)) {
    v1 <- n.inds[i]
    my.inds2[[i]] <- my.inds[[v1]]
    names(my.inds2)[i] <- names(my.inds)[v1]
  }

  ncfp <- c(
    "channel_1", "channel_2", "channel_3", "channel_4", "channel_5",
    "channel_6"
  )
  if (!is.null(my.palette)) {
    cfp <- rep(my.palette, 100)
  } else {
    cfp <- c(
      "cornflowerblue", "chartreuse4", "gold2", "red",
      "orange", "purple"
    )
  }

  col.list <- list(NA)
  att1 <- numeric()
  list.data <- list(NA)

  if (exists("list.data.covarrubias")) {
    list.data <- env$list.data.covarrubias
  } else {
    list.ladders <- lapply(my.inds2, function(x) {
      y <- x[, channel.ladder]
      return(y)
    })
    list.data <- lapply(list.ladders, Fragman::find.ladder,
      ladder = ladder,
      draw = FALSE, dev = dev,
      warn = warn, method = method, init.thresh = ladd.init.thresh
    )
  }

  list.models <- lapply(list.data, function(da) {
    y <- da[[3]]
    x <- da[[1]]
    mod <- stats::lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5),
      data = da
    )
    return(mod)
  }
  )

  list.models.inv <- lapply(list.data, function(da) {
    x <- da[[3]]
    y <- da[[1]]
    mod <- stats::lm(y ~ x, data = da)
    return(mod)
  }
  )

  xx <- lapply(my.inds2, function(x, channel) {
    1:length(x[, channel])
  }, channel = channel)

  newxx <- numeric()
  newyy <- numeric()
  new.whole.data <- list(NA)

  for (h in 1:length(xx)) {
    h1 <- n.inds[h]
    newxx <- as.vector(predict(list.models[[h1]],
                               newdata = data.frame(x = xx[[h]])))
    newyy <- my.inds2[[h]][, channel]
    new.whole.data[[h]] <- list(xx = newxx, yy = newyy)
  }

  top <- max(unlist(lapply(new.whole.data, function(x) {
    max(x$yy)
  })))
  bott <- min(unlist(lapply(new.whole.data, function(x) {
    min(x$yy)
  })))

  list.weis <- list(NA)
  lower.bounds <- numeric()
  for (k in 1:length(my.inds2)) {
    newtt <- init.thresh
    lower.bounds[k] <- newtt
    plant <- Fragman::big.peaks.col(new.whole.data[[k]]$yy, newtt)
    plant$wei <- new.whole.data[[k]]$xx[plant$pos]
    plant <- Fragman::separate(plant, shift, type = "bp")
    list.weis[[k]] <- plant
  }

  list.weis <- lapply(list.weis, function(x) {
    x$wei <- round(x$wei, digits = 4)
    return(x)
  })
  names(list.weis) <- names(my.inds2)

  if (length(panel) > 0) {
    list.weis <- lapply(list.weis, reals_rev1,
      panel = panel, shi = shift, ploidy = ploidy, left.cond = left.cond,
      right.cond = right.cond, windowL = windowL, windowR = windowR
    )
    list.weis2 <- lapply(list.weis, FUN = homo_panel_rev1, panel = panel,
      windowL = windowL, windowR = windowR
    )
  } else {
    list.weis2 <- list.weis
  }

  if (plotting == TRUE) {
    layout(matrix(1:pref, pref, 1))
    if (length(panel) > 0) {
      xm <- round(min(panel, na.rm = TRUE) - 10, digits = 0)
      xl <- round(max(panel, na.rm = TRUE) + 10, digits = 0)
    } else {
      xm <- 0
      xl <- max(ladder)
    }

    # Create output directory
    if (is.null(plotdir)) {
      stop("The 'plotdir' parameter cannot be NULL.")
    }


    message(paste0("Writing plots to directory '", plotdir, "'\n"))
    #dirplot <- paste0(getwd(), "/", plotdir, "/")
    #dir.create(dirplot, showWarnings = TRUE)

    if (!grepl("^(/|~|[A-Za-z]:)", plotdir)) {
      dirplot <- file.path(getwd(), plotdir)
    } else {
      dirplot <- normalizePath(plotdir, winslash = "/", mustWork = FALSE)
    }

    if (!dir.exists(dirplot)) {
      dir.create(dirplot, recursive = TRUE)
    }


    for (g in 1:length(n.inds)) {
      hh4 <- n.inds[g]

      if (length(which(new.whole.data[[g]]$xx > xm & new.whole.data[[g]]$xx <
        xl)) > 0) {
        mylim <- max(new.whole.data[[g]]$yy[which(new.whole.data[[g]]$xx >
          xm & new.whole.data[[g]]$xx < xl)], na.rm = TRUE) +
          100
      } else {
        mylim <- 1000
      }
      if (is.infinite(mylim)) {
        mylim <- 1000
      }

      nm <- gsub("mic_", "", mic)

      pdf_file_path <- file.path(
        dirplot,
        paste0(
          "samp", names(list.models)[hh4],
          "_ch", channel,
          "_rplot.pdf"
        )
      )

      grDevices::pdf(pdf_file_path, width = 8, height = 4)


      plot(new.whole.data[[g]]$xx, new.whole.data[[g]]$yy,
        type = "l", col = cfp[channel], xaxt = "n",
        xlim = c(xm, xl), ylim = c(-200, mylim), ylab = "Intensity",
        main = paste(nm, names(list.models)[hh4]),
        sub = paste0("windowL = ", windowL, ", windowR = ", windowR,
                     ", shift = ", shift),
        xlab = "bp",
        lwd = 2, las = 2
      )
      graphics::axis(1, at = c(xm:xl), labels = xm:xl, cex.axis = 0.8)
      rect(
        xleft = (list.weis2[[g]]$wei - windowL),
        ybottom = (bott - 200), xright = (list.weis2[[g]]$wei + windowR),
        ytop = (top + 1000), col = Fragman::transp("lightpink",0.2),border = NA
      )
      graphics::abline(
        v = list.weis[[g]]$wei, lty = 3, col = "blue",
        cex = 0.5
      )
      graphics::abline(
        v = list.weis2[[g]]$wei, lty = 3, col = "red",
        cex = 0.5
      )
      graphics::abline(
        h = lower.bounds[g], lty = 2, col = "chocolate",
        cex = 0.5
      )
      graphics::legend("topright",
        legend = c("Peak found", "Panel peak","Panel window",
                   "Minimum Detected"),
        col = c("blue", "red", Fragman::transp("lightpink", 0.3), "chocolate"),
        bty = "n", lty = c(3, 3, 1, 3), lwd = c(1, 1, 3, 1), cex = 0.75
      )

      grDevices::dev.off()
      message(paste0("Saved ", g, "/", length(n.inds), " plots"))
      utils::flush.console()
    }

    message("Cleaning up...")

    pdflist <- paste0(
      dirplot,
      dir(dirplot, "^samp")
    )

    qpdf::pdf_combine(
      input = pdflist,
      output = paste0(
        dirplot,
        paste0("all_", nm, "_scores.pdf")
      )
    )

    # # cleanup files
    lapply(pdflist, file.remove)
  }

  return(list.weis2)
}
