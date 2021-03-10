#' Area Biplot
#'
#' @description
#' Consider an (n x m) centered data matrix \eqn{X} and let \eqn{rank(X) = r}.
#' Alternatively to the ordinary NIPALS decomposition of \eqn{X}, where \eqn{X = T P'},
#' this package uses the resulting matrices from the extended version of the NIPALS
#' decomposition (\eqn{X = G H P'}) to determine \eqn{n} triangles whose areas are
#' used to visually estimate the \eqn{n} elements of a specific column of \eqn{X}
#' (a variable of interest). After a 90-degree rotation of the sample points, the
#' triangles are drawn regarding the following points:
#' 1. the origin of the axes.
#' 2. the sample points.
#' 3. the vector endpoint representing the selected variable.
#' @description
#' Just keep in mind that The extended NIPALS decomposition, \eqn{X = G H P'}, is
#' equivalent to the SVD decomposition, \eqn{X = U D V'}, being that:
#' 1. \eqn{G} is the matrix containing in its columns the normalized score vectors
#' of \eqn{X}, i.e., the normalized columns of \eqn{T}. If \eqn{t} is the i-th score
#' vector of the matrix \eqn{T}, then the i-th column of  \eqn{G} will be
#' \eqn{g = t / || t || }, which will correspond to the i-th left singular vector \eqn{u}.
#' 2. If \eqn{t} is the i-th column of \eqn{T}, then \eqn{|| t || = \sqrt(t' t)} gives
#' the i-th singular value of \eqn{X}. In addition, \eqn{H} is the diagonal matrix
#' containing these singular values in decreasing order, i.e., \eqn{H = D}.
#' 3. \eqn{P} is the loadings matrix, which is equivalent to the \eqn{V} matrix that
#' contains the right singular vectors of \eqn{X}.
#'
#'
#' @param L               A (n x 2) matrix containing normalized score vectors \eqn{g} (or left singular
#'                        vectors \eqn{u}).
#' @param S               An appropriate (2 x 2) diagonal matrix containing the corresponding singular
#'                        values in decreasing order.
#' @param R               A (m x 2) matrix containing the corresponding loading vectors (or right singular
#'                        vectors).
#' @param ord.row         The row  of \eqn{R} used as the base of the triangle, e.g.,
#'                        if 1 is provided, then the first row of \eqn{R} will be taken.
#' @param mode            a string providing the way the singular values will be allocated. The default
#'                        is "SS", i.e., the similar spread proposed by Gower et al.. Alternatively,
#'                        one can choose the "HJ" method (see more in Details).
#' @param tri.rgb         The hexadecimal color and alpha transparency code for the triangle. The
#'                        default is #19FF811A (green and 90% of transparency).
#' @param bg.col          A string providing the color of the background. The default is #001F3D (blue).
#' @param plot.title      A string providing the main title. The default is NONE.
#' @param plot.title.col  A string specifying the color of the main title text. The default is "FFFFFF"
#'                        (white).
#' @param plot.title.font An integer providing the style of the main title text. The default is 1 (normal
#'                        text).
#' @param plot.title.cex  A number indicating the amount by which the main title text should be scaled
#'                        relative to the default. 1 = default, 1.5 is 50% larger, and so on.
#' @param plot.sub        A string providing a sub-title. The default is NONE.
#' @param plot.sub.col    A string specifying the color of the sub-title text. The default is "FFFFFF"
#'                        (white).
#' @param plot.sub.font   An integer providing the style of the main title text. The default is 1 (normal
#'                        text).
#' @param plot.sub.cex    A number indicating the amount by which the sun-title text should be scaled
#'                        relative to the default. 1 = default, 1.5 is 50% larger, and so on.
#' @param plot.cex        A number indicating the expansion or contraction factor used to specify the
#'                        point size. The default is 0.6 (40% smaller).
#' @param plot.col        A string specifying the color of the points. The default is "FFFFFF" (white).
#' @param plot.pch        An integer specifying the shape of the points. The default is 21 (circle)
#' @param plot.xlab       A string specifying a label to the horizontal axis. The default is NONE.
#' @param plot.ylab       A string specifying a label to the vertical axis. The default is NONE.
#' @param plot.xlim       The limits for the x axis.
#' @param plot.ylim       The limits for the y axis.
#' @param axis.col        A string specifying the color of the axis. The default is #FFFFFF (white).
#' @param axis.cex        A number indicating the expansion or contraction factor used to specify
#'                        the tick label. The default is 0.7
#' @param axis.font       An integer providing the style of the tick label. The default is 1 (normal
#'                        text).
#' @param axis.asp        A number specifying the aspect ratio of the axes. The default is 1.
#' @param points.lab      A vector of characters containing the names of the data matrix rows.
#' @param var.lab         A string providing the variable name used as triangle base.
#' @param text.col.var    A string specifying the color of the variable label text. The default is
#'                        "FFFFFF" (white).
#' @param text.cex        A number indicating the expansion or contraction factor used to specify
#'                        the point labels. The default is 0.5.
#' @param text.font       An integer providing the style of the point labels. The default is 2 (bold).
#' @param text.col        A string specifying the color of the point labels text. The default is
#'                        "FFFFFF" (white).
#' @param text.pos        An integer providing the position of the point labels. The default is 3
#'                        (above).
#' @param arrow.lwd       A number specifying the line width of the arrow. The default is 1.
#' @param arrow.len       The length of the edges of the arrow head (in inches). The default is 0.1.
#' @param arrow.col       A string specifying the color of the arrow. The default is "FFFFFF"
#'                        (white).
#'
#'
#' @return An area biplot is produced on the current graphics device.
#'
#'
#' @author
#' Alberto Silva <albertos@ua.pt>, Adelaide Freitas <adelaide@ua.pt>
#'
#'
#' @references
#' J.C. Gower, P.J.F. Groenen, M. van de Velden (2010). Area Biplots. Journal of Computational and
#' Graphical Statistics, v.19 (1), pp. 46-61. \doi{10.1198/jcgs.2010.07134}
#'
#'
#' @examples
#' library(nipals)
#' data(uscrime)
#' Y = uscrime[, -1]
#'
#' # first case: scale is false
#' nip = nipals(Y, ncomp = 2, center = TRUE, scale = FALSE, force.na = TRUE)
#' L = nip$scores
#' R = nip$loadings
#' S = diag(nip$eig[1:2])
#' areabiplot(L, S, R, 5, points.lab = c(uscrime[, 1]),var.lab= "burglary")
#'
#' # second case: scale is true
#' nip = nipals(Y, ncomp = 2, center = TRUE, scale = TRUE, force.na = TRUE)
#' L = nip$scores
#' R = nip$loadings
#' S = diag(nip$eig[1:2])
#' areabiplot(L, S, R, 4, points.lab = c(uscrime[, 1]),var.lab= "assault")
#'
#'
#' @details
#' 1. If the variables (the columns of X) are measured in different units or
#' their variability differs considerably, one could perform a variance scaling
#' to get better visual results on the graph (see Examples). In this case, the
#' percentage of variance explained by the first principal components might decrease.
#' 2. The "HJ" mode is reserved for an application under implementation.
#'
#'
#' @export
#' @importFrom            grDevices rgb
#' @importFrom            graphics arrows par polygon axis text
#' @import                nipals
#'
#'
#-------------------------------------------------------------------------------------------------------

## area biplot function
    areabiplot <-
      function(L, S, R, ord.row, mode = NULL, tri.rgb = NULL, bg.col = NULL, plot.title = NULL,
            plot.title.col = NULL, plot.title.font = NULL, plot.title.cex = NULL, plot.sub = NULL,
            plot.sub.col = NULL, plot.sub.font = NULL, plot.sub.cex = NULL, plot.cex = NULL,
            plot.col = NULL, plot.pch = NULL, plot.xlab = NULL, plot.ylab = NULL, plot.xlim = NULL,
            plot.ylim = NULL, points.lab = NULL, var.lab = NULL, text.col.var = NULL, text.cex = NULL,
            text.font = NULL, text.col = NULL, text.pos = NULL, axis.col = NULL, axis.cex = NULL,
            axis.font = NULL, axis.asp = NULL, arrow.lwd = NULL, arrow.len = NULL, arrow.col = NULL)
        {

          ## dimensions
            n <- nrow(L)
            m <- nrow(R)

          ## scaling mode
            if (is.null( mode )) mode <- "SS"

            if ( mode == "HJ" ) {

              ## HJ scaling mode
                A <- L %*% S
                B <- R %*% S

            } else if ( mode == "SS" ) {

              ## similar spread mode
                q   <-  (n / m)^(1 / 4)
                sig <-  S^(1 / 2)
                A   <-  q * L %*% sig
                B   <-  (1 / q) * R %*% sig

            } else { warning ( "the scaling mode provided was not recognized. Try 'SS' or 'HJ'!" )}

          ## rotate sample points
            rotate  <- matrix(c(0, -1, 1, 0), 2)
            M <- A %*% rotate
            b <- c(B[ord.row, c(1, 2)])

          ## other default parameters

            if (is.null( tri.rgb ))          tri.rgb         <- "#19FF811A"
            if (is.null( bg.col ))           bg.col          <- "#001F3D"
            if (is.null( plot.title.col ))   plot.title.col  <- "#FFFFFF"
            if (is.null( plot.title.font ))  plot.title.font <- 1
            if (is.null( plot.title.cex ))   plot.title.cex  <- 1
            if (is.null( plot.sub.col ))     plot.sub.col    <- "#FFFFFF"
            if (is.null( plot.sub.font ))    plot.sub.font   <- 1
            if (is.null( plot.sub.cex ))     plot.sub.cex    <- 0.8
            if (is.null( plot.cex ))         plot.cex        <- 0.6
            if (is.null( plot.col ))         plot.col        <- "#FFFFFF"
            if (is.null( plot.pch ))         plot.pch        <- 21
            if (is.null( plot.xlab ))        plot.xlab       <- ""
            if (is.null( plot.ylab ))        plot.ylab       <- ""
            if (is.null( plot.xlim ))        plot.xlim       <- c(min( min( M[, 1]), b[1]), max( max( M[, 1] ), b[1] ))
            if (is.null( plot.ylim ))        plot.ylim       <- c(min( min( M[, 2]), b[2]), max( max( M[, 2] ), b[2] ))
            if (is.null( text.cex ))         text.cex        <- 0.5
            if (is.null( text.font ))        text.font       <- 2
            if (is.null( text.col ))         text.col        <- "#FFFFFF"
            if (is.null( text.pos ))         text.pos        <- 3
            if (is.null( axis.col ))         axis.col        <- "#FFFFFF"
            if (is.null( axis.cex ))         axis.cex        <- 0.7
            if (is.null( axis.font ))        axis.font       <- 1
            if (is.null( axis.asp ))         axis.asp        <- 1
            if (is.null( arrow.lwd ))        arrow.lwd       <- 1
            if (is.null( arrow.len ))        arrow.len       <- 0.1
            if (is.null( arrow.col ))        arrow.col       <- "#FFFFFF"
            if (is.null( text.col.var ))     text.col.var    <- "#FF0000"

          ## graphical parameter
            temppar <- par(bg = bg.col)
            on.exit(par(temppar), add = TRUE)

          ## plot sample points
            plot(M[, 1], M[, 2], pch = plot.pch, cex = plot.cex, col  = plot.col, xlim = plot.xlim,
                 ylim = plot.ylim, main = plot.title, col.main = plot.title.col, font.main = plot.title.font,
                 cex.main = plot.title.cex, sub = plot.sub, col.sub = plot.sub.col, font.sub = plot.sub.font,
                 cex.sub = plot.sub.cex, axes = FALSE, xlab = plot.xlab, ylab = plot.ylab, asp = axis.asp)
            text(M, labels = points.lab, cex = text.cex, font = text.font, col = text.col, pos = text.pos)
            text(x = b[1], y = b[2], labels = var.lab, cex = text.cex, font = text.font, col = text.col.var, pos = text.pos)
            axis(1, col = axis.col, col.axis = axis.col, col.ticks = axis.col, cex.axis = axis.cex, font = axis.font)
            axis(2, col = axis.col, col.axis = axis.col, col.ticks = axis.col, cex.axis = axis.cex, font = axis.font)

          ## draw the triangles
            for (j in 1:n) {
              v1 = c(0, M[j, 1], b[1])
              v2 = c(0, M[j, 2], b[2])
              polygon(x = v1, y = v2, col = tri.rgb, border = NA)
            }

          ## draw the base of the triangle
            arrows(0, 0, x1 = b[1], y1 = b[2], lwd = arrow.lwd, length = arrow.len, col = arrow.col)

        }
