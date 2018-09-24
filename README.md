Context-dependent Space Filling curves
======================================

-   Why: 2D gridded data are harder to work with than 1D equally spaced
    data.
-   How: A space-filling algorithm inspired by the work of [Dafner,
    Cohen-Or and
    Matias (2000)](http://theory.stanford.edu/~matias/papers/eg2000.pdf).
-   What: Unlike Hilbert and Peano curves, which are universal (i.e.,
    context-independent), the algorithm proposed here adapts to the data
    provided as an input matrix. This matrix can be rectangular, and may
    even contain missing values (under certain constraints), as shown
    below. Results are both mesmerizing (IMO…) and potentially useful to
    model 2D gridded data as if it were 1D, while preserving some
    locality (spatial autocorrelation).

Examples
========

Setting up…

    library(scurvy)

    ## Loading required package: Matrix

    ## Loading required package: ggplot2

    verbose = FALSE

Let’s start with a mock square dataset. For the algorithm to work, the
data matrix must have even numbers of rows and columns.

    r = 8
    square = mapply(1:24, FUN = function(j) mapply(1:24, FUN = function(i) {
      sin((i-1)*pi/r) + cos((j-1)*pi/r) + 8 * sin((i-1)*pi/r) * cos((j-1)*pi/r)
    }))
    image(Matrix(square))

![](README_files/figure-markdown_strict/square-1.png)

The algorithm draws a context-dependent path through the gridpoints so
that it spends as much time as possible in each “bubble” before moving
on to the next one.

    s = sfc(data = square, lat = nrow(square):1, lon = 1:ncol(square), verbose = verbose)
    plotPath(s)

![](README_files/figure-markdown_strict/squarepath-1.png)

Let’s replace one of the bubbles with NAs.

    squareNA = square
    squareNA[9:16, 13:20] = NA
    image(Matrix(squareNA))

![](README_files/figure-markdown_strict/squareNA-1.png)

We can see the impact of NAs on the path.

    sNA = sfc(data = squareNA, lat = nrow(square):1, lon = 1:ncol(square), verbose = verbose)
    plotPath(sNA)

![](README_files/figure-markdown_strict/squareNApath-1.png)

Not all datasets with NAs allow valid paths. It is better to pass a full
dataset without NAs into a function like `preprocess`, included in the
package. Here, for example, we fill with NAs most values above 0, except
those needed to make a valid path. Note that NAs come in 2x2 blocks that
end in even row and column numbers.

    squareB = preprocess(data = square, thresh = 0, verbose = verbose)
    image(Matrix(squareB))

![](README_files/figure-markdown_strict/squareB-1.png)

And shown below, the path adapts nicely to the NAs in the image.

    s = sfc(data = squareB, lat = nrow(square):1, lon = 1:ncol(square), verbose = verbose)
    plotPath(s)

![](README_files/figure-markdown_strict/squareBpath-1.png)

OK, enough mock data. Let us look at a 2 degree gridded global
topography of the world.

    cETOPO4 = getCoarseETOPO(4)
    image(Matrix(cETOPO4$data))

![](README_files/figure-markdown_strict/coarseEtopo-1.png)

We are only interested in bathymetry (depth &lt; 0 m), though, so we use
the function `preprocessBathymetry`.

    cETOPO4$data = preprocessBathymetry(cETOPO4$data, neritic = -1, verbose = verbose)
    s4 = sfc(data = cETOPO4$data, lat = cETOPO4$lat, lon = cETOPO4$lon, verbose = verbose)
    plotPath(s4)

![](README_files/figure-markdown_strict/coarsepath-1.png)

Finally, we run the algorithm on the whole ETOPO30 bathymetry data.
Warning: this may take an hour to run, because it works with a large
(64k x 64k) sparse matrix.

    data = preprocessBathymetry(etopo30, neritic = -2000, verbose = verbose)
    lat = seq(89.75, -89.75, by=-0.5)
    lon = seq(-179.75, 179.75, by=0.5)
    sEtopo = sfc(data = data, lat = lat, lon = lon, verbose = verbose)
    plotPath(sEtopo)

![](README_files/figure-markdown_strict/etopopath-1.png)
