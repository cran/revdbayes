## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)

## ----fig.width = 5, fig.align='center'----------------------------------------
library(revdbayes)
# Set a threshold at the 90% quantile
thresh <- quantile(newlyn, probs = 0.90)
postsim <- kgaps_post(newlyn, thresh, k = 1)
plot(postsim, xlab = expression(theta))

