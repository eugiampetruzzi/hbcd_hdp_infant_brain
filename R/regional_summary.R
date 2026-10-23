parcel_network_summary <- function(results) {
  thick <- results[results$metric == "Parcel cortical thickness", ]
  thick$network <- sub("^([LR]_)?(.*)_[0-9]+$", "\\2", thick$feature)
  thick$thicker <- thick$q < .05 & thick$beta > 0
  thick$thinner <- thick$q < .05 & thick$beta < 0
  aggregate(
    cbind(parcels = 1, significant = thicker, thinner) ~ network,
    data = thick,
    FUN = sum
  )
}
