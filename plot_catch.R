# Plot observed catches against modelled catches
plot_catch <- function(params, catch) {
    heights <- catch$count / catch$dl / sum(catch$count)

    # Create an empty plot with the correct x and y limits
    plot(NULL, xlim = c(min(catch$length), max(catch$length + catch$dl)),
         ylim = c(0, max(heights)), xlab = "Length [cm]", ylab = "Density",
         main = "Histogram with Areas Representing Counts")

    # Draw the bars using rect()
    rect(catch$length, 0, catch$length + catch$dl, heights,
         col = "blue", border = "blue")

    # We use lengths instead of weights in the catch data
    lengths <- (params@w / params@species_params$a)^(1/params@species_params$b)
    # Catch density as a function of weight
    model_catch <- params@initial_n * getFMort(params)
    # Normalise the density
    model_catch <- model_catch  / sum(model_catch * params@dw)
    # Convert to density as a function of length
    model_catch <- model_catch * params@species_params$b * params@w / lengths

    # Add the fitted PDF to the plot
    lines(lengths, model_catch, col = 'red', lwd = 2)
    legend('topright', legend = c('Observed', 'Modelled'),
           col = c('blue', 'red'), lwd = 2)
}
