# Plot observed catches against modelled catches
plot_catch <- function(params, catch) {

    hist_data <- data.frame(
        length = catch$length + catch$dl / 2,
        count = catch$count
    )
    plot(hist_data$length, hist_data$count, type = 'h', lwd = 2, col = 'blue',
         xlab = 'Length', ylab = 'Count', main = 'Observed Data and Fitted PDF')

    model_catch <- params@initial_n * getFMort(params) * params@dw
    # Scale the catch for visualization
    model_catch <- model_catch * max(hist_data$count) / max(model_catch)
    # We use lengths instead of weights in the catch data
    lengths <- (params@w / params@species_params$a)^(1/params@species_params$b)

    # Add the fitted PDF to the plot
    lines(lengths, model_catch, col = 'red', lwd = 2)
    legend('topright', legend = c('Observed Counts', 'Fitted PDF'),
           col = c('blue', 'red'), lwd = 2)
}
