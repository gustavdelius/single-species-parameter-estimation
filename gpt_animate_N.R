library(gganimate)

gpt_animate_N <- function(N, years, bin_start, logy = TRUE) {
    n_years <- length(years)
    n_bins <- length(bin_start)
    n_steps <- dim(N)[1] - 1

    # Create a data frame for animation
    N_df <- data.frame(
        year = rep(years, times = n_bins),
        weight = rep(bin_start, each = n_years),
        abundance = as.vector(N[seq(1, n_steps, length.out = n_years), ])
    )

    # Create the plot
    Nplot <- ggplot(N_df, aes(x = weight, y = abundance, group = year)) +
        geom_line() +
        scale_x_log10() +
        labs(x = "Weight", y = "Abundance", title = "Year: {closest_state}") +
        transition_states(year, transition_length = 2, state_length = 1) +
        ease_aes("linear")
    if (logy) {
        Nplot <- Nplot + scale_y_log10()
    }
    # Render the animation
    animate(Nplot, nframes = n_years, fps = 1)
}
