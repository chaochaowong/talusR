#' Volcano plots for TalusResults
#'
#' Valcano plots for TalusResults object
#'
#' @param res a \code{TalusResults} instance containing results returned by \code{talus_limma} or \code{talus_row_t_walch}.
#' @param alpha the significance adjust p-vluae cutoff for independent filtering. Default to 0.05.
#' @param lfc_threshold a non-negative values of significance log fold change threshold for independent filtering. Default to 0.5.
#' @param use_adjP logical to indicate wether to use 'adj.P.val' or 'P.Value'. Default to TRUE.
#' @param contrast_levels Default to NULL.
#' @param label_top_n
#' @param return_data logical to return data only with mutated `logp` and `sig` columns
#' @import ggplot2
#' @import ggrepel
#' @importFrom dplyr if_else mutate case_when
#' @importFrom purrr map_dfr
#'
#' @author Chao-Jen Wong
#' @examples
#'
#' @export
plot_volcano <- function(res, alpha = 0.05,
                         lfc_threshold = 0.5,
                         use_adjP = TRUE,
                         contrast_levels = NULL,
                         label_top_n = 10,
                         return_data = FALSE) {

  which_p <- if_else(use_adjP, 'adj.P.Val', 'P.Value')

  if (is.list(res)) {
    df <- map_dfr(res, .add_logp_and_sig,
                  which_p = which_p,
                  lfc_threshold = lfc_threshold,
                  alpha = alpha,
                  .id='contrast_name')

    if (!is.null(contrast_levels))
      df <- df %>%
        dplyr::mutate(contrast_name = factor(contrast_name,
                                             levels = contrast_levels))

    # ggplot
    gg <- .pre_plot_vocano(df, which_p = which_p,
                           lfc_threshold = lfc_threshold,
                           label_top_n = label_top_n) +
      facet_wrap( ~contrast_name, ncol=2)

  } else {
    df <- .add_logp_and_sig(res, which_p = which_p,
                            lfc_threshold = lfc_threshold,
                            alpha = alpha)

    # ggplot
    gg <- .pre_plot_vocano(df, which_p = which_p,
                           lfc_threshold = lfc_threshold,
                           label_top_n = label_top_n)
  }

  if (return_data)
    return(df)
  else
    return(gg)

}

.add_logp_and_sig <- function(x, which_p = 'adj.P.Val',
                              lfc_threshold = 0.5,
                              alpha = 0.05) {
  x %>%
    dplyr::mutate(
      logp = -log10(.data[[which_p]])) %>%
    dplyr::mutate(
      sig  = case_when(
        .data[[which_p]] < alpha & logFC >  lfc_threshold  ~ "Up",
        .data[[which_p]] < alpha & logFC < -lfc_threshold  ~ "Down",
        TRUE                        ~ "NS")
    )

}

.pre_plot_vocano <- function(df, which_p = 'adj.P.Val',
                             lfc_threshold = 0.5,
                             alpha = 0.05,
                             label_top_n = 10) {

  y_lab <- bquote(-log[10] ~ (.(as.name(which_p))))

  ggplot(df, aes(x = logFC, y = logp, color = sig)) +
    geom_point(alpha = 0.6, size = 1.5, show.legend = FALSE) +
    # 3. Add threshold lines
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    scale_color_manual(
      values = c(
        Down = "steelblue",
        NS   = "grey70",
        Up   = "firebrick"
      ),
      name = "Significance"
    ) +
    # 5. (Optional) label top hits
    geom_text_repel(
      data = df %>% filter(sig != "NS") %>%
        slice_max(order_by = logp, n = label_top_n),
      aes(label = Genes),
      size         = 2,
      box.padding  = 0.3,
      segment.alpha = 0.5,
      show.legend = FALSE
    ) +
    theme_bw() +
    labs(
      x     = expression(log[2]~Fold~Change),
      y     = y_lab
     ) +
    theme(panel.grid.minor = element_blank())
}
