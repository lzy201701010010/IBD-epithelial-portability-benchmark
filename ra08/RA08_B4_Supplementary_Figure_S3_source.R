# Single-panel RA08 evidence grid; reads only the frozen-derived display TSV.
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1]] else dirname(normalizePath(sys.frame(1)$ofile))
source_file <- file.path(out_dir, "RA08_B3_EXTERNAL_EVIDENCE_OVERVIEW_SOURCE.tsv")
d <- read.delim(source_file, sep = "\t", check.names = FALSE, fileEncoding = "UTF-8")
stopifnot(nrow(d) == 66L, length(unique(d$program)) == 11L,
          all(sort(unique(d$method)) == c("SM01", "SM02", "SM03")),
          all(sort(unique(d$cohort)) == c("Elmentaite", "PREDICT")))
stopifnot(all(d$predict_Q_lt_0_05[d$cohort == "Elmentaite"] == "NOT_APPLICABLE"),
          all(d$formal_status[d$cohort == "Elmentaite"] == "NOT_ESTIMABLE"))

d$x <- match(d$method, c("SM01", "SM02", "SM03")) +
       ifelse(d$cohort == "Elmentaite", 3, 0)
d$y <- 12 - as.integer(sub("SP", "", d$program))
d$label <- d$observed_sign
is_formal_mark <- d$cohort == "PREDICT" & d$predict_Q_lt_0_05 == "TRUE"
d$label[is_formal_mark] <- paste0(d$label[is_formal_mark], "*")
d$fill <- ifelse(d$observed_sign == "N/A", "na",
                 ifelse(d$cohort == "PREDICT", "formal", "descriptive"))
d$label_color <- ifelse(d$observed_sign == "N/A", "#62676d", "#18232f")

headers <- data.frame(x = c(2, 5), y = c(12.0, 12.0),
                      label = c("PREDICT formal secondary",
                                "Elmentaite descriptive only"))

p <- ggplot(d, aes(x = x, y = y)) +
  geom_tile(aes(fill = fill), width = 0.91, height = 0.86,
            colour = "#FFFFFF", linewidth = 0.5) +
  geom_text(aes(label = label, colour = label_color), size = 3.05,
            family = "sans", fontface = "plain", show.legend = FALSE) +
  geom_vline(xintercept = 3.5, colour = "#56616c", linewidth = 0.55) +
  geom_text(data = headers, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, family = "sans", fontface = "bold",
            size = 2.7, colour = "#26323d") +
  scale_fill_manual(values = c(formal = "#DCEAF5", descriptive = "#F5EBDD",
                               na = "#ECEFF1"), guide = "none") +
  scale_colour_identity() +
  scale_x_continuous(breaks = 1:6,
                     labels = rep(c("SM01", "SM02", "SM03"), 2),
                     limits = c(0.45, 6.55), expand = c(0, 0),
                     position = "top") +
  scale_y_continuous(breaks = 11:1,
                     labels = sprintf("SP%02d", 1:11),
                     limits = c(0.45, 12.38), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Additional external evidence across predefined epithelial programs",
       subtitle = "Observed CD-minus-comparator sign in frozen SM01-SM03 views",
       caption = "PREDICT: * after a sign denotes frozen q < 0.05 (one 30-test family).  Elmentaite: signs only; formal status NOT_ESTIMABLE.  SP08: frozen N/A.") +
  theme_minimal(base_family = "sans", base_size = 8) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.x = element_text(size = 7.4, face = "bold", colour = "#28323c", margin = margin(b = 6)),
        axis.text.y = element_text(size = 7.2, colour = "#28323c", margin = margin(r = 5)),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 10, face = "bold", colour = "#18232f", margin = margin(b = 3)),
        plot.subtitle = element_text(size = 7.4, colour = "#4a5661", margin = margin(b = 9)),
        plot.caption = element_text(size = 6.5, colour = "#3d4750", hjust = 0,
                                    margin = margin(t = 9)),
        plot.margin = margin(9, 10, 8, 8))

base <- file.path(out_dir, "RA08_B3_EXTERNAL_EVIDENCE_OVERVIEW")
ggsave(paste0(base, ".svg"), p, width = 183, height = 112, units = "mm",
       device = svglite, bg = "white")
ggsave(paste0(base, ".pdf"), p, width = 183, height = 112, units = "mm",
       device = grDevices::pdf, bg = "white", useDingbats = FALSE)
ggsave(paste0(base, ".png"), p, width = 183, height = 112, units = "mm",
       dpi = 300, bg = "white")
cat("Wrote single-panel SVG/PDF/PNG figure from", nrow(d), "cohort-method display rows\n")
