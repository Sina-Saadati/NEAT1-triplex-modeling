suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork); library(ragg); library(ggrepel)})
DATA_DIR <- "data"; OUT_DIR <- "figures"; dir.create(OUT_DIR, showWarnings = FALSE)
FONT <- "Times New Roman"   # change if unavailable on your system
U <- c(14, 23); A <- c(82, 101)
col <- c(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00", grey="#6E6E6E", purple="#CC79A7")
th <- theme_classic(base_size = 18, base_family = FONT) +
  theme(axis.title = element_text(face = "bold"), axis.line = element_line(linewidth = 0.7))
minus <- function(x) gsub("-", "\u2212", x)

find_file <- function(name) {
  for (ext in c(".xvg", ".txt", ".md", ".fasta", ".ct", ".csv")) {
    p <- file.path(DATA_DIR, paste0(name, ext)); if (file.exists(p)) return(p)
  }
  stop("File not found: ", name)
}
save_fig <- function(p, name, w, h) {
  for (ext in c("png", "svg")) {
    f <- file.path(OUT_DIR, paste0(name, ".", ext))
    if (ext == "png") ggsave(f, p, device = agg_png, width = w, height = h, units = "in", dpi = 600)
    else ggsave(f, p, device = svglite::svglite, width = w, height = h)
  }
  cat("Saved:", name, "\n")
}
read_xvg <- function(path, ncols, cols) {
  ln <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  ln <- ln[nzchar(ln) & !grepl("^[#@]", ln)]
  df <- read.table(text = ln, header = FALSE)[, seq_len(ncols), drop = FALSE]
  names(df) <- cols; df
}
read_fasta_md <- function(path) {
  ln <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  ln <- ln[nzchar(ln) & !grepl("^(#|\u0060\u0060\u0060)", ln)]
  nm <- character(); sq <- character()
  for (l in ln) {
    if (startsWith(l, ">")) { nm <- c(nm, sub("^>", "", l)); sq <- c(sq, "") }
    else sq[length(sq)] <- paste0(sq[length(sq)], gsub("T", "U", toupper(l)))
  }
  names(sq) <- nm; sq
}
parse_info <- function(path) {
  ln <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  idx <- grep("^Pocket [0-9]+ :", ln); out <- list()
  for (k in seq_along(idx)) {
    blk <- ln[idx[k]:(if (k < length(idx)) idx[k + 1] - 1 else length(ln))]
    getv <- function(key) as.numeric(sub(".*:\\s*", "", blk[grepl(paste0("^", key), blk)][1]))
    out[[k]] <- data.frame(pocket = as.integer(sub("Pocket ([0-9]+) :.*", "\\1", ln[idx[k]])),
                           score = getv("Score"), drug = getv("Druggability Score"), volume = getv("Volume :"))
  }
  bind_rows(out)
}