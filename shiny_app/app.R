# ============================================================================
# BayeSPeer — Bias-Aware Aggregation of Peer & Self Rankings
# Shiny companion app to:
#   Ruan, Wang, Lim & Wang, "Bias Meets Bayes: Bayesian Aggregation for
#   Peer and Self Ranking"
#
# Upload an n x n peer + self ranking matrix (columns = rankers, rows =
# individuals, diagonal = self-ranks, 1 = best). The app returns aggregated
# rankings from BayeSPeer and the 13 baseline methods compared in the paper,
# plus the real-data analyses of Section 4 (rank probability profiles,
# posterior distributions, convergence diagnostics, and the bias-vs-
# performance correlation analysis of Fig. 9).
#
# Run with:  shiny::runApp("path/to/this/folder")
# Requires:  R >= 4.0 and the packages shiny, ggplot2, GGally (the figures
#            reuse the paper's ggplot2/GGally code), truncnorm, TopKLists and
#            RobustRankAggreg (the paper's sampler and baseline methods).
# ============================================================================

local({
  need <- c("shiny", "ggplot2", "GGally", "truncnorm", "TopKLists", "RobustRankAggreg")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Please install the missing packages first: install.packages(c(",
         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE)
})
library(shiny)
library(ggplot2)
library(GGally)
source("methods.R", local = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# --------------------------------------------- paper styling (Figs 7 to 9)
# Method colors and order exactly as in the paper's figure code
# (05a_RWD_2020.R). The .selfrm variants are not shown in the paper; they use
# a lighter shade of their method's color and sit right after it.
PAPER_METHOD_COLORS <- c(
  BayeSPeer = "gray30",  BIRRA   = "#A3B7C2", RRA     = "#A9C1A3",
  Stuart    = "#9CBAA0", min     = "#C9D3B5", MC1     = "#E6D3A3",
  MC2       = "#D8C49C", MC3     = "#C1B987", CEMC.s  = "#D8A48F",
  CEMC.k    = "#C38D8D", Mean    = "#B9A0C6", Geomean = "#B8A5B3",
  L2norm    = "#C0B7D1", Median  = "#A0A7B8")
# density line types of the paper's Fig. 8, recycled for larger groups
PAPER_LINETYPES <- c("solid", "dotted", "solid", "longdash", "dotted",
                     "solid", "longdash", "dotted", "solid", "longdash")
.shade <- function(hex, f) {           # f in [0,1]: mix toward white
  m <- col2rgb(hex)
  rgb(t(m + (255 - m) * f), maxColorValue = 255)
}
method_levels <- function(methods) {
  m0  <- names(PAPER_METHOD_COLORS)
  lev <- as.vector(rbind(m0, paste0(m0, ".selfrm")))
  lev[lev %in% methods]
}
method_palette <- function(methods)
  vapply(setNames(methods, methods), function(m) {
    col <- PAPER_METHOD_COLORS[[sub("\\.selfrm$", "", m)]]
    if (grepl("\\.selfrm$", m)) .shade(col, 0.5) else col
  }, character(1))
SEQ_RAMP <- c("#FFFFFF", "#BDD7EE", "#0072B2")   # rank-probability heatmap

# ------------------------------------------------------------ size system
# Every text size in the app (page text, tables, figures) comes from TXT, a
# modular scale anchored at BASE_PX: edit BASE_PX to make everything larger or
# smaller together. Page sizes are CSS px; figure sizes are points, and since
# renderPlot draws at 72 dpi (Shiny rescales for high-DPI screens), 1 pt
# renders as 1 CSS px, so figure text matches the page text beside it. The
# figures keep the paper's design; only their text sizes come from TXT.
BASE_PX <- 14
size_step <- function(k) round(BASE_PX * 1.125^k)     # ratio 9:8 per step
TXT <- c(
  note   = size_step(-1),  # 12  secondary notes; the smallest text anywhere
  body   = size_step(0),   # 14  body text and table cells
  head   = size_step(2),   # 18  section headings (h4)
  tick   = size_step(-1),  # 12  axis tick labels
  annot  = size_step(-1),  # 12  values printed on bars and heatmap cells
  label  = size_step(0),   # 14  axis titles and variable names
  legend = size_step(0),   # 14  legend entries
  title  = size_step(1),   # 16  panel titles
  emph   = size_step(3)    # 20  headline numbers (correlations)
)

# Figure geometry: the paper's proportions, capped so each figure fits on a
# laptop screen. Maximum widths in px: max_w for the rank probability heatmap,
# dens_w for each posterior density panel (Fig. 8), pairs_w for the pairs plot
# (Fig. 9); Fig 7 rows are a strip plus a panel of row_min to row_max px.
FIG <- list(max_w = 820, dens_w = 640, pairs_w = 600,
            strip = 30, row_min = 100, row_max = 140)

# page and table CSS, generated from the same scale; tables are ruled like
# the paper's (rule above, below the header, and at the bottom; no stripes)
size_css <- sprintf(paste(
  "body { font-size:%1$dpx; }",
  ".small-note { font-size:%2$dpx; }",
  "h4 { font-size:%3$dpx; }",
  ".shiny-table { font-size:%1$dpx; border-top:1px solid #333;",
  "  border-bottom:1px solid #333; }",
  ".shiny-table > thead > tr > th { border-bottom:1px solid #333 !important; }",
  ".shiny-table > tbody > tr > td { border-top:none !important; }",
  ".table-wrap { overflow-x:auto; }",
  ".grid-head { font-size:%2$dpx; color:#6B6B6B; text-align:center; }",
  sep = "\n"), TXT[["body"]], TXT[["note"]], TXT[["head"]])
size_css <- paste0(size_css, sprintf(
  "\n.fig-narrow { max-width:%dpx; }\n.fig-dens { max-width:%dpx; }\n.fig-pairs { max-width:%dpx; }",
  FIG$max_w, FIG$dens_w, FIG$pairs_w))

# ggplot2: element_text() sizes are points, geom_text() sizes are mm
txt_mm <- function(role) TXT[[role]] / ggplot2::.pt

# Aggregated ranks are integers or exact halves (average tie-breaking)
rank_label <- function(v) {
  h <- round(v * 2) / 2
  ifelse(h == round(h), sprintf("%d", as.integer(round(h))),
         paste0(floor(h), "\u00bd"))
}

# Aggregated-rankings layout (paper Fig. 7): as many facets per row as fit
# without bar value labels touching; the height leaves room for the rotated
# method names and for the legend.
fig7_geom <- function(width, tab) {
  if (is.null(width) || !is.finite(width) || width < 100) width <- 1000
  n <- ncol(tab); nb <- nrow(tab)
  w <- function(role, k) k * 0.6 * TXT[[role]]            # approx text width, px
  # bar value labels: digits are ~0.56 em, a half sign up to ~0.84 em wide;
  # neighbouring bars need their centres at least half of both labels apart
  lab_px <- label_px(tab)
  slot <- max((lab_px[-1, , drop = FALSE] + lab_px[-nb, , drop = FALSE]) / 2) + 1
  lcols <- ceiling(nb / 20); lrows <- ceiling(nb / lcols)  # ggplot2 legend wrap
  legend_w <- max(1.05 * w("legend", 11),
                  lcols * (w("note", max(nchar(rownames(tab)))) + 35)) + 20
  room <- width - legend_w - 48
  fit <- max(1, min(4, n, floor(room / (slot * (nb + 1.2) + 8))))
  # prefer full rows (a ragged last row pushes method names into the row below);
  # keep a ragged layout only when full rows would cost more than one extra row
  ncol <- max(Filter(function(k) n %% k == 0, seq_len(fit)))
  if (ceiling(n / ncol) > ceiling(n / fit) + 1) ncol <- fit
  nrow <- ceiling(n / ncol); ragged <- n %% ncol != 0
  names_h <- w("tick", max(nchar(rownames(tab)))) * sin(55 * pi / 180) + 10
  # facets about twice as wide as tall, as in the paper; flatter when there
  # are many rows, so large groups stay as short as the labels allow
  body <- if (nrow > 4) FIG$row_min else min(FIG$row_max, max(FIG$row_min, room / ncol / 2))
  row_h <- FIG$strip + body +
           if (ragged) names_h else 0                # ragged: names under every row
  legend_h <- (lrows + 2.5) * 1.3 * TXT[["legend"]] + 30
  list(ncol = ncol, nrow = nrow, ragged = ragged, body = body,
       bar_px = (room / ncol - 8) / (nb + 1.2),            # px between bar centres
       height = max(nrow * row_h + (if (ragged) 0 else names_h) + 25, legend_h))
}

# approximate px width of bar value labels: digits ~0.56 em, a half sign up to
# ~0.84 em (Helvetica/Arial)
label_px <- function(v) (0.56 * nchar(floor(v)) + 0.84 * (v %% 1 != 0)) * TXT[["annot"]]

# Made-up exam scores so users can try the external-score analysis: simulated
# from each person's average peer rank plus noise. NOT real data.
demo_scores <- function(R, seed = 2026) {
  if (exists(".Random.seed", envir = globalenv())) {   # leave the session RNG alone
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()))
  }
  set.seed(seed)
  n <- nrow(R)
  peer <- vapply(seq_len(n), function(i) mean(R[i, -i]), numeric(1))
  z <- if (sd(peer) > 0) -(peer - mean(peer)) / sd(peer) else rep(0, n)
  round(pmin(100, pmax(40, 72 + 9 * z + rnorm(n, 0, 7))))
}

# ------------------------------------------------------------------ UI ------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { background:#FAFAF8; }
    .well { background:#FFFFFF; border:1px solid #E5E2DC; border-radius:10px;
            padding:12px 16px; }
    h4 { margin-top:18px; font-weight:600; }
    .well h4 { margin:2px 0 8px; }
    .well h4 ~ h4 { margin-top:12px; }
    .well .form-group, .well .checkbox { margin-top:0; margin-bottom:6px; }
    .well .shiny-file-input-progress { display:none; }
    .well .small-note { margin:2px 0 8px; line-height:1.35; }
    .well .row + .small-note { margin-top:0; }
    .well .run-btn { font-weight:600; }
    .btn-primary { background:#B4551F; border-color:#B4551F; }
    .btn-primary:hover { background:#93441A; border-color:#93441A; }
    .nav-tabs > li.active > a { font-weight:600; }
    .small-note { color:#6B6B6B; }
    table.data { font-variant-numeric: tabular-nums; }
  ")), tags$style(HTML(size_css))),
  titlePanel(div("BayeSPeer: Bias-Aware Peer & Self Ranking Aggregation",
                 div(class = "small-note",
                     "Companion app to “Bias Meets Bayes: Bayesian Aggregation for Peer and Self Ranking” (Ruan, Wang, Lim & Wang)"))),
  sidebarLayout(
    sidebarPanel(width = 4,
      h4("1 · Data"),
      fileInput("file", "Upload ranking matrix (CSV)",
                accept = c(".csv", ".txt", ".tsv")),
      div(class = "small-note",
          "n×n matrix: columns = rankers, rows = individuals, diagonal = self-ranks, 1 = best. Header row / ID column are auto-detected."),
      fluidRow(
        column(6, actionButton("ex1", "Class I (8)", width = "100%")),
        column(6, actionButton("ex2", "Class II (16)", width = "100%"))),
      div(class = "small-note", style = "margin-top:4px;", "Example datasets in the paper."),
      fluidRow(
        column(6, numericInput("n_manual", NULL, 8, min = 3, max = 20, step = 1)),
        column(6, actionButton("mk_grid", "Type it in", width = "100%"))),
      div(class = "small-note",
          "Or enter the matrix by hand: set the group size n (3–20) and click ‘Type it in’ to get an editable grid on the Data tab."),
      h4("2 · Settings"),
      checkboxInput("selfrm", "Add baselines without self-ranks (.selfrm)", FALSE),
      selectInput("prior", "Variance prior", width = "100%", selected = "IG_0.1",
        choices = list(
          "Inverse gamma" = c("IG(0.1, 0.1), the paper's default" = "IG_0.1",
                              "IG(0.001, 0.001)" = "IG_0.001", "IG(0.01, 0.01)" = "IG_0.01",
                              "IG(1, 1)" = "IG_1"),
          "Half-Cauchy on the standard deviations" = c(
                              "Half-Cauchy(0, 0.5)" = "HC_0.5", "Half-Cauchy(0, 1)" = "HC_1",
                              "Half-Cauchy(0, 2)" = "HC_2", "Half-Cauchy(0, 4)" = "HC_4"))),
      div(class = "small-note", "Used for all variances in both steps; step 2 sets its parameters from the step 1 draws. Other options: the sensitivity analysis priors (Section S2)."),
      fluidRow(
        column(7, numericInput("iters", "Iterations per step", 4000, min = 2000, step = 1000)),
        column(5, numericInput("seed", "Random seed", 123, min = 1, step = 1))),
      div(class = "small-note", "Each of the two steps runs this many iterations, the first half as burn-in. With 2,000 iterations, a run takes about 1 second for 5 individuals and 17 seconds for 25."),
      uiOutput("status"),
      actionButton("run", "Run analysis", class = "btn-primary run-btn", width = "100%")
    ),
    mainPanel(width = 8,
      tabsetPanel(id = "tabs",
        tabPanel("Data",
          br(), uiOutput("manual_ui"), uiOutput("data_msg"), div(class = "table-wrap", tableOutput("data_tab"))),
        tabPanel("Aggregated rankings",
          br(),
          p("Final aggregated rank for every individual under each method (1 = best), as in Fig. 7 of the paper. BayeSPeer uses the full data including self-ranks and corrects for self-evaluation bias; baselines treat all ranks as unbiased."),
          plotOutput("fig7", height = "auto"),
          br(), div(class = "table-wrap", tableOutput("rank_tab")),
          downloadButton("dl_ranks", "Download rankings (CSV)")),
        tabPanel("Rank probabilities",
          br(),
          p("BayeSPeer posterior probability that each individual occupies each rank position (Table 4 of the paper): the share of post burn-in draws in which the individual's μᵢ has that rank among all μ's."),
          div(class = "fig-narrow", plotOutput("heat", height = "auto")),
          br(), div(class = "table-wrap", tableOutput("prob_tab")),
          downloadButton("dl_prob", "Download probabilities (CSV)")),
        tabPanel("Posteriors & bias",
          br(),
          p("Posterior distributions of true performance μᵢ and self-evaluation bias βⱼ (Fig. 8 of the paper), and the posterior summary. βⱼ is the shift in ranker j's latent score for themselves (ωⱼⱼ = μⱼ + βⱼ + εⱼⱼ), so βⱼ > 0 means ranking oneself higher than one's performance μⱼ implies."),
          div(class = "fig-dens", plotOutput("dens_mu", height = "auto")),
          div(class = "fig-dens", plotOutput("dens_beta", height = "auto")),
          br(), div(class = "table-wrap", tableOutput("summ_tab")),
          downloadButton("dl_summ", "Download posterior summary (CSV)")),
        tabPanel("Bias vs performance",
          br(),
          p("Pairwise correlations as in Fig. 9 of the paper, among an optional external measure, estimated performance μ̂ᵢ (posterior mean), self-bias β̂ⱼ (posterior mean), the observed rank difference (rank from the average peer rank, self excluded, minus the self-rank), and ranker variance σ̂²ⱼ (posterior median)."),
          uiOutput("ext_ui"),
          div(class = "fig-pairs", plotOutput("pairs", height = "auto")),
          br(), div(class = "table-wrap", tableOutput("cor_tab")),
          uiOutput("interp")),
        tabPanel("About",
          br(), uiOutput("about"))
      )
    )
  )
)

# --------------------------------------------------------------- server -----
server <- function(input, output, session) {

  rv <- reactiveValues(R = NULL, src = NULL, res = NULL, manual_n = NULL,
                       ext = NULL, ext_label = "Exam Grade")

  # ---------------- data input
  observeEvent(input$ex1, { rv$R <- example_class1(); rv$src <- "Example: Class I (8 teams, paper Table 3)"; rv$res <- NULL; rv$ext <- NULL })
  observeEvent(input$ex2, { rv$R <- example_class2(); rv$src <- "Example: Class II (16 students, paper Table 5)"; rv$res <- NULL; rv$ext <- NULL })

  observeEvent(input$file, {
    tryCatch({
      path <- input$file$datapath
      sep <- if (grepl("\\.tsv$", input$file$name, ignore.case = TRUE)) "\t" else ","
      first <- readLines(path, n = 1)
      toks <- strsplit(first, sep, fixed = TRUE)[[1]]
      has_header <- any(is.na(suppressWarnings(as.numeric(toks[-1]))))
      d <- utils::read.table(path, sep = sep, header = has_header,
                             stringsAsFactors = FALSE, check.names = FALSE)
      if (!is.numeric(d[[1]])) { rn <- as.character(d[[1]]); d <- d[, -1, drop = FALSE]; rownames(d) <- rn }
      M <- as.matrix(d); storage.mode(M) <- "double"
      if (is.null(rownames(M)) || all(rownames(M) == as.character(seq_len(nrow(M)))))
        rownames(M) <- paste0("Item", seq_len(nrow(M)))
      if (nrow(M) == ncol(M)) colnames(M) <- paste0("Ranker", seq_len(ncol(M)))
      rv$R <- validate_rank_matrix(M)
      rv$src <- paste0("Uploaded: ", input$file$name)
      rv$res <- NULL
      rv$ext <- NULL
    }, error = function(e) {
      rv$R <- NULL; rv$src <- NULL; rv$res <- NULL
      showModal(modalDialog(title = "Could not use this file",
                            tags$p(conditionMessage(e)),
                            tags$p(class = "small-note",
                                   "Expected a square CSV: one column per ranker, one row per individual, each column a complete permutation of 1..n, diagonal = self-ranks."),
                            easyClose = TRUE))
    })
  })

  # ---------------- manual matrix entry
  observeEvent(input$mk_grid, {
    n <- suppressWarnings(as.integer(input$n_manual))
    if (is.na(n) || n < 3 || n > 20) {
      showModal(modalDialog(title = "Group size out of range",
                            "Manual entry supports n between 3 and 20 (use file upload for larger groups).",
                            easyClose = TRUE))
      return()
    }
    rv$manual_n <- n
    updateTabsetPanel(session, "tabs", selected = "Data")
  })

  output$manual_ui <- renderUI({
    n <- rv$manual_n
    if (is.null(n)) return(NULL)
    # prefill from the currently loaded matrix when sizes match (lets you edit an
    # example or an uploaded dataset); otherwise start blank
    pre <- if (!is.null(rv$R) && nrow(rv$R) == n) rv$R else NULL
    cell <- function(i, j) {
      style <- paste0("width:46px;text-align:center;border:1px solid #D8D4CC;",
                      "border-radius:4px;padding:3px 1px;",
                      if (i == j) "background:#FDEBD3;font-weight:600;" else "background:#FFF;")
      tags$input(id = sprintf("m_%d_%d", i, j), type = "text",
                 inputmode = "numeric", value = if (is.null(pre)) "" else pre[i, j],
                 style = style)
    }
    ids <- if (!is.null(pre) && !is.null(rownames(pre))) rownames(pre)
           else paste0("Item", seq_len(n))
    tagList(
      h4(sprintf("Manual entry (%d × %d)", n, n), style = "margin-top:4px;"),
      p(class = "small-note",
        "Each column is one ranker's complete ranking (a permutation of 1..n, 1 = best). Shaded diagonal cells are self-ranks."),
      textInput("manual_ids", "IDs (optional, comma-separated, row order)",
                value = paste(ids, collapse = ", "), width = "100%"),
      div(style = "overflow-x:auto;",
        tags$table(style = "border-collapse:separate;border-spacing:3px;",
          tags$tr(tags$th(""),
                  lapply(seq_len(n), function(j)
                    tags$th(paste0("R", j), class = "grid-head"))),
          lapply(seq_len(n), function(i)
            tags$tr(tags$td(tags$span(class = "small-note", paste0("row ", i)),
                            style = "padding-right:6px;white-space:nowrap;"),
                    lapply(seq_len(n), function(j) tags$td(cell(i, j))))))),
      br(),
      actionButton("use_manual", "Use this matrix", class = "btn-primary"),
      actionButton("hide_manual", "Close"),
      hr())
  })

  observeEvent(input$hide_manual, rv$manual_n <- NULL)

  observeEvent(input$use_manual, {
    n <- req(rv$manual_n)
    M <- matrix(NA_real_, n, n)
    for (i in seq_len(n)) for (j in seq_len(n)) {
      v <- input[[sprintf("m_%d_%d", i, j)]]
      M[i, j] <- suppressWarnings(as.numeric(trimws(v %||% "")))
    }
    if (any(is.na(M))) {
      bad <- which(is.na(M), arr.ind = TRUE)[1, ]
      showModal(modalDialog(title = "Incomplete matrix",
                            sprintf("Cell (row %d, ranker %d) is empty or not a number.", bad[1], bad[2]),
                            easyClose = TRUE))
      return()
    }
    ids <- trimws(strsplit(input$manual_ids %||% "", ",")[[1]])
    rownames(M) <- if (length(ids) == n && all(nzchar(ids))) ids
                   else paste0("Item", seq_len(n))
    colnames(M) <- paste0("Ranker", seq_len(n))
    tryCatch({
      rv$R <- validate_rank_matrix(M)
      rv$src <- "Manual entry"
      rv$res <- NULL
      rv$ext <- NULL
      rv$manual_n <- NULL
    }, error = function(e)
      showModal(modalDialog(title = "Invalid ranking matrix",
                            tags$p(conditionMessage(e)), easyClose = TRUE)))
  })

  # ---------------- external scores (entered on the Bias vs performance tab)
  output$ext_ui <- renderUI({
    req(rv$res)
    ids <- rv$res$fit$ids; n <- length(ids)
    pre <- if (!is.null(rv$ext) && length(rv$ext) == n) rv$ext else rep("", n)
    boxes <- lapply(seq_len(n), function(i)
      div(style = "display:inline-block;margin:0 10px 8px 0;text-align:center;",
          div(class = "small-note", ids[i]),
          tags$input(id = sprintf("ext_%d", i), type = "text", inputmode = "decimal",
                     value = pre[i],
                     style = "width:64px;text-align:center;border:1px solid #D8D4CC;border-radius:4px;padding:3px 2px;")))
    wellPanel(style = "background:#FFFFFF;padding:12px 16px;",
      h4("External scores (optional)", style = "margin-top:2px;"),
      p(class = "small-note",
        "Enter each individual's exam grade or other independent performance score to add it to the correlation analysis. Leave all boxes empty to analyze without it. To try the analysis out, use the made-up example scores."),
      div(style = "max-width:220px;",
          textInput("ext_label", "Score label", value = rv$ext_label)),
      div(boxes),
      actionButton("use_ext", "Apply scores", class = "btn-primary"),
      actionButton("demo_ext", "Use made-up example scores"),
      actionButton("clear_ext", "Clear"),
      if (ext_is_demo())
        p(class = "small-note", style = "color:#B4551F;margin:10px 0 0;",
          "These scores are made up for illustration (simulated from each person's average peer rank plus random noise). They are not real exam grades, so correlations involving them only demonstrate the analysis."))
  })

  ext_is_demo <- reactive(!is.null(rv$ext) && !is.null(rv$R) &&
                            isTRUE(all.equal(rv$ext, demo_scores(rv$R))))

  observeEvent(input$demo_ext, {
    req(rv$res)
    rv$ext <- demo_scores(rv$R)
    rv$ext_label <- "Exam Grade (made up)"
  })

  observeEvent(input$use_ext, {
    n <- rv$res$n
    v <- vapply(seq_len(n), function(i)
      suppressWarnings(as.numeric(trimws(input[[sprintf("ext_%d", i)]] %||% ""))),
      numeric(1))
    if (all(is.na(v))) { rv$ext <- NULL; return() }
    if (any(is.na(v))) {
      showModal(modalDialog(title = "Incomplete scores",
                            sprintf("Missing or non-numeric score for: %s. Enter a value for every individual (or clear all).",
                                    paste(rv$res$fit$ids[is.na(v)], collapse = ", ")),
                            easyClose = TRUE))
      return()
    }
    rv$ext <- v
    lab <- if (nzchar(trimws(input$ext_label %||% ""))) trimws(input$ext_label) else "External score"
    # a "(made up)" tag only stays on while the made-up scores are in use
    if (!isTRUE(all.equal(v, demo_scores(rv$R)))) lab <- sub(" \\(made up\\)$", "", lab)
    rv$ext_label <- lab
  })

  observeEvent(input$clear_ext, { rv$ext <- NULL; rv$ext_label <- "Exam Grade" })

  output$data_msg <- renderUI({
    if (is.null(rv$R))
      return(tagList(
        p("Upload a ranking matrix or load an example dataset to begin."),
        p(class = "small-note",
          "Format: entry (i, j) is the rank that ranker j gives individual i (1 = best). Every column must be a complete permutation of 1..n; the diagonal holds self-ranks.")))
    n <- nrow(rv$R)
    p(strong(rv$src), sprintf(" — %d individuals × %d rankers. Diagonal (self-ranks) shown in bold.", n, n))
  })

  output$data_tab <- renderTable({
    req(rv$R)
    D <- as.data.frame(rv$R)
    for (j in seq_len(ncol(D))) D[[j]] <- as.character(D[[j]])
    for (j in seq_len(ncol(D))) D[j, j] <- paste0("<b>", D[j, j], "</b>")
    cbind(" " = rownames(rv$R), D)
  }, sanitize.text.function = identity, striped = FALSE, spacing = "xs", align = "r")

  # ---------------- run
  observeEvent(input$run, {
    req(rv$R)
    M <- 2 * max(1000, round((input$iters %||% 4000) / 2))   # even, so burn-in = M/2
    withProgress(message = "BayeSPeer Gibbs sampler", value = 0, {
      prog <- function(f) setProgress(value = f,
        detail = if (f < 0.5) "Step 1: vague priors" else "Step 2: empirical-Bayes priors")
      pr <- strsplit(input$prior %||% "IG_0.1", "_", fixed = TRUE)[[1]]
      fit <- bayespeer_fit(rv$R, M = M, burnin = M / 2, seed = input$seed,
                           prior = pr[1], prior_value = as.numeric(pr[2]), progress = prog)
      setProgress(value = 1, detail = "Baseline methods…")
      agg <- aggregate_all(rv$R, fit = fit, include_selfrm = isTRUE(input$selfrm),
                           seed = input$seed)
      rv$res <- list(fit = fit, tab = agg$table, n = nrow(rv$R))
    })
    updateTabsetPanel(session, "tabs", selected = "Aggregated rankings")
  })

  output$status <- renderUI({
    if (is.null(rv$R)) return(p(class = "small-note", "No data loaded."))
    if (is.null(rv$res)) return(p(class = "small-note", "Data ready — press Run."))
    s <- rv$res$fit$settings
    p(class = "small-note", sprintf(
      "Done: %s iterations per step (burn-in %s), seed %d, %s prior, %d methods.",
      format(s$M, big.mark = ","), format(s$burnin, big.mark = ","),
      s$seed, s$prior_label, nrow(rv$res$tab)))
  })

  # ---------------- Fig 7: aggregated rankings
  fig7_geom_now <- function() {
    req(rv$res)
    fig7_geom(session$clientData$output_fig7_width, rv$res$tab)
  }
  output$fig7 <- renderPlot({
    req(rv$res)
    tab <- rv$res$tab
    g <- fig7_geom_now()
    lev <- method_levels(rownames(tab))
    d <- data.frame(
      method  = factor(rep(rownames(tab), ncol(tab)), levels = lev),
      team    = factor(rep(colnames(tab), each = nrow(tab)), levels = colnames(tab)),
      ranking = as.vector(tab))
    d$label <- rank_label(d$ranking)
    # where two neighbouring value labels would touch, lift the right one a line
    n <- rv$res$n
    U0  <- n * (1 + (TXT[["annot"]] + 4) / g$body)          # y range without lifts
    upx <- g$body / U0                                      # px per rank unit
    lift <- (TXT[["annot"]] + 2) / upx
    d$y_lab <- d$ranking
    for (tm in levels(d$team)) {
      i <- which(d$team == tm); i <- i[order(d$method[i])]
      lifted <- FALSE
      for (k in seq_along(i)[-1]) {
        a <- i[k - 1]; b <- i[k]
        touch <- (label_px(d$ranking[a]) + label_px(d$ranking[b])) / 2 + 2 > g$bar_px &&
                 abs(d$y_lab[b] - d$y_lab[a]) * upx < TXT[["annot"]] + 2
        if (touch && !lifted) { d$y_lab[b] <- d$y_lab[b] + lift; lifted <- TRUE }
        else lifted <- FALSE
      }
    }
    top <- max(U0, max(d$y_lab) + (TXT[["annot"]] + 4) / upx)
    # paper Fig. 7 code (05a_RWD_2020.R); text sizes from the size system
    ggplot(d, aes(x = method, y = ranking, fill = method)) +
      geom_bar(stat = "identity", position = "dodge", linewidth = 0.5) +
      labs(title = NULL, x = NULL, y = "Aggregated Ranking") +
      guides(fill = guide_legend(title = "Aggregation\nMethod")) +
      geom_text(aes(y = y_lab, label = label), vjust = 0, size = txt_mm("annot")) +
      scale_fill_manual(values = method_palette(lev)) +
      facet_wrap(~team, ncol = g$ncol, scales = if (g$ragged) "free_x" else "fixed") +
      scale_y_continuous(limits = c(0, top),
                         breaks = Filter(function(b) b <= n, pretty(c(0, n)))) +
      theme_bw() +
      theme(strip.text.x = element_text(size = TXT[["title"]], face = "bold"),
            axis.title.y = element_text(size = TXT[["label"]], face = "bold"),
            axis.text.y  = element_text(size = TXT[["tick"]]),
            legend.title = element_text(size = TXT[["legend"]], face = "bold"),
            legend.text  = element_text(size = TXT[["note"]]),
            axis.text.x  = element_text(size = TXT[["tick"]], angle = 55,
                                        vjust = 1, hjust = 1))
  }, height = function() fig7_geom_now()$height)

  output$rank_tab <- renderTable({
    req(rv$res)
    tab <- rv$res$tab[method_levels(rownames(rv$res$tab)), , drop = FALSE]
    cbind(Method = rownames(tab), as.data.frame(tab))
  }, striped = FALSE, spacing = "xs", digits = 1)

  output$dl_ranks <- downloadHandler(
    filename = function() "aggregated_rankings.csv",
    content = function(f) write.csv(rv$res$tab, f))

  # ---------------- rank probability heatmap (Table 4)
  rank_rows <- function(P) sub("^Rank", "Rank ", rownames(P))   # as in Table 4
  output$heat <- renderPlot({
    req(rv$res)
    P <- rv$res$fit$rank_prob
    d <- data.frame(
      rank = factor(rep(rank_rows(P), ncol(P)), levels = rev(rank_rows(P))),
      item = factor(rep(colnames(P), each = nrow(P)), levels = colnames(P)),
      p    = as.vector(P))
    d$lab <- ifelse(d$p >= 0.005, sprintf("%.2f", d$p), "")
    ggplot(d, aes(x = item, y = rank, fill = p)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = lab, colour = p > 0.6), size = txt_mm("annot"),
                show.legend = FALSE) +
      scale_colour_manual(values = c(`FALSE` = "#333333", `TRUE` = "white")) +
      scale_fill_gradientn(colours = SEQ_RAMP, limits = c(0, 1),
                           name = "Posterior\nprobability") +
      labs(x = NULL, y = NULL) +
      theme_minimal() +
      theme(axis.text.y  = element_text(size = TXT[["tick"]]),
            axis.text.x  = element_text(size = TXT[["tick"]], angle = 55,
                                        vjust = 1, hjust = 1),
            legend.title = element_text(size = TXT[["legend"]], face = "bold"),
            legend.text  = element_text(size = TXT[["tick"]]),
            panel.grid   = element_blank())
  }, height = function() { req(rv$res); max(300, 30 * rv$res$n + 120) })

  output$prob_tab <- renderTable({
    req(rv$res)
    P <- rv$res$fit$rank_prob
    cbind(" " = rank_rows(P), as.data.frame(round(P, 3)))
  }, striped = FALSE, spacing = "xs", digits = 3, align = "r")

  output$dl_prob <- downloadHandler(
    filename = function() "rank_probabilities.csv",
    content = function(f) write.csv(rv$res$fit$rank_prob, f))

  # ---------------- posteriors (Fig 8) + summary
  # paper Fig. 8 code (05a_RWD_2020.R); text sizes from the size system
  dens_plot <- function(M, labels, xlab, tag, lwd) {
    d <- data.frame(value = as.vector(M),
                    who   = factor(rep(labels, each = nrow(M)), levels = labels))
    ggplot(d, aes(x = value, fill = who, linetype = who, colour = who)) +
      geom_density(alpha = 0.4, linewidth = lwd) +
      scale_linetype_manual(values = rep_len(PAPER_LINETYPES, length(labels))) +
      guides(fill     = guide_legend(ncol = dens_cols(length(labels))),
             colour   = guide_legend(ncol = dens_cols(length(labels))),
             linetype = guide_legend(ncol = dens_cols(length(labels)))) +
      theme_minimal() +
      labs(x = xlab, y = "Density", tag = tag) +
      theme(axis.title.x = element_text(size = TXT[["label"]], face = "bold"),
            axis.title.y = element_text(size = TXT[["label"]], face = "bold"),
            axis.text.x  = element_text(size = TXT[["tick"]]),
            axis.text.y  = element_text(size = TXT[["tick"]]),
            legend.text  = element_text(size = TXT[["legend"]]),
            legend.title = element_blank(),
            plot.tag     = element_text(size = TXT[["head"]], face = "bold"))
  }
  dens_cols <- function(n) ceiling(n / 12)                # legend columns
  dens_height <- function() {             # the paper's panel aspect, legend fit
    req(rv$res)
    w <- session$clientData$output_dens_mu_width %||% FIG$dens_w
    n <- rv$res$n
    max(220, 1.4 * TXT[["legend"]] * ceiling(n / dens_cols(n)) + 80, w / 2.6)
  }
  output$dens_mu <- renderPlot({
    req(rv$res)
    dens_plot(rv$res$fit$draws$mu, rv$res$fit$ids, expression(~ mu[i]), "A", 1.2)
  }, height = dens_height)
  output$dens_beta <- renderPlot({
    req(rv$res)
    dens_plot(rv$res$fit$draws$beta, rv$res$fit$ids, expression(~ beta[j]), "B", 1)
  }, height = dens_height)

  output$summ_tab <- renderTable({
    req(rv$res)
    s <- rv$res$fit$summary
    data.frame(
      ID = s$id, `Final rank` = s$final_rank,
      `mu (mean)` = sprintf("%.2f", s$mu_mean),
      `mu 95% CI` = sprintf("[%.2f, %.2f]", s$mu_lo, s$mu_hi),
      `beta (mean)` = sprintf("%.2f", s$beta_mean),
      `beta 95% CI` = sprintf("[%.2f, %.2f]", s$beta_lo, s$beta_hi),
      `sigma2 (median)` = sprintf("%.2f", s$sig2_med),
      check.names = FALSE)
  }, striped = FALSE, spacing = "xs")

  output$dl_summ <- downloadHandler(
    filename = function() "bayespeer_posterior_summary.csv",
    content = function(f) write.csv(rv$res$fit$summary, f, row.names = FALSE))

  # ---------------- Fig 9: bias vs performance
  pairs_df <- reactive({
    req(rv$res)
    ext <- if (!is.null(rv$ext) && length(rv$ext) == rv$res$n) rv$ext else NULL
    analysis_frame(rv$R, rv$res$fit, external = ext)
  })

  output$pairs <- renderPlot({
    df <- pairs_df()
    X <- df[, -1, drop = FALSE]
    # paper Fig. 9 code (05C_RWD_2018.R): plotmath column labels, GGally
    # ggpairs with correlations above the diagonal; text sizes from TXT
    # (each label on two lines so it fits the strips at the app's figure size)
    two_lines <- function(x) {             # plotmath for a free-text label
      x <- gsub("'", "", trimws(x)); sp <- gregexpr(" ", x)[[1]]
      if (nchar(x) <= 12 || sp[1] < 0) return(sprintf("'%s'", x))
      k <- sp[which.min(abs(sp - nchar(x) / 2))]
      sprintf("atop('%s', '%s')", substr(x, 1, k - 1), substr(x, k + 1, nchar(x)))
    }
    lab <- c(External       = two_lines(rv$ext_label),
             Performance_mu = "atop('Performance', (hat(mu)[i]))",
             SelfBias_beta  = "atop('Self-bias', (hat(beta)[j]))",
             RankDifference = "atop('Rank', 'Difference')",
             Variability_s2 = "atop('Variance', (hat(sigma)[j]^2))")
    labs_k <- unname(lab[colnames(X)])
    colnames(X) <- paste0("v", seq_len(ncol(X)))
    ggpairs(X, columnLabels = labs_k, labeller = label_parsed,
            upper = list(continuous = wrap("cor", size = txt_mm("emph")))) +
      theme_bw() +
      theme(strip.text = element_text(size = TXT[["note"]]),
            axis.text  = element_text(size = TXT[["tick"]]))
  }, height = function() {
    # the paper's Fig. 9 aspect (11 x 7.5), but each row at least tall enough
    # for its rotated strip label on the right
    k <- ncol(pairs_df()) - 1
    max((session$clientData$output_pairs_width %||% FIG$pairs_w) / 1.47,
        k * 6.5 * TXT[["note"]] + 75)
  })

  output$cor_tab <- renderTable({
    df <- pairs_df()
    C <- cor(df[, -1, drop = FALSE])
    cbind(" " = rownames(C), as.data.frame(round(C, 3)))
  }, striped = FALSE, spacing = "xs", digits = 3)

  output$interp <- renderUI({
    df <- pairs_df(); req(nrow(df) > 3)
    p(class = "small-note",
      sprintf("Pearson correlations. Marks after a value: p ≤ 0.1 (.), 0.05 (*), 0.01 (**), 0.001 (***). With n = %d, these correlations have large sampling variability.%s",
              rv$res$n, if (ext_is_demo()) " The exam scores in this analysis are made up." else ""))
  })

  # ---------------- about
  output$about <- renderUI({
    tagList(
      h4("What this app does"),
      p("BayeSPeer models each observed rank rᵢⱼ through a latent importance score ωᵢⱼ = μᵢ + βⱼ·𝕀(i = j) + εᵢⱼ, where μᵢ is true performance, βⱼ is ranker j's self-evaluation bias (applied only to the diagonal self-rank), and εᵢⱼ ~ N(0, σ²ⱼ) captures ranker-specific noise. It is fitted with the paper's two-step empirical-Bayes Gibbs sampler: step 1 uses vague priors (μ_β ~ U(−10, 10) and, for the variances, the prior chosen under Settings, by default IG(0.1, 0.1) as in the paper), step 2 uses the same prior family with its parameters set from the step-1 draws, and each step runs the chosen number of iterations with the first half as burn-in. Results come from the step-2 draws after burn-in: the aggregated ranking (by posterior mean of μ), rank probability profiles, credible intervals, and bias estimates."),
      h4("Methods compared (as in the paper)"),
      tags$ul(
        tags$li(strong("BayeSPeer"), " — this paper's bias-aware Bayesian hierarchical model."),
        tags$li(strong("Bayesian:"), " BIRRA (Badgeley et al., 2015)."),
        tags$li(strong("Distribution-based:"), " RRA (Kolde et al., 2012), Stuart (Stuart et al., 2003), min, from the RobustRankAggreg package."),
        tags$li(strong("Markov chain:"), " MC1, MC2, MC3 (Dwork et al., 2001), from the TopKLists package."),
        tags$li(strong("Optimization:"), " CEMC.s, CEMC.k: Cross-Entropy Monte Carlo with Spearman / Kendall distance (Lin & Ding, 2009), from the TopKLists package."),
        tags$li(strong("Summary statistics:"), " arithmetic mean, geometric mean, L2 norm, median of ranks (Borda variants; Lin, 2010), from the TopKLists package.")),
      h4("Implementation notes"),
      tags$ul(
        tags$li("Complete rankings without ties are required (as assumed in the paper)."),
        tags$li("The BayeSPeer sampler follows BayeSPeer.pre() and bayespeer_2step() in the paper's code repository step by step, so the same seed gives the same draws."),
        tags$li("The baselines are called as in the paper's analysis scripts. The .selfrm variants set the self-ranks to missing and pass the shortened lists to the same functions."),
        tags$li("CEMC is a stochastic search, so its ordering of nearly tied individuals can change with the random seed.")),
      h4("Reference"),
      p("Ruan, H., Wang, K., Lim, J., and Wang, X. “Bias Meets Bayes: Bayesian Aggregation for Peer and Self Ranking.”"))
  })
}

shinyApp(ui, server)
