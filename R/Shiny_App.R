#' @title Launch the Interactive BioSeqAligner App
#'
#' @description
#' Starts a Shiny web application that ties the package together: paste or upload
#' sequences, run pairwise or multiple sequence alignment, and view colour-coded
#' alignment plots, dot plots, summary metrics, and a neighbour-joining tree —
#' all in one place, with downloadable results.
#'
#' @details Requires the \pkg{shiny} package (and \pkg{ape} for the tree tab).
#'   Install with \code{install.packages(c("shiny", "ape"))}.
#'
#' @return Called for its side effect of launching the app. Does not return.
#'
#' @examples
#' \dontrun{
#' launch_aligner()
#' }
#'
#' @export
launch_aligner <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install it with install.packages('shiny').")
  }

  ui <- shiny::fluidPage(
    shiny::titlePanel("BioSeqAligner"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::radioButtons("mode", "Mode",
          c("Pairwise" = "pair", "Multiple (MSA)" = "msa")
        ),
        shiny::conditionalPanel(
          "input.mode == 'pair'",
          shiny::textInput("seq1", "Sequence 1", "ACGTGGATCGA"),
          shiny::textInput("seq2", "Sequence 2", "ACGTGCATCGA"),
          shiny::selectInput("method", "Method",
            c("Global (Needleman-Wunsch)" = "global",
              "Local (Smith-Waterman)" = "local")
          )
        ),
        shiny::conditionalPanel(
          "input.mode == 'msa'",
          shiny::textAreaInput("fasta", "FASTA input",
            ">s1\nACGTGGAA\n>s2\nACGTGCAA\n>s3\nACGTGGATA",
            rows = 8
          )
        ),
        shiny::selectInput("alphabet", "Scoring",
          c("Nucleotide" = "nucleotide", "Protein (BLOSUM62)" = "BLOSUM62")
        ),
        shiny::numericInput("gap_open", "Gap open", 10, min = 0),
        shiny::numericInput("gap_extend", "Gap extend", 1, min = 0),
        shiny::actionButton("run", "Align", class = "btn-primary")
      ),
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Alignment", shiny::verbatimTextOutput("aln_text")),
          shiny::tabPanel("Plot", shiny::plotOutput("aln_plot", height = "400px")),
          shiny::tabPanel("Metrics", shiny::tableOutput("metrics")),
          shiny::tabPanel("Tree", shiny::plotOutput("tree", height = "400px"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    submat <- shiny::reactive(scoring_matrix(input$alphabet))

    result <- shiny::eventReactive(input$run, {
      if (input$mode == "pair") {
        list(type = "pair", obj = align(input$seq1, input$seq2,
          method = input$method, submat = submat(),
          gap_open = input$gap_open, gap_extend = input$gap_extend
        ))
      } else {
        seqs <- read_fasta(strsplit(input$fasta, "\n")[[1]])
        list(type = "msa", obj = msa_align(seqs,
          submat = submat(),
          gap_open = input$gap_open, gap_extend = input$gap_extend
        ))
      }
    })

    output$aln_text <- shiny::renderPrint({
      shiny::req(result())
      print(result()$obj)
    })

    output$aln_plot <- shiny::renderPlot({
      shiny::req(result())
      plot(result()$obj)
    })

    output$metrics <- shiny::renderTable({
      shiny::req(result())
      if (result()$type == "pair") {
        alignment_stats(result()$obj)
      } else {
        data.frame(
          column = seq_along(conservation_scores(result()$obj)),
          conservation = round(conservation_scores(result()$obj), 3)
        )
      }
    })

    output$tree <- shiny::renderPlot({
      shiny::req(result())
      if (result()$type == "msa" && requireNamespace("ape", quietly = TRUE)) {
        nj_tree(result()$obj)
      } else {
        plot.new()
        text(0.5, 0.5, "Tree available for MSA mode (needs the 'ape' package).")
      }
    })
  }

  shiny::shinyApp(ui, server)
}
