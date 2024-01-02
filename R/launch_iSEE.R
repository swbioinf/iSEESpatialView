

categorical_color_fun <- function(n){
  if (n <= 37) {
    # Less than 37 colours, use something from colour brewer
    # (joining a bunch of palettes, best colours up front)
    multiset <-  c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"),
                   RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(8, "Dark2"))
    return(multiset[1:n])
  }
  else{
    # More that 37, well at least it looks pretty
    return(rainbow(n))
  }
}

get_ecm <- function() {
  ecm <- ExperimentColorMap(

    # The default is viridis::viridis
    # https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html#the-color-scales
    # Setting continous is entirely a matter of taste. I find magma easier to read than viridis.
    all_continuous = list(
      assays  = viridis::magma,
      colData = viridis::magma,
      rowData = viridis::magma

    ),
    all_discrete = list(
      colData = categorical_color_fun,
      rowData = categorical_color_fun
    )
  )
  return(ecm)
}


get_initial_plots <- function(plot_mode) {

  initial_plots <- NULL # default is nULL

  if (plot_mode == "spatial") {

    # These options are all sce-contents agnostic.
    initial_plots <- c(
      # Show umap with clusters by default
      ReducedDimensionPlot(
        DataBoxOpen=TRUE,
        ColorBy="Column data",
        VisualBoxOpen=TRUE,
        PanelWidth=6L),
      # Show gene expression plot separated (and coloured) by cluster, by default.
      FeatureAssayPlot(XAxis = "Column data",
                       DataBoxOpen=TRUE,
                       VisualBoxOpen=TRUE,
                       ColorBy="Column data",
                       PanelWidth=6L
      ),
      # Gene list is better wide
      RowDataTable(PanelWidth=12L),

      # For cell level observations (QC.)
      ColumnDataPlot(PanelWidth=6L,
                     DataBoxOpen=TRUE,
                     VisualBoxOpen=TRUE )

      # Gene levels observations (RowDataPlot)are sometimes handy in QC steps, (total gene expression) but not so much for later browsing.
      # The defalts for this tend to plot a big blob of gene names too - slow and ugly.
      # leave it for people to add manually when wanted.
    )
  }


  return(initial_plots)

}



#' Convert Seurat to viewable SCE
#'
#' Take a processed seurat object and convert to a SCE object with some features
#' for spatial browsing in iSEE. Mainly - two extra dimensional reductions.
#'
#'    * local_spatial  : X,Y Coordinate of each cell within an FOV. Scaled to the 0-1 range for each fov.
#'    * global : Overall X,Y coorindates of each cell on the slide.
#'
#' If not subsetted to fov / slide in iSEE, these will be overplotted. (No
#' concept of different slides)
#'
#' @param sce single cell experiment object to display (Suggest generating with convert_seurat_to_viewable_sce() )
#' @param plot_mode How to setup plot panels in interface. Currently only 'spatial' or NULL (iSEE defaults).
#' @param ecm Custom colourmap used in display (handy for sample colours. See iSEE docs.)
#' @param max_cat More than this many categoricals, and data will be treated as continuous. (default=40)
#' @return None. Launches iSEE interface.
#' @examples
#'
#' \dontrun{
#' launch_iSEE_with_spatial_view(sce)
#' }
#' @export
launch_iSEE_with_spatial_view <- function(sce, plot_mode="spatial", ecm=get_ecm(), max_cat=40) {

  initial_plots=get_initial_plots(plot_mode = plot_mode)

  registerAppOptions(sce, color.maxlevels=max_cat)

  app <- iSEE(sce,
              colormap = ecm,
              initial  = initial_plots
  )


  shiny::runApp(app)

}





