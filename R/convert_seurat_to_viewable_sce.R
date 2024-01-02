


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
#' @param so Seurat object (already processed, annotated). Must be Seurat V5+!
#' @param assay Assay within Seurat object. For cosmx that's 'Nanostring'.
#'        Expected to have counts and data layers.
#' @param fov_id Column in metadata with unique fov identifiers (Must be unique across all slides)
#' @param reductions Vector of reductions in Seurat object to copy across. If none specified, will take all.
#' @return SingleCellExperiment object
#' @examples
#'
#' \dontrun{
#' sce <- convert_seurat_to_viewable_sce(so, fov_id = 'fov_name');
#' }
#' @export
convert_seurat_to_viewable_sce <- function(so, assay="Nanostring", fov_id = 'fov_name', reductions=c()) {

  counts_matrix <- Seurat::GetAssayData(so, assay = assay, layer = 'counts')
  norm_matrix   <- Seurat::GetAssayData(so, assay = assay, layer = 'data')

  # gene info
  gene_table <- so@assays[[assay]]@meta.data
  rownames(gene_table) <- rownames(so@assays[[assay]])

  # cell info
  anno_table <- so@meta.data

  # Add x/y coords

  # One image (slide) at a time
  images <- names(so@images)
  #coords_with_fov <- function(so, image) {
  #  coords <- Seurat::GetTissueCoordinates(object=so, image = image)
  #  coords$image <-
  #}
  coords_list      <- lapply(FUN=Seurat::GetTissueCoordinates, X=images, object=so)
  coords           <- dplyr::bind_rows(coords_list)
  rownames(coords) <- coords$cell

  # overall coords
  coords.global <- coords[,1:2]
  coords.global <- coords.global[colnames(so),]


  # Create coordinates in the 0..1 space
  # because I can't zoom
  # THese make sense when subset to a single fov only, else they'll be overplotted.

  fov_lookup <- stats::setNames(so@meta.data[,fov_id], nm=colnames(so))

  coords.local <- coords
  coords.local$fov_name <- fov_lookup[coords.local$cell]
  coords.local <-
    coords.local  |> dplyr::group_by_('fov_name')  |>
    dplyr::mutate(xmin = min(x),
           ymin = min(y),
           xmax = max(x),
           ymax = max(y))  |>
    dplyr::ungroup()  |>
    dplyr::mutate(
      x = (x-xmin) * 1/(xmax-xmin),
      y = (y-ymin) * 1/(ymax-ymin)
    )  |>
    dplyr::select(x,y,cell)  |>
    as.data.frame()
  rownames(coords.local) <- coords.local$cell
  coords.local <- coords.local[,1:2]



  # Build a simple SCE
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=counts_matrix,
                                                         norm=norm_matrix),
                              colData=S4Vectors::DataFrame(anno_table),
                              rowData=S4Vectors::DataFrame(gene_table)
  )

  # Jam in reductions
  # without recreating them from scratch.
  # https://rdrr.io/bioc/SingleCellExperiment/man/reducedDims.html
    # Check for specified reductions, or get all present.
  if (length(reductions)==0) {
    reductions <- names(so@reductions)
  }
  if( !all(reductions %in% names(so@reductions)) ) {
      stop(paste("Can't find all reductions in seurat object: see only ",
                 paste(names(so@reductions))))
  }
  # Add them one by one
  for (reduction in reductions) {
    SingleCellExperiment::reducedDim(sce, type=reduction ) <- so@reductions[[reduction]]@cell.embeddings
  }


  # Then add spatial coordinates as another dimension
  SingleCellExperiment::reducedDim(sce, type="local_spatial" )  <-  coords.local[colnames(sce),]
  SingleCellExperiment::reducedDim(sce, type="global_spatial" ) <-  coords.global[colnames(sce),]

  return(sce)

}
