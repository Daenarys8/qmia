#' Generate a Quarto-based report for microbiome analysis
#'
#' Generates a complete Quarto site for microbiome analysis based on input parameters.
#' The function saves parameters in an `.Rds` file and renders the Quarto site.
#'
#' @param tse.path \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}} 
#' object containing taxonomic abundance data.
#'
#' @param indices \code{Character scalar}. The index or parameter used for the analysis. 
#' (Default: c("shannon"))
#'
#' @param group.var \code{Character scalar}. Name of the grouping variable in \code{colData(tse)}. 
#' (Default: "diet")
#'
#' @param id.var \code{Character scalar}. Name of the ID column in \code{colData(tse)}. 
#' (Default: "id")
#'
#' @param time.var \code{Character scalar}. Name of the column in \code{colData(tse)} that contains timepoints. 
#' (Default: "timepoint")
#'
#' @param tp1 \code{Character scalar}. The baseline timepoint to be used for analysis. 
#' (Default: "before")
#'
#' @param tp2 \code{Character scalar}. The follow-up timepoint to be used for analysis. 
#' (Default: "after")
#'
#' @param comp.group \code{Character vector}. Comparisons to be made before the timepoint.
#' (Default: \code{c()})
#'
#' @param comp.paired \code{Character vector}. Paired comparisons between timepoints.
#' (Default: \code{c()})
#'
#' @param lmm_formula \code{Character scalar}. The formula to be used for linear mixed model analysis.
#' (Default: "y ~ diet * timepoint + (1 | id)")
#'
#' @param adjust.method \code{Character scalar}. The adjustment method to use for multiple comparisons
#' (Default: "fdr")
#'
#' @param other_fields \code{Character scalar}. Other fields to be included in the analysis.
#' (Default: "default")
#' @param analysis \code{Character vector}. Types of analyses to perform. 
#' Options: c("alpha", "beta", "ratio", "daa"). (Default: "alpha")
#'
#' @return None. A Quarto site is generated and rendered.
#' @name generateMia
#'
#' @examples
#' \dontrun{
#' generateMia(
#'   tse = my_tse_object, 
#'   index = "default", 
#'   group.var = "diet", 
#'   id.var = "id", 
#'   time.var = "timepoint", 
#'   tp1 = "before", 
#'   tp2 = "after", 
#'   comp.group = c("Group1", "Group2"),
#'   comp.paired = c("Group1", "Group2"),
#'   lmm_formula = "y ~ diet * timepoint + (1 | id)", 
#'   adjust.method = "fdr", 
#'   other_fields = "default"
#' )
#' }
NULL

#' @importFrom quarto quarto_render
#' @importFrom readr write_rds
#' @importFrom dplyr bind_rows
#' @importFrom fs file_exists
#' @export
setGeneric("generateMia", function(
        tse.path, 
        indices = c("shannon"), 
        group.var = "diet", 
        id.var = "id", 
        time.var = "timepoint", 
        tp1 = "before", 
        tp2 = "after", 
        comp.group = NULL, 
        comp.paired = NULL, 
        lmm_formula = "y ~ diet * timepoint + (1 | id)", 
        adjust.method = "fdr", 
        other_fields = "default",
        analysis = "alpha"
) {
    standardGeneric("generateMia")
})

#' @rdname generateMia
#' @export
setMethod("generateMia", signature = "ANY", function(
        tse.path, 
        indices = c("shannon"), 
        group.var = "diet", 
        id.var = "id", 
        time.var = "timepoint", 
        tp1 = "before", 
        tp2 = "after", 
        comp.group = NULL, 
        comp.paired = NULL, 
        lmm_formula = "y ~ diet * timepoint + (1 | id)", 
        adjust.method = "fdr", 
        other_fields = "default",
        analysis = "alpha"
) {
    
    valid_analyses <- c("alpha", "beta", "ratio", "daa")
    
    # Ensure all requested analyses are valid
    invalid <- setdiff(analysis, valid_analyses)
    if (length(invalid) > 0) {
        stop("Invalid analysis type(s): ", paste(invalid, collapse = ", "))
    }
    
    # Render the Quarto site
    for (a in analysis) {
    tryCatch({
        if (a == "alpha") {
        lapply(indices, function(index) {
            orig_dir <- dirname("alpha/alpha.qmd")
            temp_qmd <- file.path(orig_dir, paste0("alpha_", index, ".qmd"))
            file.copy("alpha/alpha.qmd", temp_qmd)
            
            quarto::quarto_render(
                input = temp_qmd,
                execute_params = list(
                    tse.path = tse.path,
                    index = index,
                    group.var = group.var,
                    id.var = id.var,
                    time.var = time.var,
                    tp1 = tp1,
                    tp2 = tp2,
                    comp.group = comp.group, 
                    comp.paired = comp.paired, 
                    lmm_formula = lmm_formula,
                    adjust.method = adjust.method,
                    other_fields = other_fields
                )
            )
            file.remove(temp_qmd)
        })
        } else if (a == "ratio") {
            quarto::quarto_render(
                input = "ratio/ratio.qmd",
                execute_params = list(
                    tse.path = tse.path,
                    group.var = group.var,
                    id.var = id.var,
                    time.var = time.var,
                    tp1 = tp1,
                    tp2 = tp2,
                    comp.group = comp.group,
                    comp.paired = comp.paired,
                    lmm_formula = lmm_formula,
                    adjust.method = adjust.method,
                    other_fields = other_fields
                )
            )
            
        } else if (a == "beta") {
            quarto::quarto_render(
                input = "beta/beta.qmd",
                execute_params = list(
                    tse.path = tse.path,
                    group.var = group.var,
                    id.var = id.var,
                )
            )
            
        } else if (a == "daa") {
            quarto::quarto_render(
                input = "daa/daa.qmd",
                execute_params = list(
                    tse.path = tse.path,
                    group.var = group.var,
                )
            )
            
        } else {
            stop("Unknown analysis type: ", a)
        } 
        message("Successfully rendered")
    }, error = function(e) {
        stop("Render failed: ", e$message)
    })
    }
})
