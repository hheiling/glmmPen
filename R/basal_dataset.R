#' @title Basal dataset: A composition of cancer datasets with top scoring pairs
#' (TSPs) as covariates and binary response indicating if the subject's cancer
#' subtype was basal-like.
#' 
#' A dataset composed of four datasets combined from studies that contain 
#' gene expression data from subjects with several types of cancer.
#' Two of these datasets contain gene expression data for subjects with 
#' Pancreatic Ductal Adenocarcinoma (PDAC), one dataset contains data for 
#' subjects with Breast Cancer, and the fourth dataset contains data for subjects 
#' with Bladder Cancer. The response of interest is whether or not the subject's
#' cancer subtype was the basal-like subtype.
#' See articles Rashid et al. (2020) "Modeling Between-Study Heterogeneity for
#' Improved Replicability in Gene Signature Selection and Clinical Prediction"
#' and Moffitt et al. (2015) "Virtual microdissection identifies distinct tumor- and 
#' stroma-specific subtypes of pancreatic ductal adenocarcinoma" 
#' for further details on these four datasets. 
#' 
#' @usage data("basal")
#' 
#' @format A list containing the following elements:
#' \describe{
#'   \item{y}{binary response vector; 1 indicates that the subject's cancer
#'   was of the basal-like subtype, 0 otherwise}
#'   \item{X}{matrix of 50 top scoring pair (TSP) covariates}
#'   \item{group}{factor indicating which cancer study the observation belongs to,
#'   which are given the following descriptions:
#'   UNC PDAC, TCGA PDAC, TCGA Bladder Cancer, and UNC Breast Cancer}
#'   \item{Z}{model matrix for random effects; organized first by variable, then
#'   by group (i.e. by cancer study)}
#' }
"basal"