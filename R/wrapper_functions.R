#' @useDynLib glmmPen
#' @export
Znew_gen2.Rfunction = function(U, Z, group, cols, n, q, d, pBigMat, J){
  Znew_gen2(U, Z, group, cols, n, q, d, pBigMat, J)
}