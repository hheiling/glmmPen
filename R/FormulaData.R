# Convert formula and data options into y, X, Z, and group

#' Extracting Useful Vectors and Matrices from Formula and Data Information
#' 
#' Takes the model \code{formula} and an optional data frame and converts them into y, X, Z, and group output.
#' 
#' @param formula a two-sided linear formula object describing both the fixed-effects and 
#' random-effects part of the model, with the response on the left of a ~ operator and the terms, 
#' sepearated by + operators, on the right. Random-effects terms are distinguished by vertical bars 
#' ("|") separating expression for design matrices from grouping factors. (\code{formula} should be 
#' of the same format needed for \code{\link[lme4]{glmer}} in package \pkg{lme4}.) Only one grouping factor 
#' will be recognized. The random effects covariates need to be a subset of the fixed effects covariates.
#' @param data an optional data frame containing the variables named in \code{formula}. Although 
#' \code{data} is optional, the package authors \emph{strongly} recommend its use. If \code{data} is 
#' omitted, variables will be taken from the environment of \code{formula} (if specified as a formula).
#' @param na.action a function that indicates what should happen when the data contain NAs. Only the 
#' options \code{na.omit} and \code{na.pass} are recognized by this function.
#' 
#' @return a list with the following elements:
#' \item{y}{a numeric vector of the response variable}
#' \item{X}{a model matrix with the fixed covariates}
#' \item{Z}{a sparse model matrix for the random effects}
#' \item{group}{a factor vector of the grouping variable}
#' \item{cnms}{a vector of column names of the random effects}
#' \item{group_name}{character name of the group variable}
#' \item{flist}{a list of grouping factors using inf the random-effects terms}
#' \item{frame}{a model frame including all fixed and random covariates, the response, and the 
#' grouping variable}
#' 
#' @importFrom lme4 mkReTrms nobars subbars findbars
#' @export
formulaData = function(formula, data = NULL, na.action){
  # Note: need to deal with potential offsets and weights
  
  ## substitute | for +
  formula_full = subbars(formula)
  ## Find fixed effects part of formula
  formula_nobars = nobars(formula)
  
  if(!is.null(data)){
    if(class(data) != "data.frame"){
      cat(sprintf("class(data) = %s", class(data)))
      stop("data must be of class 'data.frame'")
    }else{
      # Proceed as normal
      
      # Deal with NAs
      if(identical(na.action, na.omit) | identical(na.action, "na.omit")){ 
        data = na.omit(data[,colnames(model.frame(formula_full, data = data))])
      }else if(identical(na.action, na.pass) | identical(na.action, "na.pass")){
          data = data
      }else{
        print(na.action)
        stop("na.action must be na.omit or na.pass for formulaData function")
      }
      
      ## Full data frame with fixed effects, random effects, and group var
      frame_full = model.frame(formula_full, data = data)
      ## Get fixed effects X matrix
      ## mm = model.matrix 
      mm = model.matrix(formula_nobars, data)
      
      ## Creating response variable Y
      Y = model.response(frame_full)
    }
    
  }else{ # if is.null(data)
    ## Full data frame with both fixed and random effects
    mf = model.frame(formula_full)
    
    # Check here to see if any character vars
    if(any(is.character(mf))){
      warning("If data = NULL, all variables - including the group variable - must be numeric \n")
      stop("Convert character variables to numeric factors")
    }
    
    ## Creating response variable Y
    Y = model.response(mf)
    
    frame_full = data.frame(Y=Y)
    for(j in 2:ncol(mf)){
      if(class(mf[,j]) == "matrix"){
        frame_full = data.frame(frame_full, mf[,j])
      }else{ # class(mf[,j]) == "numeric"
        # nm = colnames(mf)[j]
        frame_full = data.frame(frame_full, mf[,j])
        colnames(frame_full) = c(colnames(frame_full)[-ncol(frame_full)], names(mf)[j])
      }
    }
    
    ## Get fixed effects X matrix
    ## mm = model.matrix
    mm = model.matrix(formula_nobars)
    # if(length(unique(attr(mm, "assign"))) != length(attr(mm, "assign"))){
    #   # Find duplicates
    #   dup = attr(mm,"assign")[duplicated(attr(mm,"assign"))]
    #   trms = attr(terms(f_full), "term.labels")[dup]
    #   # Remove prefixes to colnames
    #   for(t in trms){
    #     colnames(mm) = gsub(t, "", colnames(mm))
    #   }
    # }
    
    # Deal with NAs
    if(identical(na.action, na.omit) | identical(na.action, "na.omit")){  
      na_rows = rowSums(is.na(frame_full))
      frame_full = na.omit(frame_full)
      Y = Y[-na_rows]
      mm = mm[-na_rows,]
    }else{
      warning("This function not equipted to deal with NA values. \n
              Please check that data does not contain NA values. \n", immediate. = T)
    }
    
    
  } # End if/else
  
  ## Identify random effects
  ### If no | (no random effects listed) then stop - mkBlist called by mkReTrms gives this error
  reExprs = findbars(formula)
  reTrms = mkReTrms(reExprs, frame_full)
  
  # t(Zt) from mkReTrms: columns organized by group level, then vars within group level
  Zt = reTrms$Zt
  
  # Change group condition below?
  if(length(reTrms$flist) > 1){
    stop("procedure can only handle one group")
  }else{
    group = reTrms$flist[[1]]
    group_name = names(reTrms$flist)
  }
  
  d = nlevels(group)
  numvars = nrow(Zt)/d
  Z = Matrix(0, nrow = ncol(Zt), ncol = nrow(Zt), sparse = T)
  # mkReTrms Zt rows: organized first by level of group, then vars 
  # Want Z columns organized first by vars, then levels of group within vars
  for(lev in 1:numvars){
    Z[,(d*(lev-1)+1):(d*lev)] = Matrix::t(Zt[seq(lev, by = numvars, length.out = d),])
  }
  
  # Problem if variable specified in formula is the intercept / a constant column 
  # and -1 not included in formula
  constant_cols = mm[,apply(mm, 2, var, na.rm=TRUE) == 0]
  
  if(class(constant_cols) == "matrix"){ 
    # If true, more than one column with zero variance
    # If only one column, class(constant_cols) == "numeric"
    stop("Variable(s) in formula has zero variance (constant column).
         Either remove this variable from formula or specify -1 in formula")
  } 
  
  ## Make sure colnames random effects subset of colnames mm
  cnms = reTrms$cnms[[1]] # Assume only one group
  if(sum(!(cnms %in% colnames(mm))) > 0){
    stop("random effects must be a subset of fixed effects")
  }
  
  return(list(y = Y, X = mm, Z = Z, group = group, cnms = cnms, 
              group_name = group_name, flist = reTrms$flist, frame = frame_full))
}
