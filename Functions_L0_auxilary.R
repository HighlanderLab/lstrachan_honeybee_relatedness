#Level 0 Auxilary Fuctions----

# Number of colonies----
nColonies <- function(colonies) {
  if (!"Colonies" %in% class(colonies)) {
    stop("colonies has to be a class Colonies")
  } 
  return(length(colonies@colonies))
}

# nWorkers----
nWorkers <- function(colony) {
  n = ifelse(!is.null(colony@workers), colony@workers@nInd, 0)
  return(n)
}


# nDrones----
nDrones <- function(colony) {
  n = ifelse(!is.null(colony@drones), colony@drones@nInd, 0)
  return(n)
}


# nVirginQueens----
nVirginQueens <- function(colony) {
  n = ifelse(!is.null(colony@virgin_queens), colony@virgin_queens@nInd, 0)
  return(n)
}

# nFathers----
nFathers <- function(colony) {
  if (is.null(colony@queen)) {
    n = 0
  } else {
    n = ifelse(!is.null(colony@queen@misc$fathers), colony@queen@misc$fathers@nInd, 0)
  }
  return(n)
}

# isQueenMated----
isQueenMated <- function(x) {
  if ("Pop" %in% class(x)) {
    return(!is.null(x@misc$fathers))
  } else if ("Colony" %in% class(x)) {
    if (!is.null(x@queen)) {
      return(!is.null(x@queen@misc$fathers))
    } else {
      return(FALSE)
    }
  }
}

# Extract the year of birth of the queen----
#' @rdname extractQueenYOB
#' @method extractQueenYOB
#' @title Extract the queen's year of birth
#' @usage \method{extractQueenYOB}(colony)
#' @description Extract the year of birth of the queen \code{colony@queen@misc$yearOfBirth} slot
#' @param colony AlphaSimR population object
#' 
#' @example 
#' #'#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' setQueenAge(colony, 1)
#' extractQueenYOB(colony)
#' 
#' @return Integer, the year of birth of the queen.
#' 
#' @export
#' 
extractQueenYOB <- function(colony) {
  return(colony@queen@misc$yearOfBirth)
}



# Compute the age of the queen----

#' @rdname computeQueenAge
#' @method computeQueenAge
#' @title COmputer the queen's age in years
#' @usage \method{computeQueenAge}(colony, year)
#' @description Compute the age of the queen from the \code{colony@queen@misc$yearOfBirth} slot
#' @param colony AlphaSimR population object
#' @param currentYear Integer, current year
#' 
#' @example 
#' #'#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' setQueenAge(colony, year = 1)
#' extractQueenYOB(colony)
#' computerQueenAge(colony, currentYear = 5)
#' 
#' @return Integer, the age of the queen.
#' 
#' @export
#' 
computeQueenAge <- function(x, currentYear) {
  if ("Pop" %in% class(x)) {
    return(currentYear - x@misc$yearOfBirth)
  } else if ("Colony" %in% class(x)) {
    if (!is.null(x@queen)) {
      return(currentYear - x@queen@misc$yearOfBirth)
    }
  }
}




# Get colony IDs from the colonies----
#' @rdname getIDs
#' @method getIDs
#' @title Get the colonies IDs from the colonies
#' @usage \method{getIDs}(colonies)
#' @description Get the colony IDs from the colonies
#' @param colonies 
#'
#' @example 
#' @return Colony IDs
#' @export
#' 
getIDs <- function(colonies, ID) {
  return(sapply(colonies@colonies, FUN = function(x) x@id))
}

