#Level 3 Colonies Functions----


#Create  Colonies ----
#' @title Create Colonies
#'
#' @description
#' Creates a new \code{\link{Colonies-class}} from one or more
#' \code{\link{Colony-class}} and/or \code{\link{Colonies-class}}
#' objects.
#'
#' @param ... one or more \code{\link{Colony-class}} and/or
#' \code{\link{Colonies-class}} objects.
#'
#' @return Returns an empty object of \code{\link{Colonies-class}}
#'
#' @examples
#' @export

createColonies = function(..., n = NULL){
  if (is.null(n)) {
    input = list(...)
    class = sapply(input, "class")
    stopifnot(all(class=="Colony" | class=="Colonies") | all(class=="NULL"))
    output = new("Colonies", colonies=input)    
  } else {
    output = new("Colonies", colonies=vector(mode = "list", length = n))
  }
  
  return(output)
}

#Add colony to Colonies----
#' @title Add colony to the Colonies
#'
#' @description
#'
#' @param
#'
#' @return
#'
#' @examples
#' @export

addColonyToTheColonies= function(colony, Colonies){
  if (class(colony) != "Colony") {
    message("The colony parameter is not a Colony object.")
  }
  Colonies@colonies = append(Colonies@colonies, list(colony))
  return(Colonies)
}

# Select colonies ----

#' @rdname selectColonies
#' @method selectColonies
#' @title Select the colonies from the colony list based on IDs.
#' @usage \method{selectColonies}(colony_list, colony_ids)
#' @description Select the colonies from the list of all colonies based
#' on colony IDs and return a list of selected colonies.
#' @param colony_list 
#' @param colony_ids
#'
#' @example inst/examples/examples_selectColonies.R
#' @return A list of selected colonies.
#' @export
#' 
selectColonies <- function(colonies, ID = NULL, p = NULL) {
  if (!is.null(ID)) {
    ret = colonies[sapply(colonies@colonies, FUN = function(x) x@id %in% ID)]
  } else if (!is.null(p)) {
    lPull <- as.logical(rbinom(n = nColonies(colonies), size = 2, p = p))
    if (any(lPull)) {
      ret <- colonies[lPull]
    } else {
      ret = NULL
    }
  } else {
    stop("Provide either ID or p!")
  }
  return(ret)
}



# Pull colonies ----

#' @rdname pullColonies
#' @method pullColonies
#' @title Pull the colonies from the colony list based on IDs.
#' @usage \method{pullColonies}(colonies, colonyIds)
#' @description Pull the colonies from the list of all colonies based
#' on colony IDs and return two lists: a list of selected colonies and 
#' updated original colonies
#' @param colonies 
#' @param colonyIDs
#'
#' @example 
#' @return Two lists: a list of selected colonies and an updated inpute colonies
#' @export
#' 
pullColonies <- function(colonies, ID = NULL, p = NULL) {
  
  if (!is.null(ID)) {
    pulledColonies <- selectColonies(colonies, ID)
    remainingColonies <- removeColonies(colonies, ID)
  } else if (!is.null(p)) {
    lPull <- as.logical(rbinom(n = nColonies(colonies), size = 2, p = p))
    message(paste0("Pulling out ", sum(lPull), " colonies."))
    if (any(lPull)) {
      ids = getIDs(colonies)
      pulledColonies <- selectColonies(colonies, ids[lPull])
      remainingColonies <- removeColonies(colonies, ids[lPull])
    } else {
      pulledColonies = createColonies()
      remainingColonies = colonies
    }
  } else {
    stop("Provide either ID or p!")
  }
  return(list(pulledColonies = pulledColonies, remainingColonies = remainingColonies))
}


# Remove colonies ----

#' @rdname removeColonies
#' @method removeColonies
#' @title Remove the colonies from the colony list based on IDs.
#' @usage \method{removeColonies}(colony_list, colony_ids)
#' @description Remove the colonies from the list of all colonies based
#' on colony IDs and return a list of remaining colonies.
#' @param colony_list 
#' @param colony_ids
#'
#' @example 
#' @return A list of remaining colonies.
#' @export
#' 
removeColonies <- function(colonies, ID) {
  return(colonies[!sapply(colonies@colonies, FUN = function(x) x@id %in% ID)])
}


# Create multiple Virgin colonies ----

createMultipleVirginColonies = function(founderPop, nColonies) {
  ret = createColonies(n = nColonies)
  virginQueens = selectInd(founderPop, nInd = nColonies, use = "rand")
  for (colony in 1:nColonies) {
    ret@colonies[[colony]] = createColony(virgin_queens = virginQueens[colony])
  }
  return(ret)
}


# Create multiple mated colonies ----

createMultipleMatedColonies = function(founderPop, nColonies, nAvgFathers) {
  ret = createColonies(n = nColonies)
  queensID = sample(founderPop@id, size = nColonies, replace = FALSE)
  queenMatch = founderPop@id[founderPop@id %in% queensID]
  queens = founderPop[queenMatch]
  DPQMatch = founderPop@id[!founderPop@id %in% queensID]
  DPQs = founderPop[DPQMatch]
  DCA =  createFounderDrones(queenPop = DPQs, nDronesPerQueen = 10)
  fatherPackages = pullDronePackagesFromDCA(DCA, n = nColonies,nAvgFathers = nAvgFathers)
  for (colony in 1:nColonies) {
    ret@colonies[[colony]] = createColony(queen = queens[colony], 
                                          fathers = fatherPackages[[colony]])
  }
  
  return(ret)
}


# 
# # Build up colonies (add workers and drones)----
# 
buildUpColonies = function(colonies, nWorkers, nDrones) {
  nCol = nColonies(colonies)
  for (colony in 1:nCol) {
    colonies@colonies[[colony]] = buildUpColony(colony = colonies[[colony]], 
                                                nWorkers = nWorkers,
                                                nDrones = nDrones) 
  }
  return(colonies)
}

# 
# # reQueen Colonies----
# 
reQueenColonies = function(colonies, queens) {
  nCol = nColonies(colonies)
  if (nInd(queens) < nCol) {
    stop("Not enough queens!")
  }
  for (colony in 1:nCol) {
    colonies@colonies[[colony]] = reQueenColony(colony = colonies[[colony]], 
                                                queen = queens[colony]) 
  }
  return(colonies)
}




# 
# # Collapse the colonies----
# 
# collapseColonies <- function(colonies, ID) {
#   return(removeColonies(colonies, ID))
# }




# Supersede the colonies----

supersedeColonies <- function(colonies
                              #, crossVirginQueen = FALSE, fathers = NULL, pWorkers = 1, pDrones = 1
) {
  nCol = nColonies(colonies)
  if (nCol == 0) {
    return(createColonies())
  }
  for (colony in 1:nCol) {
    colonies@colonies[[colony]] = supersedeColony(colonies[[colony]])
  }
  return(colonies)
}




# Swarm the colonies----

swarmColonies <- function(colonies
                          #, crossVirginQueen = FALSE, fathers = NULL, pWorkers = 1, pDrones = 1
) {
  nCol = nColonies(colonies)
  if (nCol == 0) {
    return(list(swarms = createColonies(n = 0), remnants = createColonies(n = 0)))
  }
  # if (length(fathers) < nCol) {
  #   stop("Not enought fathers)
  # }
  ret = list(swarms = createColonies(n = nCol), remnants = createColonies(n = nCol))
  
  for (colony in 1:nCol) {
    tmp = swarmColony(colonies[[colony]])
    ret$swarms@colonies[[colony]] = tmp$swarm
    ret$remnants@colonies[[colony]] = tmp$remnant
  }
  return(ret)
  
}


# Split the colonies----

splitColonies <- function(colonies
                          #, crossVirginQueen = FALSE, fathers = NULL, pWorkers = 1, pDrones = 1
) {
  nCol = nColonies(colonies)
  if (nCol == 0) {
    return(list(splits = createColonies(n = 0), remnants = createColonies(n = 0)))
  }
  # if (length(fathers) < nCol) {
  #   stop("Not enought fathers)
  # }
  ret = list(splits = createColonies(n = nCol), remnants = createColonies(n = nCol))
  
  for (colony in 1:nCol) {
    tmp = splitColony(colonies[[colony]])
    ret$splits@colonies[[colony]] = tmp$split
    ret$remnants@colonies[[colony]] = tmp$remnant
  }
  return(ret)
  
}

# Cross the colonies----

crossColonies <- function(colonies, DCA, nAvgFathers
                          #, crossVirginQueen = FALSE, fathers = NULL, pWorkers = 1, pDrones = 1
) {
  nCol = nColonies(colonies)
  # if (length(fathers) < nCol) {
  #   stop("Not enought fathers)
  # }
  ret = createColonies(n = nCol)
  nFathers = rpois(n = nCol, lambda = nAvgFathers)
  fatherPackages = pullDronePackagesFromDCA(DCA, n = nCol, nAvgFathers = nAvgFathers)
  
  for (colony in 1:nCol) {
    ret@colonies[[colony]] <- crossColony(colonies[[colony]],
                                          fathers = fatherPackages[[colony]]
    )
  }
  return(ret)
  
}

# Set Pheno colonies ----
setPhenoColonies <- function(colonies, FUN = NULL, ...) {
  nCol = nColonies(colonies)
  for (colony in 1:nCol) {
    colonies@colonies[[colony]] <- setPhenoColony(colonies[[colony]],
                                                  FUN = FUN, ...)
  }
  return(colonies)
}