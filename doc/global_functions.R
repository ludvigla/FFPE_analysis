library(wholebrain)
data("atlasIndex")
data("EPSatlas")
data("ontology")

get_allen_brain_annotations <- function (
  regi,
  xy,
  image.size.micron = 8705,
  image.size.pixel = 2000,
  spot.radius = 27.5
) {
  # Define pixel size and radius
  pixels.per.um <- (image.size.pixel/image.size.micron)
  pixelsize <- 1/pixels.per.um # micrometer
  # make polygon dataframe
  x <- xy$x; y <- xy$y
  spots.contours <- do.call("rbind", lapply(seq_along(x), function(i){circle.perimeter(x[i], y[i], i, spot.radius, pixelsize)}))
  # create a soma list object to contain the spot info in
  spotSoma <- list(x = x,
       y = y,
       intensity <- rep(100, length(x)),
       area = rep(pi*(spot.radius^2)*pixelsize, length(x)),
       contour.x = spots.contours$x,
       contour.y = spots.contours$y,
       contour.ID = spots.contours$id)
  # create segmentation output list object 
  seg.Spots <- list(filter = NULL, soma = spotSoma)
  # get identity of features in atlas and in both coordinate systems
  datasetSpots <- wholebrain::get.cell.ids(regi, seg.Spots, forward.warp = TRUE)[, c("x", "y", "id", "color", "acronym", "name")]
  rownames(datasetSpots) <- rownames(xy)
  return(datasetSpots)
}
  
# function to get polygon for spot
circle.perimeter <- function ( 
  x, y, 
  id, 
  spot.radius, 
  pixelsize
){
  points.per.circle <- 40
  theta = seq(0, 2*pi,length = points.per.circle)
  x = x + (spot.radius*cos(theta))/pixelsize
  y = y + (spot.radius*sin(theta))/pixelsize
  return(data.frame(x, y, id))
}


find_parents <- function (
  children_acronyms,
  target_acronyms,
  ann
) {
  stopifnot(children_acronyms %in% ann$acronym, target_acronyms %in% ann$acronym)
  conv_acr_to_id <- setNames(ann$id, ann$acronym)
  conv_id_to_acr <- setNames(as.character(ann$acronym), ann$id)
  children_ids <- conv_acr_to_id[children_acronyms]
  target_ids <- conv_acr_to_id[target_acronyms]
  rownames(ann) <- ann$id
  # Build hierarchy
  tree <- setNames(lapply(ann$id, function(id) {
    parents <- id
    cur.id <- id
    while (!997 == cur.id) {
      cur.parent <- ann[paste0(cur.id), "parent"]
      parents <- c(parents, cur.parent)
      cur.id <- cur.parent
    }
    return(parents)
  }), nm = ann$id)
  
  # Check 
  selected_branches <- tree[paste0(children_ids)]
  checks <- setNames(data.frame(do.call(rbind, lapply(selected_branches, function(br) {
    target_ids %in% br
  }))), target_acronyms)
  checks <- checks[!rowSums(checks) > 2, ]
  results <- data.frame(acronym = unlist(apply(checks, 1, function(x) {
    res <- target_acronyms[x]
    if (length(res) == 0) {
      res <- NA
    }
    if (length(res) == 2) {
      if (res[2] %in% tree[[conv_acr_to_id[res[1]]]]) {
        res <- res[1]
      } else {
        res <- res[2]
      }
    }
    return(res)
  })), row.names = rownames(checks))
  results <- na.omit(results)
  final_res <- data.frame(children_acronym = conv_id_to_acr[rownames(results)], 
                          target_acronym = results$acronym, 
                          row.names = conv_id_to_acr[rownames(results)])
  return(final_res)
}


get.region.coordinates <- function (
  acronym, 
  registration
) {
  coordinate <- registration$coordinate
  ont <- ontology
  id_to_acr <- setNames(as.character(ont$acronym), nm = as.character(ont$id))
  ont$grandparent <- wholebrain::id.from.acronym(wholebrain::get.acronym.parent(id_to_acr[ont$parent]))
  ont$ggparent <- wholebrain::id.from.acronym(wholebrain::get.acronym.parent(id_to_acr[ont$grandparent]))
  ont$gggparent <- wholebrain::id.from.acronym(wholebrain::get.acronym.parent(id_to_acr[ont$ggparent]))
  k <- which(abs(coordinate - atlasIndex$mm.from.bregma) == min(abs(coordinate - atlasIndex$mm.from.bregma)))
  plate.info <- EPSatlas$plate.info[[k]]
  get.outline <- function(acronym, registration) {
    id <- id.from.acronym(acronym)
    index1 <- index2 <- index3 <- index4 <- which(plate.info$structure_id == id)
    #if (length(index1) == 0) {
      if (length(index1) == 0) {
        id <- ont$id[which(ont$parent %in% id)]
        index1 <- which(plate.info$structure_id %in% id)
      }
      if (length(index2) == 0) {
        id <- ont$id[which(ont$parent %in% id)]
        index2 <- which(plate.info$structure_id %in% id)
      }
      if (length(index3) == 0) {
        id <- ont$id[which(ont$parent %in% id)]
        index3 <- which(plate.info$structure_id %in% id)
      }
      if (length(index4) == 0) {
        id <- ont$id[which(ont$parent %in% id)]
        index4 <- which(plate.info$structure_id %in% id)
      }
      region <- registration$atlas$outlines[c(index1, index2, index3, index4)]
      region <- lapply(1:length(region), function(x) {
        data.frame(xT = c(region[[x]]$xlT, region[[x]]$xrT), 
                   yT = c(region[[x]]$ylT, region[[x]]$yrT), 
                   x = c(region[[x]]$xl, region[[x]]$xr), 
                   y = c(region[[x]]$yl, region[[x]]$yr), 
                   right.hemisphere = c(rep(FALSE, length(region[[x]]$xl)), 
                                        rep(TRUE, length(region[[x]]$xr))), name = acronym.from.id(id[x]))
      })
      region <- do.call("rbind", region)
    return(region)
  }
  if (length(acronym) == 1) {
    region <- get.outline(acronym, registration)
  } else {
    region <- lapply(acronym, get.outline)
    region <- do.call("rbind", region)
  }
  #if (length(region) > 1) 
  #  scale.factor <- mean(dim(registration$transformationgrid$mx)/c(registration$transformationgrid$height, 
  #                                                                 registration$transformationgrid$width))
  #region[, 1:4] <- region[, 1:4] * (1/scale.factor)
  region$query <- acronym
  return(region)
}

.select.ids <- function(x, side=c("top", "bottom"), mid=NULL) {
  side <- match.arg(side)
  if (side == "top")  {
    x <- sort(x, decreasing=TRUE)
    if (is.null(mid))
      return(names(x))
    else 
      return(names(x)[which(x > mid)])
  } else if (side == "bottom") {
    x <- sort(x, decreasing=FALSE)
    if (is.null(mid)) 
      return(names(x))
    else 
      return(names(x)[which(x < mid)])
  }
}

.rbo.ext <- function(x, y, p, k, uneven.lengths = TRUE) {
  if (length(x) <= length(y)) {
    S <- x
    L <- y
  } else {
    S <- y
    L <- x
  }
  l <- min(k, length(L))
  s <- min(k, length(S))
  
  if (uneven.lengths) {
    Xd <- sapply(1:l, function(i) length(intersect(S[1:i], L[1:i])))
    ((1-p) / p) *
      ((sum(Xd[seq(1, l)] / seq(1, l) * p^seq(1, l))) +
         (sum(Xd[s] * (seq(s+1, l) - s) / (s * seq(s+1, l)) * p^seq(s+1, l)))) +
      ((Xd[l] - Xd[s]) / l + (Xd[s] / s)) * p^l  
  } else {
    #stopifnot(l == s)
    k <- min(s, k)
    Xd <- sapply(1:k, function(i) length(intersect(x[1:i], y[1:i])))
    Xk <- Xd[k]
    (Xk / k) * p^k + (((1-p)/p) * sum((Xd / seq(1,k)) * p^seq(1,k)))
  }
}

rbo <- function (
  s, 
  t, 
  p, 
  k=floor(max(length(s), length(t))/2), 
  side=c("top", "bottom"), 
  mid=NULL, 
  uneven.lengths = TRUE
) {
  side <- match.arg(side)
  if (!is.numeric(s) | !is.numeric(t))
    stop("Input vectors are not numeric.")
  if (is.null(names(s)) | is.null(names(t)))
    stop("Input vectors are not named.")
  ids <- switch(side,
                "top"=list(s=.select.ids(s, "top", mid), t=.select.ids(t, "top", mid)),
                "bottom"=list(s=.select.ids(s, "bottom", mid), t=.select.ids(t, "bottom", mid))
  )
  min(1, .rbo.ext(ids$s, ids$t, p, k, uneven.lengths = uneven.lengths))
}

enrichment_score <- function(X, Y) {
  ((Matrix::rowMeans(X) + 1) / (Matrix::rowMeans(Y) + 1)) * ((Matrix::rowMeans(X > 0) + 1) / (Matrix::rowMeans(Y > 0) + 1))
}


pairwise.scores <- function(object, label = "target_acronym", slot = "data", assay = "SCT") {
  data <- GetAssayData(object, slot = slot, assay = assay)
  grp <- object[[ ]][, label]
  lvls <- na.omit(unique(as.character(grp)))
  d <- do.call(cbind, lapply(lvls, function(lvl) {
    X <- data[, grp %in% lvl]
    Y <- data[, !grp %in% lvl]
    enrichment_score(X, Y)
  }))
  colnames(d) <- lvls
  return(d)
}
