setchar <- function(x, val, text, list = FALSE, name = NULL,
                    stop.at.error = TRUE, addtext = "",
                    return.NULL = TRUE, nchar.equal = FALSE,
                    setNA = FALSE) {
  val <- unique(val)
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  nval <- length(val)
  ##
  if (is.numeric(x)) {
    numeric.x <- TRUE
    idx <- x
    idx[idx < 1] <- NA
    idx[idx >= nval + 1] <- NA
  }
  else {
    numeric.x <- FALSE
    ##
    if (length(unique(tolower(x))) != length(unique(x)) |
        length(unique(tolower(val))) != length(unique(val)))
      idx <- charmatch(x, val, nomatch = NA)
    else
      idx <- charmatch(tolower(x), tolower(val), nomatch = NA)
  }
  ##
  if ((anyNA(idx) || any(idx == 0)) && !setNA) {
    if (list)
      first <- "List element '"
    else
      first <- "Argument '"
    ##
    if (missing(text)) {
      if (numeric.x) {
        if (nval == 1)
          vlist <- "1"
        else if (nval == 2)
          vlist <- "1 or 2"
        else
          vlist <- paste("between 1 and", nval)
      }
      else {
        if (nval == 1)
          vlist <- paste0('"', val, '"')
        else if (nval == 2)
          vlist <- paste0('"', val, '"', collapse = " or ")
        else
          vlist <- paste0(paste0('"', val[-nval], '"', collapse = ", "),
                          ', or ', '"', val[nval], '"')
      }
      ##
      if (stop.at.error)
        stop(first, name, "' must be ", vlist, addtext, ".", call. = FALSE)
      else {
        if (return.NULL)
          return(NULL)
        else
          return(x)
      }
    }
    else {
      if (stop.at.error)
        stop(first, name, "' ", text, ".", call. = FALSE)
      else {
        if (return.NULL)
          return(NULL)
        else
          return(x)
      }
    }
  }
  ##
  if (is.null(x))
    return(NULL)
  else
    res <- val[idx]
  ##
  if (nchar.equal && nchar(res) != nchar(x))
    res <- x
  ##
  res
}

setsv <- function(x) {
  if (is.null(x))
    res <- "desirable"
  else {
    res <- setchar(x, c("good", "bad"), stop.at.error = FALSE)
    ##
    if (!is.null(res))
      res <- switch(res, good = "desirable", bad = "undesirable")
    else
      res <- x
  }
  ##
  setchar(res, c("desirable", "undesirable"))
}

chkclass <- function(x, class, name = NULL) {
  ##
  ## Check class of R object
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  n.class <- length(class)
  if (n.class == 1)
    text.class <- paste0('"', class, '"')
  else if (n.class == 2)
    text.class <- paste0('"', class, '"', collapse = " or ")
  else
    text.class <- paste0(paste0('"', class[-n.class], '"', collapse = ", "),
                         ', or ', '"', class[n.class], '"')
  ##
  if (!inherits(x, class))
    stop("Argument '", name,
         "' must be an object of class \"",
         text.class, "\".", call. = FALSE)
  ##
  invisible(NULL)
}

chknumeric <- function(x, min, max, zero = FALSE, length = 0,
                       name = NULL, single = FALSE, integer = FALSE) {
  if (!missing(single) && single)
    length <- 1
  ##
  ## Check numeric variable
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  x <- x[!is.na(x)]
  if (length(x) == 0)
    return(NULL)
  ##
  if (!is.numeric(x))
    stop("Non-numeric value for argument '", name, "'.",
         call. = FALSE)
  ##
  if (length && length(x) != length)
    stop("Argument '", name, "' must be a numeric of length ", length, ".",
         call. = FALSE)
  ##
  if (!missing(min) & missing(max)) {
    if (zero & min == 0 & any(x <= min, na.rm = TRUE))
      stop("Argument '", name, "' must be positive.",
           call. = FALSE)
    else if (any(x < min, na.rm = TRUE))
      stop("Argument '", name, "' must be larger equal ",
           min, ".", call. = FALSE)
  }
  ##
  if (missing(min) & !missing(max)) {
    if (zero & max == 0 & any(x >= max, na.rm = TRUE))
      stop("Argument '", name, "' must be negative.",
           call. = FALSE)
    else if (any(x > max, na.rm = TRUE))
      stop("Argument '", name, "' must be smaller equal ",
           min, ".", call. = FALSE)
  }
  ##
  if ((!missing(min) & !missing(max)) &&
      (any(x < min, na.rm = TRUE) | any(x > max, na.rm = TRUE)))
    stop("Argument '", name, "' must be between ",
         min, " and ", max, ".", call. = FALSE)
  ##
  if (integer && any(!is_wholenumber(x))) {
    if (length(x) == 1)
      stop("Argument '", name, "' must be an integer.",
           call. = FALSE)
    else
      stop("Argument '", name, "' may only contain integers.",
           call. = FALSE)
  }
  ##
  invisible(NULL)
}

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.numeric(x))
    res <- abs(x - round(x)) < tol
  else
    res <- NA
  ##
  res
}
