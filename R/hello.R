#' @title Hello
#'
#' @description We don't need any description here.
#'
#' @param x The first parameter
#'
#' @return Some random text
#' @export
#'
#' @examples
#' hello("John")
#' \dontrun{
#' hello("Steve")
#' }
hello <- function(x) {
  print(paste0("Hello Mr.", x, ". You have come to the right place"))
}
