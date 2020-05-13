#' @title Hello guys
#'
#' @description This should be straightforward, no ?
#'
#' @param x The first letter of your family name
#'
#' @examples
#' hello_chris("z")
#' @export
#' @import testthat

hello_chris <- function(x = "w") {

  if(x == "z") {
    print("Hello Christian !")
  } else if (x == "l") {
    print("Hello Chris !")
  } else if (x == "r"){
    print("Hi Mr Rolf, nice to meet you !")
  } else {
    print("Are you a Chris ?")
  }

}
