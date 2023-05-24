#' Probability-model linearizing transformations
#'
#' @return A scale transformer
#' @export
#'
#' @examples
#' weibull_trans()
weibull_trans <- function() {
  scales::trans_new('weibull',
    function(x) log(-log(1 - x)),
    function(x) 1 - exp(-exp(x)))}

frechet_trans <- function() {
  scales::trans_new('weibull',
    function(x) -log(-log(x)),
    function(x) exp(-exp(-x)))}

arrhenius_trans <- function() {
  scales::trans_new(
    name = 'arrhenius',
    transform = function(x) 1 / (x + 273.15),
    inverse = function(x) 1 / x - 273.15,
    domain = c(-273.15, Inf))}

arrh <- function(temp_c) 1 / (temp_c + 273.15)
