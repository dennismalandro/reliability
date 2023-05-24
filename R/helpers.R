#' Consruct fitdistcens-compatible input from {{ time }} and {{ fail }} columns
#'
#' @param data A data frame of (possibly) censored lifetimes
#' @param time A quasiquoted column of lifetimes
#' @param fail A quasiquoted column of censoring info
#'
#' @return A data frame appended with left/right
#' @export
#'
#' @examples
#' lifetimes <- data.frame(
#'   time = 1:5,
#'   cnsrd = c(TRUE, rep(FALSE, 4)))
#'
#' add_censored_df(lifetimes, fail = !cnsrd)
add_censored_df <-
  function(data, time = time, fail = fail) {
    data |>
      dplyr::mutate(
        left = {{ time }},
        right = dplyr::if_else({{ fail }}, left, Inf)) |>
      as.data.frame()}

#' Append predicted quantiles from survreg model
#'
#' @param data A tibble containing p's from which to make predict quantiles
#' @param model A survreg model
#'
#' @return A tibble appended with predicted quantiles
#' @export
#'
#' @details p's are gotten from km or ppts
#'
add_mle_qtile <- function(data, model) {
  km <- data$km

  predict(model, data,
    type = 'quantile',
    p = km) |>
    base::unique() |>
    base::as.vector() |>
    tibble::as_tibble_col(column_name = 'qtile') |>
    dplyr::bind_cols(data, .)}

#' Extract Weibull alpha from a (tidied) survreg model
#' @param tidy a tidied survreg model
#' @return numeric(1) Weibull alpha
#' @export
get_alpha <- function(tidy) {
  intercept <- purrr::pluck(tidy, 'estimate', 1)
  exp(intercept)}

#' Extract Weibull beta from a (tidied) survreg model
#' @param tidy a tidied survreg model
#' @return numeric(1) Weibull beta
#' @export
get_beta <- function(tidy) {
  log_scale <- purrr::pluck(tidy, 'estimate', 2)
  exp(-log_scale)}

#' Get ppoints (alternative to km)
#' @param time a numeric vector
#' @return numeric vector of plotting positions (the same length as time)
#' @export
get_ppts <- function(time) order_by(
  desc(time), ppoints(time))


#' Append ppts to lifetime tibble
#'
#' @param data A tibble with a time column
#'
#' @return A tibble with 1 - ppts appended
#' @export
#'
add_ppts <- function(data) {
  ppts <- get_ppts(data$time)
  data %>% dplyr::mutate(ppts = 1 - ppts)}

#' Append(?) km points for probability plotting
#'
#' @param data A tibble with time and fail columns
#'
#' @return A tibble with km's appended
#' @export
#'
add_km <- function(data) {
  fit <- survival::survfit(
    survival::Surv(time, fail) ~ 1,
    data = data) %>%
    summary()

  km <- tibble::tibble(
    time = fit$time,
    km = 1 - fit$surv) %>%
    dplyr::mutate(km0 = dplyr::lag(km, default = 0)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      km = base::mean(dplyr::c_across(c(km, km0))),
      km0 = NULL)

    data %>% dplyr::left_join(km, by = 'time')}

#' Fit a survreg model
#'
#' @param data tibble containing lifetimes and censoring info
#'
#' @return A survreg object
#' @export
#'
fit_surv <- function(data) survival::survreg(
  survival::Surv(time, fail) ~ 1,
  data = data)
