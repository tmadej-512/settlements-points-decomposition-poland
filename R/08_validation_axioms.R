# ============================================================
# R/08_validation_axioms.R
# ------------------------------------------------------------
# Formal validation that the implemented Shapley decomposition
# satisfies the axioms it is supposed to satisfy. This is the
# "does the method do what it claims to do, mathematically"
# check referenced in the paper's validation section, distinct
# from the empirical/statistical validation in R/09 and R/10.
#
# Axioms checked, given the coalition value function v(S) =
# scenario_from_set(..., active_set = S):
#
#   1. EFFICIENCY:    sum_j phi_j  ==  v(full) - v(empty)  ==  total_change
#                      (checked to within floating-point tolerance)
#
#   2. NULL PLAYER:   if component j is a "dummy" (switching it on
#                      changes v(S) by exactly the same amount for
#                      every coalition S, including by exactly the
#                      coalition-value contribution when toggled alone
#                      -- the degenerate case being v(S u {j}) = v(S)
#                      for ALL S, i.e. j has literally no effect),
#                      then phi_j == 0. We test this constructively by
#                      injecting a synthetic dummy component into the
#                      coalition-value function (see test_null_player()).
#
#   3. SYMMETRY:      if components j and k are interchangeable (i.e.
#                      v(S u {j}) - v(S) == v(S u {k}) - v(S) for every
#                      coalition S not containing j or k), then
#                      phi_j == phi_k. Tested constructively the same way.
#
#   4. ADDITIVITY (of the value function, not the allocation): not
#      separately unit-tested here, because it is a property of how
#      v(S) itself would behave under summed games, not something to
#      check against a single empirical game; efficiency (1) is the
#      operationally relevant guarantee for this paper's purposes.
#
# All four checks are run against TOY coalition-value functions
# with known structure (so the "correct" Shapley value is known
# analytically), NOT against the fitted NFZ model -- that would only
# tell us the model's predictions are internally consistent, not that
# our Shapley implementation is computing the Shapley value correctly.
# The empirical check against the real model (does
# sum(N,PI,LAMBDA,KAPPA effects) == total_change for the model-implied
# fitted model) is done separately by check_efficiency_empirical().
# ============================================================

suppressPackageStartupMessages(library(dplyr))

#' Generic Shapley-value calculator for an arbitrary coalition-value function
#'
#' Used only by the axiom unit tests below, on toy `v_fun`s with known
#' analytic Shapley values -- NOT used in the main NFZ pipeline (which
#' calls decompose_shapley_2y() directly for efficiency, since that
#' function also needs the model-specific newdata machinery).
#'
#' @param players Character vector of player names.
#' @param v_fun Function taking a character vector (subset of `players`)
#'   and returning a single numeric (the coalition's value). v_fun(character(0))
#'   should return the empty-coalition value (usually 0).
#' @return Named numeric vector of Shapley values, one per player.
shapley_value_generic <- function(players, v_fun) {
  k <- length(players)
  perms <- if (k <= 1) list(players) else combinat_permutations_n(players)
  contrib <- matrix(0, nrow = length(perms), ncol = k, dimnames = list(NULL, players))

  for (i in seq_along(perms)) {
    perm <- perms[[i]]
    active <- character(0)
    v_prev <- v_fun(active)
    for (p in perm) {
      active <- c(active, p)
      v_curr <- v_fun(active)
      contrib[i, p] <- v_curr - v_prev
      v_prev <- v_curr
    }
  }
  colMeans(contrib)
}

#' All permutations of an arbitrary-length vector (small n only; used for
#' axiom tests with 2-4 toy players, not performance-critical)
combinat_permutations_n <- function(v) {
  if (length(v) <= 1) return(list(v))
  res <- list()
  for (i in seq_along(v)) {
    rest <- v[-i]
    for (sub in combinat_permutations_n(rest)) {
      res <- c(res, list(c(v[i], sub)))
    }
  }
  res
}

#' AXIOM 1 (constructive check): efficiency on a toy additive game
#'
#' v(S) = sum of weights of players in S. The Shapley value of an
#' additive game is exactly each player's own weight, and the sum
#' across players must equal v(full).
test_efficiency_toy <- function(tol = 1e-8) {
  weights <- c(N = 3.5, PI = -1.2, LAMBDA = 7.0, KAPPA = 0.4)
  v_fun <- function(S) sum(weights[S])
  phi <- shapley_value_generic(names(weights), v_fun)

  list(
    passed = all(abs(phi - weights) < tol) && abs(sum(phi) - sum(weights)) < tol,
    shapley_values = phi,
    expected = weights,
    max_abs_error = max(abs(phi - weights))
  )
}

#' AXIOM 2 (constructive check): null / dummy player
#'
#' Add a 5th "dummy" component whose presence never changes the
#' coalition value (v(S u {D}) = v(S) for all S). Its Shapley value
#' must be exactly zero, and the other four players' Shapley values
#' must be unaffected by its presence.
test_null_player <- function(tol = 1e-8) {
  weights <- c(N = 3.5, PI = -1.2, LAMBDA = 7.0, KAPPA = 0.4)
  v_fun_with_dummy <- function(S) sum(weights[intersect(S, names(weights))])  # "DUMMY" contributes nothing
  players <- c(names(weights), "DUMMY")
  phi <- shapley_value_generic(players, v_fun_with_dummy)

  list(
    passed = abs(phi["DUMMY"]) < tol && all(abs(phi[names(weights)] - weights) < tol),
    shapley_values = phi,
    dummy_value = phi["DUMMY"]
  )
}

#' AXIOM 3 (constructive check): symmetry
#'
#' Two players with identical marginal contribution profiles must
#' receive identical Shapley values. Construct a game where players
#' "A" and "B" are interchangeable (symmetric weight) and confirm
#' phi_A == phi_B.
test_symmetry <- function(tol = 1e-8) {
  weights <- c(A = 2.0, B = 2.0, C = 5.0, D = -1.0)  # A and B tied
  v_fun <- function(S) sum(weights[S]) + 0.5 * (("A" %in% S) * ("C" %in% S)) +
    0.5 * (("B" %in% S) * ("C" %in% S))  # symmetric interaction: A-C and B-C identical
  phi <- shapley_value_generic(names(weights), v_fun)

  list(
    passed = abs(phi["A"] - phi["B"]) < tol,
    shapley_values = phi,
    diff_A_B = phi["A"] - phi["B"]
  )
}

#' EMPIRICAL EFFICIENCY CHECK against the real fitted NFZ model
#'
#' Confirms that, for the actual model chain and actual disease data,
#' the four Shapley component effects sum to the model-implied total_change
#' to within numerical tolerance. This is the check that belongs in
#' the paper's results/validation section as direct evidence the
#' implementation is correct on real data, complementing the toy-game
#' axiom checks above (which confirm the *algorithm* is correct in the
#' abstract).
#'
#' @param shapley_result Output of decompose_shapley_2y() (a 1-row tibble).
#' @param tol Absolute tolerance for the residual.
#' @return A 1-row tibble summarising the check, suitable for a paper table.
check_efficiency_empirical <- function(shapley_result, tol = 1e-6) {
  sum_effects <- with(shapley_result,
                       N_effect + PI_effect + LAMBDA_effect + KAPPA_effect)
  tibble::tibble(
    baseline_year   = shapley_result$baseline_year,
    comparison_year = shapley_result$comparison_year,
    total_change    = shapley_result$total_change,
    sum_of_effects  = sum_effects,
    residual        = shapley_result$residual,
    passed          = abs(shapley_result$residual) < tol
  )
}

#' Run the full axiom validation suite and print a readable summary
run_axiom_validation_suite <- function() {
  res_eff <- test_efficiency_toy()
  res_null <- test_null_player()
  res_sym <- test_symmetry()

  results <- tibble::tibble(
    axiom = c("Efficiency", "Null player", "Symmetry"),
    passed = c(res_eff$passed, res_null$passed, res_sym$passed),
    detail = c(
      sprintf("max |phi - true weight| = %.2e", res_eff$max_abs_error),
      sprintf("dummy player Shapley value = %.2e", res_null$dummy_value),
      sprintf("phi_A - phi_B = %.2e", res_sym$diff_A_B)
    )
  )

  message("==================================================")
  message("SHAPLEY AXIOM VALIDATION (toy games, known ground truth)")
  message("==================================================")
  for (i in seq_len(nrow(results))) {
    status <- if (results$passed[i]) "PASS" else "FAIL"
    message(sprintf("[%s] %-12s %s", status, results$axiom[i], results$detail[i]))
  }
  message("==================================================")

  invisible(results)
}
