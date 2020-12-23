pvalue_change <-
  function(insertion_len,
           num_motifs,
           motif_lib,
           indel_info,
           motif_scores,
           prior,
           trans_mat,
           sample_size) {
    x <- which(lapply(indel_info, `[[`, 2) == insertion_len)
    indel_info_selected <- indel_info[x]
    num_indels <- length(x)
    selected_motif_scores <- list()
    selected_motif_scores$match_pos_short <-
      motif_scores$match_pos_short[x, , drop = FALSE]
    selected_motif_scores$match_pos_long <-
      motif_scores$match_pos_long[x, , drop = FALSE]
    selected_motif_scores$log_lik_ratio <-
      motif_scores$log_lik_ratio[x, , drop = FALSE]
    selected_motif_scores$log_lik_short <-
      motif_scores$log_lik_short[x, , drop = FALSE]
    selected_motif_scores$log_lik_long <-
      motif_scores$log_lik_long[x, , drop = FALSE]
    selected_motif_scores$motif <- motif_scores$motif
    selected_motif_scores$k <- motif_scores$k
    results2 <- list()
    result_id <- 1

    for (indel_id in seq_along(indel_info_selected)) {
      for (motif_id in seq_along(motif_lib)) {
        # Step 2. Compute single allele p-values
        this_indel_info <- indel_info_selected[[indel_id]]
        pwm <- motif_lib[[motif_id]]
        if (length(indel_info_selected) == 1) {
          scores <-
            t(
              c(
                selected_motif_scores$log_lik_long[motif_id],
                selected_motif_scores$log_lik_short[motif_id]
              )
            )
        } else{
          scores <-
            t(
              c(
                selected_motif_scores$log_lik_long[indel_id, motif_id],
                selected_motif_scores$log_lik_short[indel_id, motif_id]
              )
            )
        }

        p_value_affinity <- rep(0, 2)
        for (j in seq(2)) {
          if (j == 1) {
            # for long sequence
            sample_seq_len <-
              2 * nrow(pwm) - 2 + this_indel_info$insertion_len
            if (length(indel_info_selected) == 1) {
              reference_score <-
                selected_motif_scores$log_lik_long[motif_id]
            } else{
              reference_score <-
                selected_motif_scores$log_lik_long[indel_id, motif_id]
            }

          } else  if (j == 2) {
            # j=2 for short sequence
            sample_seq_len <- 2 * nrow(pwm) - 2
            if (length(indel_info_selected) == 1) {
              reference_score <-
                selected_motif_scores$log_lik_short[motif_id]
            } else{
              reference_score <-
                selected_motif_scores$log_lik_short[indel_id, motif_id]
            }

          }
          # Compute theta parameter in importance sampling distribution
          theta <- .Call(
            "test_find_theta",
            pwm,
            prior,
            trans_mat,
            # Importance sample scores will have
            # average value of reference score.
            reference_score,
            sample_seq_len,
            package = "atIndel"
          )
          p_value_affinity[j] <- pval_with_less_var(
            .Call(
              "compute_p_values",
              # PWM
              pwm,
              # MC stationary distribution
              prior,
              # transition matrix
              trans_mat,
              scores[, j],
              # theta parameter in importance sampling
              theta,
              sample_size,
              # The sequence length
              sample_seq_len,
              # Use 1 for mean log lik scores
              loglik_type = 0,
              package = "atIndel"
            )[, seq(4)] # the first 4 columns are p-values
            # The last 4 columns are conditional p-values and are not used.
          )[, 1]
        }


        mat_d <-
          comp_indel_mat_d(pwm, prior, this_indel_info$insertion_len)
        score_diff <- c(scores[, 1] - scores[, 2])
        # reference_score is used to compute the theta parameter in importance
        # sampling.
        reference_score <-
          selected_motif_scores$log_lik_long[indel_id, motif_id] - selected_motif_scores$log_lik_short[indel_id, motif_id]
        p_value_change <-
          .Call(
            "p_value_change_indel",
            # Markov Chain transition matrix
            trans_mat,
            # Markov Chain stationary distribution
            prior,
            # The D matrix used to induce binding affinity change
            mat_d,
            # Insertion length
            this_indel_info$insertion_len,
            # PWM
            pwm,
            # Adjusted PWM
            (pwm + 0.25) / 2,
            score_diff,
            c(log(p_value_affinity[1]) - log(p_value_affinity[2])),
            # This is used to compute the theta parameter in importance
            # sampling.
            reference_score,
            sample_size,
            loglik_type = 0,
            package = "atIndel"
          )
        results2[[result_id]] <- list(
          motif_scores = scores,
          p_value_change = list(
            rank = pval_with_less_var(p_value_change$rank)[, 1],
            score = pval_with_less_var(p_value_change$score)[, 1]
          ),
          p_value_affinity1 = p_value_affinity[1],
          p_value_affinity2 = p_value_affinity[2]
        )
        result_id <- result_id + 1
      }
    }
    r <- plyr::ldply (results2, data.frame)
    r <-
      data.frame(id <-
                   rep(names(indel_info_selected), each = length(motif_lib)),
                 motif = rep(names(motif_lib), num_indels),
                 r)
    colnames(r)[1] <- "id"
    r
  }

mac_pvalue_change <- function(i,
                              s,
                              k,
                              num_cores,
                              num_motifs,
                              motif_lib,
                              indel_info,
                              motif_scores,
                              prior,
                              trans_mat,
                              sample_size) {
  a <- list()
  a <- BiocParallel::bpmapply(
    function(x)
      mac_pvalue_change.par(
        i = i,
        insertion_len = x,
        k = k,
        num_cores = num_cores,
        num_motifs = num_motifs,
        motif_lib = motif_lib,
        indel_info = indel_info,
        motif_scores = motif_scores,
        prior = prior,
        trans_mat = trans_mat,
        sample_size = sample_size
      ),
    s,
    BPPARAM = MulticoreParam(workers = num_cores),
    SIMPLIFY = FALSE
  )
  do.call(rbind.data.frame, a)
}

mac_pvalue_change.par <-
  function(i,
           insertion_len,
           k,
           num_cores,
           num_motifs,
           motif_lib,
           indel_info,
           motif_scores,
           prior,
           trans_mat,
           sample_size) {
    x <- which(lapply(indel_info, `[[`, 2) == insertion_len)
    indel_info_selected <- indel_info[x]
    sequence_len <- length(x)
    selected_motif_scores <- list()
    selected_motif_scores$match_pos_short <-
      motif_scores$match_pos_short[x,]
    selected_motif_scores$match_pos_long <-
      motif_scores$match_pos_long[x,]
    selected_motif_scores$log_lik_ratio <-
      motif_scores$log_lik_ratio[x,]
    selected_motif_scores$log_lik_short <-
      motif_scores$log_lik_short[x,]
    selected_motif_scores$log_lik_long <-
      motif_scores$log_lik_long[x,]
    selected_motif_scores$motif <- motif_scores$motif
    selected_motif_scores$k <- motif_scores$k
    prior <- prior
    trans_mat <- trans_mat
    results2 <- list()
    result_id <- 1

    if (num_motifs >= k * num_cores + i) {
      nm <- c(0:k) * num_cores + i
    } else{
      nm <- c(0:(k - 1)) * num_cores + i
    }

    for (indel_id in seq_along(indel_info_selected)) {
      for (motif_id in c(nm)) {
        # Step 2. Compute single allele p-values
        indel_info <- indel_info_selected[[indel_id]]
        pwm <- motif_lib[[motif_id]]
        if (length(indel_info_selected) == 1) {
          scores <-
            t(
              c(
                selected_motif_scores$log_lik_long[motif_id],
                selected_motif_scores$log_lik_short[motif_id]
              )
            )
        } else{
          scores <-
            t(
              c(
                selected_motif_scores$log_lik_long[indel_id, motif_id],
                selected_motif_scores$log_lik_short[indel_id, motif_id]
              )
            )
        }

        p_value_affinity <- rep(0, 2)
        for (j in seq(2)) {
          if (j == 1) {
            # for long sequence
            sample_seq_len <-
              2 * nrow(pwm) - 2 + indel_info$insertion_len
            if (length(indel_info_selected) == 1) {
              reference_score <-
                selected_motif_scores$log_lik_long[motif_id]
            } else{
              reference_score <-
                selected_motif_scores$log_lik_long[indel_id, motif_id]
            }

          } else  if (j == 2) {
            # j=2 for short sequence
            sample_seq_len <- 2 * nrow(pwm) - 2
            if (length(indel_info_selected) == 1) {
              reference_score <-
                selected_motif_scores$log_lik_short[motif_id]
            } else{
              reference_score <-
                selected_motif_scores$log_lik_short[indel_id, motif_id]
            }

          }
          # Compute theta parameter in importance sampling distribution
          theta <- .Call(
            "test_find_theta",
            pwm,
            prior,
            trans_mat,
            # Importance sample scores will have
            # average value of reference score.
            reference_score,
            sample_seq_len,
            package = "atIndel"
          )
          p_value_affinity[j] <- pval_with_less_var(
            .Call(
              "compute_p_values",
              # PWM
              pwm,
              # MC stationary distribution
              prior,
              # transition matrix
              trans_mat,
              scores[, j],
              # theta parameter in importance sampling
              theta,
              sample_size,
              # The sequence length
              sample_seq_len,
              # Use 1 for mean log lik scores
              loglik_type = 0,
              package = "atIndel"
            )[, seq(4)] # The first 4 columns are p-values
            # The last 4 columns are conditional p-values and are not used.
          )[, 1]
        }


        mat_d <-
          comp_indel_mat_d(pwm, prior, indel_info$insertion_len)
        score_diff <- c(scores[, 1] - scores[, 2])
        # reference_score is used to compute the theta parameter in importance
        # sampling
        reference_score <-
          selected_motif_scores$log_lik_long[indel_id, motif_id] - selected_motif_scores$log_lik_short[indel_id, motif_id]
        p_value_change <-
          .Call(
            "p_value_change_indel",
            # Markov Chain transition matrix
            trans_mat,
            # Markov Chain stationary distribution
            prior,
            # The D matrix used to induce binding affinity change
            mat_d,
            # Insertion length
            indel_info$insertion_len,
            # PWM
            pwm,
            # Adjusted PWM
            (pwm + 0.25) / 2,
            score_diff,
            c(log(p_value_affinity[1]) - log(p_value_affinity[2])),
            # This is used to compute the theta parameter in importance
            # sampling.
            reference_score,
            sample_size,
            loglik_type = 0
          )
        results2[[result_id]] <- list(
          motif_scores = scores,
          p_value_change = list(
            rank = pval_with_less_var(p_value_change$rank)[, 1],
            score = pval_with_less_var(p_value_change$score)[, 1]
          ),
          p_value_affinity1 = p_value_affinity[1],
          p_value_affinity2 = p_value_affinity[2]
        )
        result_id <- result_id + 1
      }
    }
    r <- ldply(results2, data.frame)
    r <-
      data.frame(id <-
                   rep(names(indel_info_selected), each = length(nm)),
                 motif = rep(names(motif_lib[nm]), sequence_len),
                 r)
    colnames(r)[1] <- "id"
    r
  }

make_insertion_tbl <- function(a, long_insertion) {
  a <- merge(a, long_insertion, by = "id", all.x = TRUE)
  a$ref.score <-
    ifelse(a$insertion, a$motif_scores.2, a$motif_scores.1)
  a$mutation.score <-
    ifelse(a$insertion, a$motif_scores.1, a$motif_scores.2)
  a$ref.pval <-
    ifelse(a$insertion, a$p_value_affinity2, a$p_value_affinity1)
  a$mutation.pval <-
    ifelse(a$insertion, a$p_value_affinity1, a$p_value_affinity2)
  a2 <- a[, c(1, 2, 9, 10, 12, 13, 14, 15, 5, 6, 11)]
  a2
}

make_motifscore_insertion_tbl <- function(a, long_insertion) {
  a <- merge(a, long_insertion, by = "id", all.x = TRUE)
  a$ref.score <-
    ifelse(a$insertion, a$log_lik_short, a$log_lik_long)
  a$mutation.score <-
    ifelse(a$insertion, a$log_lik_long, a$log_lik_short)
  a2 <- a[, c(1, 2, 6, 7, 9, 10, 5, 8)]
  a2
}


#' @name indel_motif_scores
#' @title Compute the motif scores given a motif library and a list of indels.
#' @description Compute the motif scores given a motif library and a list of
# indels.
#' @param motif_lib A list of the position weight matrices for the motifs.
#' @param indel_info A list object. Each element corresponds to an indel and
#'' includes the following components:
#' \tabular{ll}{
#' inserted_sequence \tab An integer vector representing the deroxyribose
#' sequence around each indel location for the longer sequence. For an
#' insertion, this is the inserted sequence. For an deletion, this is the
#' reference allele sequence. Nucleotides are coded as "A"-1, "C"-2, "G"-3,
#' "T"-4.\cr
#' insertion_len \tab The length of the inserted / deleted sequence.\cr
#' insertion \tab A boolean variable indicating whether the indel is an
#' insertion (TRUE) or a deletion (FALSE).\cr
#' ref \tab A character vector TODO\cr
#' alt \tab A character vector TODO\cr
#' }
#' @param num_cores An integer for the number of parallel processes used for
#' parallel computation.
#' @details TODO.
#' @return A list object of position weight matrices.
#' @author Qinyi Zhou \email{qinyi.zhou@utdallas.edu},
#' Chandler Zuo \email{chandler.c.zuo@@gmail.com},
#' Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}
#' @examples
#' data(example)
#' indel_motif_scores(motif_lib, indel_info, num_cores=1)
#' @useDynLib atIndel
#' @import BiocParallel Rcpp
#' @export
indel_motif_scores <-
  function(motif_lib, indel_info, num_cores = 1) {
    motifs <- names(motif_lib)
    ids <- names(indel_info)
    num_motifs <- length(motif_lib)
    sequence_len <- length(indel_info)
    num_cores <- min(num_cores, num_motifs)
    k <- as.integer(num_motifs / num_cores)
    insertion <- unlist(lapply(indel_info, `[[`, 3))
    ref <- unlist(lapply(indel_info, `[[`, 4))
    alt <- unlist(lapply(indel_info, `[[`, 5))
    long_insertion <-
      data.frame(
        id = ids,
        ref = ref,
        alt = alt,
        insertion = insertion
      )
    if (num_cores > 1) {
      if (Sys.info()[["sysname"]] == "Windows") {
        snow <- BiocParallel::SnowParam(workers = num_cores, type = "SOCK")
        motif_score_par <- function(i,
                                    k,
                                    num_cores,
                                    num_motifs,
                                    sequence_len,
                                    motif_lib,
                                    indel_info) {
          k <- list()
          if (num_motifs >= ((k - 1) * num_cores + i)) {
            nm <- c(0:(k - 1)) * num_cores + i
          } else{
            nm <- c(0:(k - 2)) * num_cores + i
          }
          motif_scores <- .Call(
            "comp_indel_motif_scores",
            motif_lib[nm],
            indel_info,
            # select the log-lik type here
            loglik_type = 0,
            package = "atIndel"
          )
          k$num <- nm
          k$motif_scores <- motif_scores
          k
        }

        motif_score_par_list <-
          BiocParallel::bpmapply(
            function(x)
              motif_score_par(
                i = x,
                k = k,
                num_cores = num_cores,
                num_motifs = num_motifs,
                sequence_len = sequence_len,
                motif_lib = motif_lib,
                indel_info = indel_info
              ),
            seq(num_cores),
            BPPARAM = snow,
            SIMPLIFY = FALSE
          )
        nm <- unlist(lapply(motif_score_par_list, `[[`, 1))
        id <- rep(ids, num_motifs)
        ins <- rep(insertion, num_motifs)
        motif <- rep(motifs[nm], each = sequence_len)
        ms <- lapply(lapply(motif_score_par_list, `[[`, 2), `[[`, 1)
        ml <- lapply(lapply(motif_score_par_list, `[[`, 2), `[[`, 2)
        short <-
          lapply(lapply(motif_score_par_list, `[[`, 2), `[[`, 4)
        long <- lapply(lapply(motif_score_par_list, `[[`, 2), `[[`, 5)
        ratio <-
          lapply(lapply(motif_score_par_list, `[[`, 2), `[[`, 3)
        match_pos_short <- ms[[1]]
        match_pos_long <- ml[[1]]
        log_lik_short <- short[[1]]
        log_lik_long <- long[[1]]
        log_lik_ratio <- ratio[[1]]
        for (i in (2:num_cores)) {
          match_pos_short <- cbind(match_pos_short, ms[[i]])
          match_pos_long <- cbind(match_pos_long, ml[[i]])
          log_lik_short <- cbind(log_lik_short, short[[i]])
          log_lik_long <- cbind(log_lik_long, long[[i]])
          log_lik_ratio <- cbind(log_lik_ratio, ratio[[i]])
        }
        motif_scores <- list()
        motif_scores$match_pos_short <- match_pos_short
        motif_scores$match_pos_long <- match_pos_long
        motif_scores$log_lik_ratio <- log_lik_ratio
        motif_scores$log_lik_short <- log_lik_short
        motif_scores$log_lik_long <- log_lik_long
        result <- list()
        x <-
          data.frame(
            id,
            motif,
            log_lik_short = as.vector(log_lik_short),
            log_lik_long = as.vector(log_lik_long),
            log_lik_ratio = as.vector(log_lik_ratio)
          )
        result$table <-
          make_motifscore_insertion_tbl(x, long_insertion)
        result$list <- motif_scores
        result$list$motif <- motifs
        result$list$k <- nm
        result$list$insertion <- insertion
        result$list$ref <- ref
        result$list$alt <- alt
        result
      } else{
        mac_motif_score_par <-
          function(i,
                   k,
                   num_cores,
                   num_motifs,
                   sequence_len,
                   motif_lib,
                   indel_info) {
            k <- list()
            if (num_motifs >= k * num_cores + i) {
              nm <- c(0:k) * num_cores + i
            } else{
              nm <- c(0:(k - 1)) * num_cores + i
            }
            motif_scores <- .Call(
              "comp_indel_motif_scores",
              motif_lib[nm],
              indel_info,
              # select the log-lik type here
              loglik_type = 0,
              package = "atIndel"
            )
            k$num <- nm
            k$motif_scores <- motif_scores
            k
          }
        motif_score_par_list <-
          BiocParallel::bpmapply(
            function(x)
              mac_motif_score_par(
                i = x,
                k = k,
                num_cores = num_cores,
                num_motifs = num_motifs,
                sequence_len = sequence_len,
                motif_lib = motif_lib,
                indel_info = indel_info
              ),
            seq(num_cores),
            BPPARAM = MulticoreParam(workers = num_cores),
            SIMPLIFY = FALSE
          )
        nm <- unlist(lapply(motif_score_par_list, `[[`, 1))
        id <- rep(ids, num_motifs)
        motif <- rep(motifs[nm], each = sequence_len)
        ms <- lapply(lapply(motif_score_par_list, `[[`, 2), `[[`, 1)
        ml <- lapply(lapply(motif_score_par_list, `[[`, 2), `[[`, 2)
        short <-
          lapply(lapply(motif_score_par_list, `[[`, 2), `[[`, 4)
        long <- lapply(lapply(motif_score_par_list, `[[`, 2), `[[`, 5)
        ratio <-
          lapply(lapply(motif_score_par_list, `[[`, 2), `[[`, 3)
        match_pos_short <- ms[[1]]
        match_pos_long <- ml[[1]]
        log_lik_short <- short[[1]]
        log_lik_long <- long[[1]]
        log_lik_ratio <- ratio[[1]]
        for (i in (2:num_cores)) {
          match_pos_short <- cbind(match_pos_short, ms[[i]])
          match_pos_long <- cbind(match_pos_long, ml[[i]])
          log_lik_short <- cbind(log_lik_short, short[[i]])
          log_lik_long <- cbind(log_lik_long, long[[i]])
          log_lik_ratio <- cbind(log_lik_ratio, ratio[[i]])
        }
        motif_scores <- list()
        motif_scores$match_pos_short <- match_pos_short
        motif_scores$match_pos_long <- match_pos_long
        motif_scores$log_lik_ratio <- log_lik_ratio
        motif_scores$log_lik_short <- log_lik_short
        motif_scores$log_lik_long <- log_lik_long
        result <- list()
        x <-
          data.frame(
            id,
            motif,
            log_lik_short = as.vector(log_lik_short),
            log_lik_long = as.vector(log_lik_long),
            log_lik_ratio = as.vector(log_lik_ratio)
          )
        result$table <-
          make_motifscore_insertion_tbl(x, long_insertion)
        result$list <- motif_scores
        result$list$motif <- motifs
        result$list$k <- nm
        result$list$insertion <- insertion
        result$list$ref <- ref
        result$list$alt <- alt
      }
    } else{
      motif_scores <- .Call(
        "comp_indel_motif_scores",
        motif_lib,
        indel_info,
        # select the log-lik type here
        loglik_type = 0,
        package = "atIndel"
      )

      id <- rep(ids, num_motifs)
      motif <- rep(motifs, each = sequence_len)
      log_lik_short <- as.vector(motif_scores$log_lik_short)
      log_lik_long <- as.vector(motif_scores$log_lik_long)
      log_lik_ratio <- as.vector(motif_scores$log_lik_ratio)

      result <- list()
      x <-
        data.frame(
          id,
          motif,
          log_lik_short = as.vector(log_lik_short),
          log_lik_long = as.vector(log_lik_long),
          log_lik_ratio = as.vector(log_lik_ratio)
        )
      result$table <- make_motifscore_insertion_tbl(x, long_insertion)
      result$list <- motif_scores
      result$list$motif <- motifs
      result$list$k <- seq(num_motifs)
      result$list$insertion <- insertion
      result$list$ref <- ref
      result$list$alt <- alt
    }
    return(result)
  }


#' @name indel_p_values
#' @title Compute the motif scores given a motif library and a list of indels.
#' @description Compute the motif scores given a motif library and a list of
# indels.
#' @param motif_lib A list of the position weight matrices for the motifs.
#' @param indel_info A list object. Each element corresponds to an indel and
#'' includes the following components:
#' \tabular{ll}{
#' inserted_sequence \tab An integer vector representing the deroxyribose
#' sequence around each indel location for the longer sequence. For an
#' insertion, this is the inserted sequence. For an deletion, this is the
#' reference allele sequence. Nucleotides are coded as "A"-1, "C"-2, "G"-3,
#' "T"-4.\cr
#' insertion_len \tab The length of the inserted / deleted sequence.\cr
#' insertion \tab A boolean variable indicating whether the indel is an
#' insertion (TRUE) or a deletion (FALSE).\cr
#' ref \tab A character vector TODO\cr
#' alt \tab A character vector TODO\cr
#' }
#' @param motif_scores A data frame in the same format as the output of
#' \code{\link{indel_motif_scores}}.
#' @param prior A numeric vector for the prior distribution parameters of
# the Markov Chain model for background sequences.
#' @param trans_mat A numeric matrix for the transition matrix parameters
#  of the Markov Chain model for background sequences.
#' @param sample_size An integer for the importance sampling sample size.
#' @param num_cores An integer for the number of parallel processes.
#' @details TODO.
#' @return A list object of position weight matrices.
#' @author Qinyi Zhou \email{qinyi.zhou@utdallas.edu},
#' Chandler Zuo \email{chandler.c.zuo@@gmail.com},
#' Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}
#' @examples
#' data(example)
#' motif_scores <- indel_motif_scores(motif_lib, indel_info)$list
#' indel_p_values(
#'   motif_lib=motif_lib,
#'   indel_info=indel_info,
#'   motif_scores=motif_scores,
#'   prior=prior,
#'   trans_mat=trans_mat,
#'   sample_size=100,
#'   num_cores=1
#' )
#' @useDynLib atIndel
#' @import BiocParallel Rcpp plyr
#' @export
indel_p_values <-
  function(motif_lib,
           indel_info,
           motif_scores,
           prior,
           trans_mat,
           sample_size,
           num_cores = 1) {
    indel_info_sorted <- rlist::list.sort(indel_info, insertion_len)
    m_list <- unlist(lapply(indel_info_sorted, `[[`, 2))
    s <- unique(m_list)
    motif_lib <- motif_lib[motif_scores$k]
    ids <- names(indel_info_sorted)
    num_motifs <- length(motif_lib)
    num_cores <- min(num_cores, num_motifs)
    k <- as.integer(num_motifs / num_cores)
    motif_scores <- motif_scores
    prior <- prior
    insertion <- unlist(lapply(indel_info_sorted, `[[`, 3))
    ref <- unlist(lapply(indel_info_sorted, `[[`, 4))
    alt <- unlist(lapply(indel_info_sorted, `[[`, 5))
    long_insertion <-
      data.frame(
        id = ids,
        ref = ref,
        alt = alt,
        insertion = insertion
      )

    if (num_cores > 1) {
      if (Sys.info()[["sysname"]] == "Windows") {
        results2 <- list()
        result_id <- 1
        snow <- SnowParam(workers = num_cores, type = "SOCK")
        motif_lib <- motif_lib
        motif_scores <- motif_scores
        win_pvalue_change.par <-
          function(i,
                   insertion_len,
                   k,
                   num_cores,
                   num_motifs,
                   motif_lib,
                   indel_info,
                   motif_scores,
                   prior,
                   trans_mat,
                   sample_size) {
            x <- which(lapply(indel_info, `[[`, 2) == insertion_len)
            indel_info_selected <- indel_info[x]
            sequence_len <- length(x)
            motif_scores <- list()
            motif_scores$match_pos_short <-
              motif_scores$match_pos_short[x,]
            motif_scores$match_pos_long <-
              motif_scores$match_pos_long[x,]
            motif_scores$log_lik_ratio <-
              motif_scores$log_lik_ratio[x,]
            motif_scores$log_lik_short <-
              motif_scores$log_lik_short[x,]
            motif_scores$log_lik_long <-
              motif_scores$log_lik_long[x,]
            motif_scores$motif <- motif_scores$motif
            motif_scores$k <- motif_scores$k
            results2 <- list()
            result_id <- 1
            if (num_motifs >= k * num_cores + i) {
              nm <- c(seq_len(k)) * num_cores + i
            } else{
              nm <- c(seq_len(k - 1)) * num_cores + i
            }

            for (indel_id in seq_along(indel_info_selected)) {
              for (motif_id in c(nm)) {
                # Step 2. Compute single allele p-values
                indel_info <- indel_info_selected[[indel_id]]
                pwm <- motif_lib[[motif_id]]
                if (length(indel_info_selected) == 1) {
                  scores <-
                    t(c(
                      motif_scores$log_lik_long[motif_id],
                      motif_scores$log_lik_short[motif_id]
                    ))
                } else{
                  scores <-
                    t(c(
                      motif_scores$log_lik_long[indel_id, motif_id],
                      motif_scores$log_lik_short[indel_id, motif_id]
                    ))
                }

                p_value_affinity <- rep(0, 2)
                for (j in seq(2)) {
                  if (j == 1) {
                    # for long sequence
                    sample_seq_len <-
                      2 * nrow(pwm) - 2 + indel_info$insertion_len
                    if (length(indel_info_selected) == 1) {
                      reference_score <-
                        motif_scores$log_lik_long[motif_id]
                    } else{
                      reference_score <-
                        motif_scores$log_lik_long[indel_id, motif_id]
                    }

                  } else  if (j == 2) {
                    # j=2 for short sequence
                    sample_seq_len <- 2 * nrow(pwm) - 2
                    if (length(indel_info_selected) == 1) {
                      reference_score <-
                        motif_scores$log_lik_short[motif_id]
                    } else{
                      reference_score <-
                        motif_scores$log_lik_short[indel_id, motif_id]
                    }

                  }
                  # Compute theta parameter in importance sampling distribution
                  theta <- .Call(
                    "test_find_theta",
                    pwm,
                    prior,
                    trans_mat,
                    # Importance sample scores will have
                    # average value of reference score.
                    reference_score,
                    sample_seq_len,
                    package = "atIndel"
                  )
                  p_value_affinity[j] <- pval_with_less_var(
                    .Call(
                      "compute_p_values",
                      # PWM
                      pwm,
                      # MC stationary distribution
                      prior,
                      # transition matrix
                      trans_mat,
                      scores[, j],
                      # theta parameter in importance sampling
                      theta,
                      sample_size,
                      # The sequence length
                      sample_seq_len,
                      # Use 1 for mean log lik scores
                      loglik_type = 0,
                      package = "atIndel"
                    )[, seq(4)] # the first 4 columns are p-values
                    # the last 4 columns are conditional p-values and are not useful here
                  )[, 1]
                }


                mat_d <-
                  comp_indel_mat_d(pwm, prior, indel_info$insertion_len)
                score_diff <- c(scores[, 1] - scores[, 2])
                # reference_score is used to compute the theta parameter in importance sampling
                reference_score <-
                  motif_scores$log_lik_long[indel_id, motif_id] - motif_scores$log_lik_short[indel_id, motif_id]
                p_value_change <-
                  .Call(
                    "p_value_change_indel",
                    # Markov Chain transition matrix
                    trans_mat,
                    # Markov Chain stationary distribution
                    prior,
                    # The D matrix used to induce binding affinity change
                    mat_d,
                    # Insertion length
                    indel_info$insertion_len,
                    # PWM
                    pwm,
                    # Adjusted PWM
                    (pwm + 0.25) / 2,
                    score_diff,
                    c(log(p_value_affinity[1]) - log(p_value_affinity[2])),
                    # This is used to compute the theta parameter in importance
                    # sampling.
                    reference_score,
                    sample_size,
                    loglik_type = 0,
                    package = "atIndel"
                  )
                results2[[result_id]] <- list(
                  motif_scores = scores,
                  p_value_change = list(
                    rank = pval_with_less_var(p_value_change$rank)[, 1],
                    score = pval_with_less_var(p_value_change$score)[, 1]
                  ),
                  p_value_affinity1 = p_value_affinity[1],
                  p_value_affinity2 = p_value_affinity[2]
                )
                result_id <- result_id + 1
              }
            }
            r <- ldply (results2, data.frame)
            r <-
              data.frame(id <-
                           rep(names(indel_info_selected), each = length(nm)),
                         motif = rep(names(motif_lib[nm]), sequence_len),
                         r)
            colnames(r)[1] <- "id"
            r
          }

        win_pvalue_change <- function(i,
                                      s,
                                      k,
                                      num_cores,
                                      num_motifs,
                                      motif_lib,
                                      indel_info,
                                      motif_scores,
                                      prior,
                                      trans_mat,
                                      sample_size) {
          a <- list()
          a <- BiocParallel::bpmapply(
            function(x)
              win_pvalue_change.par(
                i = i,
                insertion_len = x,
                k = k,
                num_cores = num_cores,
                num_motifs = num_motifs,
                motif_lib = motif_lib,
                indel_info = indel_info,
                motif_scores = motif_scores,
                prior = prior,
                trans_mat = trans_mat,
                sample_size = sample_size
              ),
            s,
            BPPARAM = MulticoreParam(workers = num_cores),
            SIMPLIFY = FALSE
          )
          do.call(rbind.data.frame, a)
          # b[order(id),]
        }

        # pval_par_list<-data.frame()
        #
        # for(i in seq(s)){
        #   l1<-BiocParallel::bpmapply(
        #     function(x,y)
        #       win_pvalue_change(
        #         i = x,
        #         insertion_len = y,
        #         k = k,
        #         num_cores = num_cores,
        #         num_motifs = num_motifs,
        #         motif_lib = motif_lib,
        #         indel_info = indel_info,
        #         motif_scores=motif_scores,
        #         prior=prior,
        #         trans_mat=trans_mat,
        #         sample_size=sample_size
        #       ),
        #     seq(num_cores),s[i],
        #     BPPARAM = snow,
        #     SIMPLIFY =FALSE
        #   )
        #   pval_par_list<-rbind(pval_par_list,do.call(rbind.data.frame,l1))
        # }
        # a<-pval_par_list[order(pval_par_list$id),]
        # rownames(a) = seq_len(nrow(a))
        # a
        pval_par_list <- list()
        pval_par_list <-
          BiocParallel::bpmapply(
            function(x)
              win_pvalue_change(
                i = x,
                s = s,
                k = k,
                num_cores = num_cores,
                num_motifs = num_motifs,
                motif_lib = motif_lib,
                indel_info = indel_info,
                motif_scores = motif_scores,
                prior = prior,
                trans_mat = trans_mat,
                sample_size = sample_size
              ),
            seq(num_cores),
            BPPARAM = snow,
            SIMPLIFY = FALSE
          )
        b <- do.call(rbind.data.frame,  pval_par_list)
        b <- b[order(b$id),]
        rownames(b) <- seq_len(nrow(b))
        b <- make_insertion_tbl(b, long_insertion)
        b
      } else
      {
        pval_par_list <-
          BiocParallel::bpmapply(
            function(x)
              mac_pvalue_change(
                i = x,
                s = s,
                k = k,
                num_cores = num_cores,
                num_motifs = num_motifs,
                motif_lib = motif_lib,
                indel_info = indel_info,
                motif_scores = motif_scores,
                prior = prior,
                trans_mat = trans_mat,
                sample_size = sample_size
              ),
            seq(num_cores),
            BPPARAM = MulticoreParam(workers = num_cores),
            SIMPLIFY = FALSE
          )

        b <- do.call(rbind.data.frame,  pval_par_list)
        b <- b[order(b$id),]
        rownames(b) <- seq_len(nrow(b))
        b <- make_insertion_tbl(b, long_insertion)
        b
      }

      # pval_par_list <-
      #   BiocParallel::bpmapply(
      #     function(x,y)
      #       mac_pvalue_change.par(
      #         i = x,
      #         insertion_len = y,
      #         k = k,
      #         num_cores = num_cores,
      #         num_motifs = num_motifs,
      #         motif_lib = motif_lib,
      #         indel_info = indel_info,
      #         motif_scores=motif_scores,
      #         prior=prior,
      #         trans_mat=trans_mat,
      #         sample_size=sample_size
      #       ),
      #     seq(num_cores),s,
      #     BPPARAM = MulticoreParam(workers =num_cores),
      #     SIMPLIFY =FALSE
      #   )
      #do.call(rbind.data.frame,  pval_par_list)
      # a<-pval_par_list[order(pval_par_list$id),]
      # rownames(a) = seq_len(nrow(a))
      # a
    } else
    {
      #pval_par_list<-data.frame()

      # for(i in seq(s)){
      #   pval_par_list<-rbind(pval_par_list,pvalue_change(insertion_len =s[i],
      #                                   num_motifs = num_motifs,
      #                                   motif_lib = motif_lib,
      #                                   indel_info = indel_info,
      #                                   motif_scores=motif_scores,
      #                                   prior=prior,
      #                                   trans_mat=trans_mat,
      #                                   sample_size=sample_size))
      # }
      pval_par_list <-
        mapply(
          function(x)
            pvalue_change(
              insertion_len = x,
              num_motifs = num_motifs,
              motif_lib = motif_lib,
              indel_info = indel_info,
              motif_scores = motif_scores,
              prior = prior,
              trans_mat = trans_mat,
              sample_size = sample_size
            ),
          s,
          SIMPLIFY = FALSE
        )

      b <- do.call(rbind.data.frame,  pval_par_list)
      b <- b[order(b$id),]
      rownames(b) <- seq_len(nrow(b))
      b <- make_insertion_tbl(b, long_insertion)
      b
    }

  }
