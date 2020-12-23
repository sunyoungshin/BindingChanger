#' @name motif_lib
#' @title A sample motif library.
#' @description A list of the position weight matrices corresponding to motifs,
#' loaded by 'data(example)'.
#' @docType data
#' @format A list object.
#' @author Qinyi Zhou \email{qinyi.zhou@utdallas.edu},
#' Chandler Zuo \email{chandler.c.zuo@@gmail.com},
#' Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}
NULL

#' @name indel_info
#' @title A data set for indel information.
#' @description This is a list object loaded by 'data(example)'. Each element
#' corresponds to an indel and includes the following components:
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
#' @docType data
#' @format A list object.
#' @author Qinyi Zhou \email{qinyi.zhou@utdallas.edu},
#' Chandler Zuo \email{chandler.c.zuo@@gmail.com},
#' Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}
NULL

#' @name prior
#' @title Default stationary distribution for nucleotide sequences in the
#' reference genome.
#' @description This parameter is fitted using 61bp windowns around the SNPs in
#' the NHGRI catalog. Loaded by 'data(default_par)'.
#' @docType data
#' @format A numeric vector.
#' @author Qinyi Zhou \email{qinyi.zhou@utdallas.edu},
#' Chandler Zuo \email{chandler.c.zuo@@gmail.com},
#' Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}
NULL

#' @name trans_mat
#' @title Default transition probability matrix for nucleotide sequences in the
#' reference genome.
#' @description ?? This parameter is fitted using 61bp windowns around the
#' SNPs in the NHGRI catalog. Loaded by 'data(default_par)'.
#' @docType data
#' @format A 4 by 4 numeric matrix.
#' @author Qinyi Zhou \email{qinyi.zhou@utdallas.edu},
#' Chandler Zuo \email{chandler.c.zuo@@gmail.com},
#' Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}
NULL
