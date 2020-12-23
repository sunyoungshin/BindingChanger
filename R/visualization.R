revert.columns <- function(mat) {
  mat[, rev(seq(ncol(mat)))]
}

rep.col <- function(x, n) {
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

makevector <- function(x, strand_sign) {
  if (strand_sign == 1) {
    if (x == 1) {
      a <- c(1, 0, 0, 0)
    } else if (x == 2) {
      a <- c(0, 1, 0, 0)
    } else if (x == 3) {
      a <- c(0, 0, 1, 0)
    } else{
      a <- c(0, 0, 0, 1)
    }
  } else{
    if (x == 1) {
      a <- c(0, 0, 0, 1)
    } else if (x == 2) {
      a <- c(0, 0, 1, 0)
    } else if (x == 3) {
      a <- c(0, 1, 0, 0)
    } else{
      a <- c(1, 0, 0, 0)
    }
  }
  a
}

makepwm <- function(x, strand_sign) {
  a <- mapply(makevector, x, strand_sign)
  a <- unlist(a)
  b <- matrix(a, ncol = length(x))
  b
}

makemotif <- function(x, strand_sign) {
  p <- t(x[[1]])
  if (strand_sign == 1) {
    rownames(p) <- c("A", "C", "G", "T")
  } else{
    p <- revert.columns(p)
    rownames(p) <- c("A", "C", "G", "T")
  }
  p
}

makepfull <- function(p, left, right) {
  if (left > 0 & right > 0) {
    p1 <-
      cbind(rep.col(c(0.25, 0.25, 0.25, 0.25), left), p, rep.col(c(0.25, 0.25, 0.25, 0.25), right))
    p1
  } else if (left > 0 & right == 0) {
    p1 <- cbind(rep.col(c(0.25, 0.25, 0.25, 0.25), left), p)
    p1
  } else{
    p1 <- cbind(p, rep.col(c(0.25, 0.25, 0.25, 0.25), right))
    p1
  }
}

plotbinding <- function(seq_pool, seq, score, motif_pool, motif) {
  j1 <- which(names(seq_pool) == names(seq))
  j2 <- which(score$motif == names(motif))
  insertion <- seq[[1]]$insertion
  m <- seq[[1]]$insertion_len
  #long seq
  long_best_match <- score$match_pos_long[j1, j2]
  long_strand_sign <- sign(long_best_match)
  seque <- seq[[1]]$inserted_sequence
  n <- length(seque)
  n1 <- (n - m) / 2
  lm <- nrow(motif[[1]])
  s1 <- n1 - lm + 2
  e1 <- n1 + m + lm - 1
  longer <- seque[s1:e1]
  nlong <- length(longer)
  #short seq
  short_best_match <- score$match_pos_short[j1, j2]
  short_strand_sign <- sign(short_best_match)
  s2 <- n1
  e2 <- n1 + m + 1
  shorter <- seque[c(s1:s2, e2:e1)]
  nshort <- length(shorter)
  # plot_l<-05*(nl+2)

  plot.new()
  if (insertion) {
    if (long_strand_sign == 1) {
      left <- long_best_match - s1
      right <- nlong - left - lm
      a <- makepwm(longer, long_strand_sign)
      rownames(a) <- c("A", "C", "G", "T")
      p <- makemotif(motif, long_strand_sign)
      p_full <- makepfull(p, left, right)
      # pushViewport(viewport(
      #   y = unit(.5, "npc") - unit(2, "lines"),
      #   height = unit(1, "npc") - unit(3, "lines")
      # ))
      pushViewport(viewport(y = .125, height = .25))
      plotMotifLogo(
        p_full,

        yaxis = FALSE,
        xaxis = FALSE,
        xlab = "",
        ylab = "",
        newpage = FALSE,
        margins = c(2, 2.5, 1.5, 2.5)
      )
      #add 5'-3' notation
      grid.lines(
        x = c(
          convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
          1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
        ),
        y = unit(1, "lines"),
        gp = gpar(
          col = "blue",
          lwd = 1.5,
          xpd = NA
        ),
        arrow = arrow(
          length = unit(0.1, "inches"),
          angle = 15,
          ends = "last"
        )
      )
      grid.text(
        "3'",
        x = unit(1, "npc") - unit(1, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(.5, "lines")
      )
      grid.text(
        "5'",
        x = unit(2, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(.5, "lines")
      )
      # plotMotifLogo(a)
      popViewport()
      pushViewport(viewport(y = .325, height = .25))
      #par(mar = c(4, 3, 1.5, 2))
      plotMotifLogo(
        a,
        #paste("insertion=",insertion),
        xaxis = FALSE,
        yaxis = FALSE,
        xlab = "",
        ylab = "(+)",
        newpage = FALSE,
        margins = c(2, 2, 1.5, 2)
      )
      grid.text(
        "3'",
        x = unit(1, "npc") - unit(1, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(.5, "lines")
      )
      grid.text(
        "5'",
        x = unit(2, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(.5, "lines")
      )
      grid.rect(
        just = "center",
        width = (m) / (nlong + 2),
        height = unit(0.6, "npc"),
        gp = gpar(
          col = "blue",
          lty = 3,
          lwd = 2,
          fill = NA
        )
      )

      if (short_strand_sign == 1) {
        a <- makepwm(shorter, short_strand_sign)
        rownames(a) <- c("A", "C", "G", "T")
        l1 <- floor(m / 2)
        r1 <- m - l1
        left <- short_best_match - s1 + l1
        right <- nlong - left - lm
        if (l1 == 0) {
          a_full <- cbind(a, rep.col(c(0.25, 0.25, 0.25, 0.25), r1))
        } else{
          a_full <-
            cbind(rep.col(c(0.25, 0.25, 0.25, 0.25), l1), a, rep.col(c(0.25, 0.25, 0.25, 0.25), r1))
        }
        p <- makemotif(motif, short_strand_sign)
        p_full <- makepfull(p, left, right)
        popViewport()
        pushViewport(viewport(y = .575, height = .25))
        plotMotifLogo(
          a_full,
          xaxis = FALSE,
          yaxis = FALSE,
          xlab = "",
          ylab = "(+)",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "3'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        grid.text(
          "5'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        #grid.rect(x=((nlong-m)/2-0.25)/nlong, width = m/nlong, gp=gpar(col="blue", lty=3, lwd=2, fill=NA))
        popViewport()
        pushViewport(viewport(y = .825, height = .25))
        plotMotifLogo(
          p_full,
          paste(names(motif), "binding change by InDel", names(seq)),
          font = "Helvetica-Bold",
          ncex = 1.2,
          yaxis = FALSE,
          xaxis = FALSE,
          xlab = "",
          ylab = "",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.lines(
          x = c(
            convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
            1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
          ),
          y = unit(1, "lines"),
          gp = gpar(
            col = "blue",
            lwd = 1.5,
            xpd = NA
          ),
          arrow = arrow(
            length = unit(0.1, "inches"),
            angle = 15,
            ends = "last"
          )
        )
        grid.text(
          "3'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        grid.text(
          "5'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        popViewport()
        pushViewport(viewport(y = .375, height = .25))
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm) / (nlong + 2)
          ), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm + m) / (nlong + 2)
          ), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm + l1) / (nlong + 2)
          ), "npc"),
          y = unit(c(2, 1.2), "npc"),
          gp = gpar(fill = "black")
        )
      } else{
        a <- makepwm(shorter, short_strand_sign)
        rownames(a) <- c("A", "C", "G", "T")
        l1 <- floor(m / 2)
        r1 <- m - l1
        right <- (nlong - m) / 2 - (short_best_match + n1 + 1) + m - r1
        left <- nlong - lm - right
        if (l1 == 0) {
          a_full <- cbind(a, rep.col(c(0.25, 0.25, 0.25, 0.25), r1))
        } else{
          a_full <-
            cbind(rep.col(c(0.25, 0.25, 0.25, 0.25), l1), a, rep.col(c(0.25, 0.25, 0.25, 0.25), r1))
        }
        p <- makemotif(motif, short_strand_sign)
        p_full <- makepfull(p, left, right)
        title <-
          ifelse(
            insertion,
            "Best match to the reference sequence m=",
            "Best match to the mutation sequence m="
          )
        popViewport()
        pushViewport(viewport(y = .575, height = .25))
        plotMotifLogo(
          a_full,
          xaxis = FALSE,
          yaxis = FALSE,
          xlab = "",
          ylab = "(-)",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "5'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        grid.text(
          "3'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        popViewport()
        pushViewport(viewport(y = .825, height = .25))
        plotMotifLogo(
          p_full,
          paste(names(motif), "binding change by InDel", names(seq)),
          font = "Helvetica-Bold",
          ncex = 1.2,
          yaxis = FALSE,
          xaxis = FALSE,
          xlab = "",
          ylab = "",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "5'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        grid.text(
          "3'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )

        grid.lines(
          x = c(
            convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
            1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
          ),
          y = unit(1, "lines"),
          gp = gpar(
            col = "blue",
            lwd = 1.5,
            xpd = NA
          ),
          arrow = arrow(
            length = unit(0.1, "inches"),
            angle = 15,
            ends = "first"
          )
        )

        popViewport()
        pushViewport(viewport(y = .325, height = .25))
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm) / (nlong + 2)
          ), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm + m) / (nlong + 2)
          ), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm + l1) / (nlong + 2)
          ), "npc"),
          y = unit(c(2, 1.2), "npc"),
          gp = gpar(fill = "black")
        )
      }
    } else if (long_strand_sign == -1) {
      #best match is negative
      right <- (nlong - m) / 2 - (long_best_match + n1 + 1)
      left <- nlong - lm - right
      a <- makepwm(longer, long_strand_sign)
      rownames(a) <- c("A", "C", "G", "T")
      p <- makemotif(motif, long_strand_sign)
      p_full <- makepfull(p, left, right)
      title <-
        ifelse(
          insertion,
          "Best match to the mutation sequence m=",
          "Best match to the reference sequence m="
        )
      # pushViewport(viewport(
      #   y = unit(.5, "npc") - unit(2, "lines"),
      #   height = unit(1, "npc") - unit(3, "lines")
      # ))
      pushViewport(viewport(y = .125, height = .25))
      plotMotifLogo(
        p_full,
        yaxis = FALSE,
        xaxis = FALSE,
        xlab = "",
        ylab = "",
        newpage = FALSE,
        margins = c(2, 2.5, 1.5, 2.5)
      )
      grid.text(
        "5'",
        x = unit(1, "npc") - unit(1, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(2.5, "lines")
      )
      grid.text(
        "3'",
        x = unit(2, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(2.5, "lines")
      )
      grid.lines(
        x = c(
          convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
          1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
        ),
        y = unit(1, "lines"),
        gp = gpar(
          col = "blue",
          lwd = 1.5,
          xpd = NA
        ),
        arrow = arrow(
          length = unit(0.1, "inches"),
          angle = 15,
          ends = "first"
        )
      )
      # plotMotifLogo(a)
      popViewport()
      pushViewport(viewport(y = .325, height = .25))
      #par(mar = c(4, 3, 1.5, 2))
      plotMotifLogo(
        a,
        xaxis = FALSE,
        yaxis = FALSE,
        xlab = "",
        ylab = "(-)",
        newpage = FALSE,
        margins = c(2, 2.5, 1.5, 2.5)
      )
      grid.text(
        "5'",
        x = unit(1, "npc") - unit(1, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(2.5, "lines")
      )
      grid.text(
        "3'",
        x = unit(2, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(2.5, "lines")
      )
      grid.rect(
        just = "center",
        width = (m) / (nlong + 2),
        height = unit(0.6, "npc"),
        gp = gpar(
          col = "blue",
          lty = 3,
          lwd = 2,
          fill = NA
        )
      )
      #grid.rect(x=(0.25+0.5*m+(nlong-m)/2)/nlong, width = m/nlong, gp=gpar(col="blue", lty=3, lwd=2, fill=NA))
      if (short_strand_sign == 1) {
        a <- makepwm(shorter, short_strand_sign)
        rownames(a) <- c("A", "C", "G", "T")
        l1 <- floor(m / 2)
        r1 <- m - l1
        left <- short_best_match - s1 + l1
        right <- nlong - left - lm
        if (l1 == 0) {
          a_full <- cbind(a, rep.col(c(0.25, 0.25, 0.25, 0.25), r1))
        } else{
          a_full <-
            cbind(rep.col(c(0.25, 0.25, 0.25, 0.25), l1), a, rep.col(c(0.25, 0.25, 0.25, 0.25), r1))
        }
        p <- makemotif(motif, short_strand_sign)
        p_full <- makepfull(p, left, right)
        popViewport()
        pushViewport(viewport(y = .575, height = .25))
        plotMotifLogo(
          a_full,
          xaxis = FALSE,
          yaxis = FALSE,
          xlab = "",
          ylab = "(+)",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "3'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        grid.text(
          "5'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        popViewport()
        pushViewport(viewport(y = .825, height = .25))
        plotMotifLogo(
          p_full,
          paste(names(motif), "binding change by InDel", names(seq)),
          font = "Helvetica-Bold",
          ncex = 1.2,
          yaxis = FALSE,
          xaxis = FALSE,
          xlab = "",
          ylab = "",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.lines(
          x = c(
            convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
            1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
          ),
          y = unit(1, "lines"),
          gp = gpar(
            col = "blue",
            lwd = 1.5,
            xpd = NA
          ),
          arrow = arrow(
            length = unit(0.1, "inches"),
            angle = 15,
            ends = "last"
          )
        )
        grid.text(
          "3'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        grid.text(
          "5'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )

        popViewport()
        pushViewport(viewport(y = .325, height = .25))
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm) / (nlong + 2)
          ), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm + m) / (nlong + 2)
          ), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm + l1) / (nlong + 2)
          ), "npc"),
          y = unit(c(2, 1.2), "npc"),
          gp = gpar(fill = "black")
        )

      } else{
        a <- makepwm(shorter, short_strand_sign)
        rownames(a) <- c("A", "C", "G", "T")
        l1 <- floor(m / 2)
        r1 <- m - l1
        right <- (nlong - m) / 2 - (short_best_match + n1 + 1) + m - r1
        left <- nlong - lm - right
        if (l1 == 0) {
          a_full <- cbind(a, rep.col(c(0.25, 0.25, 0.25, 0.25), r1))
        } else{
          a_full <-
            cbind(rep.col(c(0.25, 0.25, 0.25, 0.25), l1), a, rep.col(c(0.25, 0.25, 0.25, 0.25), r1))
        }
        p <- makemotif(motif, short_strand_sign)
        p_full <- makepfull(p, left, right)
        title <-
          ifelse(
            insertion,
            "Best match to the reference sequence m=",
            "Best match to the mutation sequence m="
          )
        popViewport()
        pushViewport(viewport(y = .575, height = .25))
        plotMotifLogo(
          a_full,
          xaxis = FALSE,
          yaxis = FALSE,
          xlab = "",
          ylab = "(-)",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "5'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        grid.text(
          "3'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        popViewport()
        pushViewport(viewport(y = .825, height = .25))
        #par(mar = c(4, 3, 1.5, 2))
        plotMotifLogo(
          p_full,
          paste(names(motif), "binding change by InDel", names(seq)),
          font = "Helvetica-Bold",
          ncex = 1.2,
          yaxis = FALSE,
          xaxis = FALSE,
          xlab = "",
          ylab = "",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "5'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        grid.text(
          "3'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        grid.lines(
          x = c(
            convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
            1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
          ),
          y = unit(1, "lines"),
          gp = gpar(
            col = "blue",
            lwd = 1.5,
            xpd = NA
          ),
          arrow = arrow(
            length = unit(0.1, "inches"),
            angle = 15,
            ends = "first"
          )
        )
        popViewport()
        pushViewport(viewport(y = .325, height = .25))
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm) / (nlong + 2)
          ), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm + m) / (nlong + 2)
          ), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((l1 + lm) / (nlong + 2), (lm + l1) / (nlong + 2)
          ), "npc"),
          y = unit(c(2, 1.2), "npc"),
          gp = gpar(fill = "black")
        )

      }
    }
  } else if (insertion == 0) {
    if (long_strand_sign == 1) {
      left <- long_best_match - s1
      right <- nlong - left - lm
      a <- makepwm(longer, long_strand_sign)
      rownames(a) <- c("A", "C", "G", "T")
      p <- makemotif(motif, long_strand_sign)
      p_full <- makepfull(p, left, right)
      ##
      pushViewport(viewport(y = .825, height = .25))
      plotMotifLogo(
        p_full,
        paste(names(motif), "binding change by InDel", names(seq)),
        font = "Helvetica-Bold",
        ncex = 1.2,
        yaxis = FALSE,
        xaxis = FALSE,
        xlab = "",
        ylab = "",
        newpage = FALSE,
        margins = c(2, 2.5, 1.5, 2.5)
      )
      #add 5'-3' notation
      grid.lines(
        x = c(
          convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
          1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
        ),
        y = unit(1, "lines"),
        gp = gpar(
          col = "blue",
          lwd = 1.5,
          xpd = NA
        ),
        arrow = arrow(
          length = unit(0.1, "inches"),
          angle = 15,
          ends = "last"
        )
      )
      grid.text(
        "3'",
        x = unit(1, "npc") - unit(1, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(.5, "lines")
      )
      grid.text(
        "5'",
        x = unit(2, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(.5, "lines")
      )
      # plotMotifLogo(a)
      popViewport()
      pushViewport(viewport(y = .575, height = .25))
      #par(mar = c(4, 3, 1.5, 2))
      plotMotifLogo(
        a,
        #paste("insertion=",insertion),
        xaxis = FALSE,
        yaxis = FALSE,
        xlab = "",
        ylab = "(+)",
        newpage = FALSE,
        margins = c(2, 2, 1.5, 2)
      )
      grid.text(
        "3'",
        x = unit(1, "npc") - unit(1, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(.5, "lines")
      )
      grid.text(
        "5'",
        x = unit(2, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(.5, "lines")
      )
      grid.rect(
        just = "center",
        width = (m) / (nlong + 2),
        height = unit(0.6, "npc"),
        gp = gpar(
          col = "blue",
          lty = 3,
          lwd = 2,
          fill = NA
        )
      )

      #grid.rect(x=(0.5*(m)+(nlong-m)/2)/nlong, width = (m)/nlong, gp=gpar(col="blue", lty=3, lwd=2, fill=NA))
      if (short_strand_sign == 1) {
        left <- short_best_match - s1
        right <- nlong - left - lm
        a <- makepwm(shorter, short_strand_sign)
        rownames(a) <- c("A", "C", "G", "T")
        a_full <-
          cbind(a[, 1:(lm - 1)], rep.col(c(0.25, 0.25, 0.25, 0.25), m), a[, lm:(2 *
                                                                                  lm - 2)])
        #a_full<-cbind(a,rep.col(c(0.25,0.25,0.25,0.25),m))
        p <- makemotif(motif, short_strand_sign)
        p_full <- makepfull(p, left, right)
        #title<-ifelse(insertion==1,"Best match to the reference sequence m=","Best match to the mutation sequence m=")
        popViewport()
        pushViewport(viewport(y = .325, height = .25))
        plotMotifLogo(
          a_full,
          xaxis = FALSE,
          yaxis = FALSE,
          xlab = "",
          ylab = "(+)",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "3'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        grid.text(
          "5'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        #grid.rect(x=((nlong-m)/2-0.25)/nlong, width = m/nlong, gp=gpar(col="blue", lty=3, lwd=2, fill=NA))
        popViewport()
        pushViewport(viewport(y = .125, height = .25))
        plotMotifLogo(
          p_full,
          yaxis = FALSE,
          xaxis = FALSE,
          xlab = "",
          ylab = "",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.lines(
          x = c(
            convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
            1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
          ),
          y = unit(1, "lines"),
          gp = gpar(
            col = "blue",
            lwd = 1.5,
            xpd = NA
          ),
          arrow = arrow(
            length = unit(0.1, "inches"),
            angle = 15,
            ends = "last"
          )
        )
        grid.text(
          "3'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        grid.text(
          "5'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        popViewport()
        pushViewport(viewport(y = .125, height = .25))
        grid.lines(
          x = unit(c(lm / (nlong + 2), lm / (nlong + 2)), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((lm + m) / (nlong + 2), lm / (nlong + 2)), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )

      } else{
        right <- (nlong - m) / 2 - (short_best_match + n1 + 1) + m
        left <- nlong - lm - right

        a <- makepwm(shorter, short_strand_sign)
        rownames(a) <- c("A", "C", "G", "T")
        a_full <-
          cbind(a[, 1:(lm - 1)], rep.col(c(0.25, 0.25, 0.25, 0.25), m), a[, lm:(2 *
                                                                                  lm - 2)])
        p <- makemotif(motif, short_strand_sign)
        p_full <- makepfull(p, left, right)
        title <-
          ifelse(
            insertion,
            "Best match to the reference sequence m=",
            "Best match to the mutation sequence m="
          )
        popViewport()
        pushViewport(viewport(y = .325, height = .25))
        plotMotifLogo(
          a_full,
          xaxis = FALSE,
          yaxis = FALSE,
          xlab = "",
          ylab = "(-)",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "5'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        grid.text(
          "3'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        popViewport()
        pushViewport(viewport(y = .125, height = .25))
        #par(mar = c(4, 3, 1.5, 2))
        plotMotifLogo(
          p_full,
          yaxis = FALSE,
          xaxis = FALSE,
          xlab = "",
          ylab = "",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "5'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        grid.text(
          "3'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )

        grid.lines(
          x = c(
            convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
            1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
          ),
          y = unit(1, "lines"),
          gp = gpar(
            col = "blue",
            lwd = 1.5,
            xpd = NA
          ),
          arrow = arrow(
            length = unit(0.1, "inches"),
            angle = 15,
            ends = "first"
          )
        )

        popViewport()
        pushViewport(viewport(y = .125, height = .25))
        grid.lines(
          x = unit(c(lm / (nlong + 2), lm / (nlong + 2)), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((lm + m) / (nlong + 2), lm / (nlong + 2)), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
      }
    } else if (long_strand_sign == -1) {
      #best match is negative
      right <- (nlong - m) / 2 - (long_best_match + n1 + 1)
      left <- nlong - lm - right
      a <- makepwm(longer, long_strand_sign)
      rownames(a) <- c("A", "C", "G", "T")
      p <- makemotif(motif, long_strand_sign)
      p_full <- makepfull(p, left, right)
      title <-
        ifelse(
          insertion,
          "Best match to the mutation sequence m=",
          "Best match to the reference sequence m="
        )
      pushViewport(viewport(y = .825, height = .25))
      plotMotifLogo(
        p_full,
        paste(names(motif), "binding change by InDel", names(seq)),
        font = "Helvetica-Bold",
        ncex = 1.2,
        yaxis = FALSE,
        xaxis = FALSE,
        xlab = "",
        ylab = "",
        newpage = FALSE,
        margins = c(2, 2.5, 1.5, 2.5)
      )
      grid.text(
        "5'",
        x = unit(1, "npc") - unit(1, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(2.5, "lines")
      )
      grid.text(
        "3'",
        x = unit(2, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(2.5, "lines")
      )
      grid.lines(
        x = c(
          convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
          1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
        ),
        y = unit(1, "lines"),
        gp = gpar(
          col = "blue",
          lwd = 1.5,
          xpd = NA
        ),
        arrow = arrow(
          length = unit(0.1, "inches"),
          angle = 15,
          ends = "first"
        )
      )
      # plotMotifLogo(a)
      popViewport()
      pushViewport(viewport(y = .575, height = .25))
      #par(mar = c(4, 3, 1.5, 2))
      plotMotifLogo(
        a,
        font = "mono,Courier",
        xaxis = FALSE,
        yaxis = FALSE,
        xlab = "",
        ylab = "(-)",
        newpage = FALSE,
        margins = c(2, 2.5, 1.5, 2.5)
      )
      grid.text(
        "5'",
        x = unit(1, "npc") - unit(1, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(2.5, "lines")
      )
      grid.text(
        "3'",
        x = unit(2, "lines"),
        gp = gpar(col = "blue", cex = 1),
        y = unit(2.5, "lines")
      )
      grid.rect(
        just = "center",
        width = (m) / (nlong + 2),
        height = unit(0.6, "npc"),
        gp = gpar(
          col = "blue",
          lty = 3,
          lwd = 2,
          fill = NA
        )
      )
      #grid.rect(x=(0.25+0.5*m+(nlong-m)/2)/nlong, width = m/nlong, gp=gpar(col="blue", lty=3, lwd=2, fill=NA))
      if (short_strand_sign == 1) {
        left <- short_best_match - s1
        right <- nlong - left - lm
        a <- makepwm(shorter, short_strand_sign)
        rownames(a) <- c("A", "C", "G", "T")
        a_full <-
          cbind(a[, 1:(lm - 1)], rep.col(c(0.25, 0.25, 0.25, 0.25), m), a[, lm:(2 *
                                                                                  lm - 2)])
        p <- makemotif(motif, short_strand_sign)
        p_full <-
          cbind(rep.col(c(0.25, 0.25, 0.25, 0.25), left), p, rep.col(c(0.25, 0.25, 0.25, 0.25), right))
        title <-
          ifelse(
            insertion,
            "Best match to the reference sequence m=",
            "Best match to the mutation sequence m="
          )
        popViewport()
        pushViewport(viewport(y = .325, height = .25))
        plotMotifLogo(
          a_full,
          font = "mono,Courier",
          xaxis = FALSE,
          yaxis = FALSE,
          xlab = "",
          ylab = "(+)",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "3'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        grid.text(
          "5'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        # grid.rect(x=((nlong-m)/2)/nlong, width = m/nlong, gp=gpar(col="blue", lty=3, lwd=2, fill=NA))
        popViewport()
        pushViewport(viewport(y = .125, height = .25))
        plotMotifLogo(
          p_full,
          yaxis = FALSE,
          xaxis = FALSE,
          xlab = "",
          ylab = "",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.lines(
          x = c(
            convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
            1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
          ),
          y = unit(1, "lines"),
          gp = gpar(
            col = "blue",
            lwd = 1.5,
            xpd = NA
          ),
          arrow = arrow(
            length = unit(0.1, "inches"),
            angle = 15,
            ends = "last"
          )
        )
        grid.text(
          "3'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )
        grid.text(
          "5'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(.5, "lines")
        )

        popViewport()
        pushViewport(viewport(y = .125, height = .25))
        grid.lines(
          x = unit(c(lm / (nlong + 2), lm / (nlong + 2)), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((lm + m) / (nlong + 2), lm / (nlong + 2)), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )

      } else{
        right <- (nlong - m) / 2 - (short_best_match + n1 + 1) + m
        left <- nlong - lm - right

        a <- makepwm(shorter, short_strand_sign)
        rownames(a) <- c("A", "C", "G", "T")
        a_full <-
          cbind(a[, 1:(lm - 1)], rep.col(c(0.25, 0.25, 0.25, 0.25), m), a[, lm:(2 *
                                                                                  lm - 2)])
        p <- makemotif(motif, short_strand_sign)
        p_full <- makepfull(p, left, right)
        title <-
          ifelse(
            insertion,
            "Best match to the reference sequence m=",
            "Best match to the mutation sequence m="
          )
        popViewport()
        pushViewport(viewport(y = .375, height = .25))
        plotMotifLogo(
          a_full,
          xaxis = FALSE,
          yaxis = FALSE,
          xlab = "",
          ylab = "(-)",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "5'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        grid.text(
          "3'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        popViewport()
        pushViewport(viewport(y = .125, height = .25))
        #par(mar = c(4, 3, 1.5, 2))
        plotMotifLogo(
          p_full,
          yaxis = FALSE,
          xaxis = FALSE,
          xlab = "",
          ylab = "",
          newpage = FALSE,
          margins = c(2, 2.5, 1.5, 2.5)
        )
        grid.text(
          "5'",
          x = unit(1, "npc") - unit(1, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        grid.text(
          "3'",
          x = unit(2, "lines"),
          gp = gpar(col = "blue", cex = 1),
          y = unit(2.5, "lines")
        )
        grid.lines(
          x = c(
            convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE),
            1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)
          ),
          y = unit(1, "lines"),
          gp = gpar(
            col = "blue",
            lwd = 1.5,
            xpd = NA
          ),
          arrow = arrow(
            length = unit(0.1, "inches"),
            angle = 15,
            ends = "first"
          )
        )
        popViewport()
        pushViewport(viewport(y = .125, height = .25))
        grid.lines(
          x = unit(c(lm / (nlong + 2), lm / (nlong + 2)), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
        grid.lines(
          x = unit(c((lm + m) / (nlong + 2), lm / (nlong + 2)), "npc"),
          y = unit(c(1.2, 0.75), "npc"),
          gp = gpar(fill = "black"),
          arrow = arrow(
            length = unit(0.05, "inches"),
            ends = "last",
            type = "closed"
          )
        )
      }
    }
  }

}
