# ============================================================================
# Hand-vectorised right-hand side of models/mrna_vaccine_dasti2025.cpp
#
# The model file writes the 17 B-cell affinity classes out one by one, so its
# translated R code runs a few thousand scalar operations per call. This is
# the same $ODE with the classes as the rows of 17-row matrices and all index
# and name lookups done once, which is several times faster in the browser.
# tests/test_engine_vs_mrgsolve.R checks it against the translated model file
# at random states and the solution against mrgsolve.
#
# Derivative of the COSBI "Multiscale QSP model for mRNA vaccines";
# COSBI-SSLA licence (non-commercial), see ../LICENSE. Modified 2026-09-28.
# ============================================================================

mrna_rhs <- function(model) {
  s <- stats::setNames(seq_along(model$cmt), model$cmt)
  cls <- function(p) unname(s[paste0(p, 1:17)])
  i <- list(NB = cls("NB"), AB_N = cls("AB_N"), MB = cls("MB"), AB_M = cls("AB_M"), SP = cls("SP"),
            LP = cls("LP"), SP_BL = cls("SP_BL"), LP_BL = cls("LP_BL"), Ab_BL = cls("Ab_BL"), GCB = cls("GCB"))
  ka <- 1e-6 * 2^((1:17) - 9)
  nb0 <- stats::dnorm(-3 + 6 / 16 * (0:16))
  nb0 <- nb0 / sum(nb0) * 1.3e9 * 0.000004
  dc_is <- c("DC_IS", "DCL_IS", "DCAgLon_IS", "DCAgMon_IS", "DCAgH_IS", "DCAgMoff_IS", "DCAgLoff_IS", "DCoff_IS")
  dc_ln <- c("DCL_LN", "DCAgLon_LN", "DCAgMon_LN", "DCAgH_LN", "DCAgMoff_LN", "DCAgLoff_LN", "DC_LN")
  # neutrophils and monocytes: indices and parameter names
  inn <- lapply(c(NP = "NP", MN = "MN"), function(cc) list(
    is = s[[paste0(cc, "_IS")]], l = s[[paste0(cc, "L_IS")]], ag = s[[paste0(cc, "Ag_IS")]],
    lln = s[[paste0(cc, "L_LN")]], agln = s[[paste0(cc, "Ag_LN")]], bl = s[[paste0(cc, "_BL")]],
    beta = paste0("Beta", cc), gamma = paste0("Gamma", cc), delta = paste0("Delta", cc),
    alpha = paste0("Alpha", cc, "_IS"), alpha_bl = paste0("Alpha", cc, "_BL"), eta = paste0("Eta", cc),
    omega = paste0("Omega", cc), mu = paste0("Mu", cc), xi = paste0("Xi", cc)))
  # myeloid and plasmacytoid dendritic cells
  dcx <- lapply(c(m = "m", p = "p"), function(x0) {
    X <- paste0(x0, "DC")
    list(v = unname(s[paste0(x0, dc_is)]), u = unname(s[paste0(x0, dc_ln)]), bl = s[[paste0(X, "_BL")]],
         L = unname(s[paste0(x0, c("DCAgLon_LN", "DCAgLoff_LN"))]),
         M = unname(s[paste0(x0, c("DCAgMon_LN", "DCAgMoff_LN"))]), H = s[[paste0(x0, "DCAgH_LN")]],
         nmax = paste0("NmaxMHC_", X), beta = paste0("Beta", X), xi = paste0("Xi", X),
         delta = paste0("Delta", X), gamma = paste0("Gamma", X), mu = paste0("Mu", X), eta = paste0("Eta", X),
         omega = paste0("Omega", X), alpha = paste0("Alpha", x0, "IDC_IS"), alpha_bl = paste0("Alpha", x0, "IDC_BL"),
         betai = paste0("Beta", x0, "IDC"), trM = paste0("k_trM_", X), trH = paste0("k_trH_", X),
         atrH = paste0("k_atrH_", X), atrM = paste0("k_atrM_", X), atrL = paste0("k_atrL_", X))
  })
  i_m <- s[["mRNA"]]; i_nt <- s[["NT"]]; i_atn <- s[["AT_N"]]; i_mt <- s[["MT"]]; i_atm <- s[["AT_M"]]
  i_ft <- s[["FT"]]; i_npis <- s[["NP_IS"]]; i_mnis <- s[["MN_IS"]]
  i_mdcis <- s[["mDC_IS"]]; i_pdcis <- s[["pDC_IS"]]
  w <- function(a, b, c) (abs(a - b) / 2 + a) * c

  function(t, Y, V) {
    p <- V
    m <- ncol(Y)
    out <- matrix(0, nrow(Y), m)
    mRNA <- Y[i_m, ]
    NT <- Y[i_nt, ]; AT_N <- Y[i_atn, ]; MT <- Y[i_mt, ]; AT_M <- Y[i_atm, ]; FT <- Y[i_ft, ]

    # dendritic cells in the lymph node: maturation and T-cell stimulation
    dcs <- lapply(dcx, function(d) {
      L <- Y[d$L[1], ] + Y[d$L[2], ]; M <- Y[d$M[1], ] + Y[d$M[2], ]; H <- Y[d$H, ]
      mat0 <- L * p$p_L + M * p$p_M + H * p$p_H
      tot <- L + M + H
      mat <- ifelse(mat0 < p$p_L, 0, mat0)
      nmed <- ifelse(mat < p$p_L, 0, mat * p[[d$nmax]] / tot)
      st <- tot / (tot + NT + AT_N + AT_M + MT)
      list(L = L, M = M, H = H, tot = tot,
           DN = st * nmed / (nmed + p$KNT), DM = st * nmed / (nmed + p$KMT),
           EN = st * (nmed - p$KNT) / (nmed + p$KNT), EM = st * (nmed - p$KMT) / (nmed + p$KMT))
    })
    dm <- dcs$m; dp <- dcs$p
    dsum <- dm$tot + dp$tot
    wm <- ifelse(dsum > 0, dm$tot / dsum, 0)
    wp <- ifelse(dsum > 0, dp$tot / dsum, 0)
    D_N <- dm$DN * wm + dp$DN * wp
    D_M <- dm$DM * wm + dp$DM * wp
    E_N <- dm$EN * wm + dp$EN * wp
    E_M <- dm$EM * wm + dp$EM * wp
    tlast <- ifelse(t >= p$TD3, p$TD3, ifelse(t >= p$TD2, p$TD2, 0))
    early <- t - tlast < 0.1
    E_N <- ifelse(early, pmax(E_N, 0), E_N)
    E_M <- ifelse(early, pmax(E_M, 0), E_M)

    # B-cell receptor occupancy by free antigen
    NB <- Y[i$NB, , drop = FALSE]; AB_N <- Y[i$AB_N, , drop = FALSE]; MB <- Y[i$MB, , drop = FALSE]
    AB_M <- Y[i$AB_M, , drop = FALSE]; SP <- Y[i$SP, , drop = FALSE]; LP <- Y[i$LP, , drop = FALSE]
    GCB <- Y[i$GCB, , drop = FALSE]
    Bcells <- NB + AB_N + GCB + AB_M + MB
    Ag <- p$expAg * (p$maxMHC_mDC * (w(p$p_L, p$p_M, dm$L) + w(p$p_M, p$p_H, dm$M) + w(p$p_H, 1, dm$H)) +
                     p$maxMHC_pDC * (w(p$p_L, p$p_M, dp$L) + w(p$p_M, p$p_H, dp$M) + w(p$p_H, 1, dp$H)))
    Agc <- Ag / p$V_LN
    cB <- ka * Bcells * rep(p$BRN / p$NAV * 1e12 / p$V_LN, each = 17)
    x <- Agc / (1 + colSums(cB))
    for (k in 1:8) {
      q <- 1 + outer(ka, x)
      x <- x - (x * (1 + colSums(cB / q)) - Agc) / (1 + colSums(cB / q^2))
    }
    kx <- outer(ka, x)
    ro <- kx / (1 + kx)
    R <- ro * rep(p$BRN, each = 17)
    Fb <- R / (p$K_R + R)
    G <- (1 - ro) * Fb
    Hb <- (R - p$K_R) / (R + p$K_R)
    Bsum <- colSums(Bcells)
    P_N <- p$CC_N * FT / (p$CC_N * FT + Bsum)
    P_M <- p$CC_M * FT / (p$CC_M * FT + Bsum)
    upt <- mRNA / (p$K_mRNA + mRNA)
    sfm <- p$SF_Gamma * mRNA

    # injection site: mRNA
    out[i_m, ] <- -p$GammapDC * mRNA * Y[i_pdcis, ] - p$GammamDC * mRNA * Y[i_mdcis, ] -
      p$GammaMN * mRNA * Y[i_mnis, ] - p$GammaNP * mRNA * Y[i_npis, ] - p$kdmrna * mRNA - p$kdeg * mRNA -
      p$ksat * mRNA * 0.5 * (1 + tanh((mRNA - p$mRNA_max) / p$k_slope))
    # neutrophils and monocytes: injection site, lymph node, blood
    for (d in inn) {
      beta <- p[[d$beta]]; g <- p[[d$gamma]]; de <- p[[d$delta]]; mu <- p[[d$mu]]; xi <- p[[d$xi]]
      eta <- p[[d$eta]]; om <- p[[d$omega]]
      is <- Y[d$is, ]; l <- Y[d$l, ]; ag <- Y[d$ag, ]; lln <- Y[d$lln, ]; bl <- Y[d$bl, ]
      out[d$is, ] <- p[[d$alpha]] + eta * upt * bl - g * sfm * is - beta * is - om * is
      out[d$l, ] <- g * sfm * is - de * l - beta * l - mu * l
      out[d$ag, ] <- de * l - beta * ag - xi * ag
      out[d$lln, ] <- mu * l - de * lln - beta * lln
      out[d$agln, ] <- xi * ag + de * lln - beta * Y[d$agln, ]
      out[d$bl, ] <- p[[d$alpha_bl]] - eta * upt * bl - beta * bl + om * is
    }
    # dendritic cells: injection site, lymph node, blood
    for (d in dcx) {
      beta <- p[[d$beta]]; xi <- p[[d$xi]]; de <- p[[d$delta]]; mu <- p[[d$mu]]; g <- p[[d$gamma]]
      trM <- p[[d$trM]]; trH <- p[[d$trH]]; atrH <- p[[d$atrH]]; atrM <- p[[d$atrM]]; atrL <- p[[d$atrL]]
      eta <- p[[d$eta]]; om <- p[[d$omega]]; bi <- p[[d$betai]]
      v1 <- Y[d$v[1], ]; v2 <- Y[d$v[2], ]; v3 <- Y[d$v[3], ]; v4 <- Y[d$v[4], ]
      v5 <- Y[d$v[5], ]; v6 <- Y[d$v[6], ]; v7 <- Y[d$v[7], ]; v8 <- Y[d$v[8], ]
      u1 <- Y[d$u[1], ]; u2 <- Y[d$u[2], ]; u3 <- Y[d$u[3], ]; u4 <- Y[d$u[4], ]
      u5 <- Y[d$u[5], ]; u6 <- Y[d$u[6], ]; u7 <- Y[d$u[7], ]; bl <- Y[d$bl, ]
      out[d$v[1], ] <- p[[d$alpha]] + eta * upt * bl - g * sfm * v1 - bi * v1 - om * v1
      out[d$v[2], ] <- g * sfm * v1 - de * v2 - beta * v2 - mu * v2
      out[d$v[3], ] <- de * v2 - trM * v3 - beta * v3 - xi * v3
      out[d$v[4], ] <- trM * v3 - trH * v4 - beta * v4 - xi * v4
      out[d$v[5], ] <- trH * v4 - atrH * v5 - beta * v5 - xi * v5
      out[d$v[6], ] <- atrH * v5 - atrM * v6 - beta * v6 - xi * v6
      out[d$v[7], ] <- atrM * v6 - atrL * v7 - beta * v7 - xi * v7
      out[d$v[8], ] <- atrL * v7 - beta * v8 - xi * v8
      out[d$u[1], ] <- mu * v2 - de * u1 - beta * u1
      out[d$u[2], ] <- xi * v3 + de * u1 - trM * u2 - beta * u2
      out[d$u[3], ] <- xi * v4 + trM * u2 - trH * u3 - beta * u3
      out[d$u[4], ] <- xi * v5 + trH * u3 - atrH * u4 - beta * u4
      out[d$u[5], ] <- xi * v6 + atrH * u4 - atrM * u5 - beta * u5
      out[d$u[6], ] <- xi * v7 + atrM * u5 - atrL * u6 - beta * u6
      out[d$u[7], ] <- xi * v8 + atrL * u6 - beta * u7
      out[d$bl, ] <- p[[d$alpha_bl]] - eta * upt * bl - bi * bl + om * v1
    }

    # helper T cells
    out[i_nt, ] <- p$BetaNT * (p$NT0 - NT) - p$SigmaNT * D_N * NT
    out[i_atn, ] <- p$SigmaNT * D_N * NT + p$RhoAT * E_N * AT_N - p$BetaAT * AT_N
    out[i_mt, ] <- p$RhoAT * (1 - E_N) * p$f1 * AT_N + p$RhoAT * (1 - E_M) * p$f1 * AT_M - p$SigmaMT * D_M * MT -
      p$BetaMT * MT
    out[i_atm, ] <- p$SigmaMT * D_M * MT + p$RhoAT * E_M * AT_M - p$BetaAT * AT_M
    out[i_ft, ] <- p$RhoAT * (1 - E_N) * (1 - p$f1) * AT_N + p$RhoAT * (1 - E_M) * (1 - p$f1) * AT_M - p$BetaFT * FT

    # B cells, 17 affinity classes (rows)
    PN <- rep(P_N, each = 17); PM <- rep(P_M, each = 17)
    r <- function(v) rep(v, each = 17)
    newc <- r(p$RhoAB_N) * (1 - Hb) * PN * GCB + r(p$RhoAB_M) * (1 - Hb) * PM * AB_M
    out[i$NB, ] <- r(p$BetaNB) * (nb0 - NB) - r(p$SigmaNB) * Fb * PN * NB
    out[i$AB_N, ] <- r(p$SigmaNB) * G * PN * NB - r(p$delayB) * AB_N - r(p$BetaAB) * AB_N
    out[i$GCB, ] <- r(p$delayB) * AB_N + r(p$RhoAB_N) * Hb * PN * GCB - r(p$BetaAB) * GCB
    out[i$MB, ] <- newc * r(p$g1) - r(p$SigmaMB) * Fb * PM * MB - r(p$BetaMB) * MB
    out[i$AB_M, ] <- r(p$SigmaMB) * G * PM * MB + r(p$RhoAB_M) * Hb * PM * AB_M - r(p$BetaAB) * AB_M
    out[i$SP, ] <- newc * r(p$g2) - r(p$BetaSP + p$LambdaSP) * SP
    out[i$LP, ] <- newc * r(1 - p$g1 - p$g2) - r(p$BetaLP + p$LambdaLP) * LP

    # plasma cells and antibody in blood
    spb <- Y[i$SP_BL, , drop = FALSE]; lpb <- Y[i$LP_BL, , drop = FALSE]
    out[i$SP_BL, ] <- -r(p$BetaSP) * spb + r(p$LambdaSP) * SP
    out[i$LP_BL, ] <- -r(p$BetaLP) * lpb + r(p$LambdaLP) * LP
    out[i$Ab_BL, ] <- -r(p$BetaA) * Y[i$Ab_BL, , drop = FALSE] + r(p$AlphaA / p$NAV * 1e12) * (spb + lpb)
    out
  }
}
