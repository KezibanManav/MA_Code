# Hilfsfunkionen:

library(pracma)

# BV-Funktion mit Hysterese und +min_vol
compute_BV_hysteresis <- function(v, c_norm,
                                  th = 0.05,
                                  n_confirm = 30,
                                  vol_confirm = NA_real_,
                                  min_vol = 0) {
  ok <- is.finite(v) & is.finite(c_norm)
  v <- v[ok]; c_norm <- c_norm[ok]
  if (length(v) < 2) return(list(BV = NA_real_, reason = "INSUFFICIENT_POINTS"))
  
  o <- order(v, method = "radix")
  v <- v[o]; c_norm <- c_norm[o]
  keep <- !duplicated(v)
  v <- v[keep]; c_norm <- c_norm[keep]
  if (length(v) < 2) return(list(BV = NA_real_, reason = "INSUFFICIENT_POINTS"))
  
  v0 <- min(v, na.rm = TRUE)
  mx <- max(c_norm, na.rm = TRUE); mn <- min(c_norm, na.rm = TRUE)
  if (mx < th) return(list(BV = NA_real_, reason = "NEVER_REACHED"))
  if (mn >= th) return(list(BV = 0 + min_vol, reason = "STARTS_ABOVE"))
  
  ge <- c_norm >= th
  runs <- rle(ge)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  
  i_run <- NA_integer_
  for (r in seq_along(runs$lengths)) {
    if (!runs$values[r]) next
    s <- starts[r]; e <- ends[r]
    ok_len <- runs$lengths[r] >= n_confirm
    ok_span <- TRUE
    if (is.finite(vol_confirm)) {
      ok_span <- (v[min(e, length(v))] - v[s]) >= vol_confirm
    }
    if (ok_len && ok_span) { i_run <- r; break }
  }
  
  if (is.na(i_run)) return(list(BV = NA_real_, reason = "NO_STABLE_RUN"))
  
  s <- starts[i_run]
  j <- max(1, s - 1)
  if (s == 1) return(list(BV = 0 + min_vol, reason = "STABLE_FROM_START"))
  
  x0 <- c_norm[j]; x1 <- c_norm[s]
  v0_seg <- v[j]; v1_seg <- v[s]
  if (!is.finite(x0) || !is.finite(x1) || !is.finite(v0_seg) || !is.finite(v1_seg) ||
      v1_seg == v0_seg || x1 == x0) {
    i_first <- which(ge)[1]
    BV <- if (is.na(i_first)) NA_real_ else max(0, v[i_first] - v0) + min_vol
    return(list(BV = BV, reason = "DEGENERATE_SEGMENT_FALLBACK"))
  }
  
  alpha <- (th - x0) / (x1 - x0)
  v_cross <- v0_seg + alpha * (v1_seg - v0_seg)
  BV <- max(0, v_cross - v0) + min_vol
  list(BV = BV, reason = "OK")
}

# Kennzahlenberechnung nach Larson
compute_metrics_larson <- function(v, c, transition_type = "up",
                                   nMA = 81, L = 250000, early_window = 0,
                                   th_BV = 0.05, n_confirm = 30, min_vol = 0) {
  # Optional Baseline-Cleaning
  if (early_window > 0) {
    early_idx <- which(v <= (min(v, na.rm = TRUE) + early_window))
    base_med <- median(c[early_idx], na.rm = TRUE)
    base_mad <- mad(c[early_idx], constant = 1.4826, na.rm = TRUE)
    low_thresh <- base_med - 3 * base_mad
    
    hampel_filter <- function(x, k = 5, t0 = 3) {
      n <- length(x); y <- x
      for (i in seq_len(n)) {
        lo <- max(1, i - k)
        hi <- min(n, i + k)
        med <- median(x[lo:hi], na.rm = TRUE)
        s <- 1.4826 * median(abs(x[lo:hi] - med), na.rm = TRUE)
        if (is.finite(s) && abs(x[i] - med) > t0 * s) {
          y[i] <- med
        }
      }
      y
    }
    
    c <- hampel_filter(c, k = 5, t0 = 3)
    c <- ifelse(c < low_thresh, base_med, c)
  }
  
  # Invertieren für step-down
  if (transition_type == "down") {
    c <- max(c, na.rm = TRUE) - c
  }
  
  # Normierung nur für C (Larson-konform); V bleibt im Originalmaßstab
  c_norm <- (c - min(c)) / (max(c) - min(c))
  
  # dC/dV nach ORIGINAL-VOLUMEN v (nicht v_norm!)
  dCdV_raw <- c(diff(c_norm) / diff(v), 0)
  dCdV_raw[!is.finite(dCdV_raw)] <- 0
  
  # Moving Average (Glättung)
  dCdV <- stats::filter(dCdV_raw, rep(1/nMA, nMA), sides = 2)
  dCdV[is.na(dCdV)] <- 0
  dCdV <- as.numeric(dCdV)
  
  # absd <- abs(madCdV)
  
  # NoIP <- sum(absd > 0.1 * max(absd, na.rm = TRUE), na.rm = TRUE)
  # MRoC <- max(absd, na.rm = TRUE)
  
  thr <- 0.1 * max(dCdV, na.rm = TRUE)
  
  if (length(dCdV) < 3 || !is.finite(thr)) {
    NoIP <- NA_integer_
  } else {
    is_peak <- rep(FALSE, length(dCdV))
    for (i in 2:(length(dCdV) - 1)) {
      is_peak[i] <- (dCdV[i] > dCdV[i - 1]) &&
        (dCdV[i] > dCdV[i + 1]) &&
        (dCdV[i] >= thr)
    }
    NoIP <- sum(is_peak)
  }
  
  MRoC <- max(dCdV, na.rm = TRUE)
  
  res_bv <- compute_BV_hysteresis(v, c_norm, th = th_BV, n_confirm = n_confirm, min_vol = min_vol)
  BV <- res_bv$BV
  
  # --- Helper: finde Crossing-Position V(c = target) per linearer Interpolation ---
  find_crossing_v <- function(v, c_norm, target) {
    # sucht erstes Crossing von <target zu >=target
    idx <- which(c_norm[-length(c_norm)] < target & c_norm[-1] >= target)[1]
    if (is.na(idx)) return(NA_real_)
    # lineare Interpolation zwischen idx und idx+1
    v1 <- v[idx];   v2 <- v[idx + 1]
    c1 <- c_norm[idx]; c2 <- c_norm[idx + 1]
    if (!is.finite(v1) || !is.finite(v2) || !is.finite(c1) || !is.finite(c2) || c2 == c1) return(NA_real_)
    v1 + (target - c1) * (v2 - v1) / (c2 - c1)
  }
  
  # --- CE nach Larson (praktische Implementierung) ---
  CE <- NA_real_
  
  # V[0] und V[1] als Randpunkte im betrachteten Fenster
  V0 <- min(v, na.rm = TRUE)
  V1 <- max(v, na.rm = TRUE)
  
  V05 <- find_crossing_v(v, c_norm, 0.5)
  
  if (is.finite(V0) && is.finite(V05) && is.finite(V1) && (V0 < V05) && (V05 < V1)) {
    
    # Segment [V0, V1]
    keep <- v >= V0 & v <= V1
    v_seg <- v[keep]
    c_seg <- c_norm[keep]
    
    # Punkt bei V05 einfügen (falls nicht vorhanden)
    insert_point <- function(vv, cc, vx) {
      if (any(abs(vv - vx) < 1e-12)) return(list(v=vv, c=cc))
      j <- max(which(vv < vx))
      if (!is.finite(j) || j >= length(vv)) return(list(v=vv, c=cc))
      c_x <- cc[j] + (cc[j+1] - cc[j]) * (vx - vv[j]) / (vv[j+1] - vv[j])
      vv2 <- append(vv, vx, after = j)
      cc2 <- append(cc, c_x, after = j)
      list(v=vv2, c=cc2)
    }
    
    tmp <- insert_point(v_seg, c_seg, V05)
    v_seg <- tmp$v; c_seg <- tmp$c
    
    # sortieren
    o <- order(v_seg)
    v_seg <- v_seg[o]; c_seg <- c_seg[o]
    
    # Index von V05
    i05 <- which.min(abs(v_seg - V05))
    
    # Ii: Fläche unter C von V0 bis V05
    Ii <- pracma::trapz(v_seg[1:i05], c_seg[1:i05])
    
    # If: Fläche über C (1-C) von V05 bis V1
    If <- pracma::trapz(v_seg[i05:length(v_seg)], 1 - c_seg[i05:length(v_seg)])
    
    CE <- (Ii + If) / (V1 - V0)
  }
  
  # i_max <- which.max(absd)
  # ten <- 0.1 * max(absd, na.rm = TRUE)
  # i_left <- which.min(abs(absd[1:i_max] - ten))
  # i_right <- which.min(abs(absd[i_max:length(absd)] - ten)) + i_max - 1
  # A <- if (is.finite(i_left) && is.finite(i_right)) {
  #   (v[i_right] - v[i_max]) / (v[i_max] - v[i_left])
  # } else NA_real_
  
  # Asymmetry (A) nach Larson
  A <- NA_real_

  if (length(dCdV) >= 3 && any(is.finite(dCdV))) {

    i_max <- which.max(dCdV)
    Dmax  <- dCdV[i_max]

    if (is.finite(Dmax) && Dmax > 0) {

      ten <- 0.1 * Dmax

      # lineare Interpolation für Crossing
      interp_cross <- function(x, y, i1, i2, y0) {
        x1 <- x[i1]; x2 <- x[i2]
        y1 <- y[i1]; y2 <- y[i2]
        if (!is.finite(x1) || !is.finite(x2) || !is.finite(y1) || !is.finite(y2) || y2 == y1) return(NA_real_)
        x1 + (y0 - y1) * (x2 - x1) / (y2 - y1)
      }

      # --- Va: letztes Crossing von <ten zu >=ten vor dem Maximum ---
      Va <- NA_real_
      if (i_max > 1) {
        left_candidates <- which(dCdV[1:(i_max-1)] < ten)
        if (length(left_candidates)) {
          iL <- max(left_candidates)  # letzter Punkt < ten
          if (iL < i_max) Va <- interp_cross(v, dCdV, iL, iL + 1, ten)
        } else {
          # falls schon von Start an >= ten, setze Va auf ersten v-Punkt
          # (optional; je nach Daten kann man auch NA lassen)
          Va <- v[1]
        }
      }

      # --- Vb: erstes Crossing von >=ten zu <ten nach dem Maximum ---
      Vb <- NA_real_
      if (i_max < length(dCdV)) {
        right_candidates <- which(dCdV[(i_max+1):length(dCdV)] < ten)
        if (length(right_candidates)) {
          iR <- i_max + min(right_candidates)      # erster Punkt < ten (globaler Index)
          # optional: Bracketing-Check, um Extrapolation zu vermeiden
          if (dCdV[iR - 1] >= ten && dCdV[iR] < ten) {
            Vb <- interp_cross(v, dCdV, iR - 1, iR, ten)
          } else {
            Vb <- NA_real_
          }
        } else {
          # falls nach dem Peak nie unter ten fällt, setze Vb auf letzten v-Punkt
          # (optional; je nach Daten kann man auch NA lassen)
          Vb <- v[length(v)]
        }
      }

      # --- A berechnen ---
      Vmax <- v[i_max]
      if (is.finite(Va) && is.finite(Vb) && (Vmax > Va) && (Vb > Vmax)) {
        A <- (Vb - Vmax) / (Vmax - Va)
      }
    }
  }
  
  
  # -----------------------
  # NGHETP (Larson): Momente aus dC/dV
  # M_k = ∫ V^k * (dC/dV) dV  (volume-based)
  M0 <- pracma::trapz(v, dCdV)
  M1 <- pracma::trapz(v, v  * dCdV)
  M2 <- pracma::trapz(v, v^2 * dCdV)
  
  mu <- if (is.finite(M0) && M0 != 0) (M1 / M0) else NA_real_
  sigma2 <- if (is.finite(M0) && M0 != 0) (M2 / M0) - (M1 / M0)^2 else NA_real_
  
  NGHETP <- if (is.finite(mu) && is.finite(sigma2) && mu != 0 && sigma2 > 0) {
    L * sigma2 / (mu^2)
  } else NA_real_
  
  # -----------------------
  # GHETP (Larson):  Vc,Vd bei 50% von (dC/dV)_max um den Peak herum
  half <- 0.5 * max(dCdV, na.rm = TRUE)
  
  # Hilfsfunktion: lineare Interpolation für Crossing
  interp_cross <- function(x, y, i1, i2, y0) {
    # lineares y(x) zwischen i1 und i2, gibt x bei y=y0
    x1 <- x[i1]; x2 <- x[i2]
    y1 <- y[i1]; y2 <- y[i2]
    if (!is.finite(x1) || !is.finite(x2) || !is.finite(y1) || !is.finite(y2) || y2 == y1) return(NA_real_)
    x1 + (y0 - y1) * (x2 - x1) / (y2 - y1)
  }
  
  # links: letzte Stelle < half vor i_max und erste Stelle >= half danach
  left_idx <- which(dCdV[1:i_max] < half)
  iL <- if (length(left_idx)) max(left_idx) else NA_integer_
  Vc <- if (!is.na(iL) && iL < i_max) interp_cross(v, dCdV, iL, iL + 1, half) else NA_real_
  
  # rechts: erste Stelle < half nach i_max (also crossing von >= half zu < half)
  right_idx <- which(dCdV[i_max:length(dCdV)] < half)
  iR <- if (length(right_idx)) (i_max + min(right_idx) - 2) else NA_integer_
  Vd <- if (!is.na(iR) && iR >= i_max && (iR + 1) <= length(dCdV)) interp_cross(v, dCdV, iR, iR + 1, half) else NA_real_
  
  GHETP <- if (is.finite(Vc) && is.finite(Vd) && (Vd - Vc) > 0) {
    Vmax <- v[i_max]
    N_G  <- 5.54 * (Vmax / (Vd - Vc))^2
    L / N_G
  } else NA_real_
  
  
  # M0 <- trapz(v, rep(1, length(v)))
  # M1 <- trapz(v, v)
  # M2 <- trapz(v, v^2)
  # variance <- (M2/M0) - (M1/M0)^2
  # NGHETP <- if (variance > 0) L / variance else NA_real_
  # 
  # half <- 0.5 * max(absd, na.rm = TRUE)
  # i_left50 <- which.min(abs(absd[1:i_max] - half))
  # i_right50 <- which.min(abs(absd[i_max:length(absd)] - half)) + i_max - 1
  # GHETP <- if (is.finite(i_left50) && is.finite(i_right50)) {
  #   Vc <- v[i_left50]; Vd <- v[i_right50]
  #   FWHM <- Vd - Vc
  #   if (FWHM > 0) L / (5.54 * ((v[i_max] / FWHM)^2)) else NA_real_
  # } else NA_real_
  
  return(list(NoIP = NoIP, MRoC = MRoC, BV = BV, CE = CE, A = A, NGHETP = NGHETP, GHETP = GHETP))
}

# OIV-Berechnung 
compute_OIV <- function(metrics, limits, weights) { 
  pass <- sapply(names(metrics), function(k) { 
    val <- metrics[[k]] 
    lim <- limits[[k]] 
    if (is.na(val)) return(0) 
    if (k %in% c("MRoC", "BV")) { 
      if (val < lim) if (val < lim * 0.9) 0 else 0.5 else 1 
    } else { 
        if (val > lim) if (val > lim * 1.1) 0 else 0.5 else 1 
    } 
  }) 
  sum(pass * weights[names(metrics)]) / sum(weights) 
  }

compute_limits_3sigma <- function(normal_df, metric_names) {
  out <- lapply(metric_names, function(k) {
    x <- normal_df[[k]]
    mu <- mean(x, na.rm = TRUE)
    sd <- stats::sd(x, na.rm = TRUE)
    list(
      KIL = mu - 3 * sd,
      KIU = mu + 3 * sd
    )
  })
  names(out) <- metric_names
  out
}

compute_OIV_larson <- function(metrics, limits, weights) {
  pass <- sapply(names(metrics), function(k) {
    val <- metrics[[k]]
    if (!is.finite(val)) return(0)

    KIL <- limits[[k]]$KIL
    KIU <- limits[[k]]$KIU
    if (!is.finite(KIL) || !is.finite(KIU)) return(0)

    # innerhalb Intervall
    if (val >= KIL && val <= KIU) return(1)

    # deutlich außerhalb (+/-10% außerhalb der Grenzen)
    if (val < 0.9 * KIL || val > 1.1 * KIU) return(0)

    # Übergangsbereich
    return(0.5)
  })

  sum(pass * weights[names(metrics)]) / sum(weights)
}

