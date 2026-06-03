########################################################################################
# Berechnung der Larson-Parameter für Modul 5, Prozessschritt PostWash
########################################################################################

########################################################################################
# 1. Pakete und Hilfsfunktionen laden
########################################################################################

library(dplyr)
library(ggplot2)
library(tibble)

source("code_TA/00_helper_functions.R")

########################################################################################
# 2. Daten einlesen
########################################################################################

overview_m5_resin <- readRDS("raw_data/overview_m5_resin.rds")

# res_m5_with_resin und dat_m5 müssen aus dem vorherigen Aufbereitungsskript vorhanden sein

########################################################################################
# 3. Einstellungen für die Larson-Berechnung
########################################################################################

weights <- c(
  NoIP = 0.6,
  MRoC = 0.7,
  BV = 0.7,
  CE = 0.8,
  NGHETP = 0.8,
  A = 0.1,
  GHETP = 0.1
)

min_vol <- 0
max_vol <- 550
nMA <- 31
th_BV <- 0.05
n_confirm <- 30
early_window <- 0
L <- 140

########################################################################################
# 4. Beispielhafte Aufbereitung einer Charge
########################################################################################

first_bn <- overview_m5_resin$BatchNumber[1]

df_example <- res_m5_with_resin$data %>%
  filter(BatchNumber == first_bn, Process_Step == "PostWash") %>%
  arrange(Volumen_kumuliert) %>%
  filter(
    Volumen_kumuliert >= min_vol,
    Volumen_kumuliert <= max_vol
  )

v <- df_example$Volumen_kumuliert
c <- df_example$C004

ok <- is.finite(v) & is.finite(c)
v <- v[ok]
c <- c[ok]

keep <- c(TRUE, diff(v) > 1)
v <- v[keep]
c <- c[keep]

v_norm <- (v - min(v)) / (max(v) - min(v))
c_norm <- (c - min(c)) / (max(c) - min(c))

dv <- c(diff(c_norm) / diff(v_norm), 0)
dv[!is.finite(dv)] <- 0

madCdV <- stats::filter(dv, rep(1 / nMA, nMA), sides = 2)
madCdV[is.na(madCdV)] <- 0

# Optional: Beispielplots zur Kontrolle
# ggplot(tibble(Volumen = v, Leitfaehigkeit = c), aes(Volumen, Leitfaehigkeit)) +
#   geom_point(alpha = 0.7) +
#   labs(title = paste0(first_bn, " - Rohdaten"), x = "Volumen [L]", y = "Leitfähigkeit") +
#   theme_minimal()
#
# ggplot(tibble(Volumen = v_norm, Leitfaehigkeit = c_norm), aes(Volumen, Leitfaehigkeit)) +
#   geom_point(alpha = 0.7) +
#   labs(title = paste0(first_bn, " - normierte Rohdaten"), x = "normiertes Volumen", y = "normierte Leitfähigkeit") +
#   theme_minimal()
#
# ggplot(tibble(Volumen = v_norm, dCdV = dv), aes(Volumen, dCdV)) +
#   geom_point(alpha = 0.7) +
#   labs(title = paste0(first_bn, " - Ableitung"), x = "normiertes Volumen", y = "dC/dV") +
#   theme_minimal()
#
# ggplot(tibble(Volumen = v_norm, madCdV = madCdV), aes(Volumen, madCdV)) +
#   geom_point(alpha = 0.7) +
#   labs(title = paste0(first_bn, " - geglättete Ableitung"), x = "normiertes Volumen", y = "geglättete dC/dV") +
#   theme_minimal()

########################################################################################
# 5. Larson-Parameter für alle PostWash-Läufe berechnen
########################################################################################

results_postwash <- list()

for (bn in overview_m5_resin$BatchNumber) {
  
  df_batch <- res_m5_with_resin$data %>%
    filter(BatchNumber == bn, Process_Step == "PostWash") %>%
    arrange(Volumen_kumuliert) %>%
    filter(
      Volumen_kumuliert >= min_vol,
      Volumen_kumuliert <= max_vol,
      C004 >= 0
    )
  
  v <- df_batch$Volumen_kumuliert
  c <- df_batch$C004
  
  ok <- is.finite(v) & is.finite(c)
  v <- v[ok]
  c <- c[ok]
  
  keep <- c(TRUE, diff(v) > 1)
  v <- v[keep]
  c <- c[keep]
  
  if (length(v) < 20) {
    next
  }
  
  metrics <- compute_metrics_larson(
    v,
    c,
    transition_type = "up",
    nMA = nMA,
    L = L,
    th_BV = th_BV,
    n_confirm = n_confirm,
    min_vol = min_vol,
    early_window = early_window
  )
  
  batch_info <- overview_m5_resin %>%
    filter(BatchNumber == bn) %>%
    slice(1)
  
  results_postwash[[bn]] <- tibble(
    BatchNumber = bn,
    Label = batch_info$Label,
    ResinID = batch_info$ResinID,
    MaterialProduct = batch_info$MaterialProduct,
    FacilityID = batch_info$FacilityID,
    NoIP = metrics$NoIP,
    MRoC = metrics$MRoC,
    BV = metrics$BV,
    CE = metrics$CE,
    A = metrics$A,
    NGHETP = metrics$NGHETP,
    GHETP = metrics$GHETP
  )
}

metrics_postwash <- bind_rows(results_postwash)

########################################################################################
# 6. Zusatzinformationen ergänzen
########################################################################################

cycle_info <- overview_m5_resin %>%
  select(
    BatchNumber,
    earliest_Dat,
    CmCyclesResin,
    MaterialNumber,
    IgM_Wert,
    IgA_Wert
  )

metrics_postwash_labeled <- metrics_postwash %>%
  left_join(cycle_info, by = "BatchNumber") %>%
  mutate(Process_Step = "PostWash")

# Start-pH im PostWash
pH <- dat_m5 %>%
  select(2, 10, 11, 18)

pH_start <- pH %>%
  arrange(BatchNumber, ProcStep, Dat) %>%
  group_by(BatchNumber, ProcStep) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(BatchNumber, ProcStep, pH02) %>%
  rename(start_pH = pH02) %>%
  filter(ProcStep == "PostWash") %>%
  select(BatchNumber, start_pH)

metrics_postwash_labeled <- metrics_postwash_labeled %>%
  left_join(pH_start, by = "BatchNumber")

# Minimaler Flow im PostWash
flow <- dat_m5 %>%
  select(2, 10, 11, 13)

min_flow <- flow %>%
  filter(ProcStep == "PostWash") %>%
  filter(flow02 > 1500) %>%
  group_by(BatchNumber) %>%
  summarise(
    min_flow = min(flow02, na.rm = TRUE),
    .groups = "drop"
  )

metrics_postwash_labeled <- metrics_postwash_labeled %>%
  left_join(min_flow, by = "BatchNumber")

########################################################################################
# 7. Pseudonymisierung und Export
########################################################################################

dir.create("data_for_modeling", showWarnings = FALSE)
dir.create("pseudo_key", showWarnings = FALSE)

postwash_cols <- c(
  "BatchNumber", "MaterialNumber", "FacilityID", "MaterialProduct",
  "Process_Step", "Label", "IgM_Wert", "IgA_Wert",
  "earliest_Dat", "CmCyclesResin", "ResinID", "start_pH", "min_flow",
  "NoIP", "MRoC", "BV", "CE", "A", "NGHETP", "GHETP"
)

df_sorted <- metrics_postwash_labeled %>%
  arrange(earliest_Dat)

pseudonyms <- paste0("M5_Batch", seq_len(nrow(df_sorted)))

key_df <- data.frame(
  Pseudonym = pseudonyms,
  OriginalBatchNumber = df_sorted$BatchNumber,
  stringsAsFactors = FALSE
)

df_sorted$BatchNumber <- pseudonyms

df_export <- df_sorted[, postwash_cols]

colnames(df_export) <- c(
  "BatchNumber", "MaterialNumber", "FacilityID", "Product_Concentration",
  "Process_Step", "Label", "IgM_Wert", "IgA_Wert",
  "CHR_Start_Date", "Resin_Cycle_No", "ResinID", "start_pH", "min_flow",
  "NoIP", "MRoC", "BV", "CE", "A", "NGHETP", "GHETP"
)

write.csv(
  df_export,
  file = "data_for_modeling/m5_postwash_larson.csv",
  row.names = FALSE
)

write.csv(
  key_df,
  file = "pseudo_key/m5_key.csv",
  row.names = FALSE
)
