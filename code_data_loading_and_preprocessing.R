########################################################################################
# Datenaufbereitung für die Masterarbeit
# Einlesen der Rohdaten, LIMS-Verknüpfung und Resin-Zuordnung
########################################################################################

########################################################################################
# 1. Pakete laden
########################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(tibble)
library(rlang)
library(readxl)

########################################################################################
# 2. Daten einlesen
########################################################################################

dat_m5 <- readRDS("raw_data/dat_m5.rds")
dat_m6 <- readRDS("raw_data/dat_m6.rds")

dat_lims <- read_excel("raw_data/CM_LIMS_for_TA.xlsx") %>%
  filter(component %in% c("IgA", "IgM"), !is.na(RESULT))

dat_pH <- readRDS("raw_data/pH_Info.rds")
colnames(dat_pH) <- c("BatchNumber", "pH_Value", "Process")

########################################################################################
# 3. Hilfsfunktionen
########################################################################################

parse_datetime_smart <- function(x, tz_local = "UTC") {
  if (inherits(x, "POSIXt")) {
    return(x)
  }
  suppressWarnings(lubridate::ymd_hms(x, tz = tz_local, quiet = TRUE))
}

cummax_posixct <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  
  tzx <- attr(x, "tzone")
  if (is.null(tzx) || !nzchar(tzx)) {
    tzx <- "UTC"
  }
  
  as.POSIXct(cummax(as.numeric(x)), origin = "1970-01-01", tz = tzx)
}

first_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & trimws(x) != ""]
  
  if (length(x) > 0) {
    x[1]
  } else {
    NA_character_
  }
}

as_num_robust <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s", "", x)
  x <- sub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

########################################################################################
# 4. Modul aufbereiten und mit LIMS-Informationen verknüpfen
########################################################################################

process_module <- function(dat_mod, lims_all, index_map = c(2, 3, 4, 11, 13, 15, 10), tz_local = "UTC") {
  
  target_names <- c( "BatchNumber", "MaterialNumber", "FacilityID", "Dat", "C002", "C004", "Process_Step")
  
  if (all(c("BatchNumber", "Process_Step") %in% names(dat_mod))) {
    dat_raw_min <- dat_mod %>% select(BatchNumber, Process_Step)
    
    module_batches <- dat_mod %>%distinct(BatchNumber)
  } else {
    dat_mod <- dat_mod %>% select(all_of(index_map))
    
    names(dat_mod) <- target_names
    
    dat_raw_min <- dat_mod %>% select(BatchNumber, Process_Step)
    
    module_batches <- dat_mod %>% distinct(BatchNumber)
  }
  
  dat_mod <- dat_mod %>%
    group_by(BatchNumber) %>%
    filter(any(Process_Step == "Produktion")) %>%
    ungroup()
  
  dat_mod <- dat_mod %>%
    group_by(BatchNumber) %>%
    mutate(
      EndProd = max(Dat[Process_Step == "Produktion"], na.rm = TRUE),
      WindowStart = EndProd - minutes(5),
      Process_Step = if_else(
        Process_Step == "Produktion" & Dat >= WindowStart & Dat <= EndProd,
        "PostWash",
        Process_Step
      )
    ) %>%
    ungroup()
  
  n_raw_equil_postwash <- dat_raw_min %>%
    group_by(BatchNumber) %>%
    summarise(
      has_equil = any(Process_Step == "Equilibrieren", na.rm = TRUE),
      has_postwash = any(Process_Step == "PostWash", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(has_equil & has_postwash) %>%
    nrow()
  
  lims_for_module <- lims_all %>%
    rename(
      BatchNumber = batch,
      ResultStatus = Result_Status,
      IncreasedIgM_F = Increased_IgM
    ) %>%
    semi_join(module_batches, by = "BatchNumber") %>%
    mutate(
      component = tolower(trimws(as.character(component))),
      CmCyclesResin = as.character(CmCyclesResin),
      MaterialProduct = as.character(MaterialProduct),
      RESULT_raw = as.character(RESULT),
      RESULT_num = as_num_robust(RESULT)
    )
  
  lims_for_module <- lims_for_module %>%
    mutate(
      is_abnormal = case_when(
        MaterialProduct == "IgPro20" & component == "igm" & RESULT_num > 2.71 ~ TRUE,
        MaterialProduct == "IgPro20" & component == "iga" & RESULT_num > 16.968 ~ TRUE,
        MaterialProduct == "IgPro10" & component == "igm" & RESULT_num > 1.23 ~ TRUE,
        MaterialProduct == "IgPro10" & component == "iga" & RESULT_num > 10.648 ~ TRUE,
        TRUE ~ FALSE
      )
    ) # Die Grenzen wurden von Prozessexperten basierend auf historischen Daten geschätzt
  
  lims_batch <- lims_for_module %>%
    group_by(BatchNumber) %>%
    summarise(
      Label = if_else(any(is_abnormal, na.rm = TRUE), "abnormal", "normal"),
      MaterialProduct = first_non_missing(MaterialProduct),
      CmCyclesResin = first_non_missing(CmCyclesResin),
      n_rows_lims = n(),
      any_OOS = any(has_OOS_IgM, na.rm = TRUE),
      any_incIgM = any(has_inc_IgM, na.rm = TRUE),
      IgM_Wert = first(RESULT_raw[component == "igm" & !is.na(RESULT_raw)]),
      IgA_Wert = first(RESULT_raw[component == "iga" & !is.na(RESULT_raw)]),
      .groups = "drop"
    )
  
  dat_labeled <- dat_mod %>%
    semi_join(lims_batch, by = "BatchNumber") %>%
    inner_join(lims_batch, by = "BatchNumber")
  
  dat_final <- dat_labeled %>%
    mutate(Dat = parse_datetime_smart(Dat, tz_local = tz_local)) %>%
    arrange(BatchNumber, Dat) %>%
    group_by(BatchNumber) %>%
    mutate(Dat = cummax_posixct(Dat)) %>%
    filter(any(Process_Step == "Equilibrieren", na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(BatchNumber, Dat) %>%
    group_by(BatchNumber) %>%
    mutate(
      Zeitdiff_sec = as.numeric(difftime(lead(Dat), Dat, units = "secs")),
      Volumen_L = (C002 / 3600) * Zeitdiff_sec
    ) %>%
    group_by(BatchNumber, Process_Step) %>%
    mutate(Volumen_kumuliert = cumsum(coalesce(Volumen_L, 0))) %>%
    ungroup()
  
  list(
    data = dat_final,
    n_with_label_from_module = nrow(lims_batch),
    n_raw_equil_postwash = n_raw_equil_postwash,
    n_total_labeled = dat_final %>% distinct(BatchNumber) %>% nrow(),
    n_after_equil = dat_final %>% distinct(BatchNumber) %>% nrow(),
    class_counts = dat_final %>%
      distinct(BatchNumber, Label) %>%
      count(Label, name = "n_batches")
  )
}

########################################################################################
# 5. Batch-Übersicht und Resin-Zuordnung
########################################################################################

make_batch_overview_base <- function(res_mod) {
  df <- res_mod$data
  
  if (!"CmCyclesResin" %in% names(df)) {
    stop("Spalte 'CmCyclesResin' fehlt.")
  }
  
  df <- df %>%
    mutate(CmCyclesResin = suppressWarnings(as.numeric(CmCyclesResin)))
  
  df %>%
    group_by(BatchNumber, MaterialNumber, FacilityID, MaterialProduct, Label, IgM_Wert, IgA_Wert) %>%
    summarise(
      earliest_Dat = suppressWarnings(min(Dat, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    left_join(
      df %>%
        group_by(BatchNumber) %>%
        filter(Dat == min(Dat, na.rm = TRUE)) %>%
        summarise(
          CmCyclesResin = suppressWarnings(min(CmCyclesResin, na.rm = TRUE)),
          .groups = "drop"
        ),
      by = "BatchNumber"
    ) %>%
    mutate(CmCyclesResin = suppressWarnings(as.numeric(CmCyclesResin))) %>%
    arrange(FacilityID, earliest_Dat, BatchNumber)
}

tag_resin_on_overview <- function(overview_df,
                                  module_name,
                                  include_facility_in_id = TRUE,
                                  id_sep = "|") {
  
  overview_df %>%
    mutate(CmCyclesResin = suppressWarnings(as.numeric(CmCyclesResin))) %>%
    arrange(FacilityID, earliest_Dat, BatchNumber) %>%
    group_by(FacilityID) %>%
    mutate(
      cycle_prev = lag(CmCyclesResin),
      is_reset = case_when(
        is.na(cycle_prev) ~ TRUE,
        is.na(CmCyclesResin) ~ FALSE,
        CmCyclesResin < cycle_prev ~ TRUE,
        TRUE ~ FALSE
      ),
      resin_idx = cumsum(is_reset)
    ) %>%
    group_by(FacilityID, resin_idx) %>%
    mutate(
      ResinStart = min(earliest_Dat, na.rm = TRUE),
      ResinStartCycle = suppressWarnings(min(CmCyclesResin, na.rm = TRUE)),
      ResinStartDateTag = format(ResinStart, "%Y%m%d"),
      ResinOrdinal = dense_rank(ResinStart),
      ResinID = if (include_facility_in_id) {
        paste(module_name, FacilityID, ResinStartDateTag, sep = id_sep)
      } else {
        paste(module_name, ResinStartDateTag, sep = id_sep)
      }
    ) %>%
    ungroup() %>%
    select(-cycle_prev, -is_reset)
}

attach_resin_id_to_resdata <- function(res_mod, overview_with_resin) {
  resin_info <- overview_with_resin %>%
    distinct(
      BatchNumber,
      ResinID,
      ResinStart,
      ResinStartDateTag,
      ResinStartCycle,
      ResinOrdinal
    )
  
  res_mod$data <- res_mod$data %>%
    left_join(resin_info, by = "BatchNumber")
  
  res_mod
}

########################################################################################
# 6. Module verarbeiten
########################################################################################

res_m5 <- process_module(dat_m5, dat_lims)
res_m6 <- process_module(dat_m6, dat_lims)

########################################################################################
# 7. Resin-IDs für Modul 5
########################################################################################

overview_m5_base <- make_batch_overview_base(res_m5)

overview_m5_resin <- tag_resin_on_overview(
  overview_m5_base,
  module_name = "M5",
  include_facility_in_id = TRUE
)

# manuelle Korrektur für M5, da hier ein Resinwechsel nicht vollständig über den Zykluszähler erkannt wird
overview_m5_resin$ResinID <- ifelse(
  overview_m5_resin$resin_idx == 3 & overview_m5_resin$CmCyclesResin < 40,
  "M5|M5|20210708",
  ifelse(
    overview_m5_resin$resin_idx == 3 & overview_m5_resin$CmCyclesResin > 140,
    "M5|M5|20210930",
    overview_m5_resin$ResinID
  )
)

saveRDS(overview_m5_resin, file = "raw_data/overview_m5_resin.rds")

res_m5_with_resin <- attach_resin_id_to_resdata(
  res_m5,
  overview_m5_resin
)

########################################################################################
# 8. Resin-IDs für Modul 6
########################################################################################

overview_m6_base <- make_batch_overview_base(res_m6)

overview_m6_resin <- tag_resin_on_overview(
  overview_m6_base,
  module_name = "M6",
  include_facility_in_id = TRUE
)

saveRDS(overview_m6_resin, file = "raw_data/overview_m6_resin.rds")

res_m6_with_resin <- attach_resin_id_to_resdata(
  res_m6,
  overview_m6_resin
)

########################################################################################
# 9. Ergebnisse speichern
########################################################################################

saveRDS(res_m5_with_resin, file = "raw_data/res_m5_with_resin.rds")
saveRDS(res_m6_with_resin, file = "raw_data/res_m6_with_resin.rds")
