library(dplyr)
library(lubridate)

df1 <- read.csv("data_for_modeling/m5_postwash_larson.csv", sep = ",", header = TRUE)
df2 <- read.csv("data_for_modeling/m6_postwash_larson.csv", sep = ",", header = TRUE)
df  <- rbind(df1, df2)
df <- df %>%
  dplyr::rename(Module = FacilityID)

df$ResinID <- sub("^[^|]*\\|", "", df$ResinID)

#df$Label <- make_label(df$Label)

df <- df %>%
  mutate(
    CHR_Start_Date = ymd_hms(CHR_Start_Date, quiet = TRUE),
    year  = year(CHR_Start_Date),
    month = floor_date(CHR_Start_Date, unit = "month")
  )


overview_tbl <- df %>%
  summarise(
    Anzahl_Laeufe = n(),
    Normal = sum(Label == "normal"),
    Abnormal = sum(Label == "abnormal"),
    Anteil_Abnormal = round(100 * mean(Label == "abnormal"), 1),
    Anlagen = n_distinct(Module),
    Produktkategorien = n_distinct(Product_Concentration),
    Resin_Chargen = n_distinct(ResinID),
    Startdatum = min(CHR_Start_Date, na.rm = TRUE),
    Enddatum = max(CHR_Start_Date, na.rm = TRUE)
  )

overview_tbl



# Gemeinsamer Ausgabeordner für EDA-Grafiken

eda_out_dir <- "thesis/iu-latex-template-main"
dir.create(eda_out_dir, recursive = TRUE, showWarnings = FALSE)

overview_tbl <- tibble(
  n_runs = nrow(df),
  n_normal = sum(df$Label == "normal"),
  n_abnormal = sum(df$Label == "abnormal"),
  abnormal_rate = mean(df$Label == "abnormal"),
  n_facilities = n_distinct(df$Module),
  n_products = n_distinct(df$Product_Concentration),
  n_resins = n_distinct(df$ResinID),
  start_date = min(df$CHR_Start_Date, na.rm = TRUE),
  end_date   = max(df$CHR_Start_Date, na.rm = TRUE)
)

overview_tbl



# Analyse der Klassenverteilung der Zielvariable

library(ggplot2)

p_label <- df %>%
  count(Label) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = Label, y = n, fill = Label)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", round(100*prop,1), "%)")),
            vjust = -0.5) +
  labs(
    title = "Klassenverteilung im Datensatz",
    x = NULL,
    y = "Anzahl Läufe"
  ) +
  theme_minimal()

p_label

ggsave(
  file.path(eda_out_dir, "EDA_class_distribution.pdf"),
  plot = p_label,
  width = 6,
  height = 4
)



# Untersuchung der Label-Verteilung für unterschiedliche Produktkonzentrationen

p_product <- df %>%
  count(Product_Concentration, Label) %>%
  group_by(Product_Concentration) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(Product_Concentration, prop, fill = Label)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Label-Verteilung nach Produktkonzentration",
    x = "Produkt",
    y = "Anteil"
  ) +
  theme_minimal()

p_product

ggsave(
  file.path(eda_out_dir, "EDA_label_distribution_by_product.pdf"),
  plot = p_product,
  width = 6,
  height = 4
)

product_overview <- df %>%
  group_by(Product_Concentration) %>%
  summarise(
    Laeufe = n(),
    Abnormal = sum(Label == "abnormal"),
    Anteil_Abnormal = round(100 * mean(Label == "abnormal"), 1)
  )

product_overview

#####

product_overview <- df %>%
  group_by(Module) %>%
  summarise(
    Laeufe = n(),
    Abnormal = sum(Label == "abnormal"),
    Anteil_Abnormal = round(100 * mean(Label == "abnormal"), 1)
  )

product_overview


# Analyse der Anzahl der Läufe pro Resin-Charge

resin_overview <- df %>%
  group_by(ResinID) %>%
  summarise(
    Laeufe = n(),
    Min_Zyklus = min(Resin_Cycle_No),
    Max_Zyklus = max(Resin_Cycle_No)
  ) %>%
  arrange(desc(Laeufe))

resin_overview


resin_overview <- df %>%
  group_by(ResinID) %>%
  summarise(
    Module = first(Module),
    Laeufe = n(),
    Min_Zyklus = min(Resin_Cycle_No),
    Max_Zyklus = max(Resin_Cycle_No),
    First_Run = min(CHR_Start_Date),
    Last_Run = max(CHR_Start_Date),
    Duration_days = as.numeric(difftime(Last_Run, First_Run, units="days"))
  ) %>%
  arrange(First_Run)

resin_overview


library(dplyr)
library(ggplot2)

df_plot_time <- df %>%
  mutate(
    y = ifelse(Label == "abnormal", 1, 0)
  )

p_time_runs <- ggplot(df_plot_time,
                      aes(x = CHR_Start_Date, y = y, color = ResinID)) +
  geom_point(alpha = 0.8, size = 1.8) +
  facet_wrap(~Module, scales = "free_x") +
  labs(
    title = "Zeitliche Verteilung der Chromatografieläufe",
    x = "Zeitpunkt des Laufs",
    y = "Laufstatus (0 = normal, 1 = abnormal)",
    color = "Resin-Charge"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12, face = "bold")
  )

p_time_runs



p_resin <- df %>%
  count(ResinID, Label) %>%
  group_by(ResinID) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(ResinID = reorder(ResinID, total)) %>%
  ggplot(aes(ResinID, n, fill = Label)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Anzahl der Läufe pro Resin-Charge",
    x = NULL,
    y = "Anzahl Läufe"
  ) +
  theme_minimal()

p_resin

ggsave(
  file.path(eda_out_dir, "EDA_runs_per_resin.pdf"),
  plot = p_resin,
  width = 8,
  height = 5
)

module_dates <- df %>%
  group_by(Module) %>%
  summarise(
    Startdatum = min(CHR_Start_Date, na.rm = TRUE),
    Enddatum   = max(CHR_Start_Date, na.rm = TRUE),
    Laeufe     = n()
  ) %>%
  arrange(Module)

module_dates



# Untersuchung der zeitlichen Verteilung der Läufe

df <- df %>%
  mutate(quarter = floor_date(CHR_Start_Date, "quarter"))

p_time <- df %>%
  count(quarter, Label) %>%
  ggplot(aes(quarter, n, color = Label)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Zeitliche Verteilung der Chromatografieläufe (quartalsweise)",
    x = "Zeit",
    y = "Anzahl Läufe"
  ) +
  theme_minimal()

p_time

ggsave(
  file.path(eda_out_dir, "EDA_time_distribution_runs.pdf"),
  plot = p_time,
  width = 8,
  height = 4
)

df <- df %>%
  mutate(
    year = year(CHR_Start_Date),
    quarter = paste0("Q", quarter(CHR_Start_Date)),
    year_quarter = paste(year, quarter)
  )

p_time_quarter <- df %>%
  count(year_quarter, Label) %>%
  ggplot(aes(x = year_quarter, y = n, color = Label, group = Label)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "Zeitliche Verteilung der Chromatografieläufe (quartalsweise)",
    x = "Jahr und Quartal",
    y = "Anzahl Läufe"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_time_quarter

ggsave(
  "thesis/iu-latex-template-main/EDA_time_distribution_quarterly.pdf",
  plot = p_time_quarter,
  width = 9,
  height = 4
)




# Untersuchung der Verteilungen der wichtigsten Transitionsparameter nach Larson

library(tidyr)

larson_features <- df %>%
  select(NoIP, MRoC, BV, CE, A, NGHETP, GHETP)

larson_long <- larson_features %>%
  pivot_longer(everything(),
               names_to = "feature",
               values_to = "value")

p_larson <- ggplot(larson_long, aes(x = value)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.8) +
  facet_wrap(~feature, scales = "free") +
  labs(
    title = "Verteilungen der Larson-Transitionsparameter",
    x = NULL,
    y = "Häufigkeit"
  ) +
  theme_minimal()

p_larson

ggsave(
  file.path(eda_out_dir, "EDA_larson_distributions.pdf"),
  plot = p_larson,
  width = 9,
  height = 6
)



# Untersuchung des Zusammenhangs zwischen Resin-Alter und Auftreten abnormaler Läufe


# df <- df %>%
#   mutate(y = ifelse(Label == "abnormal", 1, 0))
# 
# p_cycles <- df %>%
#   group_by(Resin_Cycle_No) %>%
#   summarise(rate = mean(y)) %>%
#   ggplot(aes(Resin_Cycle_No, rate)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(, se=FALSE, color="blue") +
#   labs(
#     title = "Rate abnormaler Läufe in Abhängigkeit von der Resin-Zyklusnummer",
#     x = "Resin-Zyklusnummer",
#     y = "Abnormal-Rate"
#   ) +
#   theme_minimal()
# 
# p_cycles
# 
# ggsave(
#   file.path(eda_out_dir, "EDA_abnormal_rate_over_resin_cycle.pdf"),
#   plot = p_cycles,
#   width = 7,
#   height = 4.5
# )
# 
# df <- df %>%
#   mutate(y = ifelse(Label == "abnormal", 1, 0))
# 
# df_rate <- df %>%
#   group_by(Resin_Cycle_No) %>%
#   summarise(rate = mean(y))
# 
# # lineares Modell
# lm_fit <- lm(rate ~ Resin_Cycle_No, data = df_rate)
# 
# r2 <- summary(lm_fit)$r.squared
# slope <- coef(lm_fit)["Resin_Cycle_No"]
# p_val <- summary(lm_fit)$coefficients["Resin_Cycle_No", "Pr(>|t|)"]
# p_text <- ifelse(p_val < 0.001, "p < 0.001",
#                  paste0("p = ", signif(p_val, 3)))
# 
# label_text <- paste0(
#   "Slope = ", round(slope, 4), "\n",
#   "R² = ", round(r2, 3)
# )
# 
# p_cycles <- ggplot(df_rate, aes(Resin_Cycle_No, rate)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE, color = "blue") +
#   annotate(
#     "text",
#     x = 80,
#     y = 0.8,
#     label = label_text,
#     hjust = 1.1,
#     vjust = 1.5,
#     size = 4.5
#   ) +
#   labs(
#     title = "Rate abnormaler Läufe in Abhängigkeit von der Resin-Zyklusnummer",
#     x = "Resin-Zyklusnummer",
#     y = "Abnormal-Rate"
#   ) +
#   theme_minimal()
# 
# p_cycles
# 
# ggsave(
#   file.path(eda_out_dir, "EDA_abnormal_rate_over_resin_cycle_r2.pdf"),
#   plot = p_cycles,
#   width = 7,
#   height = 4.5
# )



library(dplyr)
library(ggplot2)

df_cycle <- df %>%
  mutate(
    y = ifelse(Label == "abnormal", 1, 0)
  ) %>%
  filter(!is.na(Resin_Cycle_No), !is.na(y))

# Beobachtete Abnormalrate je Zyklusnummer
df_rate <- df_cycle %>%
  group_by(Resin_Cycle_No) %>%
  summarise(
    abnormal_rate = mean(y),
    n = n(),
    .groups = "drop"
  )

# Exploratives logistisches Modell auf Einzellaufebene
glm_cycle <- glm(
  y ~ Resin_Cycle_No,
  data = df_cycle,
  family = binomial(link = "logit")
)

summary(glm_cycle)

# Vorhersagegitter
pred_grid <- tibble(
  Resin_Cycle_No = seq(
    min(df_cycle$Resin_Cycle_No, na.rm = TRUE),
    max(df_cycle$Resin_Cycle_No, na.rm = TRUE),
    by = 1
  )
)

pred_link <- predict(glm_cycle, newdata = pred_grid, type = "link", se.fit = TRUE)

pred_grid <- pred_grid %>%
  mutate(
    fit_link = pred_link$fit,
    se_link  = pred_link$se.fit,
    prob     = plogis(fit_link),
    prob_lo  = plogis(fit_link - 1.96 * se_link),
    prob_hi  = plogis(fit_link + 1.96 * se_link)
  )

# Odds Ratio pro zusätzlichem Zyklus
beta_cycle <- coef(glm_cycle)["Resin_Cycle_No"]
or_cycle   <- exp(beta_cycle)

p_val <- summary(glm_cycle)$coefficients["Resin_Cycle_No", "Pr(>|z|)"]
p_text <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", signif(p_val, 3)))

label_text <- paste0(
  "Odds Ratio pro Zyklus = ", round(or_cycle, 3), "\n",
  p_text
)

p_cycles_glm <- ggplot() +
  geom_point(
    data = df_rate,
    aes(x = Resin_Cycle_No, y = abnormal_rate),
    alpha = 0.45,
    size = 2,
    color = "grey35"
  ) +
  geom_ribbon(
    data = pred_grid,
    aes(x = Resin_Cycle_No, ymin = prob_lo, ymax = prob_hi),
    inherit.aes = FALSE,
    alpha = 0.2,
    fill = "grey50"
  ) +
  geom_line(
    data = pred_grid,
    aes(x = Resin_Cycle_No, y = prob),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  annotate(
    "text",
    x = max(df_cycle$Resin_Cycle_No, na.rm = TRUE) * 0.5,
    y = 0.88,
    label = label_text,
    hjust = 0,
    vjust = 1,
    size = 4.5
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    title = "Wahrscheinlichkeit abnormaler Läufe nach der Resin-Zyklusnummer",
    subtitle = "Punkte = beobachtete Abnormalrate je Zyklusnummer; \nLinie = modellierte Wahrscheinlichkeit; Band = 95%-Konfidenzintervall",
    x = "Resin-Zyklusnummer",
    y = "Wahrscheinlichkeit für \"abnormal\""
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 13, face = "bold"))

p_cycles_glm

ggsave(
  file.path(eda_out_dir, "EDA_abnormal_probability_over_resin_cycle.pdf"),
  plot = p_cycles_glm,
  width = 7.5,
  height = 4.8
)

# glm(y ~ Resin_Cycle_No, data=df, family="binomial")


# Analyse der linearen Zusammenhänge zwischen den wichtigsten numerischen Merkmalen

num_vars <- df %>%
  select(NoIP, MRoC, BV, CE, A, NGHETP, GHETP,
         Resin_Cycle_No, start_pH, min_flow) 

cor_mat <- cor(num_vars, use = "complete.obs")

cor_mat

library(corrplot)

pdf(file.path(eda_out_dir, "EDA_feature_correlation_matrix.pdf"), width = 7, height = 7)
corrplot(cor_mat,
         method = "color",
         type = "upper",
         tl.cex = 0.8)
dev.off()

corrplot(cor_mat,
         method = "color",
         type = "upper",
         tl.cex = 0.8)



# Vergleich der Larson-Parameter zwischen normalen und abnormalen Läufen mittels Boxplots und Rohdatenpunkten

library(dplyr)
library(tidyr)
library(ggplot2)

larson_vars <- c("NoIP", "MRoC", "BV", "CE", "A", "NGHETP", "GHETP")

df_larson_long <- df %>%
  select(Label, all_of(larson_vars)) %>%
  pivot_longer(
    cols = -Label,
    names_to = "Feature",
    values_to = "Value"
  )

p_larson_box <- ggplot(df_larson_long, aes(x = Label, y = Value, fill = Label)) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.8, color = "grey35") +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  facet_wrap(~Feature, scales = "free", ncol = 2) +
  labs(
    title = "Verteilung der Larson-Parameter nach Klasse",
    x = NULL,
    y = "Wert"
  ) +
  theme_minimal()

p_larson_box

ggsave(
  file.path(eda_out_dir, "EDA_larson_boxplots_by_label.pdf"),
  plot = p_larson_box,
  width = 10,
  height = 12
)



# Vergleich zusätzlicher metrischer Prozessvariablen zwischen normalen und abnormalen Läufen

extra_vars <- c("Resin_Cycle_No", "start_pH", "min_flow") 

df_extra_long <- df %>%
  select(Label, all_of(extra_vars)) %>%
  pivot_longer(
    cols = -Label,
    names_to = "Feature",
    values_to = "Value"
  )

p_extra_box <- ggplot(df_extra_long, aes(x = Label, y = Value, fill = Label)) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.8, color = "grey35") +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  facet_wrap(~Feature, scales = "free", ncol = 2) +
  labs(
    title = "Verteilung zusätzlicher Prozessvariablen nach Klasse",
    x = NULL,
    y = "Wert"
  ) +
  theme_minimal()

p_extra_box

ggsave(
  file.path(eda_out_dir, "EDA_metrics_boxplots_by_label.pdf"),
  plot = p_extra_box,
  width = 10,
  height = 6
)



# # Entwicklung der Abnormal-Rate über Resin-Zyklen getrennt nach Resin-Chargen
# 
# p_resin_cycles <- df %>%
#   group_by(ResinID, Resin_Cycle_No) %>%
#   summarise(rate = mean(y), .groups = "drop") %>%
#   ggplot(aes(Resin_Cycle_No, rate)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(se = FALSE) +
#   facet_wrap(~ResinID, scales = "free_x") +
#   labs(
#     title = "Abnormal-Rate über Resin-Zyklen pro Resin-Charge",
#     x = "Resin-Zyklus",
#     y = "Abnormal-Rate"
#   ) +
#   theme_minimal()
# 
# p_resin_cycles
# 
# ggsave(
#   file.path(eda_out_dir, "EDA_abnormal_rate_by_resin.pdf"),
#   plot = p_resin_cycles,
#   width = 10,
#   height = 8
# )





# Entwicklung aller metrischen Features über Resin-Zyklen, getrennt nach Resin-Charge

library(dplyr)
library(ggplot2)


# Eine Grafik pro metrischem Parameter: Verlauf über Resin-Zyklen, getrennt nach Resin-Charge


metric_vars <- c(
  "NoIP", "MRoC", "BV", "CE", "A", "NGHETP", "GHETP",
  "start_pH", "min_flow"
) 

plots_by_feature <- lapply(metric_vars, function(varname) {
  
  ggplot(df, aes(x = Resin_Cycle_No, y = .data[[varname]], color = Label)) +
    geom_point(alpha = 0.6, size = 1.2) +
    facet_wrap(~ResinID, scales = "free_x", ncol = 3) +
    labs(
      title = paste0(varname, " über Resin-Zyklen pro Resin-Charge"),
      x = "Resin-Zyklus",
      y = varname,
      color = "Klasse"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12)
    )
})

names(plots_by_feature) <- metric_vars

for (nm in names(plots_by_feature)) {
  print(plots_by_feature[[nm]])
}

for (nm in names(plots_by_feature)) {
  ggsave(
    filename = file.path(eda_out_dir, paste0("EDA_", nm, "_over_resin_cycles.pdf")),
    plot = plots_by_feature[[nm]],
    width = 10,
    height = 7
  )
}


# Univariate AUC pro Feature

library(pROC)
library(dplyr)

features <- c(
  "NoIP","MRoC","BV","CE","A",
  "NGHETP","GHETP",
  "start_pH","min_flow"
)

df$y <- ifelse(df$Label == "abnormal", 1, 0)

auc_tbl <- lapply(features, function(v){
  
  r <- roc(df$y, df[[v]], quiet = TRUE)
  
  data.frame(
    Feature = v,
    AUC = as.numeric(r$auc)
  )
  
}) %>% bind_rows()

auc_tbl <- auc_tbl %>%
  arrange(desc(AUC))

auc_tbl

auc_tbl <- lapply(features, function(v){
  
  r <- roc(df$y, df[[v]], quiet = TRUE)
  
  ci_vals <- ci.auc(r)
  
  data.frame(
    Feature = v,
    AUC = as.numeric(r$auc),
    CI_low = ci_vals[1],
    CI_high = ci_vals[3]
  )
  
}) %>% bind_rows() %>%
  arrange(desc(AUC)
  )

auc_tbl



library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

plot_vars <- c(
  "NoIP", "MRoC", "BV", "CE", "A", "NGHETP", "GHETP",
  "start_pH", "min_flow"
)

feature_levels <- c(
  "NoIP", "MRoC", "BV", "CE", "A", "NGHETP", "GHETP",
  "start_pH", "min_flow"
)

df_plot_long <- df %>%
  select(Label, all_of(plot_vars)) %>%
  pivot_longer(
    cols = -Label,
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  mutate(
    Label = factor(Label, levels = c("normal", "abnormal")),
    Feature = factor(Feature, levels = feature_levels)
  )



# Hilfsfunktion für gut lesbare p-Werte in der Grafik

format_p_for_plot <- function(p) {
  ifelse(
    is.na(p),
    "p = NA",
    ifelse(
      p < 0.001,
      "p < 0.001",
      paste0("p = ", formatC(p, format = "f", digits = 3))
    )
  )
}


# Welch-t-Tests + Cohen's d

ttest_tbl <- df_plot_long %>%
  group_by(Feature) %>%
  t_test(Value ~ Label, var.equal = FALSE) %>%
  left_join(
    df_plot_long %>%
      group_by(Feature) %>%
      cohens_d(Value ~ Label),
    by = "Feature"
  ) %>%
  mutate(
    p_label = format_p_for_plot(p),
    d_label = paste0("Cohen's d = ", sprintf("%.2f", effsize))
  ) %>%
  ungroup()

# Annotation oben: Cohen's d
ann_top <- ttest_tbl %>%
  transmute(
    Feature,
    x = 2,
    y = Inf,
    label = d_label
  )

# Annotation unten: p-Wert des Welch-t-Tests
ann_bottom <- ttest_tbl %>%
  transmute(
    Feature,
    x = 1.8,
    y = -Inf,
    label = p_label
  )


# Plot

p_metrics_violin <- ggplot(df_plot_long, aes(x = Label, y = Value, fill = Label)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(width = 0.12, alpha = 0.15, size = 0.5, color = "grey30") +
  geom_label(
    data = ann_top,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.1,
    size = 5,
    label.size = 0.2
  ) +
  geom_label(
    data = ann_bottom,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = -0.1,
    size = 5,
    label.size = 0.2
  ) +
  facet_wrap(~ Feature, scales = "free", ncol = 2) +
  labs(
    title = "Verteilung der Transitionsparameter und zusätzlicher Prozessvariablen nach Klasse",
    x = NULL,
    y = "Wert",
    fill = "Klasse"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )

p_metrics_violin

ggsave(
  file.path(eda_out_dir, "EDA_metric_boxplots_by_label_.pdf"),
  plot = p_metrics_violin,
  width = 12,
  height = 18
)


####################


library(dplyr)
library(tidyr)

summary_tbl <- df_plot_long %>%
  group_by(Feature, Label) %>%
  summarise(
    n = n(),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    q1 = quantile(Value, 0.25, na.rm = TRUE),
    q3 = quantile(Value, 0.75, na.rm = TRUE),
    IQR = IQR(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Feature, Label)

summary_tbl

summary_wide <- summary_tbl %>%
  select(Feature, Label, median) %>%
  pivot_wider(
    names_from = Label,
    values_from = median
  ) %>%
  mutate(
    median_difference = abnormal - normal
  )

summary_wide


##############################

df1 <- read.csv("data_for_modeling/m5_postwash_larson.csv", sep = ",", header = TRUE)
df2 <- read.csv("data_for_modeling/m6_postwash_larson.csv", sep = ",", header = TRUE)
df  <- rbind(df1, df2)
df <- df %>%
  dplyr::rename(Module = FacilityID)
df$ResinID <- sub("^[^|]*\\|", "", df$ResinID)

df$Label <- make_label(df$Label)

# Feature-Set
features <- c("NoIP","MRoC","BV","CE","A","NGHETP","GHETP",
              "Resin_Cycle_No", "start_pH","min_flow","Product_Concentration") 

# Nur komplette Zeilen
df_model <- df %>%
  dplyr::filter(complete.cases(dplyr::across(all_of(c(features, "Label")))))

numeric_features <- c(
  "NoIP","MRoC","BV","CE","A","NGHETP","GHETP",
  "Resin_Cycle_No","start_pH","min_flow"
)

cor_mat <- cor(
  df_model[, numeric_features],
  use = "pairwise.complete.obs",
  method = "pearson"
)

round(cor_mat, 3)

library(tidyr)
library(dplyr)

cor_table <- cor_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "Correlation") %>%
  filter(Var1 != Var2) %>%
  mutate(abs_cor = abs(Correlation)) %>%
  arrange(desc(abs_cor))

cor_table


cor_table %>%
  filter(abs_cor > 0.8)

cor_table %>%
  filter(abs_cor > 0.7)

library(ggplot2)

cor_long <- cor_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to="Var2", values_to="Correlation")

ggplot(cor_long, aes(Var1, Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       midpoint = 0, limits=c(-1,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Korrelationsmatrix der Features")



# EDA: Correlation Table + Heatmap

cor_mat <- cor(
  df_model[, numeric_features_for_corr, drop = FALSE],
  use = "pairwise.complete.obs",
  method = "pearson"
)

cor_table <- cor_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "Correlation") %>%
  dplyr::filter(Var1 != Var2) %>%
  mutate(abs_cor = abs(Correlation)) %>%
  arrange(desc(abs_cor))

high_cor_table <- cor_table %>%
  filter(abs_cor > 0.7) %>%
  arrange(desc(abs_cor))

safe_write_csv(
  as.data.frame(round(cor_mat, 3)) %>% tibble::rownames_to_column("Variable"),
  "results/correlation_matrix.csv"
)
safe_write_csv(cor_table, "results/correlation_table_long.csv")
safe_write_csv(high_cor_table, "results/high_correlations_gt_0_7.csv")

cor_long <- cor_long %>%
  mutate(
    label = sprintf("%.2f", Correlation),
    text_col = ifelse(abs(Correlation) > 0.5, "white", "black")
  )

p_cor <- ggplot(cor_long, aes(Var1, Var2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = label, color = text_col), size = 3, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1)
  ) +
  coord_fixed() +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "Feature Correlation Heatmap", x = NULL, y = NULL)

ggsave(
  file.path(eda_out_dir, "EDA_heatmap_.pdf"),
  plot = p_cor,
  width = 8,
  height = 7
)
