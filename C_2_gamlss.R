# Load necessary libraries
library(gamlss)
library(dplyr)
library(readr)
library(caret)
library(ggplot2)
library(tidyr)

Top <- "/data/projects/punim2400"

harmonized_scores_file <- file.path(Top, "derivatives", "MASTER_HARMONIZED_deltar_table_BCPHCPD.csv")
harmonized_table <- read_csv(harmonized_scores_file)

output_dir <- file.path(Top, "derivatives", "gamlss_scores_BCPHCPD")

master_table <- harmonized_table

feature_names <- c("deltar_a", "deltar_h")

confmodel = lm('deltar_a ~ MeanFD + ICV',master_table)
master_table$deltar_amyg = resid(confmodel)
confmodel = lm('deltar_h ~ MeanFD + ICV',master_table)
master_table$deltar_hipp = resid(confmodel)

# Iterative run to find optimal distribution

model_data <- master_table %>%
  mutate(
    age = as.numeric(interview_age),
    sex = as.factor(sex),
    race = as.factor(race),
    site = as.factor(site)
  )

mu_formula <- ~ fp(age,npoly=2) + sex + race + random(site)
sigma_formula <- ~fp(age,npoly=2)
nu_formula <- ~fp(age,npoly=2)
tau_formula <- ~fp(age,npoly=2)
dist_list <- c("TF", "TF2", "JSU", "SHASH", "SEP1", "SEP2", "SEP3", "SEP4", "ST1", "ST2", "ST3", "ST4", "ST5", "PE", "PE2", "SN1",
               "SN2", "exGAUS", "EGB2", "GT", "JSUo", "NET", "SHASHo", "SHASHo2", "SST")

for (feature in feature_names) {

  aic_scores <- data.frame(distribution = dist_list, aic = NA, converged = FALSE)
  y_formula <- as.formula(paste(feature, "~", as.character(mu_formula)[2]))
  
  for (i in 1:length(dist_list)) {
    dist_name <- dist_list[i]
    fam_details <- gamlss.family(dist_name)
    
    current_sigma_formula <- if (fam_details$nopar >= 2) sigma_formula else ~1
    current_nu_formula <- if (fam_details$nopar >= 3) nu_formula else ~1
    current_tau_formula <- if (fam_details$nopar >= 4) tau_formula else ~1
    
    fit <- try(gamlss(formula = y_formula, sigma.formula = current_sigma_formula, nu.formula = current_nu_formula, tau.formula = current_tau_formula, family = eval(parse(text=paste0(dist_name, "()"))), data = model_data, n.cyc = 300), silent = TRUE)
    
    if (!inherits(fit, "try-error") && isTRUE(fit$converged)) {
      aic_scores$aic[i] <- AIC(fit)
      aic_scores$converged[i] <- TRUE
    }
  }
  
  ranked_results <- aic_scores %>% 
    filter(converged == TRUE) %>% 
    arrange(aic)
  
  output_filename <- file.path(output_dir, paste0("aic_ranking_", feature, ".csv"))
  write_csv(ranked_results, output_filename)
}


# GAMLSS Modeling
# Amyg

z_scores_df <- data.frame(src_subject_id = model_data$src_subject_id)
nfolds <- 10

feature = feature_names[1]
mu_formula <- ~  fp(age,npoly=2) + sex + random(site)
sigma_formula <- ~ fp(age,npoly=2)
nu_formula <- ~ fp(age,npoly=1)
tau_formula <- ~ fp(age,npoly=1)

best_dist_name <- "SEP3"

y_formula <- as.formula(paste(feature, "~", as.character(mu_formula)[2]))
zscores_outofsample <- matrix(data = NA, nrow = nrow(model_data), ncol = 1)
folds_idx <- createFolds(model_data[[feature]], k = nfolds)

for (k in 1:nfolds) {
  test_indices <- folds_idx[[k]]
  data_train <- model_data[-test_indices, ]
  data_test <- model_data[test_indices, ]
  
  best_fam_details <- gamlss.family(best_dist_name)
  cv_sigma_formula <- if (best_fam_details$nopar >= 2) sigma_formula else ~1
  cv_nu_formula <- if (best_fam_details$nopar >= 3) nu_formula else ~1
  cv_tau_formula <- if (best_fam_details$nopar >= 4) tau_formula else ~1
  
  model <- try(
    gamlss(
      formula = y_formula, sigma.formula = cv_sigma_formula, nu.formula = cv_nu_formula, tau.formula = cv_tau_formula,
      family = eval(parse(text = paste0(best_dist_name, "()"))),
      data = data_train, na.action = na.exclude, n.cyc = 300
    )
  )
  
  if (inherits(model, "try-error") || !isTRUE(model$converged)) {
    cat(sprintf("    -> Failed to converge for feature %s on fold %d\n", feature, k))
    next
  }
  
  params <- predictAll(model, newdata = data_test, type = "response")
  p_func_name <- paste0("p", best_dist_name)
  p_func <- get(p_func_name)
  
  p_args <- list(q = data_test[[feature]], mu = params$mu)
  if("sigma" %in% names(params)) p_args$sigma <- params$sigma
  if("nu" %in% names(params)) p_args$nu <- params$nu
  if("tau" %in% names(params)) p_args$tau <- params$tau
  
  p_values <- do.call(p_func, p_args)
  z_scores_fold <- qnorm(p_values)
  z_scores_fold[is.infinite(z_scores_fold)] <- NA
  
  zscores_outofsample[test_indices, ] <- z_scores_fold
}

z_scores_df[[paste0(feature, "_zscore")]] <- zscores_outofsample[, 1]
z_scores_df[[paste0(feature, "_dist_used")]] <- best_dist_name

# Hipp
feature = feature_names[2]
mu_formula <- ~  fp(age,npoly=2) + sex + random(site)
sigma_formula <- ~ fp(age,npoly=2)
nu_formula <- ~ fp(age,npoly=1)
tau_formula <- ~ fp(age,npoly=1)

best_dist_name <- "TF2"

y_formula <- as.formula(paste(feature, "~", as.character(mu_formula)[2]))
zscores_outofsample <- matrix(data = NA, nrow = nrow(model_data), ncol = 1)
folds_idx <- createFolds(model_data[[feature]], k = nfolds)

for (k in 1:nfolds) {
  test_indices <- folds_idx[[k]]
  data_train <- model_data[-test_indices, ]
  data_test <- model_data[test_indices, ]
  
  best_fam_details <- gamlss.family(best_dist_name)
  cv_sigma_formula <- if (best_fam_details$nopar >= 2) sigma_formula else ~1
  cv_nu_formula <- if (best_fam_details$nopar >= 3) nu_formula else ~1
  cv_tau_formula <- if (best_fam_details$nopar >= 4) tau_formula else ~1
  
  model <- try(
    gamlss(
      formula = y_formula, sigma.formula = cv_sigma_formula, nu.formula = cv_nu_formula, tau.formula = cv_tau_formula,
      family = eval(parse(text = paste0(best_dist_name, "()"))),
      data = data_train, na.action = na.exclude, n.cyc = 100
    )
  )
  
  if (inherits(model, "try-error") || !isTRUE(model$converged)) {
    cat(sprintf("    -> GAMLSS failed or did not converge for feature %s on fold %d\n", feature, k))
    next
  }
  
  params <- predictAll(model, newdata = data_test, type = "response")
  p_func_name <- paste0("p", best_dist_name)
  p_func <- get(p_func_name)
  
  p_args <- list(q = data_test[[feature]], mu = params$mu)
  if("sigma" %in% names(params)) p_args$sigma <- params$sigma
  if("nu" %in% names(params)) p_args$nu <- params$nu
  if("tau" %in% names(params)) p_args$tau <- params$tau
  
  p_values <- do.call(p_func, p_args)
  z_scores_fold <- qnorm(p_values)
  z_scores_fold[is.infinite(z_scores_fold)] <- NA
  
  zscores_outofsample[test_indices, ] <- z_scores_fold
}

z_scores_df[[paste0(feature, "_zscore")]] <- zscores_outofsample[, 1]
z_scores_df[[paste0(feature, "_dist_used")]] <- best_dist_name



# Save
output_file <- file.path(output_dir, "deltar_gamlss_zscores_BCPHCPD.csv")
write_csv(z_scores_df, output_file)



# Visualization

# Amyg

get_mode <- function(x) {
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}

pred_grid <- data.frame(age = seq(min(model_data$age), max(model_data$age), length.out = 100))
pred_grid$sex <- get_mode(model_data$sex)
pred_grid$site <- get_mode(model_data$site)

centiles_to_plot <- c(0.05, 0.25, 0.50, 0.75, 0.95)

mu_formula <- ~  fp(age,npoly=2) + sex + random(site)
sigma_formula <- ~ fp(age,npoly=2)
nu_formula <- ~ fp(age,npoly=2)
tau_formula <- ~ fp(age,npoly=2)

feature = feature_names[1]
  
dist_used_for_feature <- "ST2"
y_formula <- as.formula(paste(feature, "~", as.character(mu_formula)[2]))
  
fam_details_plot <- gamlss.family(dist_used_for_feature)
plot_sigma_formula <- if (fam_details_plot$nopar >= 2) sigma_formula else ~1
plot_nu_formula <- if (fam_details_plot$nopar >= 3) nu_formula else ~1
plot_tau_formula <- if (fam_details_plot$nopar >= 4) tau_formula else ~1
  
model_full <- gamlss(
  formula = y_formula, sigma.formula = plot_sigma_formula, nu.formula = plot_nu_formula, tau.formula = plot_tau_formula,
  family = eval(parse(text = paste0(dist_used_for_feature, "()"))),
  data = model_data, n.cyc = 300
)
  
pred_params <- predictAll(model_full, newdata = pred_grid, type = "response")
q_func_name <- paste0("q", dist_used_for_feature)
q_func <- get(q_func_name)
  
predicted_centiles <- sapply(centiles_to_plot, function(p) {
  q_args <- list(p = p, mu = pred_params$mu)
  if("sigma" %in% names(pred_params)) q_args$sigma <- pred_params$sigma
  if("nu" %in% names(pred_params)) q_args$nu <- pred_params$nu
  if("tau" %in% names(pred_params)) q_args$tau <- pred_params$tau
  do.call(q_func, q_args)
})
  
plot_data <- bind_cols(pred_grid, as.data.frame(predicted_centiles))
colnames(plot_data) <- c(colnames(pred_grid), paste0("p", centiles_to_plot * 100))
  
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = starts_with("p"), names_to = "centile", values_to = "value"
  ) %>%
  mutate(is_median = ifelse(centile == "p50", "Median (50th)", "Other Centiles"))

normative_plot <- ggplot(model_data, aes(x = age, y = .data[[feature]])) +
  geom_point(alpha = 0.3, size = 2.5, color = "grey40") +
  geom_line(data = plot_data_long, aes(x = age, y = value, group = centile, color = is_median), linewidth = 1) +
  scale_color_manual(name = "Centile", values = c("Median (50th)" = "#D55E00", "Other Centiles" = "#0072B2")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom") + ylim(-0.1,0.16)
log_scaled_plot <- normative_plot + scale_x_sqrt(
  breaks = c(0.1, 0.5, 1, 2, 5, 10, 21),
  labels = c("0.1", "0.5", "1", "2", "5", "10", "21")
  )
  
plot_output_file <- file.path(output_dir, paste0("normative_plot_", feature, ".png"))
ggsave(plot_output_file, plot = normative_plot, width = 4, height = 4, dpi = 300)
  
plot_output_file_log <- file.path(output_dir, paste0("normative_plot_log_", feature, ".png"))
ggsave(plot_output_file_log, plot = log_scaled_plot, width = 4, height = 4, dpi = 300)
}
summary(model_full)
plot(model_full)



# Hipp
model_data$deltar_h = -model_data$deltar_h


mu_formula <- ~ fp(age,npoly=1) + sex + random(site)
sigma_formula <- ~ fp(age,npoly=1)
nu_formula <- ~ poly(age,npoly=2)
tau_formula <- ~ poly(age,npoly=2)

feature = feature_names[2]

dist_used_for_feature <- "PE"
y_formula <- as.formula(paste(feature, "~", as.character(mu_formula)[2]))
  
fam_details_plot <- gamlss.family(dist_used_for_feature)
plot_sigma_formula <- if (fam_details_plot$nopar >= 2) sigma_formula else ~1
plot_nu_formula <- if (fam_details_plot$nopar >= 3) nu_formula else ~1
plot_tau_formula <- if (fam_details_plot$nopar >= 4) tau_formula else ~1
  
model_full <- gamlss(
  formula = y_formula, sigma.formula = plot_sigma_formula, nu.formula = plot_nu_formula, tau.formula = plot_tau_formula,
  family = eval(parse(text = paste0(dist_used_for_feature, "()"))),
  data = model_data, n.cyc = 300
)
  
pred_params <- predictAll(model_full, newdata = pred_grid, type = "response")
q_func_name <- paste0("q", dist_used_for_feature)
q_func <- get(q_func_name)
  
predicted_centiles <- sapply(centiles_to_plot, function(p) {
  q_args <- list(p = p, mu = pred_params$mu)
  if("sigma" %in% names(pred_params)) q_args$sigma <- pred_params$sigma
  if("nu" %in% names(pred_params)) q_args$nu <- pred_params$nu
  if("tau" %in% names(pred_params)) q_args$tau <- pred_params$tau
  do.call(q_func, q_args)
})
  
plot_data <- bind_cols(pred_grid, as.data.frame(predicted_centiles))
colnames(plot_data) <- c(colnames(pred_grid), paste0("p", centiles_to_plot * 100))
  
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = starts_with("p"), names_to = "centile", values_to = "value"
  ) %>%
  mutate(is_median = ifelse(centile == "p50", "Median (50th)", "Other Centiles"))

normative_plot <- ggplot(model_data, aes(x = age, y = .data[[feature]])) +
  geom_point(alpha = 0.3, size = 2.5, color = "grey40") +
  geom_line(data = plot_data_long, aes(x = age, y = value, group = centile, color = is_median), linewidth = 1) +
  scale_color_manual(name = "Centile", values = c("Median (50th)" = "#D55E00", "Other Centiles" = "#0072B2")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom") + ylim(-0.1,0.16)
log_scaled_plot <- normative_plot + scale_x_sqrt(
  breaks = c(0.1, 0.5, 1, 2, 5, 10, 21),
  labels = c("0.1", "0.5", "1", "2", "5", "10", "21")
  )
  
plot_output_file <- file.path(output_dir, paste0("normative_plot_", feature, ".png"))
ggsave(plot_output_file, plot = normative_plot, width = 4, height = 4, dpi = 300)
  
plot_output_file_log <- file.path(output_dir, paste0("normative_plot_log_", feature, ".png"))
ggsave(plot_output_file_log, plot = log_scaled_plot, width = 4, height = 4, dpi = 300)
}
summary(model_full)
plot(model_full)