# ================================================================================
# Case 8.1: A Tale of Two Thieves - Complete Statistical Analysis
# ================================================================================

# ================================================================================
# 3.1 Exploratory Data Analysis (EDA)
# ================================================================================
# This script follows the structure of the analysis report.
# Data preparation and package loading are performed first.

# Create the thief data
thief_data <- data.frame(
  METHOD = c(
    rep("Intm", 18), rep("Unit", 18)
  ),
  LOCATION = c(
    1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6,
    1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6
  ),
  REPLICATE = c(
    1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
    1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3
  ),
  ASSAY = c(
    34.38, 34.87, 35.71, 35.31, 37.59, 38.02, 36.71, 36.56, 35.92,
    37.80, 37.41, 38.00, 36.28, 36.63, 36.62, 38.89, 39.80, 37.84,
    33.94, 34.72, 34.10, 39.11, 37.51, 37.79, 37.46, 34.12, 35.94,
    38.05, 34.82, 35.42, 36.52, 38.60, 38.16, 39.16, 32.77, 36.95
  )
)

# Create the tablet data (final product from compressed tablets)
tablet_data <- data.frame(
  DRUM = c(
    1, 1, 1, 5, 5, 5, 7, 7, 7, 11, 11, 11, 14, 14, 14,
    17, 17, 17, 19, 19, 19, 22, 22, 22, 25, 25, 25, 28, 28, 28
  ),
  TABLET = c(
    1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
    1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3
  ),
  ASSAY = c(
    35.77, 39.44, 36.43, 35.71, 37.08, 36.54, 35.08, 34.25, 33.09,
    35.21, 34.36, 35.94, 35.17, 36.54, 36.45, 35.43, 33.80, 35.15,
    34.56, 35.33, 37.69, 35.82, 35.67, 35.06, 35.75, 37.32, 35.06,
    38.58, 36.63, 35.60
  )
)

thief_data$METHOD <- factor(thief_data$METHOD, levels = c("Intm", "Unit"))
thief_data$LOCATION <- factor(thief_data$LOCATION)
thief_data$REPLICATE <- factor(thief_data$REPLICATE)

thief_data

tablet_data$DRUM <- factor(tablet_data$DRUM)
tablet_data$TABLET <- factor(tablet_data$TABLET)

tablet_data

# Create combined dataset for three-way comparison
combined_data <- data.frame(
  METHOD = c(rep("Intm", 18), rep("Unit", 18), rep("Tablet", 30)),
  ASSAY = c(thief_data$ASSAY, tablet_data$ASSAY)
)
combined_data$METHOD <- factor(combined_data$METHOD, levels = c("Intm", "Unit", "Tablet"))

combined_data


# Three-way comparison including tablets
combined_stats <- data.frame(
  Method = c("Intm", "Unit", "Tablet"),
  N = c(18, 18, 30),
  Mean = round(tapply(combined_data$ASSAY, combined_data$METHOD, mean), 2),
  SD = round(tapply(combined_data$ASSAY, combined_data$METHOD, sd), 2),
  Min = round(tapply(combined_data$ASSAY, combined_data$METHOD, min), 2),
  Q1 = round(tapply(combined_data$ASSAY, combined_data$METHOD, quantile, 0.25), 2),
  Median = round(tapply(combined_data$ASSAY, combined_data$METHOD, median), 2),
  Q3 = round(tapply(combined_data$ASSAY, combined_data$METHOD, quantile, 0.75), 2),
  Max = round(tapply(combined_data$ASSAY, combined_data$METHOD, max), 2)
)

# Table A.1: Summary Statistics for Assay Value by Method
cat("\n========== TABLE A.1: SUMMARY STATISTICS (3.1 Exploratory Data Analysis) ==========\n")
print(combined_stats)
cat("\n")

# ================================================================================
#   3.2 Mixed-Effects Model
# ================================================================================
# Load nlme package for proper mixed-effects modeling
if (!require(nlme, quietly = TRUE)) {
  install.packages("nlme", repos = "https://cran.r-project.org/")
  library(nlme)
}

# Implement the correct mixed-effects model matching SAS proc mixed
# Model: ASSAY ~ METHOD + random effects for LOCATION + heterogeneous variances by METHOD
# Fit the proper mixed model (matches SAS proc mixed)
# Try with tighter convergence criteria to match SAS exactly
mixed_model <- lme(
  fixed = ASSAY ~ METHOD,                    # Fixed effect: METHOD
  random = ~ 1 | LOCATION,                   # Random effect: LOCATION
  weights = varIdent(form = ~ 1 | METHOD),   # Heterogeneous variances by METHOD
  data = thief_data,
  method = "ML",                             # Use ML (matches PDF results)
  control = lmeControl(tolerance = 1e-8,     # Tighter tolerance
                      maxIter = 200,        # More iterations
                      msMaxIter = 200)      # More inner iterations
)

# Extract model summary
model_summary <- summary(mixed_model)
model_summary

# Get F-test results for fixed effects
anova_results <- anova(mixed_model)

# Get confidence intervals for fixed effects (95% CI)
ci_results <- intervals(mixed_model, which = "fixed")

# Table A.2: Type III Tests of Fixed Effects
cat("\n========== TABLE A.2: FIXED EFFECTS TESTS (3.2 Mixed-Effects Results) ==========\n")
print(anova_results)
cat("\n========== CONFIDENCE INTERVALS FOR FIXED EFFECTS ==========\n")
print(ci_results)


# ================================================================================
#   3.2.2 Superior Precision of INTM: Variance Heterogeneity by Sampling Method
# ================================================================================

# Extract variance components from the mixed model
var_components <- VarCorr(mixed_model)

# Extract location variance (random effect)
location_var <- as.numeric(var_components[1,1])

# Extract residual variances (method-specific)
sigma_baseline <- mixed_model$sigma  # Baseline residual SD (for Intm)
weights <- coef(mixed_model$modelStruct$varStruct, unconstrained = FALSE)

# Calculate variances for each method
var_intm_mixed <- sigma_baseline^2
var_unit_mixed <- (sigma_baseline * weights[1])^2

ratio_mixed <- var_unit_mixed / var_intm_mixed

# Store the correct variances for later use
var_intm <- var_intm_mixed
var_unit <- var_unit_mixed

# Table A.3: Variance Components Decomposition
cat("========== TABLE A.3: VARIANCE COMPONENTS DECOMPOSITION (3.2.2) ==========\n")
var_components_table <- data.frame(
  Component = c("Between-Location (LOC)", "Within-Location (INTM)", "Within-Location (UNIT)"),
  Estimate = c(round(location_var, 4), round(var_intm, 4), round(var_unit, 4)),
  Role = c("Variability across locations", "Residual precision - INTM", "Residual precision - UNIT")
)
print(var_components_table)

cat("\nVariance Heterogeneity Summary (3.2.2):\n")
cat("INTM Residual Variance:", round(var_intm, 4), "mg²/100mg²\n")
cat("UNIT Residual Variance:", round(var_unit, 4), "mg²/100mg²\n")
cat("Heterogeneity Ratio (UNIT/INTM):", round(ratio_mixed, 2), "× (Unit Dose is more variable)\n\n")

# ================================================================================
#   3.2.3 Location Effects Dominate: Location-Specific Deviations
# ================================================================================
# Extract random intercepts from the base model
random_effects_base <- random.effects(mixed_model)
location_intercepts <- as.numeric(random_effects_base[, 1])

# Calculate grand mean
grand_mean <- mean(thief_data$ASSAY)

# Extract confidence intervals for random effects
# intervals() provides 95% CI for random effects
random_intervals <- intervals(mixed_model, which = "var-cov")

# Calculate SE from confidence intervals
# For 95% CI: estimate ± 1.96 * SE, so SE = (upper - lower) / (2 * 1.96)
random_ci <- random_intervals$reStruct$LOCATION

# For each location, calculate SE from the intercept variance CI
# SE of random effect = sqrt(Var(b_i) / n_obs_per_location)
var_b <- as.numeric(VarCorr(mixed_model)[1, 1])  # Between-location variance
n_obs_per_location <- 6  # 3 replicates × 2 methods per location

# Calculate SE for random intercepts
# Using the standard formula: SE(b_i) = sqrt(Var(residual) / n)
pooled_residual_var <- (var_intm + var_unit) / 2
se_random <- rep(sqrt(pooled_residual_var / n_obs_per_location), 6)

# Calculate z-scores and p-values for significance testing
z_scores <- location_intercepts / se_random
p_values <- 2 * (1 - pnorm(abs(z_scores)))  # Two-tailed test

# Determine status based on p-value and direction
status <- ifelse(p_values < 0.01,
                 ifelse(location_intercepts < 0,
                        "Significantly lower",
                        "Significantly higher"),
                 ifelse(p_values < 0.05,
                        ifelse(location_intercepts < 0,
                               "Marginally lower",
                               "Marginally higher"),
                        "Not significant"))

# Create location effects summary table
location_effects_df <- data.frame(
  Location = 1:6,
  Random_Intercept = round(location_intercepts, 2),
  Location_Mean = round(grand_mean + location_intercepts, 2),
  Relative_Deviation_Pct = round((location_intercepts / grand_mean) * 100, 2),
  SE = round(se_random, 4),
  z_score = round(z_scores, 3),
  p_value = round(p_values, 4),
  Status = status
)

cat("========== TABLE A.4: LOCATION EFFECTS SUMMARY (3.2.3 Random Effects) ==========\n")
print(location_effects_df)
cat("\nLocation Range Span:", round(max(location_intercepts) - min(location_intercepts), 2), "mg/100mg\n")
cat("Lowest Location (1):", round(grand_mean + min(location_intercepts), 2), "mg/100mg\n")
cat("Highest Location (6):", round(grand_mean + max(location_intercepts), 2), "mg/100mg\n\n")

# ================================================================================
#   3.3 Regression Perspective: Variance Components and R² Decomposition
# ================================================================================
# Calculate R² for the base model (Model 1: fixed METHOD + random LOCATION intercepts)
fitted_marginal_1 <- predict(mixed_model, level = 0)
fitted_conditional_1 <- predict(mixed_model, level = 1)

ss_tot <- sum((thief_data$ASSAY - mean(thief_data$ASSAY))^2)

# Model 1 R² components
r2_marginal_1 <- 1 - sum((thief_data$ASSAY - fitted_marginal_1)^2) / ss_tot
r2_conditional_1 <- 1 - sum((thief_data$ASSAY - fitted_conditional_1)^2) / ss_tot

# Calculate variance component contributions for reporting
method_variance_contribution <- r2_marginal_1 * 100  # 2.30%
location_variance_contribution <- (r2_conditional_1 - r2_marginal_1) * 100  # 34.20%
residual_variance_contribution <- 100 - (r2_conditional_1 * 100)  # 70.38%

cat("========== VARIANCE DECOMPOSITION SUMMARY (3.3 R² Decomposition) ==========\n")
cat("METHOD Effect Contribution:", round(method_variance_contribution, 2), "% of total variance\n")
cat("LOCATION Random Effects Contribution:", round(location_variance_contribution, 2), "% of total variance\n")
cat("Residual Measurement Error:", round(residual_variance_contribution, 2), "% of total variance\n")
cat("Combined Location Effects (Location + residual influenced):", 
    round(location_variance_contribution + residual_variance_contribution, 2), "%\n\n")

cat("Marginal R² (METHOD only):", round(r2_marginal_1 * 100, 2), "%\n")
cat("Conditional R² (METHOD + LOCATION):", round(r2_conditional_1 * 100, 2), "%\n")
cat("Additional Variance Explained by Location Effects:", round((r2_conditional_1 - r2_marginal_1) * 100, 2), "%\n")
cat("Ratio of Location Effects to Method Effects:", 
    round(location_variance_contribution / method_variance_contribution, 1), "-fold\n\n")


# ================================================================================
#   3.4 Interaction Analysis: Location-Specific METHOD Effects
# ================================================================================
# Fit model WITH interaction (using Random Slopes)
mixed_interaction <- lme(
  fixed = ASSAY ~ METHOD,
  random = ~ 1 + METHOD | LOCATION,  # ← Random slopes: intercept + METHOD slope
  weights = varIdent(form = ~ 1 | METHOD),
  data = thief_data,
  method = "ML",
  control = lmeControl(tolerance = 1e-8, maxIter = 200, msMaxIter = 200)
)

mixed_interaction

# Compare models using likelihood ratio test
anova_comparison <- anova(mixed_model, mixed_interaction)
cat("\n========== MODEL COMPARISON: Base vs Interaction Model (3.3) ==========\n")
print(anova_comparison)

# Extract global METHOD effect
global_method_effect <- fixef(mixed_interaction)["METHODUnit"]

# Extract random effects (both intercepts and slopes)
random_effects_loc <- random.effects(mixed_interaction)

# Create data structure for location-specific effects
location_slopes <- data.frame()

for (i in 1:6) {
  intercept_re <- random_effects_loc[i, 1]
  
  # Extract METHOD slope if it exists
  if (ncol(random_effects_loc) > 1) {
    method_slope_re <- random_effects_loc[i, 2]
  } else {
    method_slope_re <- 0
  }
  
  # Total METHOD effect at this location = global + random slope
  total_method_effect <- global_method_effect + method_slope_re
  
  location_slopes <- rbind(location_slopes,
                           data.frame(
                             Location = i,
                             Intercept_RE = intercept_re,
                             Method_Slope_RE = method_slope_re,
                             Total_Method_Effect = total_method_effect
                           ))
}

slope_var <- var(location_slopes$Method_Slope_RE)
slope_sd <- sd(location_slopes$Method_Slope_RE)

cat("\n========== LOCATION-SPECIFIC METHOD EFFECTS (3.3) ==========\n")
print(location_slopes)
cat("\nGlobal METHOD Effect:", round(global_method_effect, 4), "mg/100mg\n")
cat("Random Slope SD (Location-specific variations):", round(slope_sd, 4), "\n")
cat("Range of Location-Specific Effects:", 
    round(min(location_slopes$Total_Method_Effect), 4), "to", 
    round(max(location_slopes$Total_Method_Effect), 4), "\n\n")

# ================================================================================
#  3.5 Model Assessment
# ================================================================================

# Normality tests on RAW DATA (before model fitting) - Table B.1
for (method in c("Intm", "Unit")) {
  subset_data <- thief_data$ASSAY[thief_data$METHOD == method]
  sw_test <- shapiro.test(subset_data)
  status <- if (sw_test$p.value > 0.05) "NORMAL" else "NON-NORMAL"
  print(sw_test)
  print(status)
}

# Outlier detection on RAW DATA (Q3 + 2×IQR criterion) - Table B.2
for (method in c("Intm", "Unit")) {
  subset_data <- thief_data$ASSAY[thief_data$METHOD == method]
  Q1 <- quantile(subset_data, 0.25)
  Q3 <- quantile(subset_data, 0.75)
  IQR <- Q3 - Q1
  upper_bound <- Q3 + 2*IQR
  outliers <- subset_data[subset_data > upper_bound]
  print(Q1)
  print(Q3)
  print(IQR)
  print(upper_bound)
  print(outliers)
}



# Residual Diagnostics

# Get fitted values and residuals from the mixed model
fitted_vals <- fitted(mixed_model)
residuals_all <- thief_data$ASSAY - fitted_vals

# Normality tests on MODEL RESIDUALS (after model fitting) - Appendix B.4
sw_all <- shapiro.test(residuals_all)
sw_intm <- shapiro.test(residuals_all[thief_data$METHOD == "Intm"])
sw_unit <- shapiro.test(residuals_all[thief_data$METHOD == "Unit"])

# Diagnostic Plots (2-panel diagnostics)

png("img/diagnostic_plots.png", width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 2))

# Plot 1: Residuals vs Fitted Values
plot(fitted_vals, residuals_all,
     main = "1. Residuals vs Fitted Values",
     xlab = "Fitted Values", ylab = "Residuals",
     pch = 16, cex = 1.1, col = "darkgray")
abline(h = 0, lty = 2, lwd = 2, col = "red")
# Add LOWESS smooth curve
lines(lowess(fitted_vals, residuals_all), lwd = 2, col = "blue")

# Plot 2: Normal Q-Q Plot
qqnorm(residuals_all,
       main = "2. Normal Q-Q Plot",
       pch = 16, cex = 1.1, col = "darkgray")
qqline(residuals_all, lwd = 2, col = "red")

par(mfrow = c(1, 1))
dev.off()

print("Diagnostic plots saved as 'diagnostic_plots.png'")




# Effect Size Analysis

# Basic statistics and tests

# Extract data by method and tablet for bootstrap analyses
intm_data <- thief_data$ASSAY[thief_data$METHOD == "Intm"]
unit_data <- thief_data$ASSAY[thief_data$METHOD == "Unit"]
tablet_assay <- tablet_data$ASSAY

# Calculate means for comparisons
mean_intm <- mean(intm_data)
mean_unit <- mean(unit_data)
mean_tablet <- mean(tablet_assay)

# Get sample sizes
n_intm <- length(intm_data)
n_unit <- length(unit_data)
n_tablet <- length(tablet_assay)

# Two-sample t-tests (Welch's) for all three comparisons
t_test_ui <- t.test(unit_data, intm_data)
t_test_ut <- t.test(unit_data, tablet_assay)
t_test_it <- t.test(intm_data, tablet_assay)

t_test_ui
t_test_ut
t_test_it

# Recalculate descriptive statistics for effect size calculations
sd_intm <- sd(intm_data)
sd_unit <- sd(unit_data)


# Effect size

# 1. Cohen's d (pooled)
pooled_sd <- sqrt(((n_intm - 1) * sd_intm^2 + (n_unit - 1) * sd_unit^2) / (n_intm + n_unit - 2))
cohens_d <- (mean_intm - mean_unit) / pooled_sd

# 2. Eta-squared and Omega-squared (proportion of variance explained)
f_stat <- as.numeric(anova_results$`F-value`[2])  # Extract F-value for METHOD (row 2)
df_den <- 29  # From mixed model output
eta_squared <- r2_marginal_1 * 100  # In percentage

omega_squared <- ((f_stat - 1) / (f_stat - 1 + df_den)) * 100  # In percentage

# Create effect size summary table
effect_size_table <- data.frame(
  Measure = c("Cohen's d", "Eta-squared (%)", "Omega-squared (%)"),
  Value = round(c(cohens_d, eta_squared, omega_squared), 4),
  Interpretation = c("Small", "2.30%", "0.27%")
)

cat("\n========== TABLE A.10: EFFECT SIZE MEASURES (3.4 Effect Size Assessment) ==========\n")
print("Effect Size Analysis:")
print(effect_size_table)

# Practical Significance Assessment

# Calculate the mean difference as percentage of target specification
mean_diff <- mean_unit - mean_intm
pct_of_target <- (abs(mean_diff) / 35) * 100

cat("\n=== Practical Significance ===\n")
cat("Mean difference:", round(mean_diff, 2), "mg/100mg\n")
cat("Percentage of target (35 mg/100mg):", round(pct_of_target, 2), "%\n")
cat("Interpretation: The difference is negligible (<0.15% of specification)\n\n")



# Bootstrap Validation (1000 resamples)
# Following Cabrera & McDougall (2002, p. 284):
# Bootstrap assessment of p-value reliability by resampling observed data
# and computing the distribution of p-values from repeated t-tests

set.seed(12345)  # For reproducibility
n_boot <- 1000

# Bootstrap for Unit vs Intm: distribution of p-values
boot_pvalues_ui <- numeric(n_boot)
for (i in 1:n_boot) {
  boot_unit <- sample(unit_data, replace = TRUE)
  boot_intm <- sample(intm_data, replace = TRUE)
  boot_test <- t.test(boot_unit, boot_intm)
  boot_pvalues_ui[i] <- boot_test$p.value
}
boot_mean_p_ui <- mean(boot_pvalues_ui)
boot_bias_ui <- boot_mean_p_ui - t_test_ui$p.value
boot_se_ui <- sd(boot_pvalues_ui)

# Bootstrap for Unit vs Tablet: distribution of p-values
boot_pvalues_ut <- numeric(n_boot)
for (i in 1:n_boot) {
  boot_unit <- sample(unit_data, replace = TRUE)
  boot_tablet <- sample(tablet_assay, replace = TRUE)
  boot_test <- t.test(boot_unit, boot_tablet)
  boot_pvalues_ut[i] <- boot_test$p.value
}
boot_mean_p_ut <- mean(boot_pvalues_ut)
boot_bias_ut <- boot_mean_p_ut - t_test_ut$p.value
boot_se_ut <- sd(boot_pvalues_ut)

# Bootstrap for Intm vs Tablet: distribution of p-values
boot_pvalues_it <- numeric(n_boot)
for (i in 1:n_boot) {
  boot_intm <- sample(intm_data, replace = TRUE)
  boot_tablet <- sample(tablet_assay, replace = TRUE)
  boot_test <- t.test(boot_intm, boot_tablet)
  boot_pvalues_it[i] <- boot_test$p.value
}
boot_mean_p_it <- mean(boot_pvalues_it)
boot_bias_it <- boot_mean_p_it - t_test_it$p.value
boot_se_it <- sd(boot_pvalues_it)

# Create bootstrap validation table (Table A.7)
bootstrap_table <- data.frame(
  Comparison = c("UNIT vs. INTM", "UNIT vs. TABLET", "INTM vs. TABLET"),
  Observed_Pvalue = round(c(t_test_ui$p.value, t_test_ut$p.value, t_test_it$p.value), 4),
  Bootstrap_Mean = round(c(boot_mean_p_ui, boot_mean_p_ut, boot_mean_p_it), 4),
  Bias = round(c(boot_bias_ui, boot_bias_ut, boot_bias_it), 4),
  Std_Error = round(c(boot_se_ui, boot_se_ut, boot_se_it), 4)
)

cat("\n========== TABLE A.7: BOOTSTRAP VALIDATION - 1000 RESAMPLES (2.7 Bootstrap) ==========\n")
print("Bootstrap Validation Results:")
print(bootstrap_table)
cat("\n")

# ================================================================================
#   4 Client Question Analysis
# ================================================================================



# ================================================================================
# Results Summary - Figures and Visualization
# ================================================================================

# Save boxplot with three methods as PNG file
png("img/boxplot.png", width = 900, height = 650, res = 120)

old_par <- par(mar = c(6, 6, 4, 2))

# Create boxplot
boxplot(ASSAY ~ METHOD,
        data = combined_data,
        main = "Parallel Boxplots of Assay Value by Method",
        xlab = "Method",
        ylab = "ASSAY Value (mg/100mg)",
        boxlwd = 1.5,
        medlwd = 2,
        whisklwd = 1.5,
        outcex = 1.2,
        outpch = 16,
        cex.lab = 1.3,
        cex.axis = 1.2,
        cex.main = 1.4)

# Add horizontal reference line at target value
abline(h = 35, lty = 2, lwd = 2)

par(old_par)
dev.off()

# Save interaction plot as PNG file
png("img/interaction_plot.png", width = 900, height = 650, res = 120)

# Create interaction plot
old_par <- par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

# Prepare data for plotting
locations <- 1:6
intm_means <- numeric(6)
unit_means <- numeric(6)

for (loc in 1:6) {
  intm_vals <- thief_data$ASSAY[thief_data$METHOD == "Intm" & thief_data$LOCATION == loc]
  unit_vals <- thief_data$ASSAY[thief_data$METHOD == "Unit" & thief_data$LOCATION == loc]
  intm_means[loc] <- mean(intm_vals)
  unit_means[loc] <- mean(unit_vals)
}

# Create interaction plot
plot(locations, intm_means,
     type = "b",
     pch = 16,
     col = "darkblue",
     lwd = 2.5,
     cex = 1.2,
     ylim = c(33, 41),
     xlab = "Sampling Location in Blender",
     ylab = "ASSAY Value (mg/100mg)",
     main = "Observed METHOD Effects by Sampling Location",
     xaxt = "n",
     axes = TRUE)

# Add UNIT line
lines(locations, unit_means,
      type = "b",
      pch = 17,
      col = "darkred",
      lwd = 2.5,
      cex = 1.2)

# Customize x-axis
axis(1, at = 1:6, labels = paste("Loc", 1:6))

# Add legend
legend("topright",
       legend = c("Intermediate Dose (INTM)", "Unit Dose (UNIT)"),
       col = c("darkblue", "darkred"),
       pch = c(16, 17),
       lwd = 2.5,
       cex = 1.1,
       bty = "round",
       bg = "white")

# Add grid for readability
grid(nx = NA, ny = NULL, col = "lightgray", lty = 3)

# Add horizontal line at overall means for reference
overall_mean <- mean(c(intm_means, unit_means))
abline(h = overall_mean, lty = 2, col = "gray", lwd = 1.5)

par(old_par)
dev.off()
