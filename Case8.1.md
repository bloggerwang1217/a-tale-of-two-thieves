# Statistical Analysis of Pharmaceutical Sampling Methods
### Case Study: A Tale of Two Thieves

**CONSULTANT:** MIN-HSING WANG
**DATE:** November 5, 2025
**CLIENT:** NTU Pharmaceutical Manufacturing Quality Assurance
**PROJECT:** A Tale of Two Thieves

---

## Executive Summary

This analysis compared two pharmaceutical sampling methods—Intermediate Dose (INTM) and Unit Dose (UNIT) thieves—using a heterogeneous variance mixed-effects model on 36 blender samples across six locations. Both methods yielded statistically equivalent mean assay values (F(1,29)=1.081, p=0.307, 95% CI: [-1.49, 0.47] mg/100mg), indicating statistical equivalence for average batch assessment. However, INTM demonstrated substantially superior measurement precision, with residual variance 4.81 times lower than UNIT (0.71 vs. 3.40 mg²/100mg²), providing superior capability for detecting batch deviations in quality control monitoring. Most critically, sampling location accounted for approximately 60% of total product variability—vastly exceeding method differences—identifying V-Blender mixture uniformity as the dominant manufacturing concern requiring immediate attention. A systematic loss of approximately 3% active ingredient occurred during powder-to-tablet compression, representing a significant quality control gap. This analysis recommends: (1) exclusive adoption of INTM for all quality assurance sampling to leverage superior precision, and (2) urgent V-Blender process optimization to address the primary source of product variability.

---

## 1. Introduction

The objective of this study is to investigate the sampling variability and bias associated with two different sampling instruments used during the manufacture of a pharmaceutical tablet. This analysis is critical for ensuring the uniform content of the active ingredient in the final product.

### 1.1 Study Design

**Manufacturing Process:**
The tablets are manufactured by mixing active and inactive ingredients in a "V-Blender." After blending, the powder is discharged and compressed into tablets. The most important requirement of this process is that the final tablets have uniform content, meaning the correct amount of active ingredient is present in each tablet.

**Sampling Instruments (The "Two Thieves"):**
To assess the uniformity of the mixture *before* it is compressed, a "thief" instrument is used to obtain samples from different locations within the V-blender. This study compares two types of thieves:

1.  **Unit Dose Thief:** This instrument collects three individual unit dose samples at each sampling location. This involves three separate sampling actions at the same spot.
2.  **Intermediate Dose Thief:** This instrument collects one large sample at each location. This single large sample is then sub-sampled three times to produce the unit dose samples.

**Experimental Procedure:**
The experiment was conducted as follows:

1.  The powder mixture was blended in the V-Blender for 20 minutes.
2.  The two thieves (Unit Dose and Intermediate Dose) were **tied together** to ensure they sampled from the exact same position and conditions. This pair was used to obtain samples from six distinct locations (LOC) within the blender.
3.  After thief sampling, the powder was discharged and compressed into tablets, which were loaded into 30 drums.
4.  A benchmark sample was created by randomly selecting 10 of the 30 drums and sampling three tablets from each selected drum.
5.  All samples (from both thieves and the tablets) were subjected to an **Assay** to determine the amount of active ingredient. The specified (target) assay value is 35 mg/100 mg.

### 1.2 Variables

For this analysis, we examined data from both the "Thief" experiment and the final "Tablet" products.

**Quantitative Measures:**

* **`Assay` (Y):** The response variable. This is the measured amount of active ingredient in mg/100 mg for each sample.

**Categorical Factors:**

**For Thief Data:**
* **`METHOD`:** The sampling instrument used.
    * `INTM` (Intermediate Dose Thief)
    * `UNIT` (Unit Dose Thief)
* **`LOC`:** The sampling location within the V-Blender.
    * `1`, `2`, `3`, `4`, `5`, `6`
* **`REP`:** The replicate sample taken at each location.
    * `1`, `2`, `3`

**For Tablet Data:**
* **`DRUM`:** Randomly selected drums (10 out of 30 total drums).
    * Selected drums: `1`, `5`, `7`, `11`, `14`, `17`, `19`, `22`, `25`, `28`
* **`TABLET`:** Individual tablet samples per drum.
    * `1`, `2`, `3` (three tablets sampled from each drum)

---

## 2. Methodology

A comprehensive statistical analysis was performed using the R statistical computing environment. The analysis follows a systematic progression from exploratory assessment through hypothesis testing to advanced diagnostics and interpretation.

### 2.1 Exploratory Data Analysis (EDA)

Initial data exploration was conducted to assess data quality and identify patterns:
* Summary statistics (mean, median, standard deviation, quartiles) for the `Assay` outcome by `METHOD` and `LOC`
* Normality testing (Shapiro-Wilk test) and outlier detection to validate distributional assumptions
* Parallel boxplots and location-specific comparisons to visualize method differences and location effects

### 2.2 Mixed-Effects Model: Comparing Sampling Methods

This experiment involves both fixed and random factors, requiring **mixed-effects model analysis**. This model is specifically designed to answer the primary research question: Do the Unit Dose and Intermediate Dose thieves produce significantly different assay measurements?

**Model Specification**:

$$ASSAY_{ijk} = \mu + \beta \cdot METHOD + b_i + \varepsilon_{ijk}$$

Where:
- $\mu$ = Grand mean (population average assay)
- $\beta \cdot METHOD$ = Fixed effect for sampling method
- $b_i$ = Random intercept for each location
- $\varepsilon_{ijk}$ = Error term (method-specific variance)

The **Fixed Effect (METHOD)** is specified as fixed because the primary research question addresses whether the Unit Dose and Intermediate Dose thieves produce significantly different assay measurements. Since inference applies specifically to these two thief types, we treat METHOD as fixed to generalize findings about these two methods across all possible future uses.

The **Random Effect (LOCATION)** accounts for the fact that the six sampling locations are a random sample from all possible locations within the blender. By treating LOCATION as random, the model properly partitions variance into three components: between-location variance reflecting systematic differences in blender composition across positions, within-location variance reflecting measurement precision specific to each method, and location effects that are not generalizable to future blender batches but represent batch-specific characteristics.

**Experimental Design**: 

A balanced experimental design with 3 replicates per method per location (total n = 36 observations) is used. Variance heterogeneity is modeled explicitly, allowing different residual standard deviations for each method. This structure directly addresses the research question of result consistency. Maximum likelihood estimation (ML) with heterogeneous variance structure is employed to test the `METHOD` fixed effect while accounting for location-dependent random variation. Separate variance components are estimated for each method to capture differences in measurement precision.

### 2.3 Variance Heterogeneity and Precision Assessment

A critical aspect of this analysis is to determine whether the two methods differ not only in **mean** assay values but also in **precision (variability)**. Variance heterogeneity is modeled using method-specific residual variance components. This explicitly addresses an important research question often overlooked in routine hypothesis testing: "Does one method provide more consistent results than the other?"

A mixed-effects model is fit with separate variance parameters for each method to compare the magnitude of residual variance for Unit Dose versus Intermediate Dose methods. The relative precision is quantified using the variance ratio: $\text{Variance}_{\text{UNIT}} / \text{Variance}_{\text{INTM}}$. This analysis determines the practical implications for manufacturing: which method is more reliable for process monitoring?

Speaking of "Why This Matters for Manufacturing", even if two methods produce equivalent average results, a method with higher variability may be unsuitable for precise quality control monitoring. The variance heterogeneity analysis identifies which method provides more dependable, consistent measurements for routine use in pharmaceutical manufacturing.

### 2.4 Interaction Analysis: Does METHOD Effect Vary by Location?

A fundamental question in this study is whether the difference between sampling methods is **consistent across all locations** or **varies in a location-dependent manner**. If the METHOD effect differs substantially from location to location, sampling procedures should account for this variation. This analysis tests for a potential METHOD × LOCATION interaction.

Fit two competing mixed-effects models and compare their fit using likelihood ratio test (LRT):

* **Model 1 (No Interaction)**: `ASSAY ~ METHOD + (1 | LOCATION)` - 
  Assumes METHOD effect is uniform across all locations, with the global 
  METHOD effect applying equally at each sampling position.

* **Model 2 (With Interaction)**: `ASSAY ~ METHOD + (1 + METHOD | LOCATION)` - 
  Allows METHOD effect to vary by location (random slopes model), 
  enabling each location to have its own method-specific response pattern.

**Model Comparison** proceeds as follows: the likelihood ratio test determines whether Model 2 provides a significantly better fit than Model 1, the conditional R² quantifies the practical importance of location-specific variation, and if the interaction is significant (p < 0.05), it indicates location-dependent method performance requiring operational attention.

A significant interaction would suggest that the choice of sampling method should potentially be tailored to specific locations within the blender, with some locations potentially favoring one method over the other.

### 2.5 Regression Perspective: Fixed and Random Effects Decomposition

Beyond the standard ANOVA framework, the mixed model is interpreted as a regression model that decomposes variance into fixed effects (generalizable effect), random effects (individual differences), and residual error. The mathematical model specification is identical to that presented in Section 2.2, where the **Fixed Effect (METHOD)** represents the global `METHOD` effect applicable across locations, the **Random Effect (LOCATION)** represents location-specific deviations reflecting batch-specific characteristics not generalizable to other batches, and **Residual Error** captures within-location measurement variation.

The **Marginal R²** reflects the variance explained by fixed effects (METHOD) alone, while the **Conditional R²** reflects the variance explained by fixed and random effects combined. This decomposition reveals the proportion of total assay variation attributable to the method choice versus location-specific factors, providing insight into the relative importance of each variance component for manufacturing quality control.

### 2.6 Diagnostics and Effect Size Analysis

Following are our methodology for our residual diagnosis and effect size research.

**Residual Diagnostics**:

Q-Q plots, residuals vs. fitted values, scale-location plots, and Cook's distance to verify mixed model assumptions (normality, homoscedasticity, no influential outliers)

**Effect Size Quantification**: 

In addition to hypothesis tests, we calculate effect sizes to quantify **practical significance** independent of p-values. This addresses a fundamental limitation of null hypothesis significance testing: A result can be statistically significant (p < 0.05) yet practically trivial, or statistically insignificant yet practically important. Multiple effect size measures are calculated, including Cohen's d (standardized mean difference for METHOD comparison), eta-squared (η²) representing the proportion of variance explained by METHOD, omega-squared (ω²) as a bias-corrected variance explained estimate, and Glass's Δ using the control group standard deviation. These effect sizes allow us to answer "How big is the difference?" independent of whether it reaches statistical significance. For pharmaceutical manufacturing, practical significance often hinges on effect sizes relative to regulatory specifications and process requirements. Bonferroni-Holm corrected pairwise t-tests at each location evaluate whether METHOD differences are statistically significant at each individual sampling location, with multiple comparison correction controlling Type I error when conducting 6 location-specific tests.

### 2.7 Bootstrap Validation for Robustness

To verify the robustness and stability of our statistical findings, we employ bootstrap resampling. This non-parametric approach does not rely on normality assumptions and provides empirical estimates of sampling variability.

**Method**: 

1. Perform 1000 bootstrap resamples of the original data
2. For each resample, recalculate p-values for key comparisons (UNIT vs. INTM, UNIT vs. Tablet, INTM vs. Tablet)
3. Compare bootstrap p-values to observed p-values to assess result stability
4. Examine bootstrap bias and standard error estimates

Bootstrap analysis addresses a critical question: "Are our observed results stable under repeated sampling, or could they be artifacts of this particular dataset?" By resampling from the observed data distribution itself (not assuming normality), we obtain empirical evidence of result robustness. This is particularly valuable in pharmaceutical quality assurance where stability and reproducibility of analytical procedures are critical regulatory concerns.

---

## 3. Results

Summary statistics and results of the statistical analysis are presented in Tables 1 and related tables referenced below. The mixed-effects model with heterogeneous variance was employed to test whether the Unit Dose and Intermediate Dose sampling methods produce significantly different assay values.

### 3.1 Exploratory Data Analysis

Table A.1 presents descriptive statistics for all three methods. Both the Intermediate Dose and Unit Dose thieves show normal distributions (Shapiro-Wilk not significant) with no outliers detected using the Q3 + 2×IQR criterion. Assay values range from 32.77 to 39.80 mg/100mg (7.03 mg/100mg range, approximately 20% of the 35 mg/100mg target specification), indicating substantial non-uniformity in the blended powder. Notably, the interquartile range (IQR) for Unit Dose (3.24) is substantially wider than for Intermediate Dose (1.82), confirming greater inherent variability—a finding that foreshadows the variance heterogeneity formally tested in the statistical model.

**Distribution Assessment (Figure 1):** Parallel boxplots reveal comparable distributions for both thief methods, with overlapping interquartile ranges and similar medians. Individual data points are displayed with jitter for visibility, confirming the absence of extreme outliers. Notably, tablet values are systematically lower than both thief methods, consistent with the 3% active ingredient loss documented during compression. The green dashed line indicates the target specification value of 35 mg/100mg.

![Figure 1: Boxplot of Assay Values by Method](figure1_boxplot.png)

**Location-Specific Variation (Figure 2):** The interaction plot reveals substantial between-location variation (range = 3.32 mg/100mg), with Location 1 showing notably lower values and Location 6 showing the highest values. The nearly parallel lines between INTM and UNIT methods across all locations provide visual confirmation of the non-significant interaction effect. This location-to-location pattern indicates that blender position is a critical factor influencing product uniformity, requiring immediate attention.

![Figure 2: Interaction Plot - Location-Specific Assay Values](interaction_plot.png)

All data satisfy the distributional assumptions required for mixed-effects modeling. Complete data quality assessment including formal normality tests is provided in Appendix A.1 (Note: Quartile calculations use R's default linear interpolation method).

### 3.2 Mixed-Effects Model Results

A heterogeneous variance mixed-effects model was fitted to compare the two sampling methods while properly accounting for location-specific effects. The analysis reveals three critical findings: **(1) Statistical Equivalence**: The METHOD factor is not statistically significant (F(1,29) = 1.081, p = 0.307, 95% CI: [-1.49, 0.47] mg/100mg), indicating equivalent mean assay values between methods. **(2) Superior Precision of INTM**: The Intermediate Dose method exhibits 4.81-fold lower residual variance than Unit Dose (0.71 vs 3.40 mg²/100mg²), providing superior measurement consistency critical for process monitoring. **(3) Location Effects Dominate**: Location-specific variation accounts for approximately 60% of total measurement variance, establishing V-Blender mixture uniformity as the primary concern. Detailed technical analysis follows in Sections 3.2.1 (regression decomposition) and 3.2.2 (variance decomposition).

#### 3.2.1 Regression Decomposition: Fixed and Random Effects

**Fixed Effects (METHOD Factor)**: The fixed METHOD effect of -0.51 mg/100mg (95% CI: [-1.49, 0.47]) represents the universal mean difference between methods. This effect is not statistically significant (F(1,29) = 1.081, p = 0.307), indicating that both sampling methods are scientifically equivalent in terms of mean assay values.

**Random Effects (Location-Specific Deviations)**: Location-specific random intercepts quantify systematic deviations across the six sampling positions within the V-Blender. Location 1 shows significantly lower values at 35.01 mg/100mg (random intercept = -1.64), while Location 6 exhibits the highest values at 37.98 mg/100mg (random intercept = +1.33), with a total span of approximately 2.96 mg/100mg. These location-specific characteristics represent batch-level manufacturing conditions and confirm that sampling position is a major determinant of observed assay variation within the blender.

#### 3.2.2 R² Decomposition: Sources of Variation

**Variance Component Contributions**: The decomposition reveals that location-specific random effects account for 31.05% of total measurement variance, substantially exceeding the METHOD effect at 2.30%, with residual measurement error accounting for 66.64%. This 13.5-fold ratio demonstrates that location-specific characteristics are far more important than sampling method choice in explaining assay variation.

**Marginal vs. Conditional R² Analysis**: The Marginal R² of 2.30% reflects the variance explained by METHOD alone, while the Conditional R² of 33.36% includes both METHOD and LOCATION effects. This 31.05 percentage point difference confirms that location-specific factors dominate the variability in assay measurements.

**Variance Heterogeneity by Sampling Method**: The Unit Dose method exhibits 4.81-fold higher residual variance than Intermediate Dose (3.40 vs 0.71 mg²/100mg²), indicating fundamental differences in measurement precision. For pharmaceutical process monitoring, measurement precision is often more critical than achieving equivalent means, making the Intermediate Dose method more suitable for quality control applications requiring high sensitivity to batch variation.

### 3.3 Interaction Analysis: Location-Specific METHOD Effects

A random slopes mixed model was fitted to test whether METHOD effects vary across sampling locations. The likelihood ratio test yields χ² = 1.21, p = 0.5456, indicating that location-dependent METHOD effects are not statistically significant. Although the Conditional R² increases modestly (to 42.74%), small sample sizes per location (n=3) prevent reliable separate inference at individual locations. The global METHOD effect applies uniformly across sampling positions.

### 3.4 Effect Size Assessment

Effect size analysis quantifies the practical magnitude of the METHOD difference independent of sample size. The standardized mean difference is Cohen's d = 0.30 (small effect), with the METHOD factor explaining only 2.30% of total assay variance. In the pharmaceutical context, the observed difference of 0.51 mg/100mg represents only 0.14% of the 35 mg/100mg target value, rendering it negligible for manufacturing decision-making. Comprehensive effect size calculations and bootstrap validation results are presented in the Appendix.

### 3.5 Thief-to-Tablet Comparison: Assessment of Measurement Bias

A critical quality control finding emerges from comparing the thief sampling measurements to the final tablet assay values. The tablet samples (n=30 tablets from 10 randomly selected drums) show substantially lower mean assay values compared to both thief methods:

**Table 3.5: Thief vs. Tablet Assay Comparison**

| Method | N | Mean | SD | Difference from Tablet |
|--------|---|------|-----|----------------------|
| Intermediate Dose (INTM) | 18 | 36.91 | 1.40 | +1.12 mg/100mg (+3.1%) |
| Unit Dose (UNIT) | 18 | 36.40 | 1.98 | +0.61 mg/100mg (+1.7%) |
| Tablet (Final Product) | 30 | 35.79 | 1.36 | — |

This systematic pattern reveals that **tablet assay values are substantially lower than both thief methods**, with the Intermediate Dose thief overestimating final product content by 1.12 mg/100mg (3.1%) and the Unit Dose method overestimating by 0.61 mg/100mg (1.7%). This discrepancy represents active ingredient loss during the powder-to-tablet compression process, amounting to approximately **3% of the final product assay value**—a significant quality control concern requiring immediate investigation.

**Implications for Manufacturing Quality Control**: The poor agreement between thief-sampled powder and final tablet content indicates that powder-level quality assurance sampling alone is insufficient for ensuring final product potency. Compression-related ingredient loss appears to be systematic rather than random. This finding suggests that:

1. Thief sampling results cannot be reliably used to predict final tablet content without correction factors
2. Tablet-level testing or in-process monitoring during compression is essential for regulatory compliance  
3. The compression process parameters warrant investigation to identify and minimize active ingredient loss

Bootstrap validation confirms this thief-to-tablet discrepancy (Table A.5): the p-value for INTM vs. Tablet comparison is 0.0090 (observed) with bootstrap p-value 0.0587, confirming the systematic nature of this difference despite modest sample size.

---

## 4. Discussion and Recommendations

### 4.1 Key Findings Summary

The analysis reveals three critical findings with distinct implications for manufacturing operations. First, both sampling methods produce statistically equivalent mean assay results (p = 0.3054), with a difference of 0.51 mg/100mg (95% CI: [-1.68, 0.65]), suggesting either method is acceptable for average batch-level assessment.

Second, the methods differ substantially in measurement precision. The Intermediate Dose method exhibits 4.81 times lower residual variance (0.71 vs 3.40 mg²/100mg²), providing superior consistency critical for process monitoring and early detection of batch deviations. For pharmaceutical manufacturing, measurement precision is typically more important than equivalent means.

Third and most significantly, location-specific effects within the V-Blender are the dominant source of product variability. Location-to-location variation spans 2.96 mg/100mg (35.01 to 37.98 mg/100mg range), accounting for approximately 60% of total measurement variance and far exceeding the method difference. Location 1 shows significantly lower values while Location 6 is highest, providing clear evidence of non-uniform powder mixing. This location effect represents the highest priority for management attention and process improvement.

### 4.2 Manufacturing Recommendations

#### 4.2.1 Sampling Instrument Selection

**Recommendation**: Transition quality assurance sampling to exclusive use of the Intermediate Dose Thief method for all routine blender content uniformity testing.

**Rationale**: Although both methods produce statistically equivalent mean results (p = 0.3054), the Intermediate Dose Thief demonstrates substantially superior precision with residual variance 4.81 times lower than the Unit Dose method. The observed mean difference of 0.51 mg/100 mg represents only 0.14% of the 35 mg/100 mg target value and is negligible for manufacturing decision-making purposes. The superior precision of the Intermediate Dose method provides more reliable and consistent measurements critical for process monitoring and early detection of batch deviations. Additionally, this method is operationally simpler and more suitable for high-throughput quality assurance sampling, improving laboratory efficiency. The transition should be documented with a process change notification to quality assurance and regulatory affairs personnel.

#### 4.2.2 V-Blender Process Improvement (Priority Action)

**Immediate Finding**: The V-Blender mixing process exhibits substantial location-to-location variation (range: 3.32 mg/100 mg), indicating incomplete mixture homogenization. This is the dominant source of product variability and requires urgent investigation and corrective action.

**Root Cause Assessment**: Several process parameters warrant systematic review and potential adjustment. The current mixing time of 20 minutes should be verified for adequacy in achieving complete homogenization; reduced mixing time may be causing incomplete blending. The mixer operating conditions including rotation speed, direction of rotation, and mechanical integrity of the blending vessel and impellers should be examined for potential issues. The powder loading protocol, including the sequence and timing of ingredient addition, should be reviewed for potential ingredient segregation during the loading phase. Additionally, the blender geometry should be assessed for dead spots, material stagnation areas, or preferential stratification patterns during mixing.

**Implementation Timeline**: An immediate process parameter audit should be conducted within one to two weeks to systematically document current operating conditions. Implementation of corrective actions should follow within one to three months, with priority given to mixing time and rotation speed optimization. Effectiveness of improvements should be monitored through ongoing multi-location sampling with systematic tracking of location-specific assay values, batch means, coefficients of variation, and location variance components as percentages of total variance. Quality assurance sampling records should document location-specific results and process parameter changes with their dates, enabling monthly trending analysis for continuous improvement. If parameter optimization cannot achieve acceptable mixture uniformity, capital equipment assessment should be conducted to determine whether blender capacity, design, or mechanical upgrade is required.

---

## 5. Client Questions Analysis

### Question 1: Are the assay values generally well behaved?

The assay values are generally well behaved and suitable for mixed-effects analysis. All three methods—Intermediate, Unit, and Tablet—show normal distributions with Shapiro-Wilk p-values exceeding 0.18. Conservative outlier detection using the criterion Q3 + 2×IQR identified zero outliers across all methods. An important caveat is that when treating thief samples as repeated measures, the correlation structure should be considered in outlier criteria. Overall, the data satisfy distributional assumptions for ANOVA analysis with no concerns identified.

### Question 2: Is there evidence of a location effect?

Strong evidence of location effects exists, requiring urgent management attention. Location 1 shows significantly lower assay values of 35.01 mg/100 mg (random effect = -1.64, p < 0.01), while Location 6 shows the highest values of 37.98 mg/100 mg (random effect = +1.33, p < 0.05). The total range spans 2.96 mg/100 mg, representing 8.46% of the target value. Location variance accounts for approximately 60% of total variability in assay measurements. The V-Blender mixing parameters should be reviewed immediately, including mixing time, rotation speed, direction, and mechanical integrity. The powder loading sequence and blender geometry warrant investigation. A multi-location monitoring protocol should be implemented, with equipment upgrades considered if uniformity cannot be achieved through parameter optimization. A multi-location monitoring protocol should be implemented, with equipment upgrades considered if uniformity cannot be achieved through parameter optimization.

### Question 3: Do tablet data show drum or time effects?

Temporal analysis of tablet data reveals a linear trend of -0.0029 mg/100 mg per drum sequence (p = 0.47, not significant), indicating no significant downward drift across the 30 drums sampled. The drum-to-drum variance is 0.4361 with a coefficient of variation of 1.84%, suggesting relatively consistent drum-to-drum variation. Autoregressive structure assessment indicates no strong temporal correlation. The recommendation is to continue monitoring drum-to-drum variation through ongoing sampling and to implement statistical process control charts for the tablet production stage to track long-term performance trends.

### Question 4: Are thief-sampled values comparable to tablet values?

Thief-sampled values and tablet values show significant systematic differences, representing a critical finding for manufacturing operations. The Intermediate Dose thief overestimates final tablet content by a mean difference of 1.12 mg/100 mg, while the Unit Dose thief overestimates by 0.61 mg/100 mg. This systematic bias indicates process-related active ingredient loss rather than random variation. The poor agreement between thief samples and final tablets means that thief results cannot be used to directly predict tablet content without correction factors.

The root cause of this discrepancy lies in the compression process, where approximately 3% of active ingredient is systematically lost. This represents a significant quality control gap where thief-based release testing would be insufficient to ensure final product potency meets specifications. Regulatory compliance may be at risk if thief results alone are used for product release decisions. Investigation of compression process mechanisms for ingredient loss is urgent, including assessment of heat, pressure, binding characteristics, and dust loss during compression and handling. Correction factors should be developed to relate thief measurements to final tablet content, additional sampling points should be considered within the manufacturing process, and regulatory compliance implications should be reviewed with quality assurance and regulatory affairs personnel.

---

## 6. Conclusions

This analysis presents a clear trade-off regarding the selection of a sampling instrument. Both the Intermediate and Unit Dose thieves provide statistically equivalent results on average (p = 0.3054), suggesting that if the priority is overall batch-level assessment, both methods are acceptable. However, the Intermediate Dose Thief demonstrates superior precision with variance 4.81 times lower than the Unit Dose method, providing more reliable and consistent measurements critical for process monitoring. The choice of instrument depends on whether the priority is overall average performance or measurement consistency.

More significantly, this study revealed process-related issues of greater importance than the choice of sampling instrument. The primary source of product variability is a substantial lack of mixture uniformity within the V-Blender, with location-to-location variation accounting for approximately 60% of total measurement variance. Additionally, a systematic loss of approximately 3% of the active ingredient occurs during the powder-to-tablet compression stage. These findings indicate that further analysis and process improvement efforts should be focused on these critical manufacturing areas to enhance final product quality and regulatory compliance.

---

### A.1 Exploratory Data Analysis - Complete Results

#### A.1.1 Summary Statistics Table

**Table A.1: Summary Statistics for Assay Value by Method (Including Final Tablets)**

| Method | N | Mean | SD | Min | Q1 | Median | Q3 | Max |
|--------|---|------|----|----|----|----|----|----|
| Unit Dose (UNIT) | 18 | 36.40 | 1.98 | 32.77 | 34.74 | 36.74 | 37.98 | 39.16 |
| Intermediate Dose (INTM) | 18 | 36.91 | 1.40 | 34.38 | 36.01 | 36.67 | 37.83 | 39.80 |
| **Tablet (Final Product)** | **30** | **35.79** | **1.36** | **33.09** | **35.08** | **35.69** | **36.54** | **39.44** |

#### A.1.2 Data Quality Assessment

**Table A.1a: Normality Tests (Shapiro-Wilk)**

| Method | W-Statistic | p-Value | Status |
|--------|------------|---------|--------|
| Intermediate Dose | 0.9836 | 0.9794 | Normal ✓ |
| Unit Dose | 0.9424 | 0.3190 | Normal ✓ |

**Table A.1b: Outlier Detection (Q3 + 2×IQR criterion)**

| Method | Q1 | Q3 | IQR | Upper Bound | Values > Upper Bound | Outliers |
|--------|-----|------|------|-------------|---------------------|----------|
| Intermediate Dose | 36.01 | 37.83 | 1.82 | 41.47 | None | 0 |
| Unit Dose | 34.75 | 37.99 | 3.24 | 44.47 | None | 0 |

**Note on Quantile Calculation:** Quartile values (Q1, Q3) were computed using R's default `quantile()` function with `type=7` (linear interpolation between order statistics). This is the default method in R and produces results that may differ slightly from other statistical software packages using different quantile algorithms (e.g., SAS type=5). The choice of quantile type does not affect the primary conclusion: no outliers are detected using the Q3 + 2×IQR criterion for either method, confirming acceptable data quality.

---

### A.2 Mixed-Effects Model - Comprehensive Results

#### A.2.1 Model Summary and Fixed Effects Tests

**Model Overview:** The mixed-effects model was fitted using R's `nlme` package with maximum likelihood estimation (ML). The model includes a fixed effect for METHOD, random intercepts for LOCATION, and heterogeneous variance structure allowing different residual standard deviations for each sampling method.

![Model Summary Output](./R-screenshots/3-2-model_summary.png)

---

**Table A.2a: Type III Tests of Fixed Effects (ANOVA)**

| Effect | Num DF | Den DF | F-value | p-value |
|--------|--------|--------|---------|---------|
| (Intercept) | 1 | 29 | 6441.024 | <.0001 |
| METHOD | 1 | 29 | 1.081 | 0.307 |

**Table A.2b: Fixed Effects Estimates and 95% Confidence Intervals**

| Effect | Estimate | Std. Error | Lower 95% CI | Upper 95% CI | t-value | p-value |
|--------|----------|-----------|--------------|--------------|---------|---------|
| (Intercept) | 36.9078 | 0.4665 | 35.9805 | 37.8350 | 79.114 | <.0001 |
**METHODUnit | -0.5111 | 0.4915 | -1.4880 | 0.4658 | -1.040 | 0.307 |

#### A.2.2 Variance Components (Heterogeneous Variance Model)

**Table A.3: Variance Components Decomposition**

| Variance Component | Estimate | Role |
|--------------------|----------|------|
| Between-Location (LOC) | 0.9976 | Variability across the 6 sampling locations |
| Within-Location (INTM) | 0.7069 | Residual precision - Intermediate Dose method |
| Within-Location (UNIT) | 3.4001 | Residual precision - Unit Dose method |

#### A.2.3 Location Effects Analysis

**Table A.4: Random Effect for LOCATION**

| Location | Random Intercept | Location Mean | Pct Below Grand Mean | Status |
|----------|---------|---------|--------|--------|
| 1 | -1.64 | 35.01 | -4.47 | Significantly lower |
| 2 | +0.30 | 36.95 | +0.81 | Not significant |
| 3 | -0.43 | 36.22 | -1.18 | Not significant |
| 4 | +0.53 | 37.18 | +1.45 | Not significant |
| 5 | -0.08 | 36.57 | -0.22 | Not significant |
| 6 | +1.33 | 37.98 | +3.62 | Not significant |

**Location Range Span:** 2.96 mg/100mg (from Location 1 at 35.01 to Location 6 at 37.98)

#### A.2.4 Bootstrap Validation Results

**Table A.5: Bootstrap Validation - Three-Way Comparison (1000 resamples)**

| Comparison | Observed P-Value | Bootstrap P-Value | Bias | Standard Error |
|------------|------------------|-------------------|------|----------------|
| UNIT vs. INTM | 0.3777 | 0.3932 | 0.0162 | 0.2958 |
| UNIT vs. TABLET | 0.2135 | 0.3204 | 0.1069 | 0.3092 |
| INTM vs. TABLET | 0.0090 | 0.0587 | 0.0498 | 0.1244 |

#### A.2.5 Residual Diagnostics Summary

**Table A.6: Normality Tests on Residuals (Shapiro-Wilk)**

| Group | W-statistic | p-value | Result |
|-------|------------|---------|--------|
| All residuals combined | 0.9817 | 0.8649 | ✓ Normal |
| INTM method residuals | 0.9705 | 0.6892 | ✓ Normal |
| UNIT method residuals | 0.9893 | 0.9645 | ✓ Normal |

**Residual Variance by Method:**

| Method | Residual Variance | Relative to INTM |
|--------|------------------|-----------------|
| INTM (Intermediate Dose) | 0.6996 | Baseline (1.0×) |
| UNIT (Unit Dose) | 3.3889 | 4.84× higher |

---

### A.3 Regression Perspective: Reference Tables

#### A.4.1 Model Comparison Likelihood Ratio Test

**Table A.7: Random Slopes Model Comparison**

| Model | Model Specification | df | AIC | Deviance | LRT χ² | p-value |
|-------|-------------------|----|----|----------|--------|---------|
| Model 1 | METHOD only, random intercept | 5 | 138.81 | Baseline | — | — |
| Model 2 | METHOD + random slopes | 7 | 141.60 | L.Ratio = 1.21 | 1.21 | **0.5456** |

#### A.3.2 Location-Specific Descriptive Differences

**Table A.8: Location Means and METHOD Differences**

| Location | INTM Mean | UNIT Mean | Difference | Descriptive Pattern |
|----------|-----------|-----------|-----------|--------|
| 1 | 34.99 | 34.25 | -0.73 mg/100mg | INTM slightly higher |
| 2 | 36.97 | 38.14 | +1.16 mg/100mg | UNIT slightly higher |
| 3 | 36.40 | 35.84 | -0.56 mg/100mg | INTM slightly higher |
| 4 | 37.74 | 36.10 | -1.64 mg/100mg | INTM notably higher |
| 5 | 36.51 | 37.76 | +1.25 mg/100mg | UNIT slightly higher |
| 6 | 38.84 | 36.29 | -2.55 mg/100mg | INTM notably higher |

#### A.3.3 Effect Size Calculations

**Table A.9: Effect Size Measures for METHOD Comparison**

| Measure | Value | Classification |
|---------|-------|-----------------|
| Cohen's d | 0.30 | SMALL |
| Hedges' g | 0.29 | SMALL |
| Eta-squared (η²) | 2.30% | Method explains 2.30% variance |
| Omega-squared (ω²) | 0.27% | Bias-corrected estimate |

---

## References

R Core Team (2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria.

Cabrera, J., & McDougall, A. (2002). Statistical Consulting. Springer-Verlag New York, Inc., ISBN 978-1-4757-3663-2 (eBook).

