# ==============================================================================
# BAYESIAN DISTANCE SAMPLING - KEPADATAN ORANG UTAN (DATA SARANG)
# Versi Khusus Eksekusi Server GitHub (Headless)
# Prepared by : M Rifqi, Tasya, Edy S, Siti Nur Badriyah, Alfons 
# ==============================================================================

library(readxl)
library(dplyr)
library(rjags)
library(coda)
library(ggplot2)

# STEP 1: IMPORT DATA
df_nest = read_excel("data_nest.xlsx", sheet = "ou nest")
df_clean = df_nest %>% filter(!is.na(PPD))
x_obs = df_clean$PPD / 1000 
n_obs = length(x_obs) 

df_transect = df_clean %>% 
  group_by(`ID Transek`) %>% 
  summarise(L_km = first(`Panjang Jalur_Km`)) 
L = sum(df_transect$L_km) 
W = 0.025                 
Area_Survei = 2 * W * L   
Area_Studi = 15.0         

# STEP 2: KALIBRASI PRIOR
p_nest = 0.85 
r_nest = 1.0  
t_nest = 602  
multiplier = p_nest * r_nest * t_nest 

prior_ou_lower = 0.12
prior_ou_upper = 0.25
prior_nest_lower = prior_ou_lower * multiplier
prior_nest_upper = prior_ou_upper * multiplier
prior_nest_mean  = mean(c(prior_nest_lower, prior_nest_upper))

sd_nest = (prior_nest_upper - prior_nest_lower) / 4
var_nest = sd_nest^2
shape_prior = (prior_nest_mean^2) / var_nest
rate_prior  = prior_nest_mean / var_nest

M = 300 
y = c(rep(1, n_obs), rep(0, M - n_obs))
x = c(x_obs, rep(NA, M - n_obs))

jags_data = list(
  y = y, x = x, M = M, W = W, Area_Survei = Area_Survei, Area_Studi = Area_Studi, 
  shape_prior = shape_prior, rate_prior = rate_prior, multiplier = multiplier, n_obs = n_obs
)

# STEP 3: MODEL JAGS
cat("
model {
  D_nest ~ dgamma(shape_prior, rate_prior)
  D_ou <- D_nest / multiplier
  lambda <- D_nest * Area_Survei
  psi <- lambda / M
  
  sigma ~ dunif(0.001, W)
  for (i in 1:M) {
    z[i] ~ dbern(psi)
    x[i] ~ dunif(0, W)
    p_det[i] <- exp(- (x[i]^2) / (2 * sigma^2))
    mu[i] <- z[i] * p_det[i]
    y[i] ~ dbern(mu[i])
  }
  
  N_area_survei <- sum(z[1:M])
  prob_deteksi <- n_obs / max(N_area_survei, 1)
  N_Total_Individu <- D_ou * Area_Studi 
}
", file = "model_sekung.txt")

# STEP 4: MCMC
inits = function() { list(z = rep(1, M), sigma = runif(1, 0.005, 0.02)) }
params = c("D_ou", "prob_deteksi", "N_Total_Individu", "sigma")

jm = jags.model("model_sekung.txt", data = jags_data, inits = inits, n.chains = 3, n.adapt = 2000)
update(jm, n.iter = 5000)
post_samples = coda.samples(jm, variable.names = params, n.iter = 10000, thin = 5)

# STEP 5 & 6: HASIL & DIAGNOSTIK
rhat_diag = gelman.diag(post_samples, multivariate = FALSE)
df_posterior = as.data.frame(as.matrix(post_samples))

mean_2016 = 0.185  
upper_2016 = 0.25  
prob_naik_vs_mean = sum(df_posterior$D_ou > mean_2016) / nrow(df_posterior)
prob_naik_vs_upper = sum(df_posterior$D_ou > upper_2016) / nrow(df_posterior)

# STEP 7: VISUALISASI GGPLOT2
df_chains = do.call(rbind, lapply(1:length(post_samples), function(chain_idx) {
  data.frame(
    Chain = as.factor(chain_idx),
    Iteration = 1:nrow(post_samples[[chain_idx]]),
    D_ou = as.numeric(post_samples[[chain_idx]][, "D_ou"]),
    prob_deteksi = as.numeric(post_samples[[chain_idx]][, "prob_deteksi"])
  )
}))

plot_trace = ggplot(df_chains, aes(x = Iteration, y = D_ou, color = Chain)) + geom_line(alpha = 0.6, linewidth = 0.5) + scale_color_manual(values = c("#FF5A5F", "#087E8B", "#F5A623")) + labs(title = "MCMC Traceplot: Kepadatan Orang Utan", x = "Iterasi", y = "Kepadatan (Ind/km2)") + theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold"))

plot_density = ggplot(df_posterior, aes(x = D_ou)) + geom_density(fill = "#4C8BF5", alpha = 0.6, color = "#1F59B6", linewidth = 1) + geom_vline(xintercept = c(0.12, 0.25), color = "red", linetype = "dashed", linewidth = 1) + labs(title = "Estimasi Kepadatan Orang Utan", subtitle = "Garis merah: Rentang Prior 2016 (0.12 - 0.25)", x = "Kepadatan (Individu / km2)", y = "Kepadatan Probabilitas") + theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

mean_N = mean(df_posterior$N_Total_Individu)
ci_lower = quantile(df_posterior$N_Total_Individu, 0.025)
ci_upper = quantile(df_posterior$N_Total_Individu, 0.975)
plot_abundance = ggplot(df_posterior, aes(x = N_Total_Individu)) + geom_histogram(fill = "#F28C28", alpha = 0.7, color = "white", bins = 40) + geom_vline(xintercept = mean_N, color = "#B22222", linewidth = 1.2) + annotate("text", x = mean_N, y = Inf, label = paste("Mean:", round(mean_N, 1)), vjust = 2, hjust = -0.1, fontface = "bold", color = "#B22222") + labs(title = paste("Estimasi Total Populasi (Area", Area_Studi, "km²)"), subtitle = paste("95% CI:", round(ci_lower, 1), "-", round(ci_upper, 1), "individu"), x = "Estimasi Total Individu", y = "Frekuensi Sampel Posterior") + theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

mean_sigma = mean(df_posterior$sigma)
x_seq = seq(0, W, length.out = 200)
p_det_seq = exp(- (x_seq^2) / (2 * mean_sigma^2))
df_deteksi = data.frame(Jarak_Meter = x_seq * 1000, Probabilitas = p_det_seq)
plot_detection_curve = ggplot(df_deteksi, aes(x = Jarak_Meter, y = Probabilitas)) + geom_line(color = "#2E8B57", linewidth = 1.5) + geom_area(fill = "#2E8B57", alpha = 0.2) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + labs(title = "Fungsi Probabilitas Deteksi Sarang (Half-Normal)", subtitle = paste("Rata-rata probabilitas deteksi:", round(mean(df_posterior$prob_deteksi), 2)), x = "Jarak dari Garis Transek (Meter)", y = "Probabilitas Deteksi (p)") + theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, color = "gray30"))

# STEP 8: EKSPOR FILE
sink("Ringkasan_Hasil_Bayesian.txt")
cat("=== HASIL POSTERIOR MCMC ===\n")
print(summary(post_samples))
cat("\n=== DIAGNOSTIK GELMAN-RUBIN ===\n")
print(rhat_diag)
cat(sprintf("\nProbabilitas kepadatan 2026 > Rata-rata 2016 (0.185) : %.2f%%\n", prob_naik_vs_mean * 100))
cat(sprintf("Probabilitas kepadatan 2026 > Batas Atas 2016 (0.25)  : %.2f%%\n", prob_naik_vs_upper * 100))
sink()

ggsave("1_Traceplot_Kepadatan.png", plot_trace, width = 8, height = 5, dpi = 300)
ggsave("2_Kurva_Posterior.png", plot_density, width = 8, height = 5, dpi = 300)
ggsave("3_Histogram_Total_Individu.png", plot_abundance, width = 8, height = 5, dpi = 300)
ggsave("4_Kurva_Deteksi_HalfNormal.png", plot_detection_curve, width = 8, height = 5, dpi = 300)
