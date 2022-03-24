library(tidyverse)
library(visreg)
library(ggplot2)
library(lme4)
library(lmerTest)
library(nlme)
library(PNWColors)
theme_set(theme_classic())
library(patchwork)

# Calculate NH4 for pre samples --------
# load data

#data with which bottles belong to which experimental levels
experimental_water <- read_csv("experimental_water.csv")

# three bottle per 3 treatment levels, two replicates and one matrix
# We added 200 uL of ammonium stock to the matrix sample and use the comparison between the WP sample and the WP matrix to calculate the matrix

# load standard curve data
standard_water_samples <- read_csv("experimental_water_fluorometry_standard.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

#do the calculations for the standard curve
standard_water_samples_f <- standard_water_samples %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading

#graph standard curve
# aes(x,y). dataframe, aes(x,y) then some more aesthetic stuff
ggplot(standard_water_samples_f, aes(mean_FLU, nh4_conc_final_umol_L)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model
# lm is used to fit linear models lm(y~x, data source)
sc_mod_water_samples <- lm(nh4_conc_final_umol_L ~ mean_FLU, data = standard_water_samples_f)
summary(sc_mod_water_samples)

#save coefficients
# coef is a  function which extracts model coefficients from objects returned by modeling functions
int_water_samples<- coef(sc_mod_water_samples)[1]
slope_water_samples <- coef(sc_mod_water_samples)[2]

# load fluorometry data
# mean_FLU is the name of the new variable that was calculated from the , the following is the modifying action youre doing. cbind means youre combining the columns. rowMeans means were adding a cloumn of mean values of flu (from the same row?) to the water_samples_feb2-22_fluorometry.csv
glow_water_samples <- read_csv("experimental_water_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

#calculate matrix as according to Taylor et al 2007
# % Matrix Effects = [(standard_spike - standard_zero) - (sample_spike - sample_zero)]/ (standard_spike - standard_zero) * 100

water_samples <- read_csv("experimental_water.csv")

#this is combining the fluorometry data with the sample site data. It also renamed the mean FLU as Fsm_zero
water_samples2 <- water_samples %>%
  left_join(glow_water_samples, by = c("bottle")) %>%
  mutate(Fsm_zero = mean_FLU)

#this new mutated data set just has site data and fsm spike (mean flu)
water_samples_spike <- water_samples2 %>%
  filter(sample == "matrix") %>%
  transmute(concentration = concentration,
            Fsm_spike = mean_FLU)

water_samples_zero <- water_samples2 %>%
  filter(sample != "matrix")

Fst_zero_water_samples <- standard_water_samples_f$mean_FLU[standard_water_samples_f$nh4_vol_uL == "0"]

Fst_spike_water_samples <- standard_water_samples_f$mean_FLU[standard_water_samples_f$nh4_vol_uL == "200"]

water_samples_matrix <- water_samples_zero %>% 
  left_join(water_samples_spike, by = "concentration") %>%
  mutate(
    Fst_zero = Fst_zero_water_samples,
    Fst_spike = Fst_spike_water_samples,
    matrix = 100 * ((Fst_spike_water_samples - Fst_zero_water_samples - (Fsm_spike - Fsm_zero))/(Fst_spike_water_samples - Fst_zero_water_samples)), matrix_est = 20) %>%
  select(concentration, Fsm_spike, Fst_zero, Fst_spike, matrix, matrix_est)

water_samples3 <- water_samples2 %>%
  left_join(water_samples_matrix, by = "concentration") %>%
  mutate(
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100))
  ) %>%
  mutate(int = int_water_samples, #include values for the int and slope in for every column
         slope = slope_water_samples) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc = int + slope * Fsm_cor_Taylor) %>%
  filter(sample != "matrix") %>%
  select(-sample_matrix)


