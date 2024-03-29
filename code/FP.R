library(polymorphology2)
library(scales)
library(minpack.lm) 
FP<-fread("data/FP/MSH6_FP_results_2023Aug_from_Aki/FP_results.csv", header = T)[,-20]
# Melt the data to long format
FP_long <- melt(FP, id.vars = c("Stat", "Trt", "Peptide"), variable.name = "Conc", value.name = "Value")

real_sci_format <- function(x) {
  format <- function(val) {
    if (is.na(val)) {
      return(NA)
    }
    if (val == 0) {
      return("0")
    }
    e <- floor(log10(abs(val)))
    m <- val/10^e
    if (e == 0) {
      return(paste0(formatC(m, format = "f", digits = 1)))
    }
    return(paste0(formatC(m, format = "f", digits = 1), " x 10^", e))
  }
  sapply(x, format)
}


# Split the data into mean and sd tables
FP_mean <- FP_long[Stat == "mean", .(Trt, Peptide, Conc, mean = Value)]
FP_sd <- FP_long[Stat == "sd", .(Trt, Peptide, Conc, sd = Value)]

# Merge the two tables
FP_final <- merge(FP_mean, FP_sd, by = c("Trt", "Peptide", "Conc"))

# Order the columns
setcolorder(FP_final, c("Trt", "Peptide", "Conc", "mean", "sd"))

# Create the plot

# Convert the 'Conc' column to character and then to numeric
# Convert the 'Conc' column from microM (uM) to M
FP_final[, Conc := as.numeric(as.character(Conc)) * 1e-9]
FP_final0<-FP_final
FP_final<-FP_final[Conc!=0]


# Define the sigmoid function
sigmoid_offset <- function(x, a, b, c, d) {
  d + (a / (1 + exp(-b * (x - c))))
}

# Function to fit and predict for each combination of Trt and Peptide
predict_for_combination <- function(trt, peptide) {
  sub <- FP_final[Trt == trt & Peptide == peptide]
  x <- sub$Conc
  y <- sub$mean
  
  # Initial parameter estimates
  start_list <- list(
    a = max(y) - min(y),
    b = 1,
    c = median(x),
    d = min(y)
  )
  
  # Fit the sigmoid curve to the data
  fit <- try(nlsLM(y ~ sigmoid_offset(x, a, b, c, d), start = start_list), silent = TRUE)
  
  # If the fit is successful, predict y values
  if (!inherits(fit, "try-error")) {
    new_x <- exp(seq(log(min(x)), log(max(x)), length.out = 100))
    predicted_y <- predict(fit, newdata = data.frame(x = new_x))
    return(data.table(Trt = trt, Peptide = peptide, Conc = new_x, mean = predicted_y))
  } else {
    return(NULL)
  }
}

# Loop across all combinations of Trt and Peptide
all_combinations <- CJ(unique(FP_final$Trt), unique(FP_final$Peptide))
predicted_list <- lapply(1:nrow(all_combinations), function(i) {
  predict_for_combination(all_combinations$V1[i], all_combinations$V2[i])
})

# Combine all predicted data.tables
predicted_dt <- rbindlist(predicted_list, fill = TRUE)

# Custom breaks function for log scale with intermediate ticks
log_breaks <- function(n = 10) {
  function(x) {
    main_breaks <- 10^seq(floor(log10(min(x))), ceiling(log10(max(x))), by = 1)
    sub_breaks <- unlist(lapply(main_breaks, function(b) seq(b, b*10, length.out = n+1)[-1]))
    sort(c(main_breaks, sub_breaks))
  }
}

FP_final[, Trt := factor(Trt, levels = c("WT", "3Amut"))]
predicted_dt[, Trt := factor(Trt, levels = c("WT", "3Amut"))]
# Create the plot

lines(RT, predict(fit, data = RT), col = cpl[2],xlab='',ylab='')

p <- ggplot(FP_final, aes(x = Conc, y = mean, group = Peptide, color = Peptide, shape=Peptide)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0) +
  labs(x = "Concentration (M)", y = "Flourescence asiotropy") +
  scale_x_log10(labels = real_sci_format) +  # Use the custom function here
  facet_wrap(~Trt, scales = "free_x") +
  theme_bw(base_size = 6) +
  geom_line(data=data.table(mean=predict(fit, data = RT), Conc=RT, Trt="WT", Peptide="H3K4me1"))+
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank()
  )

p

pdf("figures/FP.pdf", width = 3.5, height = 2)
print(p)
p+geom_line(data=predicted_dt)
dev.off()


max<-FP_final[Trt=="WT",.(FP=max(mean)), by=Peptide]
arat_msh6<-means_in_silico[file=="arat_msh6"]
plot(max$FP, arat_msh6$dist6)
