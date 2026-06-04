setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()



library(ggplot2)
library(dplyr)

# ---------------------------------------------------------
# 1. Load empirical data
# ---------------------------------------------------------


data <- read.table("ddm_data.txt",
                   header = TRUE,
                   sep = " ")
colnames(data) <- c("RT", "response","condition")

# ---------------------------------------------------------
# 2. Load predicted CDFs from fast-dm ML 8 parameters
# ---------------------------------------------------------

cdf_IF <- read.table("cdf_IF.lst",
                     header = FALSE)

cdf_I <- read.table("cdf_I.lst",
                    header = FALSE)

colnames(cdf_IF) <- c("time", "cdf")
colnames(cdf_I)  <- c("time", "cdf")

cdf_IF <- cdf_IF %>% filter(time >= 0)
cdf_I  <- cdf_I  %>% filter(time >= 0)

# ---------------------------------------------------------
# 3. Create empirical CDF function
# ---------------------------------------------------------

# Condition I+F
rt_IF <- data %>%
  filter(condition == "I+F") %>%
  pull(RT)

ecdf_IF <- ecdf(rt_IF)

emp_IF <- data.frame(
  time = sort(rt_IF),
  cdf  = ecdf_IF(sort(rt_IF))
)

# Condition I
rt_I <- data %>%
  filter(condition == "I") %>%
  pull(RT)

ecdf_I <- ecdf(rt_I)

emp_I <- data.frame(
  time = sort(rt_I),
  cdf  = ecdf_I(sort(rt_I))
)

# ---------------------------------------------------------
# 4. Plot empirical vs predicted CDF Condition I + F
# ---------------------------------------------------------

p1 <- ggplot() +
  
  # empirical
  geom_step(data = emp_IF,
            aes(x = time, y = cdf),
            linewidth = 1) +
  
  # predicted
  geom_line(data = cdf_IF,
            aes(x = time, y = cdf),
            linewidth = 1,
            linetype = "dashed") +
  
  labs(title = "Condition I+F",
       x = "Reaction Time (s)",
       y = "Cumulative Probability") +
  
  theme_classic()

print(p1)

# ---------------------------------------------------------
# 5. Plot for condition I
# ---------------------------------------------------------

p2 <- ggplot() +
  
  geom_step(data = emp_I,
            aes(x = time, y = cdf),
            linewidth = 1) +
  
  geom_line(data = cdf_I,
            aes(x = time, y = cdf),
            linewidth = 1,
            linetype = "dashed") +
  
  labs(title = "Condition I",
       x = "Reaction Time (s)",
       y = "Cumulative Probability") +
  
  theme_classic()

print(p2)


#-------------------------------------
## ML using full model 11 parameters
#-------------------------------------

cdf_IF_full <- read.table("cdf_IF_full.lst",
                      header = FALSE)

cdf_I_full <- read.table("cdf_I_full.lst",
                     header = FALSE)

colnames(cdf_IF_full) <- c("time", "cdf")
colnames(cdf_I_full)  <- c("time", "cdf")




# ---------------------------------------------------------
# 4. Plot empirical vs predicted CDF Condition I+F
# ---------------------------------------------------------

p1 <- ggplot() +
  
  # empirical
  geom_step(data = emp_IF,
            aes(x = time, y = cdf),
            linewidth = 1) +
  
  # predicted
  geom_line(data = cdf_IF_full,
            aes(x = time, y = cdf),
            linewidth = 1,
            linetype = "dashed") +
  
  labs(title = "Condition I+F",
       x = "Reaction Time (s)",
       y = "Cumulative Probability") +
  
  theme_classic()+

  coord_cartesian(xlim = c(0, 2.5))


print(p1)

# ---------------------------------------------------------
# 5. Plot for condition I
# ---------------------------------------------------------

p2 <- ggplot() +
  
  geom_step(data = emp_I,
            aes(x = time, y = cdf),
            linewidth = 1) +
  
  geom_line(data = cdf_I_full,
            aes(x = time, y = cdf),
            linewidth = 1,
            linetype = "dashed") +
  
  labs(title = "Condition I",
       x = "Reaction Time (s)",
       y = "Cumulative Probability") +
  
  theme_classic()+
  coord_cartesian(xlim = c(0, 2.5))

print(p2)


