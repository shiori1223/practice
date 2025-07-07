ABB411<-read.csv("ABB_Apr_Nov.csv")
ABBhs<-read.csv("ABB_h_vs_a.csv")

library(ggplot2)
library(dplyr)
library(scales)
library(ggpubr)

#●4月11月の折れ線グラフ
# sampling_date2 を日付型に変換
ABB411$sampling_date2 <- as.Date(ABB411$sampling_date2)

# ラベル表示したい日付と対応する表示形式
breaks <- as.Date(c("2023-11-22", "2024-04-10", "2024-11-20", "2025-04-10"))
labels <- c("2023-11", "2024-4", "2024-11", "2025-4")

# グラフ作成
SLC1_411 <- ggplot(ABB411, aes(x = sampling_date2, 
                               y = SLC12A5_1_methylation_level_ave,
                               group = bear_ID, 
                               color = bear_ID)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_date(breaks = breaks, labels = labels) +
  labs(
    title = "SLC12A5-1",
    x = "Sampling date",
    y = "Methylation level (%)",
    color = "bear ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# 表示
SLC1_411

# グラフ作成
SLC2_411 <- ggplot(ABB411, aes(x = sampling_date2, 
                               y = SLC12A5_2_methylation_level_ave,
                               group = bear_ID, 
                               color = bear_ID)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_date(breaks = breaks, labels = labels) +
  labs(
    title = "SLC12A5-2",
    x = "Sampling date",
    y = "Methylation level (%)",
    color = "bear ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# 表示
SLC2_411


# グラフ作成
SLC3_411 <- ggplot(ABB411, aes(x = sampling_date2, 
                               y = SLC12A5_3_methylation_level_ave,
                               group = bear_ID, 
                               color = bear_ID)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_date(breaks = breaks, labels = labels) +
  labs(
    title = "SLC12A5-3",
    x = "Sampling date",
    y = "Methylation level (%)",
    color = "bear ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# 表示
SLC3_411

# グラフ作成
SLC4_411 <- ggplot(ABB411, aes(x = sampling_date2, 
                               y = SLC12A5_4_methylation_level_ave,
                               group = bear_ID, 
                               color = bear_ID)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_date(breaks = breaks, labels = labels) +
  labs(
    title = "SLC12A5-4",
    x = "Sampling date",
    y = "Methylation level (%)",
    color = "bear ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# 表示
SLC4_411

ggarrange(SLC1_411, SLC2_411, SLC3_411, SLC4_411, ncol=2, nrow=2, common.legend = TRUE, legend="right")

#●冬眠期間・活動期間中の変化の図
library(tidyr)

hs_1 <- ABBhs %>%
  select(bear_ID, change_h_1, change_a_1) %>%
  gather(key = "period", value = "methylation_change", -bear_ID) %>%
  mutate(period = recode(period, 
                         "change_h_1" = "hibernation period", 
                         "change_a_1" = "active period")) %>%
  ggplot(aes(x = bear_ID, y = methylation_change, color = period)) +
  geom_point() +
  labs(x = "bear ID", y = " methylation change (%)", title = "SLC12A5-1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme_minimal() 

# プロットを表示
hs_1

hs_2 <- ABBhs %>%
  select(bear_ID, change_h_2, change_a_2) %>%
  gather(key = "period", value = "methylation_change", -bear_ID) %>%
  mutate(period = recode(period, 
                         "change_h_2" = "hibernation period", 
                         "change_a_2" = "active period")) %>%
  ggplot(aes(x = bear_ID, y = methylation_change, color = period)) +
  geom_point() +
  labs(x = "bear ID", y = " methylation change (%)", title = "SLC12A5-2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme_minimal() 

hs_3 <- ABBhs %>%
  select(bear_ID, change_h_3, change_a_3) %>%
  gather(key = "period", value = "methylation_change", -bear_ID) %>%
  mutate(period = recode(period, 
                         "change_h_3" = "hibernation period", 
                         "change_a_3" = "active period")) %>%
  ggplot(aes(x = bear_ID, y = methylation_change, color = period)) +
  geom_point() +
  labs(x = "bear ID", y = " methylation change (%)", title = "SLC12A5-3") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme_minimal() 

hs_4 <- ABBhs %>%
  select(bear_ID, change_h_4, change_a_4) %>%
  gather(key = "period", value = "methylation_change", -bear_ID) %>%
  mutate(period = recode(period, 
                         "change_h_4" = "hibernation period", 
                         "change_a_4" = "active period")) %>%
  ggplot(aes(x = bear_ID, y = methylation_change, color = period)) +
  geom_point() +
  labs(x = "bear ID", y = " methylation change (%)", title = "SLC12A5-4") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme_minimal() 

ggarrange(hs_1, hs_2, hs_3, hs_4, ncol=2, nrow=2, common.legend = TRUE, legend="right")

#●冬眠期間・活動期間中の1日あたりの変化の図
library(tidyr)

hsd_1 <- ABBhs %>%
  select(bear_ID, daychange_h_1, daychange_a_1) %>%
  gather(key = "period", value = "methylation_change", -bear_ID) %>%
  mutate(period = recode(period, 
                         "daychange_h_1" = "hibernation period", 
                         "daychange_a_1" = "active period")) %>%
  ggplot(aes(x = bear_ID, y = methylation_change, color = period)) +
  geom_point() +
  labs(x = "bear ID", y = " methylation change per day (%)", title = "SLC12A5-1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme_minimal() 

# プロットを表示
hsd_1

hsd_2 <- ABBhs %>%
  select(bear_ID, daychange_h_2, daychange_a_2) %>%
  gather(key = "period", value = "methylation_change", -bear_ID) %>%
  mutate(period = recode(period, 
                         "daychange_h_2" = "hibernation period", 
                         "daychange_a_2" = "active period")) %>%
  ggplot(aes(x = bear_ID, y = methylation_change, color = period)) +
  geom_point() +
  labs(x = "bear ID", y = " methylation change per day (%)", title = "SLC12A5-2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme_minimal() 

hsd_3 <- ABBhs %>%
  select(bear_ID, daychange_h_3, daychange_a_3) %>%
  gather(key = "period", value = "methylation_change", -bear_ID) %>%
  mutate(period = recode(period, 
                         "daychange_h_3" = "hibernation period", 
                         "daychange_a_3" = "active period")) %>%
  ggplot(aes(x = bear_ID, y = methylation_change, color = period)) +
  geom_point() +
  labs(x = "bear ID", y = " methylation change per day (%)", title = "SLC12A5-3") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme_minimal() 

hsd_4 <- ABBhs %>%
  select(bear_ID, daychange_h_4, daychange_a_4) %>%
  gather(key = "period", value = "methylation_change", -bear_ID) %>%
  mutate(period = recode(period, 
                         "daychange_h_4" = "hibernation period", 
                         "daychange_a_4" = "active period")) %>%
  ggplot(aes(x = bear_ID, y = methylation_change, color = period)) +
  geom_point() +
  labs(x = "bear ID", y = " methylation change per day (%)", title = "SLC12A5-4") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme_minimal() 

ggarrange(hsd_1, hsd_2, hsd_3, hsd_4, ncol=2, nrow=2, common.legend = TRUE, legend="right")

#●4月11月の推定年齢折れ線グラフ

ABB411age<-read.csv("predicted_age_result_ABB_Apr_Nov.csv")

# sampling_date2 を日付型に変換
ABB411age$sampling_date2 <- as.Date(ABB411age$sampling_date2)

# ラベル表示したい日付と対応する表示形式
breaks <- as.Date(c("2023-11-22", "2024-04-10", "2024-11-20", "2025-04-10"))
labels <- c("2023-11", "2024-4", "2024-11", "2025-4")

# グラフ作成
ABBage_411_SRM <- ggplot(ABB411age, aes(x = sampling_date2, 
                                        y = predicted_age_SRM,
                                        group = bear_ID, 
                                        color = bear_ID)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_date(breaks = breaks, labels = labels) +
  labs(
    title = "single regression model",
    
    x = "sampling date",
    y = "predicted age (year)",
    color = "bear ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# 表示
ABBage_411_SRM

ABBage_411_PCRM <- ggplot(ABB411age, aes(x = sampling_date2, 
                                         y = predicted_age_PCRM,
                                         group = bear_ID, 
                                         color = bear_ID)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_date(breaks = breaks, labels = labels) +
  labs(
    title = "principal component regression model",
    
    x = "sampling date",
    y = "predicted age (year)",
    color = "bear ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ABBage_411_ENRM <- ggplot(ABB411age, aes(x = sampling_date2, 
                                         y = predicted_age_ENM,
                                         group = bear_ID, 
                                         color = bear_ID)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_date(breaks = breaks, labels = labels) +
  labs(
    title = "elastic net regression model",
    
    x = "sampling date",
    y = "predicted age (year)",
    color = "bear ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ABBage_411_SVRM <- ggplot(ABB411age, aes(x = sampling_date2, 
                                         y = predicted_age_SVRM,
                                         group = bear_ID, 
                                         color = bear_ID)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_date(breaks = breaks, labels = labels) +
  labs(
    title = "support vector regression model",
    
    x = "sampling date",
    y = "predicted age (year)",
    color = "bear ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggarrange(ABBage_411_SRM, ABBage_411_PCRM, ABBage_411_ENRM, ABBage_411_SVRM, ncol=2, nrow=2, common.legend = TRUE, legend="right")

