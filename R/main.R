#--------
#调试（新SIF模型）----
setwd('E:/code/rTRIPLEX_CW_Flux-fixed')
source('E:/code/rTRIPLEX_CW_Flux-fixed/R/TRIPLEX_SIF_Flux.R')
Inputvariable <- read.csv('data/Inputvariable_sif_DBS_fix5_7days.csv')
Inputpara <- read.csv('data/Inputpara_sif_DBS.csv')
output_SIF <- TRIPLEX_SIF_Flux(Input_variable = Inputvariable,Input_parameter = Inputpara,overyear = FALSE)
write.csv(output_SIF,'data/output_sif_DBS_fix5_7days.csv')
print('done')



#绘图（新CW模型-GPP）----
## 计算线性模型，R方，RMSE，IA
setwd('E:/code/rTRIPLEX_CW_Flux-fixed')
output_SIF <- read.csv('data/output_sif_DBS_fix5-1.csv')
library(ggplot2)
library(ggpubr)
regre <- lm(GPP30min ~ GPP, output_SIF)
r <- formatC(summary(regre)$r.squared, format = "f", digits = 2)
rmse <- formatC(sqrt(sum(residuals(regre)^2) / (nrow(output_SIF) - 2)), format = "f", digits = 2)
ia <- formatC(1 - sum((output_SIF$GPP30min - output_SIF$GPP)^2) /
                sum((abs(output_SIF$GPP30min - mean(output_SIF$GPP)) +
                       abs(output_SIF$GPP - mean(output_SIF$GPP)))^2), format = "f", digits = 2)

# 获取线性回归方程的系数，并保留两位小数
intercept <- formatC(coef(regre)[1], format = "f", digits = 2)
slope <- formatC(coef(regre)[2], format = "f", digits = 2)

## 创建ggplot基础图
p <- ggplot(output_SIF, aes(x= GPP, y = GPP30min)) +
  geom_point(shape = 21, size = 3, color = 'blue', stroke = 1) + # 散点图，边缘蓝色，填充白色
  geom_smooth(method = 'lm', color = 'red', se = FALSE, size = 1) + # 线性拟合线
  geom_abline(intercept = 0, slope = 1, color = 'grey60', size = 1, linetype = 'dashed') + # 1:1对线线
  labs(x = expression(Observed~GPP~(g~C~m^-2~30~min^-1)),
       y = expression(Simulated~GPP~(g~C~m^-2~30~min^-1))) +  # 坐标轴标签

  theme_minimal(base_size = 15) +  # 设置主题
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 15),
        panel.grid = element_line(),
        panel.border = element_rect(fill = NA, color = 'black',size = 0.5))

## 设置坐标轴范围和刻度
min_val <- min(output_SIF$GPP, output_SIF$GPP30min) - 0.1
max_val <- max(output_SIF$GPP, output_SIF$GPP30min) * 1.1
p <- p +
  scale_x_continuous(limits = c(0, max_val), breaks = seq(0, max_val, by = 0.2)) +  # 限制坐标轴范围并设置刻度
  scale_y_continuous(limits = c(0, max_val), breaks = seq(0, max_val, by = 0.2))

## 添加r方，RMSE，和IA的注释
p <- p +
  annotate("text", x = 0.05,
           y = max_val * 0.95,
           label = paste("y =", slope, "x +", intercept),size = 5, hjust = 0) +
  annotate("text", x = 0.05,
           y = max_val * 0.85,
           label = bquote(R^2 == .(r)), size = 5, hjust = 0) +
  annotate("text", x = 0.05,
           y = max_val * 0.75,
           label = paste("RMSE =", rmse), size = 5, hjust = 0) +
  annotate("text", x = 0.05,
           y = max_val * 0.65,
           label = paste("IA =", ia), size = 5, hjust = 0) +
  annotate("text", x = 0,
           y = max_val * 0.95,
           label = "(b)", size = 5, hjust = 0)

print(p)
ggsave("./result/DBS_SIF_fix5-1.jpg", width = 8, height = 6)

#--------
# 分季度绘图（GPP）----
setwd('E:/code/rTRIPLEX_CW_Flux-yuan')
output_SIF_season <- read.csv('data/output_DBS_fix5-1.csv')
library(ggplot2)
library(ggpubr)
library(dplyr)

# 添加季节列
output_SIF_season <- output_SIF_season %>%
  mutate(Season = case_when(
    Month %in% 3:5 ~ "Spring",
    Month %in% 6:8 ~ "Summer",
    Month %in% 9:11 ~ "Autumn",
    Month == 12 ~ "Winter",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Season)) %>%
  mutate(Season = factor(Season, levels = c("Spring", "Summer", "Autumn", "Winter")))

# 创建ggplot基础图
p <- ggplot(output_SIF_season, aes(x = GPP, y = GPP30min)) +
  geom_point(shape = 21, size = 3, color = 'blue', stroke = 1) + # 散点图，边缘蓝色，填充白色
  geom_smooth(method = 'lm', color = 'red', se = FALSE, size = 1) + # 线性拟合线
  geom_abline(intercept = 0, slope = 1, color = 'grey60', size = 1, linetype = 'dashed') + # 1:1对线线
  labs(x = expression(Observed~GPP~(g~C~m^-2~30~min^-1)),
       y = expression(Simulated~GPP~(g~C~m^-2~30~min^-1))) +  # 坐标轴标签
  theme_minimal(base_size = 15) +  # 设置主题
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 15),
        panel.grid = element_line(),
        panel.border = element_rect(fill = NA, color = 'black', size = 0.5)) +
  facet_wrap(~ Season, nrow = 2, ncol = 2)  # 按季节分面

# 设置坐标轴范围和刻度
min_val <- min(output_SIF_season$GPP, output_SIF_season$GPP30min) - 0.1
max_val <- max(output_SIF_season$GPP, output_SIF_season$GPP30min) * 1.1
p <- p +
  scale_x_continuous(limits = c(0, max_val), breaks = seq(0, max_val, by = 0.2)) +  # 限制坐标轴范围并设置刻度
  scale_y_continuous(limits = c(0, max_val), breaks = seq(0, max_val, by = 0.2))


# 计算斜率和截距，并保留两位小数
season_stats <- output_SIF_season %>%
  group_by(Season) %>%
  summarise(
    slope = round(coef(lm(GPP30min ~ GPP))[2], 2),
    intercept = round(coef(lm(GPP30min ~ GPP))[1], 2),
    r_squared = round(summary(lm(GPP30min ~ GPP))$r.squared, 2),
    rmse = round(sqrt(mean(residuals(lm(GPP30min ~ GPP))^2)), 2),
    ia = round(1 - sum((GPP30min - GPP)^2) / sum((abs(GPP30min - mean(GPP)) + abs(GPP - mean(GPP)))^2), 2)
  )


# 添加线性公式、r方，RMSE，和IA的注释
season_stats <- season_stats %>%
  mutate(label = paste0("y = ", slope, "x + ", intercept,
                        "\nR² = ", r_squared,
                        "\nRMSE = ", rmse,
                        "\nIA = ", ia))

# 在ggplot中添加文本标签
p <- p +
  geom_text(data = season_stats, aes(x = 0.05, y = max_val * 0.95, label = label),
            hjust = 0, vjust = 1, size = 5, color = "black")


print(p)
ggsave("./result/DBS_Seasonal.jpg", plot = p, width = 12, height = 10)
#--------
# 日均GPP绘图----
setwd('E:/code/rTRIPLEX_CW_Flux-fixed')
output_SIF_daily <- read.csv('data/output_sif_DBS_fix5-1-daily.csv')
library(ggplot2)
library(ggpubr)
regre <- lm(GPPdaily ~ GPPob, output_SIF_daily)
r <- formatC(summary(regre)$r.squared, format = "f", digits = 2)
rmse <- formatC(sqrt(sum(residuals(regre)^2) / (nrow(output_SIF_daily) - 2)), format = "f", digits = 2)
ia <- formatC(1 - sum((output_SIF_daily$GPPdaily - output_SIF_daily$GPPob)^2) /
                sum((abs(output_SIF_daily$GPPdaily - mean(output_SIF_daily$GPPob)) +
                       abs(output_SIF_daily$GPPob - mean(output_SIF_daily$GPPob)))^2), format = "f", digits = 2)

## 获取线性回归方程的系数，并保留两位小数
intercept <- formatC(coef(regre)[1], format = "f", digits = 2)
slope <- formatC(coef(regre)[2], format = "f", digits = 2)

## 创建ggplot基础图
p <- ggplot(output_SIF_daily, aes(x= GPPob, y = GPPdaily)) +
  geom_point(shape = 21, size = 4, color = 'blue', stroke = 1) + # 散点图，边缘蓝色，填充白色
  geom_smooth(method = 'lm', color = 'red', se = FALSE, size = 1) + # 线性拟合线
  geom_abline(intercept = 0, slope = 1, color = 'grey60', size = 1, linetype = 'dashed') + # 1:1对线线
  labs(x = expression(Observed~GPP~(g~C~m^-2~day^-1)),
       y = expression(Simulated~GPP~(g~C~m^-2~day^-1))) +  # 坐标轴标签

  theme_minimal(base_size = 15) +  # 设置主题
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 15),
        panel.grid = element_line(size = 0.5),
        panel.border = element_rect(fill = NA, color = 'black',size = 0.5))

## 设置坐标轴范围和刻度
min_val <- min(output_SIF_daily$GPPob, output_SIF_daily$GPPdaily) - 0.1
max_val <- max(output_SIF_daily$GPPob, output_SIF_daily$GPPdaily) * 1
p <- p +
  scale_x_continuous(limits = c(0, max_val), breaks = seq(0, max_val, by = 2)) +  # 限制坐标轴范围并设置刻度
  scale_y_continuous(limits = c(0, max_val), breaks = seq(0, max_val, by = 2))

## 添加r方，RMSE，和IA的注释
p <- p +
  annotate("text", x = 0.05,
           y = max_val * 0.95,
           label = paste("y =", slope, "x +", intercept),size = 5, hjust = 0) +
  annotate("text", x = 0.05,
           y = max_val * 0.85,
           label = bquote(R^2 == .(r)), size = 5, hjust = 0) +
  annotate("text", x = 0.05,
           y = max_val * 0.75,
           label = paste("RMSE =", rmse), size = 5, hjust = 0) +
  annotate("text", x = 0.05,
           y = max_val * 0.65,
           label = paste("IA =", ia), size = 5, hjust = 0)

print(p)
ggsave("./result/DBS_SIF_fix5-1-daily.jpg", width = 8, height = 6)

# 日均GPP逐日绘图----
setwd('E:/code/rTRIPLEX_CW_Flux-fixed')
output_SIF_daily <- read.csv('data/output_sif_DBS_fix5-1-daily.csv')
library(ggplot2)
library(dplyr)

# 计算R方和RMSE
regre <- lm(GPPdaily ~ GPPob, output_SIF_daily)
r_squared <- formatC(summary(regre)$r.squared, format = "f", digits = 2)
rmse <- formatC(sqrt(mean(residuals(regre)^2)), format = "f", digits = 2)

# 创建唯一的月份和对应的id数据框
unique_months <- output_SIF_daily %>%
  group_by(月份) %>%
  summarise(id = first(id))

# 绘图
p <- ggplot(output_SIF_daily, aes(x = id)) +
  geom_line(aes(y = GPPob, color = "GPPob"), size = 1) +
  geom_line(aes(y = GPPdaily, color = "GPPdaily"), size = 1) +
  scale_color_manual(values = c("GPPob" = "blue", "GPPdaily" = "red")) +
  scale_x_continuous(breaks = unique_months$id, labels = unique_months$月份) +
  scale_y_continuous(breaks = seq(0, 10.5, by = 2)) +
  labs(x = "Month", y = expression(GPP~(g~C~m^-2~day^-1)), color = "图例") +
  theme_minimal(base_size = 15) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        panel.grid = element_line(size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black',size = 0.5),
        legend.position = c(0.85, 0.9),
        legend.title = element_blank()) +
  annotate("text", x = max(output_SIF_daily$id) * 0.82, y = max(output_SIF_daily$GPPdaily) * 0.82,
           label = paste("R² =", r_squared, "\nRMSE =", rmse), size = 5, hjust = 0)

print(p)
ggsave("./result/GPP_daily_plot.jpg", plot = p, width = 10, height = 6)
#--------
# 环境变量逐月绘图----
setwd('E:/code/rTRIPLEX_CW_Flux-fixed')
output_SIF_variable <- read.csv('data/output_sif_DBS_fix5-1月汇总.csv')
library(ggplot2)
library(ggpubr)

# 确保 Month 按月份大小升序排序
output_SIF_variable$Month <- factor(output_SIF_variable$Month,
                                    levels = unique(output_SIF_variable$Month))

# 绘制 Ta 的折线图（实心圆点）
plot_ta <- ggplot(output_SIF_variable, aes(x = Month, y = Ta, group = 1)) +
  geom_point(aes(shape = "Ta"), size = 3, color = "black", na.rm = TRUE) +  # 改为黑色
  geom_line(aes(linetype = "Ta"), color = "black", na.rm = TRUE) +  # 改为黑色
  scale_shape_manual(values = c("Ta" = 16)) +
  scale_linetype_manual(values = c("Ta" = "solid")) +
  scale_y_continuous(name = "Ta (℃)", limits = c(0, max(output_SIF_variable$Ta) * 1.1)) +
  scale_x_discrete(labels = function(x) gsub("月", "", x)) +  # 去掉月份中的“月”
  labs(x = "Month", shape = "Legend", linetype = "Legend") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),  # 添加边框
        legend.position = c(0.1, 0.95),  # 图例放在左上角
        legend.text = element_text(size = 15),  # 图例文本大小
        legend.title = element_blank(),  # 隐藏图例标题
        axis.text = element_text(size = 15),  # 调整坐标轴文本字号
        axis.title = element_text(size = 15),  # 调整轴标题字号
        axis.title.y.right = element_text(angle = 90))  # 右 y 轴标题从下到上排列

# 计算 RH 和 VPD 的最大值，并计算转换因子 coef_rh_vpd
max_rh <- max(output_SIF_variable$RH, na.rm = TRUE)
max_vpd <- max(output_SIF_variable$VPD, na.rm = TRUE)
coef_rh_vpd <- max_rh / max_vpd  # 计算 RH 和 VPD 的转换因子

# 绘制 RH 和 VPD 的折线图（空心圆点和实心三角），实现双 y 轴
plot_rh_vpd <- ggplot(output_SIF_variable, aes(x = Month)) +
  geom_point(aes(y = RH, shape = "RH"), size = 3, color = "black", na.rm = TRUE) +  # 改为黑色
  geom_line(aes(y = RH, linetype = "RH", group = 1), color = "black", na.rm = TRUE) +  # 改为黑色
  geom_point(aes(y = VPD * coef_rh_vpd, shape = "VPD"), size = 3, color = "black", na.rm = TRUE) +  # 改为黑色
  geom_line(aes(y = VPD * coef_rh_vpd, linetype = "VPD", group = 1), color = "black", na.rm = TRUE) +  # 改为黑色
  scale_shape_manual(values = c("RH" = 1, "VPD" = 17)) +
  scale_linetype_manual(values = c("RH" = "solid", "VPD" = "dashed")) +
  scale_y_continuous(
    name = "RH (%)",
    limits = c(0, max_rh * 1.1),
    sec.axis = sec_axis(~ . / coef_rh_vpd, name = "VPD (hPa)")  # 反向转换 VPD 的值
  ) +
  scale_x_discrete(labels = function(x) gsub("月", "", x)) +  # 去掉月份中的“月”
  labs(x = "Month", shape = "Legend", linetype = "Legend") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),  # 添加边框
        legend.position = c(0.1, 0.95),  # 图例放在左上角
        legend.title = element_blank(),  # 隐藏图例标题
        legend.text = element_text(size = 15),  # 图例文本大小
        axis.text = element_text(size = 15),  # 调整坐标轴文本字号
        axis.title = element_text(size = 15),  # 调整轴标题字号
        axis.title.y.right = element_text(angle = 90))  # 右 y 轴标题从下到上排列

# 计算 Rn 和 SIF 的最大值，并计算转换因子 coef
max_rn <- max(output_SIF_variable$Rn, na.rm = TRUE)
max_sif <- max(output_SIF_variable$SIF, na.rm = TRUE)
coef <- max_rn / max_sif

# 绘制 Rn 和 SIF 的折线图（空心矩形和*号），SIF 乘以转换因子 coef
plot_rn_sif <- ggplot(output_SIF_variable, aes(x = Month)) +
  geom_point(aes(y = Rn, shape = "Rn"), size = 3, color = "black", na.rm = TRUE) +  # 改为黑色
  geom_line(aes(y = Rn, linetype = "Rn", group = 1), color = "black", na.rm = TRUE) +  # 改为黑色
  geom_point(aes(y = SIF * coef, shape = "SIF"), size = 3, color = "black", na.rm = TRUE) +  # 改为黑色
  geom_line(aes(y = SIF * coef, linetype = "SIF", group = 1), color = "black", na.rm = TRUE) +  # 改为黑色
  scale_shape_manual(values = c("Rn" = 0, "SIF" = 8)) +
  scale_linetype_manual(values = c("Rn" = "solid", "SIF" = "dotted")) +
  scale_y_continuous(
    name = expression(Rn ~ (W~m^-2)),
    limits = c(0, max_rn * 1.1),
    sec.axis = sec_axis(~ . / coef, name = expression(SIF ~ (mW~m^-2~nm^-1~sr^-1)))
  ) +
  scale_x_discrete(labels = function(x) gsub("月", "", x)) +  # 去掉月份中的“月”
  labs(x = "Month", shape = "Legend", linetype = "Legend") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),  # 添加边框
        legend.position = c(0.1, 0.95),  # 图例放在左上角
        legend.title = element_blank(),  # 隐藏图例标题
        legend.text = element_text(size = 15),  # 图例文本大小
        axis.text = element_text(size = 15),  # 调整坐标轴文本字号
        axis.title = element_text(size = 15),  # 调整轴标题字号
        axis.title.y.right = element_text(angle = 90))  # 右 y 轴标题从下到上排列

# 保存图表
ggsave("./result/DBS_SIF_fix5-1月汇总_Ta.jpg", plot = plot_ta, width = 8, height = 6)
ggsave("./result/DBS_SIF_fix5-1月汇总_RH_VPD.jpg", plot = plot_rh_vpd, width = 8, height = 6)
ggsave("./result/DBS_SIF_fix5-1月汇总_Rn_SIF.jpg", plot = plot_rn_sif, width = 8, height = 6)

#--------
# 日内平均变化图----
setwd('E:/code/rTRIPLEX_CW_Flux-fixed')
output_SIF_variable <- read.csv('data/output_sif_DBS_fix5-1.csv')
library(ggplot2)
library(dplyr)

# 去除 time 为 5.5, 6, 6.5, 19, 19.5 的数据
filtered_data <- output_SIF_variable %>%
  filter(!time %in% c(5.5, 19.5))

# 按 time 分类汇总，计算 GPP30min、GPP 和 SIFfull 的平均值及误差条范围
daily_stats <- filtered_data %>%
  group_by(time) %>%
  summarise(
    GPP30min_avg = mean(GPP30min, na.rm = TRUE),
    GPP30min_sd = sd(GPP30min, na.rm = TRUE),  # 计算 GPP30min 的标准差
    GPP_avg = mean(GPP, na.rm = TRUE),
    GPP_sd = sd(GPP, na.rm = TRUE),  # 计算 GPP 的标准差
    SIFfull_avg = mean(SIFfull, na.rm = TRUE),
    SIFfull_sd = sd(SIFfull, na.rm = TRUE)  # 计算 SIFfull 的标准差
  ) %>%
  mutate(
    GPP30min_min = GPP30min_avg - GPP30min_sd,  # 平均值减去标准差
    GPP30min_max = GPP30min_avg + GPP30min_sd,  # 平均值加上标准差
    GPP_min = GPP_avg - GPP_sd,  # 平均值减去标准差
    GPP_max = GPP_avg + GPP_sd,  # 平均值加上标准差
    SIFfull_min = SIFfull_avg - SIFfull_sd,  # 平均值减去标准差
    SIFfull_max = SIFfull_avg + SIFfull_sd   # 平均值加上标准差
  )

# 绘制 GPP30min 的日内变化图（带误差条）
plot_gpp30min <- ggplot(daily_stats, aes(x = time, y = GPP30min_avg)) +
  geom_line(aes(color = "Simulated GPP"), size = 1) +  # 平均值折线图
  geom_point(aes(color = "Simulated GPP", shape = "Simulated GPP"), size = 2) +  # 平均值点
  geom_errorbar(aes(ymin = GPP30min_min, ymax = GPP30min_max, color = "Simulated GPP"), width = 0.2) +  # 误差条
  scale_x_continuous(
    breaks = seq(6, max(daily_stats$time, na.rm = TRUE), by = 2),  # x刻度从6开始，间隔2
    name = "Hour of day, 2023"
  ) +
   scale_color_manual(values = c("Simulated GPP" = "blue")) +
  scale_shape_manual(values = c("Simulated GPP" = 16)) +
  labs(
    x = "Hour of day, 2023",
    y = expression(Siumlated~GPP~(g~C~m^-2~day^-1)),
    color = "Legend",
    shape = "Legend"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0, 1),  # 图例放置在左上角
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    panel.border = element_rect(color = "black", fill = NA)  # 添加边框
  ) +
  annotate("text", x = max(daily_stats$time, na.rm = TRUE), y = max(daily_stats$GPP30min_max, na.rm = TRUE) * 1.1,
           label = "(a)", color = "blue", size = 5, hjust = 1)

# 绘制 GPP 的日内变化图（带误差条）
plot_gpp <- ggplot(daily_stats, aes(x = time, y = GPP_avg)) +
  geom_line(aes(color = "Observed GPP"), size = 1) +  # 平均值折线图
  geom_point(aes(color = "Observed GPP", shape = "Observed GPP"), size = 2) +  # 平均值点
  geom_errorbar(aes(ymin = GPP_min, ymax = GPP_max, color = "Observed GPP"), width = 0.2) +  # 误差条
  scale_x_continuous(
    breaks = seq(6, max(daily_stats$time, na.rm = TRUE), by = 2),  # x刻度从6开始，间隔2
    name = "Hour of day, 2023"
  ) +
  scale_color_manual(values = c("Observed GPP" = "red")) +
  scale_shape_manual(values = c("Observed GPP" = 16)) +
  labs(
    x = "Hour of day, 2023",
    y = expression(Observed~GPP~(g~C~m^-2~day^-1)),
    color = "Legend",
    shape = "Legend"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0, 1),  # 图例放置在左上角
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    panel.border = element_rect(color = "black", fill = NA)  # 添加边框
  ) +
  annotate("text", x = max(daily_stats$time, na.rm = TRUE), y = max(daily_stats$GPP_max, na.rm = TRUE) * 1.1,
           label = "(b)", color = "red", size = 5, hjust = 1)

# 绘制 SIFfull 的日内变化图（带误差条）
plot_siffull <- ggplot(daily_stats, aes(x = time, y = SIFfull_avg)) +
  geom_line(aes(color = "SIFfull"), size = 1) +  # 平均值折线图
  geom_point(aes(color = "SIFfull", shape = "SIFfull"), size = 2) +  # 平均值点
  geom_errorbar(aes(ymin = SIFfull_min, ymax = SIFfull_max, color = "SIFfull"), width = 0.2) +  # 误差条
  scale_x_continuous(
    breaks = seq(6, max(daily_stats$time, na.rm = TRUE), by = 2),  # x刻度从6开始，间隔2
    name = "Hour of day, 2023"
  ) +
  scale_color_manual(values = c("SIFfull" = "skyblue")) +
  scale_shape_manual(values = c("SIFfull" = 16)) +
  labs(
    x = "Hour of day, 2023",
    y = expression(SIFfull~(u~mol~m^-2~s^-1)),
    color = "Legend",
    shape = "Legend"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0, 1),  # 图例放置在左上角
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    panel.border = element_rect(color = "black", fill = NA)  # 添加边框
  ) +
  annotate("text", x = max(daily_stats$time, na.rm = TRUE), y = max(daily_stats$SIFfull_max, na.rm = TRUE) * 1.1,
           label = "(c)", color = "skyblue", size = 5, hjust = 1)

# 保存图表
print(plot_gpp)
ggsave("./result/DBS_SIF_fix5-1日内汇总_Siumlated_GPP_Errorbar.jpg", plot = plot_gpp30min, width = 8, height = 6)
ggsave("./result/DBS_SIF_fix5-1日内汇总_Observed_GPP_Errorbar.jpg", plot = plot_gpp, width = 8, height = 6)
ggsave("./result/DBS_SIF_fix5-1日内汇总_SIFfull_Errorbar.jpg", plot = plot_siffull, width = 8, height = 6)

# 打印图表以便检查
print(plot_siffull)
