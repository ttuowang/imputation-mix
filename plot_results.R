
# Load packages
library(showtext)
library(patchwork)
library(glue)

font_add_google("Roboto", "roboto")
font_add_google("Roboto Slab", "roboto slab")
showtext_auto()
col_palette <-  c("#386cb0","#fdb462","#7fc97f")
date <- "11_13"
level <- c("No censoring", "Imputation", glue("No imputation"))


load(file = "sim11-10-ind.RData")

avg.survest.dat.A0 <- apply(survest.dat.A0, 2, mean)
avg.survest.dat.A1 <- apply(survest.dat.A1, 2, mean)
avg.survest.imp.A0 <- apply(survest.imp.A0, 2, mean)
avg.survest.imp.A1 <- apply(survest.imp.A1, 2, mean)
avg.survest.mis.A0 <- apply(survest.mis.A0, 2, mean)
avg.survest.mis.A1 <- apply(survest.mis.A1, 2, mean)

survest.df.A0 <- data.frame(
  time = rep(t.seq, 3),
  avg.survest = c( avg.survest.dat.A0,avg.survest.imp.A0,avg.survest.mis.A0),
  cat = factor(
    rep(level, each=length(t.seq)),
    levels = level)
)

survest.df.A1 <- data.frame(
  time = rep(t.seq, 3),
  avg.survest = c(avg.survest.dat.A1,avg.survest.imp.A1, avg.survest.mis.A1),
  cat = factor(
    rep(level, each=length(t.seq)),
    levels = level)
)

p0ind <- survest.df.A0 %>%
  ggplot( aes(x=time, y=avg.survest) )+
  geom_line(aes(color=cat, linetype = cat), size=1, alpha=0.8) + 
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0, 1) +
  labs(x = "Time (month)", y = "Survival probability",
       title = "Control group (A=0) \n independent censoring", 
       color="Data", linetype="Data") +
  theme_bw() 

p1ind <- survest.df.A1 %>%
  ggplot( aes(x=time, y=avg.survest, color=cat) )+
  geom_line(aes(linetype = cat), size=1, alpha=0.8) + 
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0, 1) +
  labs(x = "Time (month)", y = "Survival probability", 
       title = "Treatment group (A=1) \n independent censoring", 
       color="Data", linetype="Data") +
  theme_bw()

load(file = "sim11-10-d.RData")

avg.survest.dat.A0 <- apply(survest.dat.A0, 2, mean)
avg.survest.dat.A1 <- apply(survest.dat.A1, 2, mean)
avg.survest.imp.A0 <- apply(survest.imp.A0, 2, mean)
avg.survest.imp.A1 <- apply(survest.imp.A1, 2, mean)
avg.survest.mis.A0 <- apply(survest.mis.A0, 2, mean)
avg.survest.mis.A1 <- apply(survest.mis.A1, 2, mean)

survest.df.A0 <- data.frame(
  time = rep(t.seq, 3),
  avg.survest = c( avg.survest.dat.A0,avg.survest.imp.A0,avg.survest.mis.A0),
  cat = factor(
    rep(level, each=length(t.seq)),
    levels = level)
)

survest.df.A1 <- data.frame(
  time = rep(t.seq, 3),
  avg.survest = c(avg.survest.dat.A1,avg.survest.imp.A1, avg.survest.mis.A1),
  cat = factor(
    rep(level, each=length(t.seq)),
    levels = level)
)

p0d <- survest.df.A0 %>%
  ggplot( aes(x=time, y=avg.survest) )+
  geom_line(aes(color=cat, linetype = cat), size=1, alpha=0.8) + 
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0, 1) +
  labs(x = "Time (month)", y = "Survival probability",
       title = "Control group (A=0) \n dependent censoring", 
       color="Data", linetype="Data") +
  theme_bw() 

p1d <- survest.df.A1 %>%
  ggplot( aes(x=time, y=avg.survest, color=cat) )+
  geom_line(aes(linetype = cat), size=1, alpha=0.8) + 
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0, 1) +
  labs(x = "Time (month)", y = "Survival probability", 
       title = "Treatment group (A=1) \n dependent censoring", 
       color="Data", linetype="Data") +
  theme_bw()

patchind <- p0ind+p1ind
patchd <- p0d+p1d

patchind / patchd +
  plot_layout(guides = "collect") &
  theme(
    plot.title = element_text(family = "roboto", size = 15, hjust = 0.5,
                              face = "bold", color="black"),
    axis.title.y = element_text(face="bold", angle=90, size=15, color="grey35"),
    axis.text  = element_text(size=13),
    axis.title.x = element_text(face="bold", size=15, color="grey35"),
    panel.grid.minor = element_blank(),
    legend.title  = element_text(face="bold", size=17),
    legend.text = element_text(size=15),
    legend.position = "bottom",
    legend.justification = "center")

ggsave(
  filename = here::here("results", date ,glue("survival_n{n}.pdf")),
  width = 10,
  height = 10,
  device = "pdf"
)


#------------------------------------------------------------------------------#

load(file = "sim11-10-ind.RData")

avg.mcfest.imp.A0 <- apply(mcfest.imp.A0, 2, mean)
avg.mcfest.imp.A1 <- apply(mcfest.imp.A1, 2, mean)
avg.mcfest.dat.A0 <- apply(mcfest.dat.A0, 2, mean)
avg.mcfest.dat.A1 <- apply(mcfest.dat.A1, 2, mean)
avg.mcfest.mis.A0 <- apply(mcfest.mis.A0, 2, mean)
avg.mcfest.mis.A1 <- apply(mcfest.mis.A1, 2, mean)

mcfest.df.A0 <- data.frame(
  time = rep(t.seq, 3),
  avg.mcf = c(avg.mcfest.dat.A0, avg.mcfest.imp.A0, avg.mcfest.mis.A0),
  cat = factor(
    rep(level, each=length(t.seq)),
    levels = level)
)

mcfest.df.A1 <- data.frame(
  time = rep(t.seq, 3),
  avg.mcf = c(avg.mcfest.dat.A1, avg.mcfest.imp.A1, avg.mcfest.mis.A1),
  cat = factor(
    rep(level, each=length(t.seq)),
    levels = level)
)

p2ind <- mcfest.df.A0 %>%
  ggplot( aes(x=time, y=avg.mcf, color=cat) )+
  geom_line(aes(linetype = cat), size=1.0, alpha=0.8) + 
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0, 2.5) +
  labs(x = "Time (month)", y = "Mean cumulative function",
       title = "Control group (A=0) \n independent censoring", 
       color="Data", linetype="Data") +
  theme_bw() 

p3ind <- mcfest.df.A1 %>%
  ggplot( aes(x=time, y=avg.mcf, color=cat) )+
  geom_line(aes(linetype = cat), size=1.0, alpha=0.8) + 
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4, 6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0, 1.8) +
  labs(x = "Time (month)", y = "Mean cumulative function", 
       title = "Treatment group (A=1) \n independent censoring", 
       color="Data", linetype="Data") +
  theme_bw()

load(file = "sim11-10-d.RData")

avg.mcfest.imp.A0 <- apply(mcfest.imp.A0, 2, mean)
avg.mcfest.imp.A1 <- apply(mcfest.imp.A1, 2, mean)
avg.mcfest.dat.A0 <- apply(mcfest.dat.A0, 2, mean)
avg.mcfest.dat.A1 <- apply(mcfest.dat.A1, 2, mean)
avg.mcfest.mis.A0 <- apply(mcfest.mis.A0, 2, mean)
avg.mcfest.mis.A1 <- apply(mcfest.mis.A1, 2, mean)

mcfest.df.A0 <- data.frame(
  time = rep(t.seq, 3),
  avg.mcf = c(avg.mcfest.dat.A0, avg.mcfest.imp.A0, avg.mcfest.mis.A0),
  cat = factor(
    rep(level, each=length(t.seq)),
    levels = level)
)

mcfest.df.A1 <- data.frame(
  time = rep(t.seq, 3),
  avg.mcf = c(avg.mcfest.dat.A1, avg.mcfest.imp.A1, avg.mcfest.mis.A1),
  cat = factor(
    rep(level, each=length(t.seq)),
    levels = level)
)

p2d <- mcfest.df.A0 %>%
  ggplot( aes(x=time, y=avg.mcf, color=cat) )+
  geom_line(aes(linetype = cat), size=1.0, alpha=0.8) + 
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0, 2.5) +
  labs(x = "Time (month)", y = "Mean cumulative function",
       title = "Control group (A=0) \n dependent censoring", 
       color="Data", linetype="Data") +
  theme_bw() 

p3d <- mcfest.df.A1 %>%
  ggplot( aes(x=time, y=avg.mcf, color=cat) )+
  geom_line(aes(linetype = cat), size=1.0, alpha=0.8) + 
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4, 6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0, 1.8) +
  labs(x = "Time (month)", y = "Mean cumulative function", 
       title = "Treatment group (A=1) \n dependent censoring", 
       color="Data", linetype="Data") +
  theme_bw()

patchind <- p2ind+p3ind
patchd <- p2d+p3d

patchind / patchd +
  plot_layout(guides = "collect") &
  theme(
    plot.title = element_text(family = "roboto", size = 15, hjust = 0.5,
                              face = "bold", color="black"),
    axis.title.y = element_text(face="bold", angle=90, size=15, color="grey35"),
    axis.text  = element_text(size=13),
    axis.title.x = element_text(face="bold", size=15, color="grey35"),
    panel.grid.minor = element_blank(),
    legend.title  = element_text(face="bold", size=17),
    legend.text = element_text(size=15),
    legend.position = "bottom",
    legend.justification = "center")

ggsave(
  filename = here::here("results", date,glue("mcf_n{n}.pdf")),
  width = 10,
  height = 10,
  device = "pdf"
)

#------------------------------------------------------------------------------#


load(file = "sim11-10-ind.RData")

Y1.df.plot <- Y1.df %>%
  mutate(time = as.numeric(time)) %>%
  mutate(A = case_when(
    A == 0 ~ "Control",
    A == 1 ~ "Treatment"
  )) %>%
  group_by(type, A, time) %>%
  summarise(
    Y1 = mean(Y1),
    .groups = 'drop'
  ) %>%
  mutate(type = factor(type, levels = level))

p4ind <- Y1.df.plot %>%
  filter(A == "Control") %>%
  ggplot(aes(x=time, y=Y1, color=type )) +
  geom_point(aes(shape=type),size=2.5,alpha=0.7) +
  geom_line(aes(linetype=type),size=1.0, alpha=0.7) +
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(-0.8, 0.1) +
  theme_bw() +
  labs(x = "Time (month)", y = expression(Y[3]), 
       title = "Control group (A=0) \n independent censoring", 
       color="Data", linetype="Data",shape="Data") 

p5ind <- Y1.df.plot %>%
  filter(A == "Treatment") %>%
  ggplot(aes(x=time, y=Y1, color=type )) +
  geom_point(aes(shape=type),size=2.5,alpha=0.7) +
  geom_line(aes(linetype=type),size=1.0, alpha=0.7) +
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(-0.8, 0.1) +
  theme_bw() +
  labs(x = "Time (month)", y = expression(Y[3]), 
       title = "Treatment group (A=1) \n independent censoring", 
       color="Data", linetype="Data",shape="Data") 

load(file = "sim11-10-d.RData")

Y1.df.plot <- Y1.df %>%
  mutate(time = as.numeric(time)) %>%
  mutate(A = case_when(
    A == 0 ~ "Control",
    A == 1 ~ "Treatment"
  )) %>%
  group_by(type, A, time) %>%
  summarise(
    Y1 = mean(Y1),
    .groups = 'drop'
  ) %>%
  mutate(type = factor(type, levels = level))

p4d <- Y1.df.plot %>%
  filter(A == "Control") %>%
  ggplot(aes(x=time, y=Y1, color=type )) +
  geom_point(aes(shape=type),size=2.5,alpha=0.7) +
  geom_line(aes(linetype=type),size=1.0, alpha=0.7) +
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(-1.0, 0.1) +
  theme_bw() +
  labs(x = "Time (month)", y = expression(Y[3]), 
       title = "Control group (A=0) \n dependent censoring", 
       color="Data", linetype="Data",shape="Data") 

p5d <- Y1.df.plot %>%
  filter(A == "Treatment") %>%
  ggplot(aes(x=time, y=Y1, color=type )) +
  geom_point(aes(shape=type),size=2.5,alpha=0.7) +
  geom_line(aes(linetype=type),size=1.0, alpha=0.7) +
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(-1.0, 0.1) +
  theme_bw() +
  labs(x = "Time (month)", y = expression(Y[3]), 
       title = "Treatment group (A=1) \n dependent censoring", 
       color="Data", linetype="Data",shape="Data") 

patchind <- p4ind+p5ind
patchd <- p4d+p5d

patchind / patchd +
  plot_layout(guides = "collect") &
  theme(
    plot.title = element_text(family = "roboto", size = 15, hjust = 0.5,
                              face = "bold", color="black"),
    axis.title.y = element_text(face="bold", angle=90, size=15, color="grey35"),
    axis.text  = element_text(size=13),
    axis.title.x = element_text(face="bold", size=15, color="grey35"),
    panel.grid.minor = element_blank(),
    legend.title  = element_text(face="bold", size=17),
    legend.text = element_text(size=15),
    legend.position = "bottom",
    legend.justification = "center")

ggsave(
  filename = here::here("results", date,glue("Y_n{n}.pdf")),
  width = 10,
  height = 10,
  device = "pdf"
)


#------------------------------------------------------------------------------#


load(file = "sim11-10-ind.RData")

Y2.df.plot <- Y2.df %>%
  mutate(time = as.numeric(time)) %>%
  mutate(A = case_when(
    A == 0 ~ "Control",
    A == 1 ~ "Treatment"
  )) %>%
  group_by(type, A, time) %>%
  summarise(
    Y2 = mean(Y2),
    .groups = 'drop'
  ) %>%
  mutate(type = factor(type, levels = level))

p6ind <- Y2.df.plot %>%
  filter(A == "Control") %>%
  ggplot(aes(x=time, y=Y2, color=type )) +
  geom_point(aes(shape=type),size=2.5,alpha=0.7) +
  geom_line(aes(linetype=type),size=1.0, alpha=0.7) +
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0.2, 0.6) +
  theme_bw() +
  labs(x = "Time (month)", y = expression(Y[4]), 
       title = "Control group (A=0) \n independent censoring", 
       color="Data", linetype="Data",shape="Data") 

p7ind <- Y2.df.plot %>%
  filter(A == "Treatment") %>%
  ggplot(aes(x=time, y=Y2, color=type )) +
  geom_point(aes(shape=type),size=2.5,alpha=0.7) +
  geom_line(aes(linetype=type),size=1.0, alpha=0.7) +
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0.2, 0.6) +
  theme_bw() +
  labs(x = "Time (month)", y = expression(Y[4]), 
       title = "Treatment group (A=1) \n independent censoring", 
       color="Data", linetype="Data",shape="Data") 

load(file = "sim11-10-d.RData")

Y2.df.plot <- Y2.df %>%
  mutate(time = as.numeric(time)) %>%
  mutate(A = case_when(
    A == 0 ~ "Control",
    A == 1 ~ "Treatment"
  )) %>%
  group_by(type, A, time) %>%
  summarise(
    Y2 = mean(Y2),
    .groups = 'drop'
  ) %>%
  mutate(type = factor(type, levels = level))

p6d <- Y2.df.plot %>%
  filter(A == "Control") %>%
  ggplot(aes(x=time, y=Y2, color=type )) +
  geom_point(aes(shape=type),size=2.5,alpha=0.7) +
  geom_line(aes(linetype=type),size=1.0, alpha=0.7) +
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0.2, 0.6) +
  theme_bw() +
  labs(x = "Time (month)", y = expression(Y[4]), 
       title = "Control group (A=0) \n dependent censoring", 
       color="Data", linetype="Data",shape="Data") 

p7d <- Y2.df.plot %>%
  filter(A == "Treatment") %>%
  ggplot(aes(x=time, y=Y2, color=type )) +
  geom_point(aes(shape=type),size=2.5,alpha=0.7) +
  geom_line(aes(linetype=type),size=1.0, alpha=0.7) +
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,4,6)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  ylim(0.2, 0.6) +
  theme_bw() +
  labs(x = "Time (month)", y = expression(Y[4]), 
       title = "Treatment group (A=1) \n dependent censoring", 
       color="Data", linetype="Data",shape="Data") 

patchind <- p6ind+p7ind
patchd <- p6d+p7d

patchind / patchd +
  plot_layout(guides = "collect") &
  theme(
    plot.title = element_text(family = "roboto", size = 15, hjust = 0.5,
                              face = "bold", color="black"),
    axis.title.y = element_text(face="bold", angle=90, size=15, color="grey35"),
    axis.text  = element_text(size=13),
    axis.title.x = element_text(face="bold", size=15, color="grey35"),
    panel.grid.minor = element_blank(),
    legend.title  = element_text(face="bold", size=17),
    legend.text = element_text(size=15),
    legend.position = "bottom",
    legend.justification = "center")

ggsave(
  filename = here::here("results", date,glue("B_n{n}.pdf")),
  width = 10,
  height = 10,
  device = "pdf"
)

