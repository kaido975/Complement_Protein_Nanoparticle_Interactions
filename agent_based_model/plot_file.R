library(dplyr)
library(ggplot2)
library(ks)

avg_step = 20
nsteps = 1400
runs = 500
n_agents <- 200
n_agents_st <- as.character(n_agents) #200 binding sites/NP 

mdf <- read.csv(paste(n_agents_st, "agents_data.csv" ,sep=""))
X <- slice(mdf[2], seq(nsteps+1, (nsteps+1)*runs , nsteps+1))

x <- n_agents - X
x <- sapply(x, as.numeric)

#Kernel Density Estimation
h <- density(x, kernel = "gaussian")$bw 
w <- 1 / (pnorm(0, mean = x, sd = h, lower.tail = FALSE) - pnorm(n_agents, mean = x, sd = h, lower.tail = FALSE))
d <- density(x, bw = h, kernel = "gaussian", weights = w / length(x), adjust = 0.5)
d$y[d$x < 0 | d$x > n_agents] <- 0

area <- sum(d$y * diff(d$x)[1])
d$y <- d$y / area
df <- data.frame(x)
p <- ggplot(df, aes(x = x)) +
  geom_histogram(aes(y =..density..),
                 breaks = seq(0, n_agents, by = 14.0),
                 fill = "#619CFF"  , alpha=.6) +
  xlab("Filled Sites") +
  ylab("Density")+
  theme_light()+
  theme(text = element_text(size = 15))

p <- p + geom_line(data = data.frame(x = d$x, y = d$y), aes(x = x, y = y), color = "#F8766D") + xlim(0, n_agents) + ggtitle(paste(n_agents_st, "Binding Sites / NP") )
print(p)
ggsave(
  filename = paste(n_agents_st, ".png" ,sep=""),
  device = "png", width = 5.5, height = 4.0)
