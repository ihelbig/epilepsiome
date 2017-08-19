# epilepsiome
```
#load and plot SCN2A studies
dxy <- read.table("SCN2A_studies.txt",header=T)
library(ggplot2)
names(dxy) <- c("x","y","ylo","yhi")
ggplot(dxy, aes(x=x, y=y, ymin=ylo, ymax=yhi))+
  geom_pointrange()+
  geom_hline(yintercept = 0.01162791, linetype=2)+
  coord_flip()+
  xlab('Study')+
  ylab('Frequency with 95% CI')+
  scale_y_continuous()
```
