library(data.table)
library(ggplot2)
library(ggthemes)
library(trend)
library(lme4)
library(lmtest)

db <- data.table(read.csv("./Data/Dati resistenza S.aureus.csv",dec=","))
db[,X:=NULL]
# Estraggo percentuali e casi e organizzo i dati in formato long
percent <- db[1:6]
setnames(percent, c("Antibiotico", "2015", "2016", "2017"))
percent <- percent[2:6]
percent <-
  melt(
    percent,
    id.vars = "Antibiotico",
    variable.name = "Anno",
    value.name = "Resistenza"
  )
absolute <- db[9:14]
setnames(absolute, c("Antibiotico", "2015", "2016", "2017"))
absolute <- absolute[2:6]
absolute <-
  melt(
    absolute,
    id.vars = "Antibiotico",
    variable.name = "Anno",
    value.name = "Resistenza"
  )
total <- db[8:9, 2:4]
setnames(total, c("2015", "2016", "2017"))
total <- total[1]
total <-
  melt(total, variable.name = "Anno", value.name = "Totale")
# Qualcosa non mi torna, la somma per anno di absolute non mi da i valori di total e
#le percentuali del 2015 superano il 100%, mentre nel 2016 e nel 2017
# Alcuni degli isolati erano multi-resistenti?
# Come sono calcolate le percentuali? Per Oxacillina vengono riportati 13, 6 e 5 casi,
#in percentuale 41.7%, 19.9% e 15.1%. I totali segnalati sono 31 25 e 19. 13/31=0.419, e ok
#ma 6/25 è 0.24 e 5/19 è 0.263, diversi dalle percentuali riportate.
# Perchè nel 2016 ho una resistenza non classificata (25 casi totali, ma nelle resistenze sono 24)?
#Manca un antibiotico?
# Andamenti separati
ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico), data = percent) +
  geom_line(show.legend = F) +
  geom_point(show.legend = F) +
  scale_y_continuous(name = "", labels = scales::percent) +
  facet_grid(Antibiotico ~ ., scales = "free_y") +
  labs(title = "Resistenze %", x = "") +
  scale_color_tableau() +
  theme_igray()
# Confronto su unico grafico
ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico), data = percent) +
  geom_point() +
  geom_line() +
  scale_y_continuous(name = "", labels = scales::percent) +
  labs(title = "Resistenze %", y = "", x = "") +
  scale_color_tableau() +
  theme_igray() +
  theme(legend.key = element_blank())
# Modelli di regressione
#l'idea è quella di calcolare la pendenza della retta interpolante e testare che sia != 0
summary(lm(Resistenza ~ as.numeric(Anno) + as.factor(Antibiotico), data =
             percent))
# Mi dice che Clinda, Levo e Oxa partono vicino al punto stimato come intercetta (poco mi importa)
# L'anno risulta significativo, quindi c'è un trend generale in diminuzione, a prescindere dal
#tipo di antibiotico
summary(lmer(Resistenza ~ as.numeric(Anno) + (1 |
                                                Antibiotico), data = percent))
# Trend generale confermato anche dal modello a effetti misti
# Modelli e correlazioni su singolo antibiotico
models <- vector("list", 0)
correlations <- vector("list", 0)
for (i in unique(percent$Antibiotico)) {
  mod <-
    lm(Resistenza ~ as.numeric(Anno), data = percent[Antibiotico == i])
  cr <-
    cor.test( ~ Resistenza + as.numeric(Anno), percent[Antibiotico == i])
  cr$Antibiotico <- i
  mod$Antibiotico <- i
  models <- c(list(mod), models)
  correlations <- c(list(cr), correlations)
}
rm(mod, cr, i)
# Spataffiatona modelli
lapply(models, summary)
# L'unico modello in cui la variabile anno è significativa è il n° 3, Levofloxacina
# Spataffiatona correlazioni
lapply(correlations, print)
# Anche qui l'unico valore interessante è quello della Levofloxacina

# Ripeto le analisi con i valori assoluti
# Assoluti
ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico
),
data = absolute) +
  geom_line(show.legend = F) +
  geom_point(show.legend = F) +
  facet_grid(Antibiotico ~ ., scales = "free_y") +
  labs(title = "Resistenze Assoluti", x = "", y = "") +
  scale_color_tableau() +
  theme_igray()

ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico
),
data = absolute) +
  geom_point() +
  geom_line() +
  labs(title = "Resistenze Assoluti", y = "", x = "") +
  scale_color_tableau() +
  theme_igray() +
  theme(legend.key = element_blank())

# Mi dice che Clinda, Levo e Oxa partono vicino al punto stimato come intercetta (poco mi importa)
summary(lm(Resistenza ~ as.numeric(Anno) + as.factor(Antibiotico), data =
             absolute))
summary(lmer(Resistenza ~ as.numeric(Anno) + (1 |
                                                Antibiotico), data = absolute))
# L'anno risulta significativo, quindi c'è un trend generale in diminuzione, a prescindere dal
#tipo di antibiotico
models_a <- vector("list", 0)
correlations_a <- vector("list", 0)
for (i in unique(absolute$Antibiotico)) {
  mod <-
    lm(Resistenza ~ as.numeric(Anno), data = absolute[Antibiotico == i])
  cr <-
    cor.test( ~ Resistenza + as.numeric(Anno), absolute[Antibiotico == i])
  cr$Antibiotico <- i
  mod$Antibiotico <- i
  models_a <- c(list(mod), models_a)
  correlations_a <- c(list(cr), correlations_a)
}
rm(mod, cr, i)
# Spataffiatona modelli
lapply(models_a, summary)
# Spataffiatona correlazioni
lapply(correlations_a, print)
# Non risulta più niente di significativoas.data.frame = T, stringsAsFactors=F))
rm(files, datapath)
# Estraggo percentuali e casi e organizzo i dati in formato long
percent <- db[1:6]
setnames(percent, c("Antibiotico", "2015", "2016", "2017"))
percent <- percent[2:6]
percent <-
  melt(
    percent,
    id.vars = "Antibiotico",
    variable.name = "Anno",
    value.name = "Resistenza"
  )
absolute <- db[9:14]
setnames(absolute, c("Antibiotico", "2015", "2016", "2017"))
absolute <- absolute[2:6]
absolute <-
  melt(
    absolute,
    id.vars = "Antibiotico",
    variable.name = "Anno",
    value.name = "Resistenza"
  )
total <- db[8:9, 2:4]
setnames(total, c("2015", "2016", "2017"))
total <- total[1]
total <-
  melt(total, variable.name = "Anno", value.name = "Totale")
# Qualcosa non mi torna, la somma per anno di absolute non mi da i valori di total e
#le percentuali del 2015 superano il 100%, mentre nel 2016 e nel 2017
# Alcuni degli isolati erano multi-resistenti?
# Come sono calcolate le percentuali? Per Oxacillina vengono riportati 13, 6 e 5 casi,
#in percentuale 41.7%, 19.9% e 15.1%. I totali segnalati sono 31 25 e 19. 13/31=0.419, e ok
#ma 6/25 è 0.24 e 5/19 è 0.263, diversi dalle percentuali riportate.
# Perchè nel 2016 ho una resistenza non classificata (25 casi totali, ma nelle resistenze sono 24)?
#Manca un antibiotico?
# Andamenti separati
ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico
),
data = percent) +
  geom_line(show.legend = F) +
  geom_point(show.legend = F) +
  scale_y_continuous(name = "", labels = scales::percent) +
  facet_grid(Antibiotico ~ ., scales = "free_y") +
  labs(title = "Resistenze %", x = "") +
  scale_color_tableau() +
  theme_igray()
# Confronto su unico grafico
ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico
),
data = percent) +
  geom_point() +
  geom_line() +
  scale_y_continuous(name = "", labels = scales::percent) +
  labs(title = "Resistenze %", y = "", x = "") +
  scale_color_tableau() +
  theme_igray() +
  theme(legend.key = element_blank())
# Modelli di regressione
#l'idea è quella di calcolare la pendenza della retta interpolante e testare che sia != 0
summary(lm(Resistenza ~ as.numeric(Anno) + as.factor(Antibiotico), data =
             percent))
# Mi dice che Clinda, Levo e Oxa partono vicino al punto stimato come intercetta (poco mi importa)
# L'anno risulta significativo, quindi c'è un trend generale in diminuzione, a prescindere dal
#tipo di antibiotico
summary(lmer(Resistenza ~ as.numeric(Anno) + (1 |
                                                Antibiotico), data = percent))
# Trend generale confermato anche dal modello a effetti misti
# Modelli e correlazioni su singolo antibiotico
models <- vector("list", 0)
correlations <- vector("list", 0)
for (i in unique(percent$Antibiotico)) {
  mod <-
    lm(Resistenza ~ as.numeric(Anno), data = percent[Antibiotico == i])
  cr <-
    cor.test( ~ Resistenza + as.numeric(Anno), percent[Antibiotico == i])
  cr$Antibiotico <- i
  mod$Antibiotico <- i
  models <- c(list(mod), models)
  correlations <- c(list(cr), correlations)
}
rm(mod, cr, i)
# Spataffiatona modelli
lapply(models, summary)
# L'unico modello in cui la variabile anno è significativa è il n° 3, Levofloxacina
# Spataffiatona correlazioni
lapply(correlations, print)
# Anche qui l'unico valore interessante è quello della Levofloxacina

# Ripeto le analisi con i valori assoluti
# Assoluti
ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico
),
data = absolute) +
  geom_line(show.legend = F) +
  geom_point(show.legend = F) +
  facet_grid(Antibiotico ~ ., scales = "free_y") +
  labs(title = "Resistenze Assoluti", x = "", y = "") +
  scale_color_tableau() +
  theme_igray()

ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico
),
data = absolute) +
  geom_point() +
  geom_line() +
  labs(title = "Resistenze Assoluti", y = "", x = "") +
  scale_color_tableau() +
  theme_igray() +
  theme(legend.key = element_blank())

# Mi dice che Clinda, Levo e Oxa partono vicino al punto stimato come intercetta (poco mi importa)
summary(lm(Resistenza ~ as.numeric(Anno) + as.factor(Antibiotico), data =
             absolute))
summary(lmer(Resistenza ~ as.numeric(Anno) + (1 |
                                                Antibiotico), data = absolute))
# L'anno risulta significativo, quindi c'è un trend generale in diminuzione, a prescindere dal
#tipo di antibiotico
models_a <- vector("list", 0)
correlations_a <- vector("list", 0)
for (i in unique(absolute$Antibiotico)) {
  mod <-
    lm(Resistenza ~ as.numeric(Anno), data = absolute[Antibiotico == i])
  cr <-
    cor.test( ~ Resistenza + as.numeric(Anno), absolute[Antibiotico == i])
  cr$Antibiotico <- i
  mod$Antibiotico <- i
  models_a <- c(list(mod), models_a)
  correlations_a <- c(list(cr), correlations_a)
}
rm(mod, cr, i)
# Spataffiatona modelli
lapply(models_a, summary)
# Spataffiatona correlazioni
lapply(correlations_a, print)
# Non risulta più niente di significativo