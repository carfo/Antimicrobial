library(data.table)
library(ggplot2)
library(ggthemes)
library(trend)
library(lme4)
library(lmtest)
library(broom)

db <-
  data.table(read.csv("./Data/Dati resistenza S.aureus.csv", dec = "."))
db[, X := NULL]
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
percent[,Resistenza := Resistenza/100]
absolute <- db[8:13]
setnames(absolute, c("Antibiotico", "2015", "2016", "2017"))
absolute <- absolute[2:6]
absolute <-
  melt(
    absolute,
    id.vars = "Antibiotico",
    variable.name = "Anno",
    value.name = "Resistenza"
  )
total <- db[7,-1]
setnames(total, c("2015", "2016", "2017"))
total <-
  melt(total, variable.name = "Anno", value.name = "Totale")
# Valori assoluti e percentuali ora tornano
# Qualcosa però continua a non tornare, i valori di total per il 2016 e 2017 sono superiori ai
#valori assoluti
# Questo non è possibile, a meno che non abbia dei casi non noti (in particolare almeno 5 nel
#2016 e 6 nel 2017) o resistenze a farmaci non considerati
# E' un buco importante a mio avviso, non so nemmeno di quanti pazienti non ho informazioni

# Andamenti separati
ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico),
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
  group = Antibiotico),
  data = percent) +
  geom_jitter(width = 0.02) +
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
correlations <- NULL
for (i in unique(percent$Antibiotico)) {
  mod <-
    lm(Resistenza ~ as.numeric(Anno), data = percent[Antibiotico == i])
  cr <-broom::tidy(cor.test( ~ Resistenza + as.numeric(Anno),
                             percent[Antibiotico == i]))
  cr$Antibiotico <- i
  mod$Antibiotico <- i
  models <- c(list(mod), models)
  correlations <- rbind.data.frame(cr, correlations)
}
rm(mod, cr, i)
# Spataffiatona modelli
lapply(models, summary)

# Plot della regressione lineare (ed SE)
ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico
),
data = percent) +
  geom_smooth(method = "lm",
              se = T,
              show.legend = F) +
  geom_point(show.legend = F) +
  scale_y_continuous(name = "", labels = scales::percent) +
  facet_grid(Antibiotico ~ ., scales = "free_y") +
  labs(title = "Resistenze %", x = "") +
  scale_color_tableau() +
  theme_igray()

ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico),
  data = percent) +
  geom_smooth(method = "lm",
              se = F,
              show.legend = F) +
  geom_jitter(width = 0.02) +
  scale_y_continuous(name = "", labels = scales::percent) +
  labs(title = "Resistenze %", y = "", x = "") +
  scale_color_tableau() +
  theme_igray() +
  theme(legend.key = element_blank())
# L'unico modello in cui la variabile anno è significativa è il n° 3, Levofloxacina
# Spataffiatona correlazioni
correlations
# Anche qui l'unico valore interessante è quello della Levofloxacina


# Ripeto le analisi con i valori assoluti
# Assoluti
ggplot(aes(
  y = Resistenza,
  x = Anno,
  col = Antibiotico,
  group = Antibiotico),
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
  group = Antibiotico),
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
correlations_a <- NULL
chi2_tests <- NULL
for (i in unique(absolute$Antibiotico)) {
  mod <-
    lm(Resistenza ~ as.numeric(Anno), data = absolute[Antibiotico == i])
  cr <- broom::tidy(cor.test( ~ Resistenza + as.numeric(Anno),
                             absolute[Antibiotico == i]))
  chi2 <-
    broom::tidy(chisq.test(x = absolute[Antibiotico == i, Resistenza],
                           simulate.p.value = T))
  chi2$Antibiotico <- i
  cr$Antibiotico <- i
  mod$Antibiotico <- i
  models_a <- c(list(mod), models_a)
  correlations_a <- rbind.data.frame(cr, correlations_a)
  chi2_tests <- rbind.data.frame(chi2, chi2_tests)
}
rm(mod, cr, chi2, i)
# Spataffiatona modelli
lapply(models_a, summary)
# Spataffiatona correlazioni
correlations_a[c(1, 2, 7)]
# Test Chi2
chi2_tests[c(2, 5)]
chi2
# Sono fuori strada? La regressione mi cerca una relazione lineare tra i dati, ma potrebbe non
#essere lineare. Il Chi2 da risultati paragonabili a quelli presentati