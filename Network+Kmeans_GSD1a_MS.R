####libraries
library(flextable)
library(GGally)
library(ggraph)
library(gutenbergr)
library(igraph)
library(Matrix)
library(network)
library(quanteda)
library(sna)
library(tidygraph)
library(tidyverse)
library(tm)
library(tibble)
library(factoextra)


##############################################################################
# repeat once for untreated and once for treated 
Data_UT_va = read.csv('D:/MiguelW12/Documents/node_list_ut.csv', header=TRUE)
Data_UT_ed = read.csv('D:/MiguelW12/Documents/edge_list_both.csv', header=TRUE)

rownames(Data_UT_va) <- Data_UT_va$Feature
ig <- igraph::graph_from_data_frame(d=Data_UT_ed, vertices=Data_UT_va, directed = FALSE)

tg <- tidygraph::as_tbl_graph(ig) %>% 
  tidygraph::activate(nodes) %>% 
  dplyr::mutate(label=name)

set.seed(12345)
# node size
tg %>%
  ggraph(layout = "fr") +
  geom_edge_arc(colour= "gray50",
                lineend = "round",
                strength = .1,
                alpha = .1) +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 point.padding = unit(0.2, "lines"), 
                 colour="gray10") +
  theme_graph(background = "white") +
  guides(edge_width = FALSE,
         edge_alpha = FALSE)

v.size <- V(tg)$Value
v.size


set.seed(12345)
###############
tg %>%
  ggraph(layout = "fr") +
  geom_edge_arc(colour= "gray50",
                lineend = "round",
                strength = .1) +
  geom_node_point(size=log(v.size)*2) +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 point.padding = unit(0.2, "lines"), 
                 size=sqrt(v.size), 
                 colour="gray10") +
  scale_edge_width(range = c(0, 2.5)) +
  scale_edge_alpha(range = c(0, .3)) +
  theme_graph(background = "white") +
  guides(edge_width = FALSE,
         edge_alpha = FALSE)



E(tg)$weight <- E(tg)$Strength
# inspect weights
head(E(tg)$weight, 10)

set.seed(12345)
# edge size 
tg %>%
  ggraph(layout = "fr") +
  geom_edge_arc(colour= "gray50",
                lineend = "round",
                strength = .1,
                aes(edge_width = weight,
                    alpha = weight)) +
  geom_node_point(size=log(v.size)*2) +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 point.padding = unit(0.2, "lines"), 
                 size=sqrt(v.size), 
                 colour="gray10") +
  scale_edge_width(range = c(2, 5)) +
  scale_edge_alpha(range = c(.5, .8)) +
  theme_graph(background = "white") +
  theme(legend.position = "top") +
  guides(edge_width = FALSE,
         edge_alpha = FALSE)

# define colors (by pathway)
mon <- c("FLUX", "TFEB")
cap <- c("SIRT1")
oth <- c("PGC1a", "COX")
met <- c("GLYCOGEN","ATP", "NAD/NADH")
# color vectors
Family <- dplyr::case_when(sapply(tg, "[")$nodes$name %in% mon ~ "FLUX",
                           sapply(tg, "[")$nodes$name %in% cap ~ "SIRT1",
                           sapply(tg, "[")$nodes$name %in% met ~ "GLYCOGEN",
                           TRUE ~ "Other")

Family
set.seed(12345)
# colors-1
tg %>%
  ggraph(layout = "eigen") +
  geom_edge_arc(colour= "black",
                lineend = "round",
                strength = .1,
                aes(edge_width = weight,
                    alpha = weight)) +
  geom_node_point(size=log(v.size)*4, 
                  aes(color=Family)) +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 point.padding = unit(0.2, "lines"), 
                 size=sqrt(v.size), 
                 colour="gray10") +
  scale_edge_width(range = c(1, 3)) +
  scale_edge_alpha(range = c(.5, .8)) +
  theme_graph(background = "white") +
  theme(legend.position = "top") +
  guides(edge_width = FALSE,
         edge_alpha = FALSE)









set.seed(12345)
# color 2
tg %>%
  ggraph(layout = "eigen") +
  geom_edge_arc(colour= "black",
                lineend = "round",
                strength = .1,
                aes(edge_width = weight,
                    alpha = weight)) +
  geom_node_point(size=(v.size)/500, 
                  aes(color=Family)) +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 point.padding = unit(0.2, "lines"), 
                 size=2, 
                 colour="gray10") +
  scale_edge_width(range = c(1, 3)) +
  scale_edge_alpha(range = c(.5, .8)) +
  theme_graph(background = "white") +
  theme(legend.position = "top") +
  guides(edge_width = FALSE,
         edge_alpha = FALSE)







set.seed(12345)
# color 3 
tg %>%
  ggraph(layout = "eigen") +
  geom_edge_arc(colour= "black",
                lineend = "round",
                strength = .1,
                aes(edge_width = weight,
                    alpha = weight)) +
  geom_node_point(size=(v.size)/500, 
                  aes(color=Family)) +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 point.padding = unit(0.2, "lines"), 
                 size=2, 
                 colour="gray10") +
  scale_edge_width(range = c(1, 3)) +
  scale_edge_alpha(range = c(.5, .8)) +
  theme_graph(background = "white") +
  theme(legend.position = "top") +
  guides(edge_width = FALSE,
         edge_alpha = FALSE)

################################
## k-means clustering 
Data = read.csv('D:/MiguelW12/Documents/K_MEANS_DATA_T.csv', header=TRUE)
rownames(Data) <- Data$Sample
Data = subset(Data, select = -c(Sample))

df <- scale(Data)
set.seed(123)
km.res <- kmeans(df, 3, nstart = 25)

fviz_cluster(km.res, data = df, show.clust.cent=TRUE,
             ellipse=TRUE, ellipse.type="confidence", ellipse.level=0.95,
)





