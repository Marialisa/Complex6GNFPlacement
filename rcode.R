#!/usr/bin/env Rscript

library(igraph)
library(foreach)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(markovchain)
library(plot3D)
library(plotly)
library(base)
library(scatterplot3d)
require(dplyr)
require(scatterplot3d)
require(fields)
require(ggplot2)
require(graphics)
require(tidyr)
require(devtools)


#population size
N=292;
#layers of multiplex network representation
M=2;


multiplex_adj <- list()  
graphs <- list()

####Data### 
#adjacency matrix from real topology#
graph_layer <- read.csv("path", header = FALSE)
#latency matrix from real topology#
graph_layer_withweights1 <- read.csv("path/250_L.csv", header = FALSE)
#interference matrix from real topology#
graph_layer_withweights2 <- read.csv("path/250_I.csv", header = FALSE)

coo <- read.csv("path/250_coos_with_rus.csv", header = FALSE)
matrix_coord<-as.matrix(coo)

colnames(graph_layer) <- c(1:N)
colnames(graph_layer_withweights1) <- c(1:N)
colnames(graph_layer_withweights2) <- c(1:N)


graph_layer_matrix<-as.matrix(graph_layer)
graph_layer_matrix_withweights1<-as.matrix(graph_layer_withweights1)
graph_layer_matrix_withweights2<-as.matrix(graph_layer_withweights2)

graph_startopology=graph.adjacency(graph_layer_matrix, mode="undirected", weighted= NULL)
grafo_withweights1=graph.adjacency(graph_layer_matrix_withweights1, mode="undirected", weighted = TRUE)
grafo_withweights2=graph.adjacency(graph_layer_matrix_withweights2, mode="undirected", weighted = TRUE)



graph_startopology<-simplify(graph_startopology, remove.multiple = TRUE, remove.loops = TRUE,
                           edge.attr.comb = igraph_opt("edge.attr.comb"))
par(mfrow=c(1,1))
plot_graph_startopology <-plot.igraph(graph_startopology, vertex.label=NA, layout=layout.fruchterman.reingold, vertex.color='orange', vertex.size=3,edge.color="grey80" )
plot(graph_startopology,directed= NULL, layout=matrix_coord,edge.arrow.size= E(graph_startopology)$weight
     , rescale=FALSE,  vertex.label=NA,
     xlim=range(matrix_coord[,1]), ylim=range(matrix_coord[,2]), vertex.label.dist=0,
     edge.color="grey", edge.width=E(graph_startopology)$weight) %>% unique(graph_startopology)

axis(1)
axis(2)


par(mfrow=c(1,2))
plot_graph_layer1 <-plot.igraph(grafo_withweights1,main="Layer 1",vertex.label=NA, layout=layout.fruchterman.reingold,vertex.color='orange', vertex.size=4, vertex.frame.color = "black", edge.width=((E(grafo_withweights1)$weight)/max(E(grafo_withweights1)$weight))*3, edge.color="grey")
plot_graph_layer2 <-plot.igraph(grafo_withweights2,main="Layer 2", vertex.label=NA, layout=layout.fruchterman.reingold,vertex.color='orange', vertex.size=4, vertex.frame.color = "black", edge.width=((E(grafo_withweights2)$weight)/max(E(grafo_withweights2)$weight))*3, edge.color="grey")
V(graph_startopology)$label.cex <- seq(0.6,5,length.out=1)  
p1_withlab<-plot.igraph(graph_startopology, layout=layout.fruchterman.reingold,vertex.color='red',vertex.label.dist=0.3, vertex.label.degree=pi/2, vertex.size=1, edge.color="black")
par(mfrow=c(1,2))
plot_graph_layer1 <-plot.igraph(grafo_withweights1,main="Layer 1",vertex.label=NA, layout=layout.fruchterman.reingold,vertex.color='lightsteelblue2', vertex.size=3, vertex.frame.color = "black", edge.width=E(grafo_withweights1)$weight, edge.color="grey50")
plot_graph_layer2 <-plot.igraph(grafo_withweights2,main="Layer 2", vertex.label=NA, layout=layout.fruchterman.reingold,vertex.color='lightsteelblue2', vertex.size=3, vertex.frame.color = "black", edge.width=E(grafo_withweights2)$weight, edge.color="grey50")
par(mfrow=c(1,2))
p2_withlab<-plot.igraph(grafo_withweights1, layout=layout.fruchterman.reingold,vertex.color='lightsteelblue2',vertex.label.cex=c(0.6),vertex.label.dist=0, vertex.label.color=c("black"),vertex.label.degree=0, vertex.size=0.5,vertex.frame.color = "black",edge.width=E(grafo_withweights1)$weight, edge.color="grey50")
p3_withlab<-plot.igraph(grafo_withweights2, layout=layout.fruchterman.reingold,vertex.color='lightsteelblue2',vertex.label.cex=c(0.6),vertex.label.dist=0, vertex.label.degree=0, vertex.label.color=c("black"), vertex.size=0.5,vertex.frame.color = "black", edge.width=E(grafo_withweights2)$weight, edge.color="grey50")



###Representation of Starting Graph into a Weighted Multiplex Network
multiplex_adj[[1]]=grafo_withweights1
multiplex_adj[[2]]=grafo_withweights2

weightedmultiplex[[1]]=grafo_withweights1
weightedmultiplex[[2]]=grafo_withweights2


strenght_node<-matrix(,nrow=N,ncol=M)
inversepart_nodo<-matrix(,nrow=N,ncol=M)
for (i in 1:M){
    for (j in 1:N){
    strenght_node[j,i]=sum(multiplex_adj[[i]][j,])
    inversepart_nodo[j,i]=sum((multiplex_adj[[i]][j,]/strenght_node[j,i])^2)
  }
}


###Multilink Computation###

graph_layer_matrix_withweights1<- ifelse(graph_layer_matrix_withweights1!=0,1,0)
graph_layer_matrix_withweights2<- ifelse(graph_layer_matrix_withweights2!=0,1,0)

multiplex_adj[[1]]<-graph_layer_matrix_withweights1
multiplex_adj[[2]]<-graph_layer_matrix_withweights2

multilink2<-multiplex_adj[[1]]+multiplex_adj[[2]]
write.table(multilink2, file="/path/matrix_multilink.txt", row.names=FALSE, col.names=FALSE)
data <- read.table("/path/matrix_multilink.txt",header = FALSE, sep = "")
crossdata <- lapply(rownames(data),function(x)sapply(colnames(data),function(y)list(x,y,data[x,y])))
crossdatatmp <- matrix(unlist(crossdata),nrow=3)
crossdatamat <- t(crossdatatmp)
colnames(crossdatamat) <- c("From","To","Value")
crossdatadf <- as.data.frame(crossdatamat,stringsAsFactors=F)
crossdatadf[,3] <- as.numeric(crossdatadf[,3])
crossdatadf
crossdatadf$To<-gsub( "V", "", as.character(crossdatadf$To))
crossdatadf
multilink2<-crossdatadf
A_multi<-c()
A_multi<- ifelse(multilink2$Value==0 | is.na(multilink2$Value), 0, multilink2$Value)
A_multi_matrix<-matrix(A_multi,ncol = N,nrow  = N, byrow = TRUE)
multilink2$amulti<-paste(A_multi)
graphmulti<-graph_from_adjacency_matrix(A_multi_matrix, mode="undirected", weighted = TRUE)
matrice_multi_adj <-get.adjacency(graphmulti)
multi_degree_nodi = degree (graphmulti, v=V(graphmulti))
id_multi_nodes<-c()
id_multi_nodes<-as_ids(V(graphmulti))
multi_data<-data.frame(id_multi_nodes,multi_degree_nodi)


plot_grafo_multidegree <-plot.igraph(graphmulti, vertex.label.dist=1,vertex.label.degree=pi/2, vertex.label.size=1, layout=layout.fruchterman.reingold,vertex.color=V(graphmulti)$color, vertex.size=multi_degree_nodi, edge.width=E(graphmulti)$weight, edge.color="black")
bet<-betweenness(graphmulti, normalized= TRUE)



## potential nodes detection ##

potential <- multi_data[ which(multi_data$multi_degree_nodi>1),] 
potential_rank <- potential[order(-potential$multi_degree_nodi),]
potential_rank_inv <- potential[order(potential$multi_degree_nodi),]

rank_cent<-data.frame(potential_rank$id_multi_nodes,bet[potential_rank$id_multi_nodes],multi_degree_nodi[potential_rank$id_multi_nodes])
rank_cent_ord<-rank_cent[order(-rank_cent$bet),]
sum_test<-bet[potential_rank$id_multi_nodes]+multi_degree_nodi[potential_rank$id_multi_nodes]
rank_cent_ord2<-data.frame(potential_rank$id_multi_nodes, sum_test)

strenght_potenziali <- ifelse(strenght_node[,1]>mean(strenght_node[,1]) & strenght_node[,2]>mean(strenght_node[,2]), 1,0)
sum_strenght<-strenght_node[,1]+strenght_node[,2]
rank_cent_ord2<-data.frame(rank_cent_ord2, sum_strenght[potential_rank$id_multi_nodes], strenght_potenziali[potential_rank$id_multi_nodes] )


N_D <- read.csv("path/250_D.csv", header = FALSE)
node_delay<-seq(0,0.025,0.050,0.075,0.100)
N_D<-N_D+node_delay*shortest.paths(graphmulti)
colnames(N_D) <- c(1:N)
N_D<-as.matrix(N_D)

Xaxis <- seq(0.001,2.5,0.001)
set <- c()
for (nd in node_delay) {
for(j in 1:2500){
  db=Xaxis[j]
  N_D[c(potential_rank$id_multi_nodes),][1,]
  check_pot<-matrix(,nrow = length(rank_cent_ord$potential_rank.id_multi_nodes), ncol = N)
  for (i in 1:length(rank_cent_ord$potential_rank.id_multi_nodes)) {
    check_pot[i,]<-unname(ifelse(N_D[c(rank_cent_ord$potential_rank.id_multi_nodes),][i,]<db,TRUE,FALSE))
  }
  check_pot=check_pot[,- (rank_cent_ord$potential_rank.id_multi_nodes)]
  check_pot
  datacheckpot<- as.data.frame(check_pot)
  datacheckpot2 <- data.frame(t(datacheckpot))
  datacheckpot2<- as.matrix(datacheckpot2)
  m <- structure(c(datacheckpot2), .Dim = c(250L,  42L), .Dimnames = list(NULL, c(1:42)))
  colnames(m)[ifelse(rowSums(m)==0, NA, max.col(m, "first"))]
  unique(colnames(m)[ifelse(rowSums(m)==0, NA, max.col(m, "first"))])
  ordineestratto<-unique(colnames(m)[ifelse(rowSums(m)==0, NA, max.col(m, "first"))])
  indices <- as.integer(ordineestratto)
  
  rank_cent_ord$potential_rank.id_multi_nodes[indices]
  nodes_in_egt<-c(rank_cent_ord2$potential_rank.id_multi_nodes)
  set[j] <- length(nodes_in_egt[indices])}
}


V(graphmulti)$color <- "white" 
V(graphmulti)$color[c(rank_cent_ord2$potential_rank.id_multi_nodes)]<-"tomato"
V(graphmulti)$label[c(rank_cent_ord2$potential_rank.id_multi_nodes)]<-rank_cent_ord2$potential_rank.id_multi_nodes
p1_withlab<-plot.igraph(graphmulti, layout=layout.fruchterman.reingold,vertex.label.dist=0.3, vertex.label.degree=pi/2, vertex.size=2, edge.color="black",edge.width=E(graphmulti)$weight)

remove_nostrenght<-subset(rank_cent_ord2, rank_cent_ord2$strenght_potenziali.potential_rank.id_multi_nodes.==0)
rank_cent_ord2<-anti_join(rank_cent_ord2, remove_nostrenght)
nodes_in_egt<-c(rank_cent_ord2$potential_rank.id_multi_nodes)

capacity<-bet[nodes_in_egt]
sizeforscore<-length(nodes_in_egt)
N_D <- read.csv("path/250_D.csv", header = FALSE)
colnames(N_D) <- c(1:N)
N_D<-as.matrix(N_D)
db=mean(N_D)
score<-c()
incentive<-c()
dataset_for_egt<-data.frame()


for (i in 1:length(nodes_in_egt)) {
  score[i]<-sum(ifelse(N_D[c(nodes_in_egt),][i,]<db,1,0))
}
incentive<-(1/(score*capacity))*10/N
dataset_for_egt<-data.frame(nodes_in_egt,capacity,score,incentive)

write.graph(graphmulti, "path/graphmulti.csv", "edgelist")
write.csv(dataset_for_egt, "path/datasetforegt.csv")

coo <- read.csv("path/250_coos_with_rus.csv", header = FALSE)
matrix_coord<-as.matrix(coo)
V(graphmulti)$color <- "white" 
V(graphmulti)$color[c(rank_cent_ord2$potential_rank.id_multi_nodes)]<-"tomato"
V(graphmulti)$label[c(rank_cent_ord2$potential_rank.id_multi_nodes)]<-rank_cent_ord2$potential_rank.id_multi_nodes
par(mfrow=c(1,1))
plot(graphmulti, layout=matrix_coord,edge.arrow.size= E(graphmulti)$weight
     , rescale=FALSE,  vertex.size=14,
     xlim=range(matrix_coord[,1]), ylim=range(matrix_coord[,2]), vertex.label.dist=0,
     edge.color="grey", edge.width=E(graphmulti)$weight) 
axis(1)
axis(2)

size_node[nodes_in_egt]<-unname(bet[nodes_in_egt])
plot(subgraph(graphmulti, nodes_in_egt),vertex.size=size_node*50,edge.width=E(graphmulti)$weight)

size_node[nodes_in_egt]<-multi_degree_nodi[nodes_in_egt]
plot(subgraph(graphmulti, nodes_in_egt),vertex.color='yellow', vertex.size=size_node,edge.width=E(graphmulti)$weight)


size_node<-c()
for(j in 1:N)
{ size_node[j]=1
}
size_node[c(rank_cent_ord2$potential_rank.id_multi_nodes)]=score

plot(graphmulti, layout=matrix_coord,edge.arrow.size= E(graphmulti)$weight
     , rescale=FALSE,  vertex.size=size_node,
     xlim=range(matrix_coord[,1]), ylim=range(matrix_coord[,2]), vertex.label.dist=0,
     edge.color="grey", edge.width=E(graphmulti)$weight) 
axis(1)
axis(2)

##########################################################################
#####EGT CODE start here #####
##########################################################################

#!/usr/bin/env Rscript
######LIB
library(igraph)
library(foreach)
library(Matrix)
require(doParallel)
registerDoParallel()  

n_step = 200
n_simulazioni = 1 
popolazione= 292
critical_mass = 3 
M= 1
dev_standard_omofilia = 1
homo_on = FALSE
eta_on = TRUE;
game = FALSE
bi = 8;
ci = 0.5;

k_fermi = 20

rimescola = TRUE

multilayer_adiacenze <- list()

grafi <- list()

layout.old  = NULL;

multilayer_omofilie <- list()

multilayer_z <- list()
elenco_nodi <- list()

lista_nodi <- list()

omega_piccolo = 0.5;
w = omega_piccolo*diag(M)  + array(c(omega_piccolo), dim=c(M,M));

RHO_MAX = 1;
RHO_MIN = 0.1;

toto_rho = c();
x1 = c();

for (i in 1:M)
{
    grafo_lay<- graphmulti
  
  matrice_adiacenze <- get.adjacency(grafo_lay, type = "both");
  matrice <- matrix(matrice_adiacenze, nrow=popolazione, ncol=popolazione)
  
  if(rimescola)
  {
    nr<-dim(matrice)[1]
    samp = sample.int(nr)
    matrice=matrice[samp,samp]
  }
  grafi[[i]] = graph.adjacency(matrice,mode="undirected",weighted=NULL)
  if(i==1)
  {
    layout.old <- layout.fruchterman.reingold(grafi[[i]],params=list(weights=E(grafi[[1]])$weight))
  }
  multilayer_adiacenze[[i]] = matrice
  elenco_nodi[[i]] = data.frame(id_nodo = 1:popolazione,  
                                color=(rnorm(popolazione, mean=0, sd=dev_standard_omofilia)),
                                is_coperative = 0,
                                payoff=0,
                                last_time_coperative = 0,
                                num_cooperazioni = 0,
                                capacity= 0,
                                score= 0,
                                incentive=1) 
  
  potential_dataset <- data.frame() 
  potential_dataset <- read.csv("path/datasetforegt.csv", header = TRUE,stringsAsFactors = FALSE, fill=FALSE, sep=",")
  potential_dataset$X <- NULL 
  
  
  for (z in 1:nrow(potential_dataset)){
    elenco_nodi[[i]]$incentive[ which(elenco_nodi[[1]]$id_nodo == potential_dataset$nodi_in_egt[z])] <- potential_dataset$incentive[z] 
    elenco_nodi[[i]]$score[ which(elenco_nodi[[1]]$id_nodo == potential_dataset$nodi_in_egt[z])] <- potential_dataset$score[z]
    elenco_nodi[[i]]$capacity[ which(elenco_nodi[[1]]$id_nodo == potential_dataset$nodi_in_egt[z])] <- potential_dataset$capacity[z]
  } 
  
  multilayer_omofilie[[i]] = matrix(, nrow = popolazione, ncol = popolazione)
  multilayer_z[[i]] = matrix(, nrow = popolazione, ncol = popolazione)
  for(j in 1:popolazione)
  {
    omofilia_nodo_j = elenco_nodi[[i]]['color'][[1]][j]
    for(k in 1:popolazione)
    {
      if(j==k){multilayer_omofilie[[i]][j,k] = 1}
      else
      {
        omofilia_nodo_k = elenco_nodi[[i]]['color'][[1]][k]
        multilayer_omofilie[[i]][j,k] = 1/(1 + abs(omofilia_nodo_j - omofilia_nodo_k))
      }
      multilayer_z[[i]][j,k] = multilayer_omofilie[[i]][j,k] * multilayer_adiacenze[[i]][j,k]
    }
  }
}

riga_matricione <- vector("list", M);  
matricione = NULL
for(i in 1:M)
{
  for(j in 1:M)
  {
    
    riga_matricione[[i]] = cbind(riga_matricione[[i]], w[i,j]* t(multilayer_z[[j]])  )
  }
  
  if(i == 1)
  {
    matricione = riga_matricione[[1]] 
  }
  else 
  { 
    matricione = rbind(matricione, riga_matricione[[i]]);
  }
}
multicentrum <- graph.adjacency(matricione, weighted=TRUE, mode="directed");
centrality <- evcent (multicentrum, directed = TRUE)

vertex_size <- rep(1, popolazione)
LAMBDA <- rep(0, popolazione)
for(i in 1:popolazione)
{
  for(k in 1:M)
  {
    LAMBDA[i] = LAMBDA[i] + centrality$vector[popolazione*(k-1) + i]
  }
  vertex_size[i] = LAMBDA[i] 
}
vertex_size = log10(1000*vertex_size) + 4


LAMBDA = LAMBDA/sum(LAMBDA)
indici_critical_mass = sort.int(LAMBDA, index.return=TRUE, decreasing = TRUE)
indici_critical_mass = indici_critical_mass$ix[1:critical_mass]

for(j in indici_critical_mass)
{
  for (i in 1:M)
  {
    elenco_nodi[[i]]["is_coperative"][j,1] = 1
  }
}


for(k in 1:M)
{
  lista_nodi[[k]] = list()
  for(i in 1:popolazione)
  {
    nodo <- list();
    
    index = 1;
    for(j in 1:popolazione)
    {
      
      if(multilayer_adiacenze[[k]][i,j]==1)
      {
        nodo[[index]]=j;  
        index=index+1;    
      } 
    }
    
    lista_nodi[[k]][[i]] = unlist(nodo);
    
  }
}


matrixList<-list(multilayer_z)
if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]

dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
finalDimension<-sum(dimensions)
finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
index<-1
for(k in 1:length(dimensions)){
  finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
  index<-index+dimensions[k]
}
Z_L = finalMatrix


omega_piccolo = 0.2;
matrice_omega = -omega_piccolo*diag(M)  + array(c(omega_piccolo), dim=c(M,M));

matrice_identica_N_X_N = diag(popolazione)
m_gotica = kronecker(matrice_omega,matrice_identica_N_X_N) + Z_L
matrice_g = expm(m_gotica)


rho <- c();
for(h in 1:n_step)
{
  rho[h] <- 0;
  elenco_nodi[["payoff"]] = 0
  
  for(lay in 1:M)
  {
    for(index in sample(popolazione))
    {
      
      adiacenze_temp = unlist(lista_nodi[[lay]][[index]]);
      coperativo = TRUE;
      if(elenco_nodi[[lay]][["is_coperative"]][index] == 0)
      {
        coperativo = FALSE
      }
      for(index_altro_player in adiacenze_temp)
      {
          if(index < index_altro_player)
        {
          antagonista_coperativo = TRUE;
          if(elenco_nodi[[lay]][["is_coperative"]][index_altro_player] == 0)
          {
            antagonista_coperativo = FALSE
          }
          if(game)
          {
            #facciamo i 4 casi del PD game
            if(coperativo == TRUE && antagonista_coperativo == TRUE)
            {
              elenco_nodi[[lay]][["payoff"]][index] = elenco_nodi[[lay]][["payoff"]][index] + bi -ci + 10*elenco_nodi[[lay]][["incentive"]][index] ;
              elenco_nodi[[lay]][["payoff"]][index_altro_player] = elenco_nodi[[lay]][["payoff"]][index_altro_player] + bi -ci;
              elenco_nodi[[lay]][["num_cooperazioni"]][index] = elenco_nodi[[lay]][["num_cooperazioni"]][index] +1;
              elenco_nodi[[lay]][["num_cooperazioni"]][index_altro_player] = elenco_nodi[[lay]][["num_cooperazioni"]][index_altro_player] +1;
            }
            else if(coperativo == TRUE && antagonista_coperativo == FALSE)
            {
              elenco_nodi[[lay]][["payoff"]][index] = elenco_nodi[[lay]][["payoff"]][index] - ci + 10*elenco_nodi[[lay]][["incentive"]][index];
              elenco_nodi[[lay]][["payoff"]][index_altro_player] = elenco_nodi[[lay]][["payoff"]][index_altro_player] + bi;
              elenco_nodi[[lay]][["num_cooperazioni"]][index] = elenco_nodi[[lay]][["num_cooperazioni"]][index] +1;
            }
            else if(coperativo == FALSE && antagonista_coperativo == TRUE)
            {
              elenco_nodi[[lay]][["payoff"]][index] = elenco_nodi[[lay]][["payoff"]][index] + bi;
              elenco_nodi[[lay]][["payoff"]][index_altro_player] = elenco_nodi[[lay]][["payoff"]][index_altro_player] - ci;
              elenco_nodi[[lay]][["num_cooperazioni"]][index_altro_player] = elenco_nodi[[lay]][["num_cooperazioni"]][index_altro_player] +1;
            }
            else if(coperativo == FALSE && antagonista_coperativo == FALSE)
            {
              elenco_nodi[[lay]][["payoff"]][index] = elenco_nodi[[lay]][["payoff"]][index] + 0;
              elenco_nodi[[lay]][["payoff"]][index_altro_player] = elenco_nodi[[lay]][["payoff"]][index_altro_player] + 0;
            }
          }
          else
          {
            #facciamo i 4 casi nell'altro gioco - SD game
            if(coperativo == TRUE && antagonista_coperativo == TRUE)
            {
              elenco_nodi[[lay]][["payoff"]][index] = elenco_nodi[[lay]][["payoff"]][index] + bi -ci*0.5 + 10*elenco_nodi[[lay]][["incentive"]][index] ;
              elenco_nodi[[lay]][["payoff"]][index_altro_player] = elenco_nodi[[lay]][["payoff"]][index_altro_player] + bi -ci*0.5;
              elenco_nodi[[lay]][["num_cooperazioni"]][index] = elenco_nodi[[lay]][["num_cooperazioni"]][index] +1;
              elenco_nodi[[lay]][["num_cooperazioni"]][index_altro_player] = elenco_nodi[[lay]][["num_cooperazioni"]][index_altro_player] +1;
            }
            else if(coperativo == TRUE && antagonista_coperativo == FALSE)
            {
              elenco_nodi[[lay]][["payoff"]][index] = elenco_nodi[[lay]][["payoff"]][index] + bi - ci +10*elenco_nodi[[lay]][["incentive"]][index];
              elenco_nodi[[lay]][["payoff"]][index_altro_player] = elenco_nodi[[lay]][["payoff"]][index_altro_player] + bi;
              elenco_nodi[[lay]][["num_cooperazioni"]][index] = elenco_nodi[[lay]][["num_cooperazioni"]][index] +1;
            }
            else if(coperativo == FALSE && antagonista_coperativo == TRUE)
            {
              elenco_nodi[[lay]][["payoff"]][index] = elenco_nodi[[lay]][["payoff"]][index] + bi;
              elenco_nodi[[lay]][["payoff"]][index_altro_player] = elenco_nodi[[lay]][["payoff"]][index_altro_player] + bi - ci;
              elenco_nodi[[lay]][["num_cooperazioni"]][index_altro_player] = elenco_nodi[[lay]][["num_cooperazioni"]][index_altro_player] +1;
            }
            else if(coperativo == FALSE && antagonista_coperativo == FALSE)
            {
              elenco_nodi[[lay]][["payoff"]][index] = elenco_nodi[[lay]][["payoff"]][index] + 0;
              elenco_nodi[[lay]][["payoff"]][index_altro_player] = elenco_nodi[[lay]][["payoff"]][index_altro_player] + 0;
            }
          }
        }
      } 
      
    }
  }  
  
  #Iteriamo su tutti i layers
  for(lay in 1:M)
  {
    for(temp_index in 1:popolazione)
    {
      elenco_nodi[[lay]][["last_time_coperative"]][temp_index] = elenco_nodi[[lay]][["is_coperative"]][temp_index]
    }
  }
  
  for(lay in 1:M)
  {
    num_coperativi = 0;
    colors <- rep(rgb(1,1,1), popolazione)
    
    for(index in sample(popolazione))
    {
      if(elenco_nodi[[lay]][["is_coperative"]][index] == 1) colors[index] = rgb(1,0,0)
      adiacenze_temp = unlist(lista_nodi[[lay]][[index]]);
       random_index = 0;
      if(length(adiacenze_temp) == 1)
      {
        random_index = adiacenze_temp[1]
      }
      else
      {
        random_index = sample(adiacenze_temp, 1);
      }
      payoff_tizio = elenco_nodi[[lay]][["payoff"]][random_index];
      payoff_mio = elenco_nodi[[lay]][["payoff"]][index]
      
      coperativo = TRUE;
      if(elenco_nodi[[lay]][["last_time_coperative"]][index] == 0)
      {
        coperativo = FALSE
      }
      else
      {
        num_coperativi = num_coperativi + 1
      }
      delta_omofilia = abs(elenco_nodi[[lay]][["color"]][index] - elenco_nodi[[lay]][["color"]][random_index]);
      prob = 1
      #DEVO COMPUTARE IL RHO_X
      RHO_X = 1;
      if(homo_on)
      {           
        if(M >1)
        {
          alpha = lay;
          
          escluso = alpha;
          max = M;
          result = 0;
          if(escluso == 1)
          {
            result = sample(2:max, 1)
          }
          else if(escluso == max)
          {
            result = sample(1:(max - 1), 1)
          }
          else
          {
            result = sample(c(1:(escluso - 1), (escluso + 1):max),1)
          }
          beta = result;
          elem_x= index;
          numeratore = 0;
          denominatore = 0;
          for(elem_y in c(index, unlist(lista_nodi[[beta]][[index]])))
          {
            comunicability = matrice_g[(alpha - 1)*popolazione + elem_x, (beta - 1)*popolazione + elem_y]
            
            
            denominatore = denominatore + comunicability;
             if(elenco_nodi[[alpha]][["last_time_coperative"]][[elem_x]] == elenco_nodi[[beta]][["last_time_coperative"]][[elem_y]] )
            {
              numeratore = numeratore + comunicability;
            }
          } 
          RHO_X = 1 - (RHO_MAX - RHO_MIN)*numeratore/denominatore;
        }
        if(!eta_on){
          RHO_X = 1;
        }	
        prob = RHO_X*(1 /( 1 + exp((payoff_mio - payoff_tizio)/(k_fermi*delta_omofilia) ) )  )
      }
      else
      {
        if(!eta_on){
          RHO_X = 1;
        }
        prob = RHO_X*(1 /( 1 + exp((payoff_mio - payoff_tizio)/k_fermi) )  )
      }
      #commutiamo
      if(prob > runif(1, 0.0, 1.0))
      {
        elenco_nodi[[lay]][["is_coperative"]][index] = elenco_nodi[[lay]][["last_time_coperative"]][random_index]
      }
      
    }
    #######
    #Stampiamo il layer
    if(h %% 10 == 1)
    {
      #layout.new <- layout.fruchterman.reingold(grafi[[lay]],params=list(niter=10,start=layout.old,weights=E(grafi[[lay]])$weight,maxdelta=1))
      
      #           file_ = paste('c:/ciao/ciao/layer_', lay, '_step_', h, '.png')
      #           png(filename=file_,width=4,height=4,units="in",res=1200)
      #           #file_ = paste('c:/ciao/ciao/layer_', lay, '_step_', h, '.pdf')
      #           #pdf(file=file_,width=4,height=4)
      #           plot.igraph(grafi[[lay]],vertex.label=NA,layout=layout.new, layout=layout.fruchterman.reingold, vertex.color=colors, vertex.size=vertex_size)
      #           dev.off()
      
      #Sys.sleep(0) 
    }
    rho[h] <- rho[h] + num_coperativi/popolazione ;
  }
  rho[h] <- rho[h] / M;
}


elenco_nodi

elenco_nodi_coop<-data.frame(elenco_nodi[[1]]$id_nodo,elenco_nodi[[1]]$num_cooperazioni)
elenco_nodi_coop_selected <- elenco_nodi_coop[elenco_nodi_coop$elenco_nodi..1...id_nodo %in% nodi_in_egt, ]
elenco_nodi_coop_selected <- elenco_nodi_coop_selected[match(nodi_in_egt, elenco_nodi_coop_selected$elenco_nodi..1...id_nodo),]


nodes_coop_selected<-elenco_nodi_coop_selected

##########################################################################
#####EGT CODE ends here #####
##########################################################################

nodes_coop_selected

last_dataset<-data.frame(nodes_coop_selected, multi_degree_nodi[nodes_coop_selected$elenco_nodi..1...id_nodo])
filtro<-last_dataset$elenco_nodi..1...num_cooperazioni+last_dataset$multi_degree_nodi.nodes_coop_selected.elenco_nodi..1...id_nodo.
last_dataset_new<-data.frame(nodes_coop_selected, filtro)
verylast_dataset<-data.frame(nodes_coop_selected$elenco_nodi..1...id_nodo,filtro)

## clustering ##
distance_matrix<-matrix(,nrow=length(nodes_in_egt),ncol=length(nodes_in_egt))
for (s in 1:length(nodes_in_egt)){
  for (j in 1:length(nodes_in_egt)){
    distance_matrix[s,j]=unname(bet[c(rank_cent_ord2$potential_rank.id_multi_nodes)]/verylast_dataset$filtro)[s]-unname(bet[c(rank_cent_ord2$potential_rank.id_multi_nodes)]/verylast_dataset$filtro)[j]
  }
}


distance_matrix<-ifelse (distance_matrix<0,0,distance_matrix)


hc1 <- hclust(as.dist(distance_matrix), method = "complete" )
plot(hc1)
clusters.list = rect.hclust(hc1 , k = 3, border = 2:6)
clusters = cutree(hc1, k = 3)
cutree(hc1,1) 

d<-data.frame(cutree(hc1,3)) 
rownames(d) <- nodes_in_egt
d_new<-d
size_node[nodes_in_egt]<-(unname(bet[c(rank_cent_ord2$potential_rank.id_multi_nodes)]/verylast_dataset$filtro))*10000
V(graphmulti)[nodes_in_egt]$color<-clusters
plot(subgraph(graphmulti, nodes_in_egt),vertex.size=size_node)


V(graphmulti)$color<-"white"
V(graphmulti)[nodes_in_egt]$color<-clusters
plot(graphmulti, layout=matrix_coord,edge.arrow.size= E(graphmulti)$weight
     , rescale=FALSE,  vertex.size=14,  vertex.label.size=1, 
     xlim=range(matrix_coord[,1]), ylim=range(matrix_coord[,2]), vertex.label.dist=0,
     edge.color="grey", edge.width=E(graphmulti)$weight) 
axis(1)
axis(2)

V(graphmulti)$color<-"white"
V(graphmulti)[nodes_in_egt]$color<-clusters
plot(graphmulti, layout=matrix_coord,edge.arrow.size= E(graphmulti)$weight
     , rescale=FALSE,  vertex.size=14, vertex.label=NA,
     xlim=range(matrix_coord[,1]), ylim=range(matrix_coord[,2]), vertex.label.dist=0,
     edge.color="grey", edge.width=E(graphmulti)$weight) 
axis(1)
axis(2)
colrs <-c("forestgreen", "steelblue1", "orange2","white")
legend(x=-1, y=11, c("NGC","CU", "DU","RU"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)


data_score_node<-data.frame(nodes_in_egt,clusters, verylast_dataset$filtro)

### placement vectors###

vectors_NGC<-data_score_node[data_score_node$clusters==3,]$nodes_in_egt
vectors_CU<-data_score_node[data_score_node$clusters==2,]$nodes_in_egt
vectors_DU<-data_score_node[data_score_node$clusters==1,]$nodes_in_egt
vectors_RU<-elenco_nodi_coop_theothers$elenco_nodi..1...id_nodo

path_RU<-list()
ngc_corrisp<-c()
for (i in 1:length(vectors_RU)){
  short_paths_forpu<-shortest_paths(graphmulti,vectors_RU[i], to=vectors_NGC, mode = c("all"), 
                                    weights = NULL, output = c("vpath"),predecessors = TRUE, 
                                    inbound.edges = FALSE)
  lenght_path_ru_ngc<-c()
  for (j in 1:length(vectors_NGC)) {
    lenght_path_ru_ngc[j]<-length(short_paths_forpu$vpath[[j]])
  }
  
  f<-which.min(lenght_path_ru_ngc)
  ngc_corrisp[i]<-short_paths_forpu$vpath[[f]][length(short_paths_forpu$vpath[[f]])]
  path_RU[[i]]<-c(unlist(short_paths_forpu$vpath[f]))
  
}

du_corrisp_1<-c()
for (i in 1:length(vectors_RU)){
  du_corrisp_1[i]<-vectors_DU[match(path_RU[[i]][2],vectors_DU)]
}



du_corrisp_2<-c()
for (i in 1:length(vectors_RU)){
   du_corrisp_2[i]<-vectors_DU[match(path_RU[[i]][3],vectors_DU)]
}


path_CU<-list()
CU_corrisp<-c()
for (i in 1:length(vectors_RU)){
  short_paths_forpu_RUCU<-shortest_paths(graphmulti,vectors_RU[i], to=vectors_CU, mode = c("all"), 
                                         weights = NULL, output = c("vpath"),predecessors = TRUE, 
                                         inbound.edges = FALSE)
  lenght_path_RU_CU<-c()
  for (j in 1:length(vectors_CU)) {
    lenght_path_RU_CU[j]<-length(short_paths_forpu_RUCU$vpath[[j]])
  }
  
  f<-which.min(lenght_path_RU_CU)
  CU_corrisp[i]<-short_paths_forpu_RUCU$vpath[[f]][length(short_paths_forpu_RUCU$vpath[[f]])]
  path_CU[[i]]<-c(unlist(short_paths_forpu_RUCU$vpath[f]))
  
}

cu_corrisp_1<-c()
for (i in 1:length(vectors_RU)){
  cu_corrisp_1[i]<-vectors_CU[match(path_RU[[i]][2],vectors_CU)]
}

cu_corrisp_2<-c()
for (i in 1:length(vectors_RU)){
  cu_corrisp_2[i]<-vectors_CU[match(path_RU[[i]][3],vectors_CU)]
}

cu_corrisp_3<-c()
for (i in 1:length(vectors_RU)){
  cu_corrisp_3[i]<-vectors_CU[match(path_RU[[i]][4],vectors_CU)]
}

RU_NGC_path<-data.frame(vectors_RU,du_corrisp_1,du_corrisp_2,cu_corrisp_1,cu_corrisp_2,cu_corrisp_3, ngc_corrisp)

check_dim1<- ifelse(length(vectors_RU)>length(vectors_DU),TRUE,FALSE)
check_dim1
check_dim2<- ifelse(length(vectors_DU)>length(vectors_CU),TRUE,FALSE)
check_dim2
check_dim3<- ifelse(length(vectors_CU)>length(vectors_NGC),TRUE,FALSE)
check_dim3

index_corrisp1<-ifelse(RU_NGC_path$cu_corrisp!=0 & RU_NGC_path$cu_corrisp2!=0,RU_NGC_path$cu_corrisp,0)
index_corrisp1
index_corrisp2<-ifelse(RU_NGC_path$cu_corrisp!=0 & RU_NGC_path$cu_corrisp2!=0,RU_NGC_path$cu_corrisp2,0)
index_corrisp2


rownames(data_score_node) <- data_score_node$nodes_in_egt
data_score_node<-data_score_node[order(data_score_node$nodes_in_egt),]

RU_NGC_path_test1 <- RU_NGC_path 

RU_NGC_path_test1<-RU_NGC_path %>% mutate(du_corrisp = coalesce(.$du_corrisp_1, .$du_corrisp_2)) %>% mutate(cu_corrisp = coalesce(.$cu_corrisp_1, .$cu_corrisp_2, .$cu_corrisp_3)) 


RU_NGC_path$du_corrisp<- RU_NGC_path_test1$du_corrisp
RU_NGC_path$cu_corrisp<- RU_NGC_path_test1$cu_corrisp


DATAFRAME_FINAL_PATH<-data.frame(RU_NGC_path$vectors_RU,RU_NGC_path$du_corrisp,RU_NGC_path$cu_corrisp,RU_NGC_path$ngc_corrisp)


write.graph(graphmulti, "path/graphmulti.csv", "edgelist")
write.csv(dataset_for_egt, "path/datasetforegt.csv")


colnames(data_score_node)[1] <- "Nodes"
colnames(data_score_node)[2]<-"Role"
colnames(data_score_node)[3]<-"Cooperation Score"

data_score_node$Role<-factor(data_score_node$Role)
levels(data_score_node$Role)

library(plyr)
data_score_node$Role<-revalue(data_score_node$Role, c("1"="DU", "2"="CU", "3"="NGC"))
colnames(DATAFRAME_FINAL_PATH)[1] <- "Ru"
colnames(DATAFRAME_FINAL_PATH)[2]<-"DU"
colnames(DATAFRAME_FINAL_PATH)[3]<-"CU"
colnames(DATAFRAME_FINAL_PATH)[4]<-"NGC"

write.csv(data_score_node,  "path/data_score_node_250.csv")
write.csv(DATAFRAME_FINAL_PATH,  "path/dataframefinal_250.csv")

sink("path/paths_250.txt")
print(path_RU)
sink()

