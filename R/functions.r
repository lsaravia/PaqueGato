# Function definitions

#' Calculates the Shannon Diversity Index using ln base e
#'
#' @param species vector of abundances of species
#'
#' @return the Shannon Diversity Index
#' 
#' @export
#'
#' @examples 
#' p <- c(100,1000,2000)
#' shannon(p)
#' 
shannon <- function(species)
{
  species <- species[species>0]
  species <- species/sum(species)
  H <- -sum(species* log(species))
  return(H)
}


#' Fix the problems of chlorophyll-a data set 
#'
#' @param chla data.frame with original data
#'
#' @return data.frame with fixed data
#' @export
#'
#' @examples 
fixClorophylData <- function(chla)
{
  actual_year <- chla$year[1]
  chla$Date <- as.Date("1900-01-01")
  
  for(i in 1:nrow(chla))
  {
    if( is.na(chla$Year[i]))
        chla$Year[i] <- actual_year
    else
        actual_year <- chla$Year[i]
    
    fecha <- ymd(paste(chla$Year[i],chla$Month[i], 1))
    if( is.na(fecha)){
      if(chla$Month[i]=="Mar")
        fecha <- ymd(paste(chla$Year[i],3, 1))
      else {
        fecha <- dmy(chla$Month[i])
        
        if(is.na(fecha)) {
          fecha <- mdy(chla$Month[i])
          
          if(is.na(fecha)) {
            temp <- strsplit(chla$Month[i]," ")[[1]][3]
            fecha <- dmy(temp)
          }
        }
      }
    }
    chla$Date[i] <- fecha

  }
  
  chla$IntegE1 <- abs(chla$IntegE1)
  chla$IntegE2 <- abs(as.numeric(chla$IntegE2))
  
  return(chla)
}


#' Read ecological networks in CSV format as edge list or adyacency matrix
#'
#' @param fileName Filename of the csv formated network
#'
#' @return an igraph object
#' @export
#'
#' @examples readEcoNetwork("econetwork.csv")
readEcoNetwork <- function(fileName){
  
  g <- lapply(fileName, function(fname){
    
    web <- read.csv(fname,  header = T,check.names = F)
    
    if( ncol(web)==2 ){
      web <- web[,c(2,1)]
      
      g <- graph_from_data_frame(web)
      
    } else {
      if( (ncol(web)-1) == nrow(web)  ) {                   # The adjacency matrix must be square
        g <- graph_from_adjacency_matrix(as.matrix(web[,2:ncol(web)]))
        
      } else {
        g <- NULL
        warning("Invalid file format: ",fileName)
      }   
    }
    
  })
  return(g)
}



#' Plot ecological network organized by trophic level, with node size determined by the node degree
#'
#' @param ig igraph object
#' @param vertexLabel logical plot vertex labels
#' @param vertexSizeFactor numeric factor to determine the size of the label with degree
#'
#' @return a plot
#' @export
#'
#' @examples
plotEcoNetworkTrophLevel <- function(g,vertexLabel=FALSE,vertexSizeFactor=5){
  
  deg <- degree(g, mode="all") # calculate the degree: the number of edges 
  # or interactions
  
  V(g)$size <- log10(deg)*vertexSizeFactor+vertexSizeFactor    # add node degrees to igraph object
  
  V(g)$frame.color <- "white"    # Specify plot options directly on the object 
  
  V(g)$color <- "orange"         #

  if(!vertexLabel)
    V(g)$label <- NA
  
  require(NetIndices)
  
  tl <- TrophInd(get.adjacency(g,sparse=F))  # Calculate the trophic level
  
  lMat <-matrix(                      # Layout matrix to specify the position of each vertix
    nrow=vcount(g),                   # Rows equal to the number of vertices (species)
    ncol=2
  )
  
  lMat[,2]<-jitter(tl$TL,0.1)              # y-axis value based on trophic level
  lMat[,1]<-runif(vcount(g))               # randomly assign along x-axis
  
  require(RColorBrewer)                # Use color brewer palettes
  
  colTL <-as.numeric(cut(tl$TL,11))   # Divide trophic levels in 11 segments
  colnet <- brewer.pal(11,"RdYlGn")   # Assign colors to trophic levels
  V(g)$color <- colnet[12-colTL]      # Add colors to the igraph object
  
  
  plot(g, edge.width=.3,edge.arrow.size=.4,
       #vertex.label=NA,
       vertex.label.color="white",
       edge.color="grey50",
       edge.curved=0.3, layout=lMat)
  
  
}



#' Calculate topological indices for ecological networks 
#'
#' @param ig igraph object
#'
#' @return a data.frame with the following fields:
#'         Size: Number of species
#'         Top:  Number of top predator species
#'         Basal: Number of basal especies 
#'         Links: number of interactions
#'         LD: linkage density
#'         Connectance: Connectance
#'         PathLength: average path length
#'         Clustering: clustering coeficient
#'         Cannib: number of cannibalistic species
#'          
#' @export
#'
#' @examples
topologicalIndicesEcoNetwork <- function(ig){
  
  df <- lapply(ig, function(g){
  
    deg <- degree(simplify(g), mode="out") # calculate the out-degree: the number of predators  
    
    V(g)$outdegree <-  deg
    
    nTop <- length(V(g)[outdegree==0]) # Top predators do not have predators
    
    deg <- degree(g, mode="in") # calculate the in-degree: the number of preys
    
    V(g)$indegree <-  deg
    
    nBasal <- length(V(g)[indegree==0]) # Basal species do not have preys 
    
    vcount(g)-nTop-nBasal
    
    size <- vcount(g)
    
    links <- ecount(g)
    
    linkDen <- links/size          # Linkage density
    
    conn <- links/size^2           # Connectance
    
    pathLength <- average.path.length(g)   # Average path length 
    
    clusCoef <- transitivity(g, type = "global") 
    
    cannib <- sum(which_loop(g))
  
    data.frame(Size=size,Top=nTop,Basal=nBasal,Links=links, LD=linkDen,Connectance=conn,PathLength=pathLength,Clustering=clusCoef, Cannib=cannib)
  })

  do.call(rbind,df)
}
  

parTopologicalIndicesEcoNetwork <- function(ig){
  
  require(doParallel)
  cn <-detectCores()
  cl <- makeCluster(cn)
  registerDoParallel(cl)
  
  df <- foreach(i=seq_along(ig), .combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% {
    
    g <- ig[[i]]
    
    deg <- degree(simplify(g), mode="out") # calculate the out-degree: the number of predators  
    
    V(g)$outdegree <-  deg
    
    nTop <- length(V(g)[outdegree==0]) # Top predators do not have predators
    
    deg <- degree(g, mode="in") # calculate the in-degree: the number of preys
    
    V(g)$indegree <-  deg
    
    nBasal <- length(V(g)[indegree==0]) # Basal species do not have preys 
    
    vcount(g)-nTop-nBasal
    
    size <- vcount(g)
    
    links <- ecount(g)
    
    linkDen <- links/size          # Linkage density
    
    conn <- links/size^2           # Connectance
    
    pathLength <- average.path.length(g)   # Average path length 
    
    clusCoef <- transitivity(g, type = "global") 
    
    cannib <- sum(which_loop(g))
    
    data.frame(Size=size,Top=nTop,Basal=nBasal,Links=links, LD=linkDen,Connectance=conn,PathLength=pathLength,Clustering=clusCoef, Cannib=cannib)
  }
  stopCluster(cl)
  
  return(df)
}
