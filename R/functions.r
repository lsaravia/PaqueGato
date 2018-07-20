# Function definitions


#' Crea un dibujo de un gato en ascii
#'
#' @param numeroDeGatos el numero de gatos a imprimir
#'
#' @return none
#' @export
#'
#' @examples
#'
#' crearGato(10)
#'
crearGato <- function(numeroDeGatos){

for(i in 1:numeroDeGatos){
  cat("  /\\___/\\
  |     |
_  *   *  _
-   /_\   -
    ---")
  }

}




#' Calculates the Shannon Diversity Index using ln base e
#'
#' @param species vector of abundances of species
#'
#' @return the Shannon Diversity Index
#'
#' @export
#'
#' @examples
#'
#' shannon(datosPruebaShannon)
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
#
#' @examples
#'
#' filename <- system.file("extdata", "Clorofila.txt", package = "PaqueGato")
#'
#' chla <- read.delim(filename)
#'
#' fixClorophylData(filename)
#'
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

    fecha <- lubridate::ymd(paste(chla$Year[i],chla$Month[i], 1))
    if( is.na(fecha)){
      if(chla$Month[i]=="Mar")
        fecha <- lubridate::ymd(paste(chla$Year[i],3, 1))
      else {
        fecha <- lubridate::dmy(chla$Month[i])

        if(is.na(fecha)) {
          fecha <- lubridate::mdy(chla$Month[i])

          if(is.na(fecha)) {
            temp <- strsplit(chla$Month[i]," ")[[1]][3]
            fecha <- lubridate::dmy(temp)
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

