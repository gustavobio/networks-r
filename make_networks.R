make.networks <- function(plantas, aves, dens = 0.2, runs = 1000) {
  networks <- vector("list", runs)
  ## Aqui assumimos que as aves consomem frutos "aleatoriamente" desde que o fruto
  ## seja menor que o tamanho do bico. Algo pra se pensar é que como o que temos são as
  ## médias, seria legal tentar assimilar a variação pelo menos nos frutos pra considerar
  ## essas relações.
  
  ## Quais links são possíveis dadas as características morfológicas das espécies?
  links.morpho <- outer(plantas$fruit, aves$gape, `<`)
  
  
  ## Estrato
  ## links.strata <- outer(plantas$strata, aves$strata, function(x, y) {
  ##  any(unlist(strsplit(x, " and |, ")), unlist(strsplit(y, " and |, ")))
  ## })
  
  links.strata <- matrix(F, nrow = length(plantas$strata), ncol = length(aves$strata))
  for (i in 1:length(plantas$strata)) {
    for (j in 1:length(aves$strata)) {
      links.strata[i, j] <- any(unlist(strsplit(plantas$strata[i], " and |, ")) %in% unlist(strsplit(aves$strata[j], " and |, ")))
    }
  }
  
  ## Ambos
  links <- ifelse(links.morpho & links.strata, 1, 0)
  
  ## Tirar as espécies sem ligações
  
  colsums <- colSums(links)
  rowsums <- rowSums(links)
  
  nomes.plantas <- tolower(unique(plantas$species)[rowsums > 0])
  # Pedir pra Jamille arrumar o nome dessa coluna
  nomes.aves <- tolower(unique(aves$specie)[colsums > 0])
  
  links <- links[rowsums > 0, colsums > 0]
  
  N.plantas <- nrow(links)
  N.aves <- ncol(links)
  
  ## Checando a densidade
  dens.simul <- sum(links)/(N.aves * N.plantas)
  
  ## Como muito provavelmente a densidade ainda é maior que a esperada, começamos
  ## a retirar alguns links baseados na probabilidade deles acontecerem dadas as
  ## abundâncias até chegarmos na densidade desejada. Aqui não estou me preocupando
  ## em fazer com que a distribuição do número de arestas de cada nível seja lognormal
  ## pq não consideramos as frequências das interações.
  
  ## Qual o número mínimo de interações de uma espécie?
  min.links <- 1
  ## Abundâncias relativas * abundâncias relativas
  mat.prob <- outer(plantas$abundance, aves$abundance, function(x, y) x/sum(plantas$abundance) * y/sum(aves$abundance))
  mat.prob <- mat.prob[rowsums > 0, colsums > 0]
  
  colsums <- colSums(links)
  rowsums <- rowSums(links)
  
  c.keep <- colsums <= min.links & colsums != 0
  r.keep <- rowsums <= min.links & rowsums != 0
  
  c.remove.n <- sum(colsums[c.keep])
  r.remove.n <- sum(rowsums[r.keep])
  ## Quantos links falta tirar?
  # links.out <- sum(links) - length(links) * dens
  # prob <- quantile(mat.prob, links.out)
  # final.mat <- ifelse(mat.prob * links >= prob, 1, 0)
  mat.prob[r.keep, c.keep] <- mat.prob[r.keep, c.keep] + links[r.keep, c.keep]
  
  for (i in seq_len(runs)) {
    final.mat <- matrix(0, ncol = N.aves, nrow = N.plantas)
    links.keep <- sample(which(links == 1), 
                         length(links) * dens,
                         prob = mat.prob[links == 1])
    final.mat[links.keep] <- 1
    (dens.simul <- sum(final.mat)/(N.aves * N.plantas))
    final.mat <- t(final.mat)
    colnames(final.mat) <- nomes.plantas
    rownames(final.mat) <- nomes.aves
    
    
    # visweb(t(mat.bin))
    # plotweb(t(mat.bin))
    # plotweb(final.mat, text.rot = 90)
    
    # networklevel(final.mat, index = "connectance")
    networks[[i]] <- final.mat
  }
  invisible(networks)
}