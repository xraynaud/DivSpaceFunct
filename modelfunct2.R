BPclip = function(intensity, owin, radii,...) {
  
  newW = affine(owin,matrix(c(1.5,0,0,1.5),ncol=2, nrow=2), c(-diff(owin$x),-diff(owin$y)))
  #  if (as.mask) {
  #    dimyx = c(diff(newW$y),diff(newW$x) )
  #  }
  X = rpoispp(intensity, win=newW)
  D = discs(X, radii=radii,...)
  return(D[owin])
}

creatematrixmask = function(x, y, array, R) {
  for (i in (-round(R-1):round(R-1))) {
    for (j in (-round(R-1):round(R-1))) {
      if (i^2 + j^2 <= round(R)^2) {
        cx = x+i
        cy = y+j
        if (cx <=0) {
          cx = cx+dim(array)[1]
        } else {
          if( cx > dim(array)[1]) {
            cx = cx-dim(array)[1]
          }
        }
        if (cy <=0) {
          cy = cy+dim(array)[2]
        } else {
          if( cy > dim(array)[2]) {
            cy = cy-dim(array)[2]
          }
        }
        array[cx,cy] = 1
        }
      }
    }
  return(array)
}
simulate = function(arena, cells,supply, rates,substratemat, productmat,  R, time, tout, showplot, verbose,torus=F) {
  arenasize = dim(arena)[1]
  nsubstrate = dim(arena)[3]
  cells_orig = cells
  
  if (torus) {
  cells = superimpose(cells,shift(cells,c(arenasize,0)),shift(cells,c(-arenasize,0)),shift(cells,c(0,arenasize)),shift(cells,c(0,-arenasize)),shift(cells,c(arenasize,arenasize)),shift(cells,c(-arenasize,arenasize)),shift(cells,c(arenasize,-arenasize)),shift(cells,c(-arenasize,-arenasize)))
  cells = cells[owin(c(0,arenasize+2*max(marks(cells)$R)),c(0,arenasize+2*max(marks(cells)$R)))]
  }
  
  zeromat =  array(0,dim=c(arenasize,arenasize))

  nspecies = length(unique(marks(cells)$sp))

  # For each individual in the model, make of binary mask of the same size of the arena containing 1 in cells that can be reach by cells, 0 otherwise.
  indmasks = lapply(1:cells$n, function(i) {creatematrixmask(round(cells[i]$x),round(cells[i]$y), matrix(0, nrow = arenasize, ncol= arenasize), marks(cells[i])$R)})
  
  # make of summing masks for all individuals of a species. This should allow to calculate substrates acquired when space is shared.  
  speciesmasks = lapply(1:nspecies,function(s) Reduce("+",lapply(as.list(which(marks(cells)$sp==s)),function(i) indmasks[[i]])))
  
  # define which matrix is the one to read(first arg), and which is the one two write (2nd arg) (see 4th arg of arena)
  switch = c(1,2)
  
  uptake= lapply(1:nspecies, function(x) rep(0,length(time)))
  if (verbose) {  
  #Some text output rto indicate give initial starting values
    print(paste("Quantity of matter in model", paste( apply(arena[,,,1],3,sum),collapse=" ")))
    print(paste("Total",sum(arena[,,,1])))
  }

  # time loop
  for (t in time) {
    if (t%%tout==0 & verbose) {
      print("")
      print(paste("Quantity of matter in model", paste( apply(arena[,,,switch[1]],3,sum),collapse=" ")))
      print(paste("Total",sum(arena[,,,switch[1]])))
    }
    # update progress bar
    
    produce =  list()
    demand = list()
    upt = list()
    ret = list()
    for (i in 1:nsubstrate) {
      produce[[i]] = zeromat
      demand[[i]] = zeromat
    }
    for (i in 1:cells$n) {
      upt[[i]] = lapply(
        1:nsubstrate,
        function(x) 
          if (x %in% which(substratemat[,marks(cells[i])$sp]!=0)) {
            substratemat[x,marks(cells[i])$sp] * rates[marks(cells[i])$sp] * indmasks[[i]] *arena[,,x,switch[1]]
          } else{
            zeromat
          }
      )
      for (s in 1:nsubstrate) {
        demand[[s]] = demand[[s]] + upt[[i]][[s]]# + ret[[i]] [[s]]
      }
    }
    for (s in 1:nsubstrate) {
      avail = ifelse(arena[,,s,switch[1]]-demand[[s]]>0,demand[[s]],arena[,,s,switch[1]])
      for (c in which(marks(cells)$sp %in% which(substratemat[s,] !=0))) {
        upt[[c]][[s]] = upt[[c]][[s]] * ifelse(demand[[s]] !=0, avail/demand[[s]], 0)
      }
      demand[[s]] = avail
    }
    
    for (c in 1:cells$n) {
      uptake[[marks(cells[c])$sp]][t] = uptake[[marks(cells[c])$sp]][t] + sum(Reduce("+",upt[[c]]))
      ret[[c]]  = lapply(as.list(1:nsubstrate),
        function(x)
          if (x %in% which(productmat[,marks(cells[c])$sp]!=0)) {
            zeromat + 
              (productmat[x,marks(cells[c])$sp]/apply(productmat,2,sum)[marks(cells[c])$sp])*indmasks[[c]]*sum(Reduce("+",upt[[c]]))/sum(indmasks[[c]])
          } else {
            zeromat
          }
      )
      if (sum(is.na(ret[[c]][[s]]))>0) {
        print(paste(c, sum(Reduce("+",upt[[c]])),sum(indmasks[[c]])))
        stop()
        }
      
      for (s in 1:nsubstrate) {
        produce[[s]] = produce[[s]] + ret[[c]][[s]]
      }
    }
    for (s in 1:nsubstrate) {
      arena[,,s,switch[2]] = arena[,,s,switch[1]] - demand[[s]] + produce[[s]]
    }
    arena[,,1,switch[2]] = arena[,,1,switch[2]] + supply
    # #show intermediate plot ? 
    # if (showplot) {
    #   g = list()
    #   for (s in 1:nsubstrate) {
    #     g[[s]] = ggplot(data=cbind(expand.grid(x=seq(1,dim(arena)[1],length.out=dim(arena)[1]),y=seq(1,dim(arena)[2],length.out=dim(arena)[2])),q=as.vector(arena[,,s,1])),aes(x=x,y=y,fill=q)) +
    #       geom_raster(show.legend=F) + 
    #       scale_fill_viridis_c(breaks=seq(0,1,length.out=32)) + 
    #       theme_bw() + 
    #       theme(aspect.ratio = 1,axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank(),axis.ticks = element_blank())
    #     #  image(seq(1,dim(arena)[1],length.out=dim(arena)[1]),seq(1,dim(arena)[2],length.out=dim(arena)[2]), arena[,,s,switch[1]],asp=1, col = heat.colors(64), breaks = seq(0,1, length.out=66)[2:66],xlab="",ylab="",xaxt="n",yaxt="n")
    #     #  points(cells[marks(cells)==s])
    #   }
    # }
    #rotate matrix switch (see arg 4 of arena)
    switch =  rev(switch)
  }
  return(list(arena = arena[,,,switch[1]], cells = cells, uptake = uptake, speciesmasks=speciesmasks, uptmat = upt))
}