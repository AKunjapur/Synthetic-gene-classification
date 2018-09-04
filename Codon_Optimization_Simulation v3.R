### Synthesis Cost Project - Codon Optimization Simulation
### Code by Neil Thompson
### Data from Aditya Kunjapur
### June 6, 2017


### Functions

Populate.Organism.Data       <- function(directory){
  ### Creates data frames for each organism, based on the data on
  ### amino acid and codon frequency from Aditya
  
  files = list.files(directory)
  
  organism.data = list()
  
  for(file in files){
    organism <- strsplit(file, split=".csv")[[1]]
    
    # import data
    organism.table <- read.csv(
                                file = paste(directory,file,sep=""),
                                header = FALSE,
                                blank.lines.skip = TRUE)
    
    # names fields
    colnames(organism.table)  <- c("codon",
                                 "amino.acid",
                                 "freq.for.amino.acid",
                                 "freq.in.all.codons",
                                 "raw.count")
    # delete empty rows
    organism.table <- organism.table[organism.table$codon != "" ,]
    
    # keep organism name with table
    attributes(organism.table)$organism.name <- organism
    
    organism.data[[organism]] <- organism.table
                                
  }
  
  return(organism.data)
}
  
Populate.Amino.Acid.And.Codon.Lists     <- function(organism.data){
  ### Fills in reference tables for the amino acids and codons
  ### for each organism, based on the raw data from Aditya
  
  amino.acid.and.codon.data <- list(amino.acid.list = list(),
                                    codon.list      = list())
    
  
  # get list of all amino acids and codons used across all organisms
  all.amino.acids = levels(unlist(lapply(organism.data, function(x) unique(x$amino.acid))))
  all.amino.acids = all.amino.acids[all.amino.acids != ""]
  
  all.codons      = levels(unlist(lapply(organism.data, function(x) unique(x$codon))))
  all.codons      = all.codons[all.codons != ""]
  
  ## build amino acid and codon list
  for(organism in names(organism.data)){
    
    ## Build data frame for amino acids
    amino.acids <- aggregate(x   = organism.data[[organism]]$freq.in.all.codons,
                             by  = list(organism.data[[organism]]$amino.acid),
                             FUN = "sum")
    
    colnames(amino.acids) <- c("name", "prob")
    
    missing.amino.acids   <- setdiff(all.amino.acids, amino.acids[[1]])
    
    if(!identical(missing.amino.acids,character(0))){
      amino.acids <- rbind(amino.acids,
                         data.frame(name = missing.amino.acids, prob = 0))
    }
    
    # normalize to address rounding issues
    amino.acids$prob <- amino.acids$prob * 1/sum(amino.acids$prob)
    
    
    # add tag to data frame
    attributes(amino.acids)$organism.name = organism
    
    # save
    amino.acid.and.codon.data$amino.acid.list[[organism]] <- amino.acids
    
    
    ## Build list for codons
    codons <- organism.data[[organism]][c("codon", "amino.acid", "freq.in.all.codons")]
    colnames(codons) <- c("codon", "name", "prob")
    
    if(!identical(missing.amino.acids,character(0))){
      codons <- rbind(codons,
                      data.frame(codon = "???", name = missing.amino.acids, prob = 1))
    }
    
    # build one data.frame for each amino acid
    codons <- split(codons, codons$name)
    codons <- codons[names(codons) != ""]
    
    # normalize to address rounding issues
    for(aa in names(codons)){
      codons[[aa]]$prob = codons[[aa]]$prob * 1/sum(codons[[aa]]$prob)
    }
    
    # add tag to data frame
    attributes(codons)$organism.name = organism
    
    # save
    amino.acid.and.codon.data$codon.list[[organism]] <- codons
    
  }
  
  return(amino.acid.and.codon.data)  
}




Generate.Amino.Acid.Sequence <- function(amino.acids, sequence.length){
  # Generates a random amino acid sequence from the frequencies given
  
  amino.acid.sequence <- sample(x       = amino.acids$name,
                                size    = sequence.length,
                                replace = TRUE,
                                prob    = amino.acids$prob)
  return(amino.acid.sequence)
}



Generate.Codon.Sequence      <- function(amino.acid.seq, codon.table){
  # Generates a random codon sequence based on an amino acid sequence based on the codon frequency of the organism
  
  # sequence
  codon.seq = rep(NA, length(amino.acid.seq))
  
  # do weighted sampling
  for(j in 1:length(amino.acid.seq)){
    amino.acid        = as.character(amino.acid.seq[j])
    new.codon         = sample(x       = as.character(codon.table[[amino.acid]]$codon),
                               size    = 1,
                               replace = TRUE,
                               prob    = as.character(codon.table[[amino.acid]]$prob))
    codon.seq[j] = new.codon
  }
  return(codon.seq)
}


Compare.Codon.Sequences      <- function(codon.seq.1, codon.seq.2){
  # Calculates the similarity of two sequences at the BASE level
  seq.1 = unlist(strsplit(codon.seq.1, split = ""))
  seq.2 = unlist(strsplit(codon.seq.2, split = ""))
  
  codon.similarity <- sum(seq.1 == seq.2)/length(seq.1)
  
  return(codon.similarity)
}


### Simulation

# Parameters and System Settings
directory       = "C:/Users/neil_t/Dropbox\ (MIT)/Synthesis\ Costs/Code/Simulation\ for\ sequence\ identity/Codon\ Usage\ CSV/"
num.repetitions = 1000
sequence.length = 1000

# Populate data on organisms
organism.data   <- Populate.Organism.Data(directory)
organisms       <- names(organism.data)
#organisms       = c("human", "e.coli", "mouse")

# construct list
amino.acid.and.codon.data <- Populate.Amino.Acid.And.Codon.Lists(organism.data)
amino.acid.list           <- amino.acid.and.codon.data[[1]]
codon.list                <- amino.acid.and.codon.data[[2]]

# Data structures
mean.similarity     = matrix(nrow     = length(organisms), 
                             ncol     = length(organisms),
                             dimnames = list(organisms, organisms))
stdev.similarity     = matrix(nrow    = length(organisms), 
                             ncol     = length(organisms),
                             dimnames = list(organisms, organisms))
all.outcomes        = list()


for(source.organism in organisms){
  for(expression.organism in organisms){
    
    simulation.outcomes   <- rep(NA, num.repetitions)
    
    for(i in 1:num.repetitions){
      
      ## Generate Sequences for comparison
      
      # Generate amino acid sequence from the source organism
      amino.acid.seq        <- Generate.Amino.Acid.Sequence(
                                  amino.acids     = amino.acid.list[[source.organism]], 
                                  sequence.length = sequence.length)
      
      # Generate a codon sequence for the amino acid sequence consistent with the source organism
      source.codon.seq      <- Generate.Codon.Sequence(
                                  amino.acid.seq  = amino.acid.seq, 
                                  codon.table     = codon.list[[source.organism]])
      
      # Generate a codon sequence for the amino acid sequence consistent with the expression organism
      expression.codon.seq <- Generate.Codon.Sequence(
                                  amino.acid.seq  = amino.acid.seq, 
                                  codon.table     = codon.list[[expression.organism]])
      
      ## Calculate the similarity between the codon sequences and save
      simulation.outcomes[i] <- Compare.Codon.Sequences(source.codon.seq, expression.codon.seq)
      
    }
    
    # Save Simulation Outcomes
    all.outcomes[[paste(source.organism, expression.organism, sep=" : ")]] <- simulation.outcomes
    mean.similarity[source.organism, expression.organism] <- mean(simulation.outcomes)
    stdev.similarity[source.organism, expression.organism] <- sd(simulation.outcomes)
    
  }
}

## Boxplot of similarity for expression in humans (plots)

# create data frame
source.organisms <- c()
identity.data <- c()

expression.organism = "H sapiens"

for(source.organism in organisms){
  source.organisms = append(source.organisms, rep(source.organism, length(all.outcomes[[paste(source.organism, expression.organism, sep=" : ")]])))
  identity.data = append(identity.data,all.outcomes[[paste(source.organism, expression.organism, sep=" : ")]])
}

boxplot.data2 <- data.frame(source.organisms, identity.data)

par(mar=c(10,5,1,1))
boxplot(boxplot.data2$identity.data ~ boxplot.data2$source.organisms, 
        col  = "grey", 
        ylim = c(0.6,0.9),
        ylab = "% Identity",
        yaxt = "n",
        las  = 2)
axis(2, at=pretty(boxplot.data2$identity.data), lab=paste0(pretty(boxplot.data2$identity.data)*100,"%"), las=TRUE, ylim = c(0.6,0.9))
abline(h=0.85, lty=3)

# Density plot for each expression organism by source organism
densityplot.data            <- data.frame(matrix(nrow=length(organisms), ncol=100, data=0))
row.names(densityplot.data) <- organisms
names(densityplot.data)     <- as.character((1:100)/100)

for(expression.organism in organisms){
  expression.data <- c()
  
  # append data from all source organisms
  for(source.organism in organisms){
    expression.data <- append(expression.data, all.outcomes[[paste(source.organism, expression.organism, sep=" : ")]])
  }
  
  # round data
  rounded.data = round(expression.data, digits=2)
  tab.data = as.data.frame(table(rounded.data), stringsAsFactors = F)
  densityplot.data[expression.organism, tab.data$rounded.data] <- tab.data$Freq
}

write.csv(densityplot.data, file = "C:/Users/neil_t/Dropbox\ (MIT)/Synthesis\ Costs/Code/MonteCarloDensities.csv")
  
save.image(paste("C:/Users/neil_t/Dropbox\ (MIT)/Synthesis\ Costs/Code/Codon\ Optimization\ Simulation\ Results", sequence.length,".RData"))
