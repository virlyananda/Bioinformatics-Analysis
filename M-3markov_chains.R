#!/usr/bin/env Rscript

# Open libraries that could help the program run successfully.
#install.packages("seqinr")
#install.packages("dplyr")
#install.packages("markovchain")
library(plyr)
library(seqinr)
library(ggplot2)
library(gridExtra)
# Create a function to plot the first order of Markov Data
plotData<-function(seq, nucleotides, title){
  # COUNT the frequency of each base (Nucleotide)
  freqNucleotide<-seqinr::count(seq,1,alphabet=nucleotides,freq=TRUE)
  # COUNT the frequency of dinucleotide
  freqDiNucleotide<-seqinr::count(seq,2,alphabet=nucleotides,freq=TRUE)
  # COUNT the frequency of trinucleotide
  freqTriNucleotide<-seqinr::count(seq,3,alphabet=nucleotides,freq=TRUE)
  # Make a data frame to draw a ggplot2 (convert each base to dataframe)
  df<-as.data.frame(freqNucleotide)
  colnames(df)<-c("Base", "Base_Proportion")
  p1<-ggplot(df, aes(x=Base, y=Base_Proportion, fill=Base)) + geom_bar(stat="identity")
  p1<-p1 + theme(legend.position="none") + ggtitle("Compositional bias of each nucleotide")
        
  df<-as.data.frame(freqDiNucleotide)
  colnames(df)<-c("Base", "Base_Proportion")
  p2<-ggplot(df, aes(x=Base, y=Base_Proportion, fill=Base)) + geom_bar(stat="identity")
  p2<-p2 + theme(legend.position="none") + ggtitle("Compositional bias of each dinucleotide")
            
  df<-as.data.frame(freqTriNucleotide)
  colnames(df)<-c("Base", "Base_Proportion")
  p3<-ggplot(df, aes(x=Base, y=Base_Proportion, fill=Base)) + geom_bar(stat="identity")
  p3<-p3 + theme(legend.position="none") + ggtitle("Compositional bias of each trinucleotide") + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=8))
  # Plot the result in 1 plot
  grid.arrange(p1, p2, p3, nrow=3, top=title)
}

# Function to generate the sequence
generateFirstOrderSeq<-function(lengthSeq, nucleotides, initialProb, firstOrderMatrix){
  # Vector for storing new sequence
  outputSeq<-character()
  firstnucleotide<-sample(nucleotides, 1, rep=TRUE, prob=initialProb)
  # Store nucleotide in the first position
  outputSeq[1]<-firstnucleotide
    
  # Create a for loop here
  for(i in 2:lengthSeq){
    prevNuc<-outputSeq[i-1]
    currentProb<-firstOrderMatrix[prevNuc,]
    #cat("\n", prevNuc, " ", currentProb, sep=" ")
    outputSeq[i]<-sample(nucleotides, 1, prob=currentProb)
  }
  return(outputSeq)
}
# Function to generate DNA sequence from the given HMM and sequence length
generateFirstOrderhmmseq<-function(lengthSeq, nucleotides, initialProb, states, transitionmatrix, emissionmatrix){
  outputSeq<-character() # Vector for storing new sequence
  mystates<-character() # Vector for storing the state of each position in the new sequence
  # State for the first position in the sequence:
  firststate<-sample(states, 1, rep=TRUE, prob=initialProb)
  # Probabilities of the current nucleotide in the "firststate"
  probabilities<-emissionmatrix[firststate,]
  # Nucleotide for the first position:
  firstnucleotide<-sample(nucleotides, 1, rep=TRUE, prob=probabilities)
  outputSeq[1]<-firstnucleotide # Store nucleotide in the first position
  mystates[1]<-firststate # Store state in the first position
        
  for(i in 2:lengthSeq){
    prevstate<-mystates[i-1] # State from previous nucleotide
    stateprobs<-transitionmatrix[prevstate,]
    # Ith position in the sequence:
    state<-sample(states, 1, rep=TRUE, prob=stateprobs)
    # Get the probabilities for the current nucleotide in the state "state":
    probabilities<-emissionmatrix[state,]
    # Nucleotides for ith position:
    nucleotide<-sample(nucleotides, 1, rep=TRUE, prob=probabilities)
    outputSeq[i]<-nucleotide # Store nucleotide for the current position in the sequence
    mystates[i]<-state # Store state for the current position in the sequence
  }
          
  for(i in 1:lengthSeq){
    nucleotide<-outputSeq[i]
    state<-mystates[i]
    print(paste("Position", i, ", State", state, ", Nucleotide = ", nucleotide))
  }
  return(outputSeq)
}

viterbi<-function(sequence, transitionmatrix, emissionmatrix){
  states<-rownames(emissionmatrix) # Get the names of the states in HMM
  v<-makeViterbimat(sequence, transitionmatrix, emissionmatrix) # Make the Viterbi matrix
  mostprobablestatepath<-apply(v, 1, function(x) which.max(x))
    
  # Print the most probable state path:
  prevnucleotide<-sequence[1]
  prevmostprobablestate<-mostprobablestatepath[1]
  prevmostprobablestatename<-states[prevmostprobablestate]
  startpos<-1
  for (i in 2:length(sequence)){
    nucleotide<-sequence[i]
    mostprobablestate<-mostprobablestatepath[i]
    mostprobablestatename<-states[mostprobablestate]
    if (mostprobablestatename != prevmostprobablestatename){
      print(paste("Positions", startpos, "-",(i-1),
                  "Most probable state = ", prevmostprobablestatename))
      startpos<-i
    }
    prevnucleotide<-nucleotide
    prevmostprobablestatename<-mostprobablestatename
  }
  print(paste("Positions", startpos, "-",i,
              "Most probable state = ", prevmostprobablestatename))
}
                                                                
makeViterbimat<-function(sequence, transitionmatrix, emissionmatrix){
  # Change sequence to uppercase
  sequence<-toupper(sequence)
  # Numbers of states in HMM
  numstates<-dim(transitionmatrix)[1]
  # Matrix with the same rows and columns as positions in the sequence and HMM:
  v<-matrix(NA, nrow=length(sequence), ncol=dim(transitionmatrix)[1])
  # Set the values in the first row inside matrix
  v[1, ]<- 0
  v[1,1]<- 1
  # Fill the matrix v:
  for (i in 2:length(sequence)){
    for (l in 1:numstates){
      statelprobnucleotidei<-emissionmatrix[l,sequence[i]]
      v[i,l]<-statelprobnucleotidei*max(v[(i-1),]*transitionmatrix[,l])
    }
  }
  return(v)
}
pdf("markov_plots.pdf")

# Define DNA alphabet to put names to objects:
nucleotides<-c("A","C","G","T")
# Vector representing probability distribution of the model
zeroOrderProbabilities<-c(0.2,0.3,0.3,0.2)
# Reference name for each base
names(zeroOrderProbabilities)<-nucleotides
# Sequence with 1000 bases for this model
zeroOrderSeq<-sample(nucleotides, 1000, rep=T, prob=zeroOrderProbabilities)
plotData(zeroOrderSeq, nucleotides, "Multinomial Model of DNA Evolution")

# Add the probability distribution per base and set the values for each
afterAprobs<-c(0.2,0.3,0.3,0.2)
afterCprobs<-c(0.1,0.41,0.39,0.1)
afterGprobs<-c(0.25,0.25,0.25,0.25)
afterTprobs<-c(0.5,0.17,0.17,0.17)
# Create 4x4 matrix that can store the probability distribution above
mytransitionmatrix<-matrix(c(afterAprobs,afterCprobs,afterGprobs,afterTprobs), 4,4, byrow = TRUE)
colnames(mytransitionmatrix)<-nucleotides
rownames(mytransitionmatrix)<-nucleotides

inProb<-c(0.4,0.1,0.1,0.4)
names(inProb)<-nucleotides

firstOrderSeq<-generateFirstOrderSeq(1000,nucleotides,inProb,mytransitionmatrix)

plotData(firstOrderSeq, nucleotides, "Markov Chain of first order")

# AT-Rich and GC-Rich State
# Define the names of the states
states<-c("AT-rich","GC-rich")
# Set probabilities of switching states, where the previous state was "AT-Rich"
ATrichprobs<-c(0.7,0.3)
# Set probabilities of switching states, where the previous state was "GC-Rich"
GCrichprobs<-c(0.1,0.9)
# Create a 2x2 Matrix
theTransitionMatrix<-matrix(c(ATrichprobs,GCrichprobs),2,2,byrow = TRUE)
rownames(theTransitionMatrix)<-states
colnames(theTransitionMatrix)<-states
# Set the probability values for AT rich state
ATrichstateprobs<-c(0.39,0.1,0.1,0.41)
GCrichstateprobs<-c(0.1,0.41,0.39,0.1)
# Create 2x4 Matrix
theEmissionMatrix<-matrix(c(ATrichstateprobs,GCrichstateprobs),2,4,byrow = TRUE)
rownames(theEmissionMatrix)<-states
colnames(theEmissionMatrix)<-nucleotides

initialProb<-c(0.5,0.5)
hmmfirstOrderSeq<-generateFirstOrderhmmseq(1000, nucleotides, initialProb, states, theTransitionMatrix, theEmissionMatrix)
plotData(hmmfirstOrderSeq, nucleotides, "Hidden Markov Model of first order")

myseq<-c("A", "A", "G", "C", "G", "T", "G", "G", "G", "G", "C", "C", "C", "C","G", "G", "C", "G", "A", "C", "A", "T", "G", "G", "G", "G", "T", "G","T", "C")
viterbi(myseq, theTransitionMatrix, theEmissionMatrix)

dev.off()

viterbi<-function(sequence, transitionmatrix, emissionmatrix){
  states<-rownames(emissionmatrix)
  # Make the Viterbi matrix v:
  v<-makeViterbimat(sequence, transitionmatrix, emissionmatrix)
  mostprobablestatepath<-apply(v, 1, function(x) which.max(x))

  # Print out the most probable state path:
  prevnucleotide<-sequence[1]
  prevmostprobablestate<-mostprobablestatepath[1]
  prevmostprobablestatename<-states[prevmostprobablestate]
  startpos<-1
  for (i in 2:length(sequence)){
      nucleotide<-sequence[i]
      mostprobablestate<-mostprobablestatepath[]
      mostprobablestatename<-states[mostprobablestate]
      if (mostprobablestatename != prevmostprobablestatename){
          print(paste("Positions",startpos,"-",(i-1),
                      "Most probable state = ", prevmostprobablestatename))
          startpos<-i
      }
      prevnucleotide<-nucleotide
      prevmostprobablestatename<-mostprobablestatename
  }
  print(paste("Positions",startpos,"-",i,
              "Most probable state = ", prevmostprobablestatename))
}
# Second function
makeViterbimat<-function(sequence, transitionmatrix, emissionmatrix){
  # Change the sequence to uppercase
  sequence<-toupper(sequence)
  # Find out how many states are in the HMM
  numstates<-dim(transitionmatrix)[1]
  # Columns as states in the HMM
  v<-matrix(NA, nrow = length(sequence), ncol=dim(transitionmatrix)[1])
  # Set the values in the first row of matrix v, representing the first position of the sequence to 0
  v[1, ]<-0
  v[1,1]<-1
  # Fill in the matrix v:
  for (i in 2:length(sequence)){
    for (l in 1:numstates){
      stateprobnucleotidei<-emissionmatrix[l,sequence[i]]
      v[i,l]<-stateprobnucleotidei*max(v[(i-1),]*transitionmatrix[,l])
    }
  }
  return(v)
}

myseq<-c("A", "A", "G", "C", "G", "T", "G", "G", "G", "G", "C", "C",
         "C", "C", "G", "G", "C", "G", "A", "C", "A", "T", "G", "G",
         "G", "G", "T", "G", "T", "C")
viterbi(myseq, theTransitionMatrix, theEmissionMatrix)
