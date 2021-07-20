#' Runs a evaluation on the E. coli genome, returns a list of accuracy value for each letter in the vocabulary
#' @param maxlen time steps to unroll for (e.g. length of semi-redundant chunks)
#' @param vocabulary char contains the vocabulary from the input char, it should be sorted.
#' @param model used model
#' @param genome used genome for the evaluation
#' @export
ecoliEvaluation <- function(model,
                            genome = "ecoli_small",
                            maxlen = 80,
                            vocabulary = c("-", "|", "a", "c", "g", "t")) {
  # load amino acid sequence 
  if (genome == "ecoli") {
    data(ecoli)
    assign("genome", ecoli)
  }
    
  if (genome == "ecoli_small") {
    data(ecoli_small)
    assign("genome", ecoli_small)
  }
  
  # get semi-redundant encoding 
  genome.encoded <- preprocessSemiRedundant(genome,
                                            vocabulary = vocabulary,
                                            maxlen = maxlen)
  ground.truth <- apply(genome.encoded$Y, 1, which.max)
  predictions <- keras::predict_classes(model, genome.encoded$X)
  
  # iterates over all char in vocabulary and calculates number of correct predictions for each character
  char.acc <- list()
  for (char in vocabulary) {
    message(paste("evaluate character",  char))
    ground.truth.char <- which(ground.truth ==  which(vocabulary == char))
    predicted.char <- predictions[ground.truth.char] # these should be predicted as "a"
    correct <- length(which(predicted.char  == which(vocabulary == char)))
    char.acc[[char]] <- as.integer(correct) / length(predicted.char) 
  } 
  return(char.acc)
}

#' Custom lambda callback that shows accuacy scores for the E. coli genome in Tensorboard
ecoliCustomScalar <-
  keras::callback_lambda(
    on_epoch_end = function(epoch, logs) {
      char.acc <- EcoliEvaluation(model)
      tensorflow::tf$summary$scalar(paste("Ecoli", vocabulary[1]), data = char.acc[[1]], step = epoch)
      tensorflow::tf$summary$scalar(paste("Ecoli", vocabulary[2]), data = char.acc[[2]], step = epoch)
      tensorflow::tf$summary$scalar(paste("Ecoli", vocabulary[3]), data = char.acc[[3]], step = epoch)
      tensorflow::tf$summary$scalar(paste("Ecoli", vocabulary[4]), data = char.acc[[4]], step = epoch)
      tensorflow::tf$summary$scalar(paste("Ecoli", vocabulary[5]), data = char.acc[[5]], step = epoch)
      tensorflow::tf$summary$scalar(paste("Ecoli", vocabulary[6]), data = char.acc[[6]], step = epoch)
    }
  )
