#' An S4 class for the prediction
#'
#'@slot next_char Next predicted char
#'@slot probability probability for the predicted char
#'@slot index index
#'@slot alternative_probability all probabilities for the possible chars
#'@slot solution sequence with the predicted char 
prediction <- setClass(Class = "prediction",
                       representation(
                         next_char = "character",
                         probability = "numeric",
                         index = "numeric",
                         alternative_probability = "matrix",
                         solution = "character"))
