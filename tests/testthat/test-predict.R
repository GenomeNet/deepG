context("predict")

test_that("Sucessful prediction", {
   
   devtools::load_all()
   library(testthat)
   
   sequence <- "AAACCNGGGTTT"
   maxlen <- 8
   filename <- tempfile(fileext = ".h5")
   
   model <- create_model_lstm_cnn(
      maxlen = maxlen,
      verbose = FALSE,
      layer_dense = 4,
      layer_lstm = 8)
   
   # test h5 output
   pred <- predict_model(layer_name = NULL, sequence = sequence,
                         filename = filename, step = 1,
                         batch_size = 1, 
                         return_states = TRUE,
                         verbose = FALSE,
                         output_type = "h5", model = model,
                         mode = "label", 
                         include_seq = TRUE)
   
   expect_true(all(pred$states >= 0))
   expect_true(all(pred$states <= 1))
   expect_equal(pred$sample_end_position, 8:12)
   
   pred_h5 <- load_prediction(filename, get_sample_position = TRUE)
   expect_equal(pred_h5$states, pred$states)
   expect_equal(pred_h5$sample_end_position, pred$sample_end_position)
   
   # test csv + padding maxlen + ...
   filename <- tempfile(fileext = ".csv")
   pred <- predict_model(layer_name = NULL, sequence = sequence,
                         filename = filename, step = 1,
                         batch_size = 2, 
                         return_states = TRUE,
                         padding_maxlen = TRUE,
                         verbose = FALSE,
                         output_type = "csv", model = model,
                         mode = "label", ambiguous_nuc = "empirical",
                         nuc_dist = c(0.1,0.4,0.4,0.1),
                         include_seq = TRUE)
   
   expect_true(all(pred$states >= 0))
   expect_true(all(pred$states <= 1))
   expect_equal(pred$sample_end_position, 0:12)
   
   pred_csv <- read.csv(filename)
   expect_equal(as.matrix(pred_csv), pred$states)
   
   # padding 
   pred <- predict_model(layer_name = NULL, sequence = "AAA",
                         filename = NULL, step = 2,
                         batch_size = 2, 
                         return_states = TRUE,
                         padding = TRUE,
                         verbose = FALSE,
                         output_type = "csv", model = model,
                         mode = "label", 
                         include_seq = TRUE)
   
   expect_true(all(pred$states >= 0))
   expect_true(all(pred$states <= 1))
   expect_equal(pred$sample_end_position, 3)
   expect_equal(nrow(pred$states), length(pred$sample_end_position))
   
   # step
   pred <- predict_model(layer_name = NULL, sequence = "AAAAACCCCC",
                         filename = NULL, step = 2,
                         batch_size = 2, 
                         return_states = TRUE,
                         padding = TRUE,
                         verbose = FALSE,
                         output_type = "csv", model = model,
                         mode = "label", 
                         include_seq = TRUE)
   
   expect_equal(pred$sample_end_position, c(8, 10))
   expect_equal(nrow(pred$states), length(pred$sample_end_position))
  
   # fasta file, by_entry
   
   devtools::load_all() ##############
   
   Sequence <- c("AAAACCCC", "TT", "AAACCCGGGTTT")
   Header <- letters[1:3]   
   df <- data.frame(Sequence, Header)
   fasta_path <- tempfile(fileext = ".fasta")
   microseq::writeFasta(df, fasta_path)
   output_path <- tempfile()
   dir.create(output_path)
   
   pred <- predict_model(layer_name = NULL, 
                         path_input = fasta_path,
                         output_format = "by_entry",
                         output_dir = output_path,
                         filename = "states.h5",
                         step = 2,
                         batch_size = 2, 
                         return_states = TRUE,
                         padding = FALSE,
                         verbose = FALSE,
                         output_type = "h5",
                         model = model,
                         mode = "label", 
                         include_seq = TRUE)
   
   expect_true(pred$sample_end_position, c(8, 10))
   expect_equal(nrow(pred$states), length(pred$sample_end_position))
   
   
})
