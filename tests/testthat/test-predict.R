context("predict")

test_that("Sucessful prediction", {
   
   #testthat::skip_if_not_installed("tensorflow")
   testthat::skip_if_not(reticulate::py_module_available("tensorflow"))
   
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
                         output_type = "h5",
                         model = model,
                         mode = "label", 
                         include_seq = TRUE)
   
   expect_true(all(pred$states >= 0))
   expect_true(all(pred$states <= 1))
   expect_equal(pred$sample_end_position, 8:12)
   
   pred_h5 <- load_prediction(filename, get_sample_position = TRUE)
   expect_equal(pred_h5$states, pred$states)
   expect_equal(pred_h5$sample_end_position, pred$sample_end_position)
   
   # batch size bigger than number of samples
   pred2 <- predict_model(layer_name = NULL, sequence = sequence,
                          filename = NULL, step = 1,
                          batch_size = 100, 
                          return_states = TRUE,
                          verbose = FALSE,
                          output_type = "h5",
                          model = model,
                          mode = "label", 
                          include_seq = TRUE)
   
   expect_true(all(abs(pred$states - pred2$states) < 1e-06))
   expect_equal(pred$sample_end_position, pred2$sample_end_position)
   
   # test csv + padding maxlen + ... (nuc_dist)
   filename <- tempfile(fileext = ".csv")
   pred <- predict_model(layer_name = NULL, sequence = sequence,
                         filename = filename, step = 1,
                         batch_size = 2, 
                         return_states = TRUE,
                         padding = "maxlen",
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
                         padding = "standard",
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
                         padding = "standard",
                         verbose = FALSE,
                         output_type = "csv", model = model,
                         mode = "label", 
                         include_seq = TRUE)
   
   expect_equal(pred$sample_end_position, c(8, 10))
   expect_equal(nrow(pred$states), length(pred$sample_end_position))
   
   # fasta file by_entry
   Sequence <- c("AAAACCCC", "TT", "AAACCCGGGTTT")
   Header <- letters[1:3]   
   df <- data.frame(Sequence, Header)
   fasta_path <- tempfile(fileext = ".fasta")
   microseq::writeFasta(df, fasta_path)
   output_path <- tempfile()
   dir.create(output_path)
   
   expect_message(
      predict_model(layer_name = NULL, 
                    path_input = fasta_path,
                    output_format = "by_entry",
                    output_dir = output_path,
                    filename = "states.h5",
                    step = 2,
                    batch_size = 2, 
                    padding = "none",
                    verbose = TRUE,
                    output_type = "h5",
                    model = model,
                    mode = "label", 
                    include_seq = TRUE)
   )
   
   h5_files <- list.files(output_path, full.names = TRUE)
   expect_true(basename(h5_files[1]) == "states_nr_1.h5")
   expect_true(basename(h5_files[2]) == "states_nr_3.h5")
   
   output_list_1 <- load_prediction(h5_files[1], get_sample_position = TRUE)
   expect_equal(output_list_1$sample_end_position, 8)
   output_list_2 <- load_prediction(h5_files[2], get_sample_position = TRUE)
   expect_equal(output_list_2$sample_end_position, c(8,10,12))
   
   # fasta file, by_entry
   h5_file <- tempfile(fileext = ".h5")
   pred <- predict_model(layer_name = NULL, 
                         path_input = fasta_path,
                         output_format = "by_entry_one_file",
                         filename = h5_file,
                         step = 2,
                         batch_size = 2, 
                         padding = "none",
                         verbose = FALSE,
                         output_type = "h5",
                         model = model,
                         mode = "label")
   
   output_list <- load_prediction(h5_file, get_sample_position = TRUE)
   expect_true(all(output_list[[1]]$states == output_list_1$states))
   expect_true(all(output_list[[1]]$sample_end_position == output_list_1$sample_end_position))
   expect_true(all(output_list[[2]]$states == output_list_2$states))
   expect_true(all(output_list[[2]]$sample_end_position == output_list_2$sample_end_position))
   
   # one pred per entry
   h5_file <- tempfile(fileext = ".h5")
   pred <- predict_model(layer_name = NULL, 
                         path_input = fasta_path,
                         output_format = "one_pred_per_entry",
                         filename = h5_file,
                         step = 2,
                         batch_size = 2, 
                         verbose = FALSE,
                         output_type = "h5",
                         model = model,
                         mode = "label")
   
   output_list <- load_prediction(h5_file)
   expect_equal(nrow(output_list$states),  nrow(df))
   
   # lm, target middle
   model <- create_model_lstm_cnn_target_middle(
      maxlen = maxlen,
      verbose = FALSE,
      layer_dense = 4,
      layer_lstm = 8)
   
   h5_file <- tempfile(fileext = ".h5")
   pred <- predict_model(layer_name = NULL, 
                         path_input = fasta_path,
                         output_format = "by_entry_one_file",
                         filename = h5_file,
                         step = 2,
                         target_len = 1,
                         batch_size = 2, 
                         padding = "standard",
                         verbose = FALSE,
                         output_type = "h5",
                         lm_format = "target_middle_lstm",
                         model = model,
                         mode = "lm")
   
   output_list <- load_prediction(h5_file, get_sample_position = TRUE)
   expect_equal(output_list[[1]]$sample_end_position, 8)
   expect_equal(output_list[[2]]$sample_end_position, 2)
   expect_equal(output_list[[3]]$sample_end_position[1], 9)
   
})
