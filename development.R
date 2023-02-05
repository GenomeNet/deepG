
devtools::load_all("deepG")


# reticulate::use_condaenv("~/projects/genomenet/install/miniconda2/envs/py37/")


devs <- tensorflow::tf$config$experimental$list_physical_devices("GPU")
devs

tensorflow::tf$config$experimental$set_memory_growth(devs[[1]], TRUE)

gpu <- tensorflow::tf$device(devs[[1]]$name)

model <- create_model_lstm_cnn(
    maxlen = 150, dropout = .49, recurrent_dropout = .7, layer_lstm = c(5, 5), layer_dense = c(4500, 2), solver = "adam",
    learning.rate = 0.0006, use.multiple.gpus = FALSE, merge.on.cpu = TRUE, gpu.num = 1, num_targets = 2,
    vocabulary.size = 4, bidirectional = TRUE, stateful = FALSE, batch.size = 2048 , compile = TRUE, kernel_size = c(83, 83),
    filters = c(47, 47), strides = c(1, 1), pool_size = c(4, 4), padding = "same", use_bias = TRUE, residual_block = FALSE,
    residual_block_length = 1, size_reduction_1Dconv = FALSE, label_input = NULL)

model
runif(100)
trained <- trainNetwork(
    train_type = "label_folder",
    model = model,
    run.name = paste0("run11", runif(1)),
    tensorboard.log = "tensorboard",
    path = c(
      "/media/int/home/aime/genomenet/bacteria_0_1/train",
      "/media/int/home/aime/genomenet/human_0_1/train"
    ),
    path.val = c(
      "/media/int/home/aime/genomenet/bacteria_0_1/validation",
      "/media/int/home/aime/genomenet/human_0_1/validation"
    ),
    validation.split = 0.4, # 0.4,
    batch.size = 2048 * 12,
    vocabulary = c("a", "c", "g", "t"),
    labelVocabulary = c("bacteria", "human"),
    epochs = 2,
    steps.per.epoch = 100L,
    shuffleFastaEntries = TRUE,
    output = list(none = FALSE, checkpoints = FALSE, tensorboard = TRUE, log = FALSE, serialize_model = FALSE, full_model = FALSE),
    output_format = "target_right",
    early_stopping_time = 3600,
    validation_only_after_training = FALSE
)


startInd <- 1:1024
numberOfSamples <- 1024
maxlen <- 150
vocabulary <- c("a", "c", "g", "t")
z <- matrix(c(1, 0, 0, 0, 0), nrow = 1173, ncol = 4)

microbenchmark::microbenchmark({
  x <- array(0, dim = c(numberOfSamples, maxlen, length(vocabulary)))
  for (i in 1:numberOfSamples) {
    start <- startInd[i]
    x[i, , ] <- z[start : (start + maxlen - 1), ]
  }
}, array(distribute(as.integer(startInd),as.integer(maxlen), as.integer(length(vocabulary)), z), c(numberOfSamples, maxlen, length(vocabulary))))



library("Rcpp")
zx <- array(distribute(1:1024, 150L, 4L, z), c(1024, 150, 4))

microbenchmark::microbenchmark({
  abind::abind(zx, zx, along = 1)
}, {
  combine(zx, zx)
})

r1 <- abind::abind(zx, zx[1:50, , ], along = 1)
r2 <- combine(zx, zx[1:50, , ])



which(is.na(abind::abind(zx, zx[1:50, , ], along = 1)))
which(is.na(combine(zx, zx[1:50, , ])))


ai

microbenchmark::microbenchmark({
  abind::abind(z, z, along = 1)
})

cx <-


keras:::as_generator.function

tools <- import_from_path("kerastools", path = system.file("python", package = "keras"))


reticulate::py_iterator


import("rpytools")$generator$RGenerator

import("rpy2")$robjects

mm <- model$fit_generator(generator = keras:::as_generator.function(gen), steps_per_epoch = 100L,
  epochs = 2L, verbose = 1L,
  callbacks = callbacks, validation_data = keras:::as_generator.function(validation_data),
  validation_steps = 40L,
  class_weight = NULL, max_queue_size = 100L,
  workers = 2L, initial_epoch = 0L, use_multiprocessing = TRUE)
