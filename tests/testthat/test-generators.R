context("generators")

test_that("Checking the generator for the Fasta files", {
  
  testpath <- file.path("fasta_2")
  vocabulary <- c("a", "c", "g", "t")
  batch_size <- 5
  maxlen <- 3
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = batch_size, maxlen = maxlen, vocabulary = vocabulary)
  
  arrays <- gen()
  
  expect_equivalent(dim(arrays[[1]])[1], batch_size)
  expect_equivalent(dim(arrays[[1]])[2], maxlen)
  expect_equivalent(dim(arrays[[1]])[3], length(vocabulary))
  expect_equivalent(dim(arrays[[2]])[1], batch_size)
  expect_equivalent(dim(arrays[[2]])[2], length(vocabulary))
  expect_equivalent(length(arrays),2)
  
  # a.fasta file starts with aaccggtt
  
  expect_equivalent(arrays[[1]][1, 1, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays[[1]][1, 2, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays[[1]][1, 3, ], c(0, 1, 0, 0)) # c
  expect_equivalent(arrays[[2]][1, ], c(0, 1, 0, 0)) # c
  
  arrays_2 <- gen()
  
  expect_equivalent(arrays_2[[1]][2, 1, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays_2[[1]][2, 2, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays_2[[1]][2, 3, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays_2[[2]][2, ], c(0, 0, 0, 1)) # t
  
  # test transition to second fasta file
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = batch_size, maxlen = maxlen, vocabulary = vocabulary)
  for (i in 1:5){
    arrays <- gen()
  }
  
  # samples start at beginning of b.fasta
  expect_equivalent(arrays[[1]][1, 1, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays[[1]][1, 2, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays[[1]][1, 3, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays[[2]][1, ], c(1, 0, 0, 0)) # a
  
  expect_equivalent(arrays[[1]][5, 1, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays[[1]][5, 2, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays[[1]][5, 3, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays[[2]][5, ], c(1, 0, 0, 0)) # a
  
  # complete one iteration (100 samples)
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = batch_size, maxlen = maxlen, vocabulary = vocabulary)
  for (i in 1:9){
    arrays <- gen()
  }
  
  # start from a.fasta again
  expect_equivalent(arrays[[1]][1, 1, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays[[1]][1, 2, ], c(1, 0, 0, 0)) # a
  expect_equivalent(arrays[[1]][1, 3, ], c(0, 1, 0, 0)) # c
  expect_equivalent(arrays[[2]][1, ], c(0, 1, 0, 0)) # c
  
  ###################
  # test for different step size
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = 4, maxlen = 3, step = 2)
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][3, 1, ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][3, 2, ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][3, 3, ], c(0, 0, 0, 1)) 
  expect_equivalent(arrays[[2]][3, ], c(0, 0, 0, 1)) 
  
  expect_equivalent(arrays[[1]][4, 1, ], c(1, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][4, 2, ], c(1, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][4, 3, ], c(1, 0, 0, 0)) 
  expect_equivalent(arrays[[2]][4, ], c(1, 0, 0, 0)) 
  
  ####
  # tests with chars outside vocabulary, vocabulary does not contain "A"
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = 5, maxlen = 3, step = 2, vocabulary = c("c", "g", "t"))
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1, ], c(0, 0, 0)) # a 
  expect_equivalent(arrays[[1]][1, 2, ], c(0, 0, 0)) # a
  expect_equivalent(arrays[[1]][1, 3, ], c(1, 0, 0)) # c
  expect_equivalent(arrays[[2]][1, ], c(1, 0, 0)) # c
  
  ####
  # test padding
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = 1, maxlen = 10, step = 4,
                            vocabulary = c("a", "c", "g", "t"))
  arrays <- gen()
  expect_equivalent(arrays[[1]][1, 1, ], c(0, 0, 0, 0))  
  expect_equivalent(arrays[[1]][1, 2, ], c(0, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][1, 3, ], c(0, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][1, 4, ], c(1, 0, 0, 0))  
  expect_equivalent(arrays[[1]][1, 5, ], c(1, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][1, 6, ], c(0, 1, 0, 0)) 
  expect_equivalent(arrays[[1]][1, 7, ], c(0, 1, 0, 0)) 
  expect_equivalent(arrays[[1]][1, 8, ], c(0, 0, 1, 0))  
  expect_equivalent(arrays[[1]][1, 9, ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][1, 10, ], c(0, 0, 0, 1)) 
  expect_equivalent(arrays[[2]][1, ], c(0, 0, 0, 1)) 
  
  # no padding
  testpath <- file.path("fasta_3")
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = 2, maxlen = 12, step = 1,
                            vocabulary = c("a", "c", "g", "t"), padding = FALSE)
  arrays <- gen()
  expect_equivalent(arrays[[1]][1, 1, ], c(1, 0, 0, 0))  
  expect_equivalent(arrays[[1]][1, 2, ], c(1, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][1, 3, ], c(0, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][1, 4, ], c(0, 1, 0, 0))  
  expect_equivalent(arrays[[1]][1, 5, ], c(0, 1, 0, 0)) 
  
  expect_equivalent(arrays[[1]][2, 1, ], c(1, 0, 0, 0))  
  expect_equivalent(arrays[[1]][2, 2, ], c(1, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][2, 3, ], c(0, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][2, 4, ], c(0, 1, 0, 0))  
  expect_equivalent(arrays[[1]][2, 5, ], c(0, 1, 0, 0)) 
  ####
  testpath <- file.path("fasta_2")
  expect_error(generator_fasta_lm())
  expect_error(generator_fasta_lm(""))
  
  expect_is(generator_fasta_lm(testpath, batch_size = batch_size, maxlen = maxlen), "function")
  expect_is(gen(), "list")
  expect_is(gen()[[1]], "array")
  expect_is(gen()[[2]], "matrix")
  
  expect_silent(generator_fasta_lm(testpath, batch_size = batch_size, maxlen = maxlen))
  
  expect_type(gen()[[1]], "double")
  expect_type(gen()[[2]], "double")
  
  ############# Test label generator (header) #############
  testpath <- file.path("fasta_2")
  
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = 5, maxlen = 3, step = 2, vocabulary = c("a", "c", "g", "t"),
                                          reverse_complement = FALSE, vocabulary_label = c("w", "x", "y"))
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1, ], c(1, 0, 0, 0)) # A  
  expect_equivalent(arrays[[1]][1, 2, ], c(1, 0, 0, 0)) # A
  expect_equivalent(arrays[[1]][1, 3, ], c(0, 1, 0, 0)) # C
  expect_equivalent(arrays[[2]][1, ], c(1, 0, 0)) # W 
  
  expect_equivalent(arrays[[1]][5, 1, ], c(1, 0, 0, 0)) # A  
  expect_equivalent(arrays[[1]][5, 2, ], c(1, 0, 0, 0)) # A
  expect_equivalent(arrays[[1]][5, 3, ], c(0, 0, 0, 1)) # T
  expect_equivalent(arrays[[2]][5, ], c(0, 1, 0)) # W 
  
  
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = 5, maxlen = 8, step = 2, vocabulary = c("a", "c", "g", "t"),
                                          reverse_complement = FALSE, vocabulary_label = c("w", "x", "y"))
  
  arrays <- gen()
  expect_equivalent(arrays[[2]][1, ], c(1, 0, 0))  
  expect_equivalent(arrays[[2]][2, ], c(0, 1, 0))  
  expect_equivalent(arrays[[2]][3, ], c(0, 0, 1))  
  expect_equivalent(arrays[[2]][4, ], c(0, 1, 0))  
  expect_equivalent(arrays[[2]][5, ], c(0, 1, 0))  
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][5, 1, ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][5, 2, ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][5, 3, ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][5, 4, ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][5, 5, ], c(0, 0, 0, 1)) 
  expect_equivalent(arrays[[1]][5, 6, ], c(0, 0, 0, 1)) 
  expect_equivalent(arrays[[1]][5, 7, ], c(0, 0, 0, 1)) 
  expect_equivalent(arrays[[1]][5, 8, ], c(0, 0, 0, 1)) 
  expect_equivalent(arrays[[2]][1, ], c(0, 0, 1))  
  expect_equivalent(arrays[[2]][2, ], c(0, 0, 1))  
  expect_equivalent(arrays[[2]][3, ], c(1, 0, 0))  
  expect_equivalent(arrays[[2]][4, ], c(0, 1, 0))  
  expect_equivalent(arrays[[2]][5, ], c(0, 0, 1))  
  
  
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = 8, maxlen = 7, step = 2, vocabulary = c("a", "c", "g", "t"),
                                          reverse_complement = FALSE, vocabulary_label = c("w", "x", "y"))
  
  arrays <- gen()
  
  # go through a/b.fasta once discard samples with target z
  expect_equivalent(arrays[[1]][8, 1, ], c(1, 0, 0, 0)) # A  
  expect_equivalent(arrays[[1]][8, 2, ], c(1, 0, 0, 0)) # A
  expect_equivalent(arrays[[1]][8, 3, ], c(0, 1, 0, 0)) # C
  expect_equivalent(arrays[[2]][8, ], c(1, 0, 0)) # W 
  
  
  ############# Test label generator (folder) #############
  directories <- c("label_folder/x", "label_folder/y", "label_folder/z")
  val <- FALSE
  gen_list <- generator_initialize(directories = directories,
                                   val = val,
                                   format = "fasta",
                                   batch_size = 6,
                                   maxlen = 2,
                                   vocabulary = c("a", "c", "g", "t"),
                                   step = 2)
  
  gen <- generator_fasta_label_folder_wrapper(val = val, path = directories, gen_list = gen_list)
  arrays <- gen()
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(0, 1, 0, 0)) 
  expect_equivalent(arrays[[1]][2, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 1, 0, 0)) 
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 1, 0, 0)) 
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 1, 0, 0)) 
  expect_equivalent(arrays[[1]][5, 1,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][5, 2,  ], c(0, 0, 0, 1)) 
  expect_equivalent(arrays[[1]][6, 1,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][6, 2,  ], c(0, 0, 0, 1)) 
  
  expect_equivalent(arrays[[2]][1,  ], c(1, 0, 0)) 
  expect_equivalent(arrays[[2]][2,  ], c(1, 0, 0))
  expect_equivalent(arrays[[2]][3,  ], c(0, 1, 0)) 
  expect_equivalent(arrays[[2]][4,  ], c(0, 1, 0))
  expect_equivalent(arrays[[2]][5,  ], c(0, 0, 1)) 
  expect_equivalent(arrays[[2]][6,  ], c(0, 0, 1))
  
  
  # test skipping file 
  for (i in 1:2){
    arrays <- gen()
  }
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(0, 1, 0, 0)) 
  expect_equivalent(arrays[[1]][2, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][5, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 2,  ], c(1, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][6, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 2,  ], c(1, 0, 0, 0)) 
  
  expect_equivalent(arrays[[2]][1,  ], c(1, 0, 0)) 
  expect_equivalent(arrays[[2]][2,  ], c(1, 0, 0))
  expect_equivalent(arrays[[2]][3,  ], c(0, 1, 0)) 
  expect_equivalent(arrays[[2]][4,  ], c(0, 1, 0))
  expect_equivalent(arrays[[2]][5,  ], c(0, 0, 1)) 
  expect_equivalent(arrays[[2]][6,  ], c(0, 0, 1))
  
  # 
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][2, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 1, 0, 0)) 
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 1, 0, 0)) 
  expect_equivalent(arrays[[1]][5, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 2,  ], c(1, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][6, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 2,  ], c(1, 0, 0, 0)) 
  
  expect_equivalent(arrays[[2]][1,  ], c(1, 0, 0)) 
  expect_equivalent(arrays[[2]][2,  ], c(1, 0, 0))
  expect_equivalent(arrays[[2]][3,  ], c(0, 1, 0)) 
  expect_equivalent(arrays[[2]][4,  ], c(0, 1, 0))
  expect_equivalent(arrays[[2]][5,  ], c(0, 0, 1)) 
  expect_equivalent(arrays[[2]][6,  ], c(0, 0, 1))
  
  ####### Test discard ambiguous nucleotides ###########
  
  testpath <- file.path("fasta_3")
  vocabulary = c("a", "c", "g", "t")
  batch_size <- 6
  maxlen <- 3
  step <- 2
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = batch_size, maxlen = maxlen, 
                            vocabulary = vocabulary, ambiguous_nuc = "discard", step = step)
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][1, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][1, ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][2, ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][3, ], c(0, 0, 0, 1))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][4, ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][5, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][5, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][5, 3,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][5, ], c(0, 0, 0, 1))
  
  expect_equivalent(arrays[[1]][6, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][6, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][6, ], c(0, 0, 0, 1))
  
  # label header
  
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = batch_size, maxlen = maxlen, 
                                          vocabulary = vocabulary, ambiguous_nuc = "discard", step = step, reverse_complement = FALSE,
                                          vocabulary_label = c("X", "Y"))
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][1, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][1, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][2, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][3, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][4, ], c(0, 1))
  
  expect_equivalent(arrays[[1]][5, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][5, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][5, 3,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][5, ], c(0, 1))
  
  expect_equivalent(arrays[[1]][6, 1,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][6, 2,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][6, 3,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][6, ], c(0, 1))
  
  # label folder
  
  directories = c("fasta_2", "fasta_3")
  gen <- get_generator(val = FALSE,
                       train_type = "label_folder",
                       path = directories,
                       format = "fasta",
                       batch_size = 6,
                       maxlen = 3,
                       ambiguous_nuc = "discard",
                       vocabulary = c("a", "c", "g", "t"),
                       reverse_complement = FALSE, 
                       step = 2)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[2]][1, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][2, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][3, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][4, ], c(0, 1))
  
  expect_equivalent(arrays[[1]][5, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][5, ], c(0, 1))
  
  expect_equivalent(arrays[[1]][6, 1,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][6, 2,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][6, 3,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][6, ], c(0, 1))
  
  ####### Test ambiguous nucleotides as 1/length(vocabulary) ###########
  
  testpath <- file.path("fasta_3")
  vocabulary = c("a", "c", "g", "t")
  batch_size <- 4
  maxlen <- 3
  step <- 2
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = batch_size, maxlen = maxlen, 
                            vocabulary = vocabulary, ambiguous_nuc = "equal", step = step)
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(1/4, 1/4, 1/4, 1/4))
  expect_equivalent(arrays[[2]][1, ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(1/4, 1/4, 1/4, 1/4))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[2]][2, ], c(1/4, 1/4, 1/4, 1/4))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(1/4, 1/4, 1/4, 1/4))
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][3, ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(1/4, 1/4, 1/4, 1/4))
  expect_equivalent(arrays[[2]][4, ], c(0, 0, 0, 1))
  
  # label header
  
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = batch_size, maxlen = maxlen,
                                          vocabulary = vocabulary, ambiguous_nuc = "equal", step = step, reverse_complement = FALSE,
                                          vocabulary_label = c("X", "Y"))
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(1/4, 1/4, 1/4, 1/4))
  expect_equivalent(arrays[[2]][1, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(1/4, 1/4, 1/4, 1/4))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[2]][2, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(1/4, 1/4, 1/4, 1/4))
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][3, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(1/4, 1/4, 1/4, 1/4))
  expect_equivalent(arrays[[2]][4, ], c(1, 0))
  
  # label folder
  
  directories = c("fasta_2", "fasta_3")
  gen <- get_generator(train_type = "label_folder",
                       val = FALSE,
                       path = directories,
                       format = "fasta",
                       batch_size = 4,
                       maxlen = 3,
                       vocabulary = c("a", "c", "g", "t"),
                       reverse_complement = FALSE, 
                       ambiguous_nuc = "equal",
                       step = 2)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[2]][1, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][2, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 3,  ], c(1/4, 1/4, 1/4, 1/4))
  expect_equivalent(arrays[[2]][3, ], c(0, 1))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(1/4, 1/4, 1/4, 1/4))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[2]][4, ], c(0, 1))
  
  ####### Test ambiguous nucleotides as "empirical" ###########
  # LM
  
  testpath <- file.path("fasta_3")
  vocabulary <- c("a", "c", "g", "t")
  batch_size <- 4
  maxlen <- 3
  step <- 2
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = batch_size, maxlen = maxlen, 
                            vocabulary = vocabulary, ambiguous_nuc = "empirical", step = step)
  arrays <- gen()
  nuc_dist <- 1/18*c(8, 2, 3, 5)
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], nuc_dist)
  expect_equivalent(arrays[[2]][1, ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], nuc_dist)
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[2]][2, ], nuc_dist)
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], nuc_dist)
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][3, ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], nuc_dist)
  expect_equivalent(arrays[[2]][4, ], c(0, 0, 0, 1))
  
  # LM second file
  
  testpath <- file.path("fasta_3")
  vocabulary <- c("a", "c", "g", "t")
  batch_size <- 4
  maxlen <- 3
  step <- 20
  gen <- generator_fasta_lm(path_corpus = testpath, batch_size = batch_size, maxlen = maxlen, 
                            vocabulary = vocabulary, ambiguous_nuc = "empirical", step = step)
  arrays <- gen()
  nuc_dist_1 <- 1/18*c(8, 2, 3, 5)
  nuc_dist_2 <- 1/17*c(3, 2, 6, 6)
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], nuc_dist)
  expect_equivalent(arrays[[2]][1, ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][2, ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][3, ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][4, ], c(0, 0, 0, 1))
  
  # label header
  
  testpath <- file.path("fasta_3")
  vocabulary <- c("a", "c", "g", "t")
  batch_size <- 4
  maxlen <- 3
  step <- 2
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = batch_size, maxlen = maxlen, 
                                          vocabulary = vocabulary, ambiguous_nuc = "empirical", step = step, reverse_complement = FALSE,
                                          vocabulary_label = c("X", "Y"))
  arrays <- gen()
  nuc_dist <- 1/18*c(8, 2, 3, 5)
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], nuc_dist)
  expect_equivalent(arrays[[2]][1, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], nuc_dist)
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[2]][2, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], nuc_dist)
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][3, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], nuc_dist)
  expect_equivalent(arrays[[2]][4, ], c(1, 0))
  
  # label folder
  
  directories = c("fasta_2", "fasta_3")
  gen <- get_generator(path = directories,
                       val = FALSE,
                       train_type = "label_folder",
                       format = "fasta",
                       batch_size = 4,
                       maxlen = 3,
                       vocabulary = c("a", "c", "g", "t"),
                       reverse_complement = FALSE, 
                       ambiguous_nuc = "empirical",
                       step = 2)
  
  arrays <- gen()
  nuc_dist <- 1/18*c(8, 2, 3, 5)
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[2]][1, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[2]][2, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 3,  ], nuc_dist)
  expect_equivalent(arrays[[2]][3, ], c(0, 1))
  
  expect_equivalent(arrays[[1]][4, 1,  ], nuc_dist)
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[2]][4, ], c(0, 1))
  
  ############# padding/amb nucleotide LM ############
  
  gen <- generator_fasta_lm(path_corpus = "fasta_3",
                            batch_size = 3,
                            maxlen = 15,
                            step = 1,
                            ambiguous_nuc = "equal")
  
  arrays <- gen()
  equal_vector <- rep(0.25, 4)
  expect_equivalent(arrays[[1]][1, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 4,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 5,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 6,  ], equal_vector)
  expect_equivalent(arrays[[1]][1, 7,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][1, 8,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][1, 9,  ], equal_vector)
  expect_equivalent(arrays[[1]][1, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][1, 11,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][1, 12,  ], equal_vector)
  expect_equivalent(arrays[[1]][1, 13,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][1, 14,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][1, 15,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][1, ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 5,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 6,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 7,  ], equal_vector)
  expect_equivalent(arrays[[1]][3, 8,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 9,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 11,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 12,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][3, 13,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][3, 14,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][3, 15,  ], equal_vector)
  expect_equivalent(arrays[[2]][3, ], c(0, 0, 0, 1))
  
  ############# padding/amb nucleotide, label_header ############
  
  gen <- generator_fasta_label_header_csv(path_corpus = "fasta_3",
                                          batch_size = 3,
                                          maxlen = 15,
                                          step = 1,
                                          vocabulary_label = c("X", "Y"),
                                          reverse_complement = FALSE,
                                          ambiguous_nuc = "empirical")
  
  nuc_dist_1 <- 1/18*c(8, 2, 3, 5)
  nuc_dist_2 <- 1/17*c(3, 2, 6, 6)
  arrays <- gen()
  expect_equivalent(arrays[[1]][1, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 4,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 5,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][1, 6,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][1, 7,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][1, 8,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][1, 9,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][1, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][1, 11,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][1, 12,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][1, 13,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][1, 14,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 14,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][1, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 5,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 6,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][3, 7,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 8,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 9,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 11,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][3, 12,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][3, 13,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][3, 14,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][3, 15,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][3, ], c(0, 1))
  
  ############# padding/amb nucleotide, label_folder ############
  
  directories = c("fasta_2", "fasta_3")
  gen <- get_generator(path = directories,
                       val = FALSE,
                       train_type = "label_folder",
                       padding = TRUE,
                       format = "fasta",
                       batch_size = 6,
                       maxlen = 15,
                       ambiguous_nuc = "equal",
                       vocabulary = c("a", "c", "g", "t"),
                       reverse_complement = FALSE, 
                       step = 1)
  
  equal_vector <- rep(0.25, 4)
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 5,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 6,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 7,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 8,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 9,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 10,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][1, 11,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][1, 12,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][1, 13,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][1, 14,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][1, 14,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][1, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][6, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 5,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 6,  ], equal_vector)
  expect_equivalent(arrays[[1]][6, 7,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][6, 8,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][6, 9,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][6, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][6, 11,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][6, 12,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][6, 13,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][6, 14,  ], equal_vector)
  expect_equivalent(arrays[[1]][6, 15,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][6, ], c(0, 1))
  
  ###### more than 2 files in one batch ######
  # LM
  
  gen <- generator_fasta_lm(path_corpus = "fasta_3",
                            batch_size = 8,
                            maxlen = 12,
                            max_iter = 10000,
                            step = 50, 
                            ambiguous_nuc = "empirical")
  
  nuc_dist_1 <- 1/18*c(8, 2, 3, 5)
  nuc_dist_2 <- 1/17*c(3, 2, 6, 6)
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][1, 4,  ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 4,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 5,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 6,  ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 4,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][3, 5,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 6,  ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 5,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 6,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 7,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][4, 8,  ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][5, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 3,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 4,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 5,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 6,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 7,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 8,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 9,  ], nuc_dist_2)
  
  expect_equivalent(arrays[[1]][6, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 3,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][6, 4,  ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][7, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 4,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 5,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 6,  ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][8, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][8, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][8, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][8, 4,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][8, 5,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][8, 6,  ], c(0, 0, 1, 0))
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 5,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 6,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 7,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][1, 8,  ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][2, 4,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][2, 5,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][2, 6,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][2, 7,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][2, 8,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 9,  ], nuc_dist_2)
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 3,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][3, 4,  ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 4,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 5,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 6,  ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][5, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 4,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 5,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][5, 6,  ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][6, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 5,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 6,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 7,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][6, 8,  ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][7, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 3,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][7, 4,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][7, 5,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][7, 6,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][7, 7,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][7, 8,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 9,  ], nuc_dist_2)
  
  expect_equivalent(arrays[[1]][8, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][8, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][8, 3,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][8, 4,  ], c(0, 1, 0, 0))
  
  # label header
  
  gen <- generator_fasta_label_header_csv(path_corpus = "fasta_3",
                                          batch_size = 8,
                                          maxlen = 12,
                                          max_iter = 10000,
                                          step = 50, 
                                          ambiguous_nuc = "empirical",
                                          reverse_complement = FALSE,
                                          vocabulary_label = c("X", "Y")
  )
  
  nuc_dist_1 <- 1/18*c(8, 2, 3, 5)
  nuc_dist_2 <- 1/17*c(3, 2, 6, 6)
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][1, 4,  ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 4,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 5,  ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 3,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][3, 4,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][3, 5,  ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 5,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 6,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][4, 7,  ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][5, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 2,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 3,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 4,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 5,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 6,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 7,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 8,  ], nuc_dist_2)
  
  expect_equivalent(arrays[[1]][6, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 3,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][6, 4,  ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][7, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 4,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 5,  ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][8, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][8, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][8, 3,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][8, 4,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][8, 5,  ], c(0, 0, 1, 0))
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 5,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][1, 6,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][1, 7,  ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][2, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 2,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][2, 3,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][2, 4,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][2, 5,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][2, 6,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][2, 7,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][2, 8,  ], nuc_dist_2)
  
  expect_equivalent(arrays[[1]][3, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][3, 3,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][3, 4,  ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 3,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 4,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 5,  ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][5, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][5, 3,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][5, 4,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][5, 5,  ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][6, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 5,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][6, 6,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][6, 7,  ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][7, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 2,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][7, 3,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][7, 4,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][7, 5,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][7, 6,  ], nuc_dist_2)
  expect_equivalent(arrays[[1]][7, 7,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 8,  ], nuc_dist_2)
  
  expect_equivalent(arrays[[1]][8, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][8, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][8, 3,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][8, 4,  ], c(0, 1, 0, 0))
  
  # label folder
  
  directories = c("fasta_2", "fasta_3")
  gen <- get_generator(path = directories,
                       train_type = "label_folder",
                       batch_size = 20,
                       maxlen = 12,
                       val = FALSE,
                       padding = TRUE,
                       ambiguous_nuc = "empirical",
                       vocabulary = c("a", "c", "g", "t"),
                       reverse_complement = FALSE, 
                       step = 1)
  
  nuc_dist_1 <- 1/18*c(8, 2, 3, 5)
  nuc_dist_2 <- 1/17*c(3, 2, 6, 6)
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][9, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][9, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][9, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][9, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][9, 5,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][9, 6,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][9, 7,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][9, 8,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][9, 9,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][9, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][9, 11,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][9, 12,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][9, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][12, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][12, 2,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][12, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][12, 4,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][12, 5,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][12, 6,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][12, 7,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][12, 8,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][12, 9,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][12, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][12, 11,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][12, 12,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][12, ], c(0, 1))
  
  expect_equivalent(arrays[[1]][18, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][18, 2,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][18, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][18, 4,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][18, 5,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][18, 6,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][18, 7,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][18, 8,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][18, 9,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][18, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][18, 11,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][18, 12,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][18, ], c(0, 1))
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][7, 1,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 2,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 3,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 4,  ], c(0, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 5,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 6,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][7, 7,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][7, 8,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][7, 9,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][7, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][7, 11,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][7, 12,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[2]][7, ], c(1, 0))
  
  expect_equivalent(arrays[[1]][14, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][14, 2,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][14, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][14, 4,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][14, 5,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][14, 6,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][14, 7,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][14, 8,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][14, 9,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][14, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][14, 11,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][14, 12,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][14, ], c(0, 1))
  
  expect_equivalent(arrays[[1]][20, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][20, 2,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][20, 3,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][20, 4,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][20, 5,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][20, 6,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][20, 7,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][20, 8,  ], nuc_dist_1)
  expect_equivalent(arrays[[1]][20, 9,  ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][20, 10,  ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][20, 11,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][20, 12,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][20, ], c(0, 1))
  
  # test quality scores LM
  
  gen <- generator_fasta_lm(path_corpus = "fastq",
                            format = "fastq",
                            batch_size = 10,
                            maxlen = 3,
                            max_iter = 10000,
                            vocabulary = c("a", "c", "g", "t"),
                            verbose = FALSE,
                            shuffle_file_order = FALSE,
                            step = 2, 
                            seed = 1234,
                            shuffle_input = FALSE,
                            file_limit = NULL,
                            path_file_log = NULL,
                            reverse_complement = FALSE,
                            output_format = "target_right",
                            ambiguous_nuc = "zeros",
                            use_quality_score = TRUE,    
                            proportion_per_seq = NULL,
                            padding = FALSE)
  
  a <- create_quality_vector(prob = quality_to_probability("J") , pos = 1, voc_length = 4)
  c <- create_quality_vector(prob = quality_to_probability("C") , pos = 2, voc_length = 4)
  g <- create_quality_vector(prob = quality_to_probability("G") , pos = 3, voc_length = 4)
  t <- create_quality_vector(prob = quality_to_probability("?") , pos = 4, voc_length = 4)
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][1, 1,  ], a)
  expect_equivalent(arrays[[1]][1, 2,  ], a)
  expect_equivalent(arrays[[1]][1, 3,  ], c)
  
  expect_equivalent(arrays[[1]][2, 1,  ], c)
  expect_equivalent(arrays[[1]][2, 2,  ], c)
  expect_equivalent(arrays[[1]][2, 3,  ], g)
  
  expect_equivalent(arrays[[1]][3, 1,  ], a)
  expect_equivalent(arrays[[1]][3, 2,  ], c)
  expect_equivalent(arrays[[1]][3, 3,  ], g)
  
  expect_equivalent(arrays[[1]][4, 1,  ], g)
  expect_equivalent(arrays[[1]][4, 2,  ], t)
  expect_equivalent(arrays[[1]][4, 3,  ], a)
  
  expect_equivalent(arrays[[1]][5, 1,  ], c)
  expect_equivalent(arrays[[1]][5, 2,  ], g)
  expect_equivalent(arrays[[1]][5, 3,  ], t)
  
  expect_equivalent(arrays[[1]][6, 1,  ], t)
  expect_equivalent(arrays[[1]][6, 2,  ], c)
  expect_equivalent(arrays[[1]][6, 3,  ], g)
  
  expect_equivalent(arrays[[1]][7, 1,  ], a)
  expect_equivalent(arrays[[1]][7, 2,  ], t)
  expect_equivalent(arrays[[1]][7, 3,  ], a)
  
  expect_equivalent(arrays[[1]][8, 1,  ], a)
  expect_equivalent(arrays[[1]][8, 2,  ], a)
  expect_equivalent(arrays[[1]][8, 3,  ], c)
  
  expect_equivalent(arrays[[2]][1, ], c)
  expect_equivalent(arrays[[2]][2, ], g)
  expect_equivalent(arrays[[2]][3, ], t)
  expect_equivalent(arrays[[2]][4, ], c)
  expect_equivalent(arrays[[2]][5, ], c)
  expect_equivalent(arrays[[2]][6, ], t)
  expect_equivalent(arrays[[2]][7, ], t)
  expect_equivalent(arrays[[2]][8, ], c)
  
  # test quality scores label
  
  gen <- generator_fasta_label_folder(path_corpus = "fastq",
                                      format = "fastq",
                                      batch_size = 10,
                                      maxlen = 3,
                                      max_iter = 10000,
                                      vocabulary = c("a", "c", "g", "t"),
                                      verbose = FALSE,
                                      shuffle_file_order = FALSE,
                                      step = 2, 
                                      seed = 1234,
                                      shuffle_input = FALSE,
                                      file_limit = NULL,
                                      path_file_log = NULL,
                                      reverse_complement = FALSE,
                                      ambiguous_nuc = "zeros",
                                      use_quality_score = TRUE,    
                                      proportion_per_seq = NULL,
                                      num_targets = 2,
                                      ones_column = 1,
                                      padding = FALSE)
  
  a <- create_quality_vector(prob = quality_to_probability("J") , pos = 1, voc_length = 4)
  c <- create_quality_vector(prob = quality_to_probability("C") , pos = 2, voc_length = 4)
  g <- create_quality_vector(prob = quality_to_probability("G") , pos = 3, voc_length = 4)
  t <- create_quality_vector(prob = quality_to_probability("?") , pos = 4, voc_length = 4)
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][1, 1,  ], a)
  expect_equivalent(arrays[[1]][1, 2,  ], a)
  expect_equivalent(arrays[[1]][1, 3,  ], c)
  
  expect_equivalent(arrays[[1]][2, 1,  ], c)
  expect_equivalent(arrays[[1]][2, 2,  ], c)
  expect_equivalent(arrays[[1]][2, 3,  ], g)
  
  expect_equivalent(arrays[[1]][3, 1,  ], a)
  expect_equivalent(arrays[[1]][3, 2,  ], c)
  expect_equivalent(arrays[[1]][3, 3,  ], g)
  
  expect_equivalent(arrays[[1]][4, 1,  ], g)
  expect_equivalent(arrays[[1]][4, 2,  ], t)
  expect_equivalent(arrays[[1]][4, 3,  ], a)
  
  expect_equivalent(arrays[[1]][5, 1,  ], a)
  expect_equivalent(arrays[[1]][5, 2,  ], c)
  expect_equivalent(arrays[[1]][5, 3,  ], g)
  
  expect_equivalent(arrays[[1]][6, 1,  ], c)
  expect_equivalent(arrays[[1]][6, 2,  ], g)
  expect_equivalent(arrays[[1]][6, 3,  ], t)
  
  expect_equivalent(arrays[[1]][7, 1,  ], t)
  expect_equivalent(arrays[[1]][7, 2,  ], c)
  expect_equivalent(arrays[[1]][7, 3,  ], g)
  
  expect_equivalent(arrays[[1]][8, 1,  ], a)
  expect_equivalent(arrays[[1]][8, 2,  ], t)
  expect_equivalent(arrays[[1]][8, 3,  ], a)
  
  expect_equivalent(arrays[[1]][9, 1,  ], a)
  expect_equivalent(arrays[[1]][9, 2,  ], t)
  expect_equivalent(arrays[[1]][9, 3,  ], a)
  
  expect_equivalent(arrays[[1]][10, 1,  ], a)
  expect_equivalent(arrays[[1]][10, 2,  ], a)
  expect_equivalent(arrays[[1]][10, 3,  ], c)
  
  expect_equivalent(arrays[[2]][1, ], c(1,0))
  expect_equivalent(arrays[[2]][10, ], c(1,0))
  
  ## test read data with quality
  
  gen <- generator_fasta_label_folder(path_corpus = "read_data",
                                      format = "fastq",
                                      batch_size = 5,
                                      maxlen = 12,
                                      max_iter = 10000,
                                      vocabulary = c("a", "c", "g", "t"),
                                      verbose = FALSE,
                                      shuffle_file_order = FALSE,
                                      step = 2, 
                                      seed = 1234,
                                      shuffle_input = FALSE,
                                      file_limit = NULL,
                                      path_file_log = NULL,
                                      read_data = TRUE,
                                      reverse_complement = FALSE,
                                      ambiguous_nuc = "zeros",
                                      use_quality_score = TRUE,    
                                      proportion_per_seq = NULL,
                                      num_targets = 2,
                                      ones_column = 1,
                                      padding = FALSE)
  
  a <- create_quality_vector(prob = quality_to_probability("J") , pos = 1, voc_length = 4)
  c <- create_quality_vector(prob = quality_to_probability("C") , pos = 2, voc_length = 4)
  g <- create_quality_vector(prob = quality_to_probability("G") , pos = 3, voc_length = 4)
  t <- create_quality_vector(prob = quality_to_probability("?") , pos = 4, voc_length = 4)
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][[1]][1, ,  ], rbind(a,a,a,c,c,c))
  expect_equivalent(arrays[[1]][[2]][1, ,  ], rbind(c,c,c,g,g,g))
  
  expect_equivalent(arrays[[1]][[1]][2, ,  ], rbind(a,c,a,c,a,c))
  expect_equivalent(arrays[[1]][[2]][2, ,  ], rbind(c,g,c,g,c,g))
  
  expect_equivalent(arrays[[1]][[1]][3, ,  ], rbind(g,g,g,t,t,t))
  expect_equivalent(arrays[[1]][[2]][3, ,  ], rbind(t,t,t,g,g,g))
  
  expect_equivalent(arrays[[1]][[1]][4, ,  ], rbind(g,t,g,t,g,t))
  expect_equivalent(arrays[[1]][[2]][4, ,  ], rbind(t,g,t,g,t,g))
  
  expect_equivalent(arrays[[1]][[1]][5, ,  ], rbind(a,a,a,c,c,c))
  expect_equivalent(arrays[[1]][[2]][5, ,  ], rbind(c,c,c,g,g,g))
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][[1]][1, ,  ], rbind(a,c,a,c,a,c))
  expect_equivalent(arrays[[1]][[2]][1, ,  ], rbind(c,g,c,g,c,g))
  
  expect_equivalent(arrays[[1]][[1]][2, ,  ], rbind(g,g,g,t,t,t))
  expect_equivalent(arrays[[1]][[2]][2, ,  ], rbind(t,t,t,g,g,g))
  
  expect_equivalent(arrays[[1]][[1]][3, ,  ], rbind(g,t,g,t,g,t))
  expect_equivalent(arrays[[1]][[2]][3, ,  ], rbind(t,g,t,g,t,g))
  
  # additional input LM
  
  gen <- generator_fasta_lm(path_corpus = "fasta_3",
                            format = "fasta",
                            batch_size = 10,
                            maxlen = 5,
                            vocabulary = c("a", "c", "g", "t"),
                            shuffle_file_order = FALSE,
                            step = 4, 
                            shuffle_input = FALSE,
                            reverse_complement = FALSE,
                            output_format = "target_right",
                            ambiguous_nuc = "zeros",
                            added_label_path = "label.csv",
                            add_input_as_seq = FALSE,
                            padding = FALSE)
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][[1]][1,], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][[1]][2,], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][[1]][3,], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][[1]][4,], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][[1]][5,], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][6,], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][7,], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][8,], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][9,], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][10,], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][[2]][10, 1, ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[2]][10, 3, ], c(0, 0, 0, 0))
  
  # additional input label_folder 
  dir <- c("label_folder/x", "label_folder/y", "label_folder/z")
  gen_list <- generator_initialize(directories = dir,
                                   format = "fasta",
                                   batch_size = 15,
                                   maxlen = 4,
                                   step = 2, 
                                   val = FALSE,
                                   padding = FALSE,
                                   added_label_path = "label.csv",
                                   add_input_as_seq = FALSE)
  
  gen <- generator_fasta_label_folder_wrapper(val = FALSE, path = dir, gen_list = gen_list) 
  arrays <- gen()
  expect_equivalent(arrays[[1]][[1]][1, ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][2, ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][3, ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][4, ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][5, ], c(1, 0, 0, 1))
  
  expect_equivalent(arrays[[1]][[1]][6, ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][[1]][7, ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][[1]][8, ], c(0, 1, 0, 1))
  expect_equivalent(arrays[[1]][[1]][9, ], c(0, 1, 0, 1))
  expect_equivalent(arrays[[1]][[1]][10, ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][[1]][11, ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][[1]][12, ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][[1]][13, ], c(0, 0, 1, 1))
  expect_equivalent(arrays[[1]][[1]][14, ], c(0, 0, 1, 1))
  expect_equivalent(arrays[[1]][[1]][15, ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][[2]][5, 1, ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[2]][5, 2, ], c(0, 0, 1, 0))
  
  expect_equivalent(arrays[[1]][[2]][10, 1, ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][[2]][10, 2, ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][[2]][15, 1, ], c(0, 0, 0, 1))
  expect_equivalent(arrays[[1]][[2]][15, 2, ], c(0, 0, 0, 1))
  
  gen <- generator_fasta_label_folder_wrapper(val = FALSE, path = dir, gen_list = gen_list) 
  arrays <- gen()
  expect_equivalent(arrays[[1]][[1]][1, ], c(1, 0, 0, 1))
  expect_equivalent(arrays[[1]][[1]][2, ], c(1, 0, 0, 1))
  expect_equivalent(arrays[[1]][[1]][3, ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][4, ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[1]][[1]][5, ], c(1, 0, 0, 0))
  
  expect_equivalent(arrays[[1]][[1]][6, ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][[1]][7, ], c(0, 1, 0, 1))
  expect_equivalent(arrays[[1]][[1]][8, ], c(0, 1, 0, 1))
  expect_equivalent(arrays[[1]][[1]][9, ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[1]][[1]][10, ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][[1]][11, ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][[1]][12, ], c(0, 0, 1, 1))
  expect_equivalent(arrays[[1]][[1]][13, ], c(0, 0, 1, 1))
  expect_equivalent(arrays[[1]][[1]][14, ], c(0, 0, 1, 0))
  expect_equivalent(arrays[[1]][[1]][15, ], c(0, 0, 1, 0))
  
  ## read data with quality and 2 classes
  
  gen <- get_generator(path = c("read_data_2/label_a", "read_data_2/label_b"),
                       train_type = "label_folder",
                       format = "fastq",
                       batch_size = 4,
                       maxlen = 12,
                       vocabulary = c("a", "c", "g", "t"),
                       verbose = FALSE,
                       shuffle_file_order = FALSE,
                       step = 1, 
                       seed = 1234,
                       shuffle_input = FALSE,
                       file_limit = NULL,
                       path_file_log = NULL,
                       reverse_complement = FALSE, 
                       val = FALSE,
                       ambiguous_nuc = "zero",
                       proportion_per_seq = NULL,
                       read_data = TRUE,
                       use_quality_score = TRUE,
                       padding = FALSE,
                       added_label_path = NULL,
                       skip_amb_nuc = NULL)
  
  arrays <- gen()
  
  a <- create_quality_vector(prob = quality_to_probability("J") , pos = 1, voc_length = 4)
  c <- create_quality_vector(prob = quality_to_probability("C") , pos = 2, voc_length = 4)
  g <- create_quality_vector(prob = quality_to_probability("G") , pos = 3, voc_length = 4)
  t <- create_quality_vector(prob = quality_to_probability("?") , pos = 4, voc_length = 4)
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][[1]][1, ,  ], rbind(a,a,a,a,a,a))
  expect_equivalent(arrays[[1]][[2]][1, ,  ], rbind(c,c,c,c,c,c))
  expect_equivalent(arrays[[1]][[1]][2, ,  ], rbind(a,a,a,a,a,a))
  expect_equivalent(arrays[[1]][[2]][2, ,  ], rbind(c,c,c,c,c,c))
  
  expect_equivalent(arrays[[1]][[1]][3, ,  ], rbind(g,g,g,g,g,g))
  expect_equivalent(arrays[[1]][[2]][3, ,  ], rbind(t,t,t,t,t,t))
  expect_equivalent(arrays[[1]][[1]][4, ,  ], rbind(g,g,g,g,g,g))
  expect_equivalent(arrays[[1]][[2]][4, ,  ], rbind(t,t,t,t,t,t))
  
  ### get output tensor from csv file + concat 
  
  testpath <- file.path("fasta_2")
  label_from_csv <- "output_label.csv"
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = 5,
                                          maxlen = 10, step = 10,
                                          vocabulary = c("a", "c", "g", "t", "Z"),
                                          reverse_complement = FALSE, 
                                          vocabulary_label = c("w", "x", "y"),
                                          format = "fasta",
                                          max_iter = 10000,
                                          verbose = FALSE,
                                          shuffle_file_order = FALSE,
                                          seed = 1234,
                                          shuffle_input = FALSE,
                                          file_limit = NULL,
                                          path_file_log = NULL,
                                          ambiguous_nuc = "zero",
                                          proportion_per_seq = NULL,
                                          read_data = FALSE,
                                          use_quality_score = FALSE,
                                          padding = TRUE,
                                          skip_amb_nuc = NULL,
                                          max_samples = NULL,
                                          concat_seq = "ZZ",
                                          added_label_path = NULL,
                                          add_input_as_seq = NULL,
                                          target_from_csv = label_from_csv)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 8, ], c(0, 0, 0, 1, 0)) 
  expect_equivalent(arrays[[1]][1, 9, ], c(0, 0, 0, 0, 1)) 
  expect_equivalent(arrays[[1]][1, 10, ], c(0, 0, 0, 0, 1))
  
  expect_equivalent(arrays[[1]][4, 1, ], c(1, 0, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][4, 2, ], c(1, 0, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][4, 3, ], c(1, 0, 0, 0, 0))
  expect_equivalent(arrays[[1]][4, 4, ], c(1, 0, 0, 0, 0))
  
  expect_equivalent(arrays[[2]][1, ], 1:4)
  expect_equivalent(arrays[[2]][2, ], 1:4)
  expect_equivalent(arrays[[2]][3, ], 1:4)
  expect_equivalent(arrays[[2]][4, ], 11:14)
  expect_equivalent(arrays[[2]][5, ], 11:14)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 8, ], c(1, 0, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][1, 9, ], c(0, 0, 0, 0, 1)) 
  expect_equivalent(arrays[[1]][1, 10, ], c(0, 0, 0, 0, 1))
  
  expect_equivalent(arrays[[1]][2, 1, ], c(1, 0, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][2, 2, ], c(1, 0, 0, 0, 0)) 
  expect_equivalent(arrays[[1]][2, 3, ], c(0, 1, 0, 0, 0))
  
  expect_equivalent(arrays[[2]][1, ], 11:14)
  expect_equivalent(arrays[[2]][2, ], 1:4)
  expect_equivalent(arrays[[2]][3, ], 1:4)
  expect_equivalent(arrays[[2]][4, ], 1:4)
  expect_equivalent(arrays[[2]][5, ], 11:14)
  
  
  ## 2 added input files LM 
  
  gen <- generator_fasta_lm(path_corpus = "fasta_3",
                            format = "fasta",
                            batch_size = 10,
                            maxlen = 5,
                            vocabulary = c("a", "c", "g", "t"),
                            shuffle_file_order = FALSE,
                            step = 4, 
                            shuffle_input = FALSE,
                            reverse_complement = FALSE,
                            output_format = "target_right",
                            ambiguous_nuc = "zeros",
                            added_label_path = c("label.csv",
                                                 "add_seq.csv"),
                            add_input_as_seq = c(FALSE, TRUE),
                            padding = FALSE)
  
  v1 <- c(0, 0, 1, 0)
  v2 <- c(1, 0, 0, 0)
  m1 <- matrix(c(1, 0, 0, 0, 
                 0, 1, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1), byrow = TRUE, ncol = 4)
  m2 <- matrix(c(0, 0, 0, 1, 
                 0, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 0, 0), byrow = TRUE, ncol = 4)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][[1]][1, ], v1)
  expect_equivalent(arrays[[1]][[1]][2, ], v1)
  expect_equivalent(arrays[[1]][[1]][3, ], v1)
  expect_equivalent(arrays[[1]][[1]][4, ], v1)
  expect_equivalent(arrays[[1]][[1]][5, ], v2)
  expect_equivalent(arrays[[1]][[1]][6, ], v2)
  expect_equivalent(arrays[[1]][[1]][7, ], v2)
  expect_equivalent(arrays[[1]][[1]][8, ], v2)
  expect_equivalent(arrays[[1]][[1]][9, ], v2)
  expect_equivalent(arrays[[1]][[1]][10, ], v1)
  
  expect_equivalent(arrays[[1]][[2]][1, , ], m1)
  expect_equivalent(arrays[[1]][[2]][2, , ], m1)
  expect_equivalent(arrays[[1]][[2]][3, , ], m1)
  expect_equivalent(arrays[[1]][[2]][4, , ], m1)
  expect_equivalent(arrays[[1]][[2]][5, , ], m2)
  expect_equivalent(arrays[[1]][[2]][6, , ], m2)
  expect_equivalent(arrays[[1]][[2]][7, , ], m2)
  expect_equivalent(arrays[[1]][[2]][8, , ], m2)
  expect_equivalent(arrays[[1]][[2]][9, , ], m2)
  expect_equivalent(arrays[[1]][[2]][10, , ], m1)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][[1]][1, ], v1)
  expect_equivalent(arrays[[1]][[1]][2, ], v1)
  expect_equivalent(arrays[[1]][[1]][3, ], v1)
  expect_equivalent(arrays[[1]][[1]][4, ], v2)
  expect_equivalent(arrays[[1]][[1]][5, ], v2)
  expect_equivalent(arrays[[1]][[1]][6, ], v2)
  expect_equivalent(arrays[[1]][[1]][7, ], v2)
  expect_equivalent(arrays[[1]][[1]][8, ], v2)
  expect_equivalent(arrays[[1]][[1]][9, ], v1)
  expect_equivalent(arrays[[1]][[1]][10, ], v1)
  
  expect_equivalent(arrays[[1]][[2]][1, , ], m1)
  expect_equivalent(arrays[[1]][[2]][2, , ], m1)
  expect_equivalent(arrays[[1]][[2]][3, , ], m1)
  expect_equivalent(arrays[[1]][[2]][4, , ], m2)
  expect_equivalent(arrays[[1]][[2]][5, , ], m2)
  expect_equivalent(arrays[[1]][[2]][6, , ], m2)
  expect_equivalent(arrays[[1]][[2]][7, , ], m2)
  expect_equivalent(arrays[[1]][[2]][8, , ], m2)
  expect_equivalent(arrays[[1]][[2]][9, , ], m1)
  expect_equivalent(arrays[[1]][[2]][10, , ], m1)
  
  ## 2 added input files, label_folder 
  
  dir <- c("label_folder/x", "label_folder/y", "label_folder/z")
  gen <- get_generator(path = dir,
                       train_type = "label_folder",
                       format = "fasta",
                       batch_size = 15,
                       maxlen = 4,
                       step = 2, 
                       val = FALSE,
                       padding = FALSE,
                       added_label_path = c("label.csv",
                                            "add_seq.csv"),
                       add_input_as_seq = c(FALSE, TRUE)
  )
  
  x1 <- c(1, 0, 0, 0)
  x2 <- c(1, 0, 0, 1)
  y1 <- c(0, 1, 0, 0)
  y2 <- c(0, 1, 0, 1)
  z1 <- c(0, 0, 1, 0)
  z2 <- c(0, 0, 1, 1)
  
  mx1 <- matrix(c(1, 0, 0, 0, 
                  1, 0, 0, 0,
                  1, 0, 0, 0,
                  1, 0, 0, 0), byrow = TRUE, ncol = 4)
  
  mx2 <- matrix(c(1, 0, 0, 0, 
                  0, 1, 0, 0,
                  1, 0, 0, 0,
                  0, 1, 0, 0), byrow = TRUE, ncol = 4)
  
  my1 <- matrix(c(0, 1, 0, 0, 
                  0, 1, 0, 0,
                  0, 1, 0, 0,
                  0, 1, 0, 0), byrow = TRUE, ncol = 4)
  
  my2 <- matrix(c(0, 1, 0, 0, 
                  0, 0, 1, 0,
                  0, 1, 0, 0,
                  0, 0, 1, 0), byrow = TRUE, ncol = 4)
  
  mz1 <- matrix(c(0, 0, 1, 0,
                  0, 0, 1, 0,
                  0, 0, 1, 0,
                  0, 0, 1, 0), byrow = TRUE, ncol = 4)
  
  mz2 <- matrix(c(0, 0, 1, 0, 
                  0, 0, 0, 1,
                  0, 0, 1, 0,
                  0, 0, 0, 1), byrow = TRUE, ncol = 4)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][[1]][1, ], x1)
  expect_equivalent(arrays[[1]][[1]][2, ], x1)
  expect_equivalent(arrays[[1]][[1]][3, ], x1)
  expect_equivalent(arrays[[1]][[1]][4, ], x1)
  expect_equivalent(arrays[[1]][[1]][5, ], x2)
  expect_equivalent(arrays[[1]][[1]][6, ], y1)
  expect_equivalent(arrays[[1]][[1]][7, ], y1)
  expect_equivalent(arrays[[1]][[1]][8, ], y2)
  expect_equivalent(arrays[[1]][[1]][9, ], y2)
  expect_equivalent(arrays[[1]][[1]][10, ], y1)
  expect_equivalent(arrays[[1]][[1]][11, ], z1)
  expect_equivalent(arrays[[1]][[1]][12, ], z1)
  expect_equivalent(arrays[[1]][[1]][13, ], z2)
  expect_equivalent(arrays[[1]][[1]][14, ], z2)
  expect_equivalent(arrays[[1]][[1]][15, ], z1)
  
  expect_equivalent(arrays[[1]][[2]][1, , ], mx1)
  expect_equivalent(arrays[[1]][[2]][2, , ], mx1)
  expect_equivalent(arrays[[1]][[2]][3, , ], mx1)
  expect_equivalent(arrays[[1]][[2]][4, , ], mx1)
  expect_equivalent(arrays[[1]][[2]][5, , ], mx2)
  expect_equivalent(arrays[[1]][[2]][6, , ], my1)
  expect_equivalent(arrays[[1]][[2]][7, , ], my1)
  expect_equivalent(arrays[[1]][[2]][8, , ], my2)
  expect_equivalent(arrays[[1]][[2]][9, , ], my2)
  expect_equivalent(arrays[[1]][[2]][10, , ], my1)
  expect_equivalent(arrays[[1]][[2]][11, , ], mz1)
  expect_equivalent(arrays[[1]][[2]][12, , ], mz1)
  expect_equivalent(arrays[[1]][[2]][13, , ], mz2)
  expect_equivalent(arrays[[1]][[2]][14, , ], mz2)
  expect_equivalent(arrays[[1]][[2]][15, , ], mz1)
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][[1]][1, ], x2)
  expect_equivalent(arrays[[1]][[1]][2, ], x2)
  expect_equivalent(arrays[[1]][[1]][3, ], x1)
  expect_equivalent(arrays[[1]][[1]][4, ], x1)
  expect_equivalent(arrays[[1]][[1]][5, ], x1)
  expect_equivalent(arrays[[1]][[1]][6, ], y1)
  expect_equivalent(arrays[[1]][[1]][7, ], y2)
  expect_equivalent(arrays[[1]][[1]][8, ], y2)
  expect_equivalent(arrays[[1]][[1]][9, ], y1)
  expect_equivalent(arrays[[1]][[1]][10, ], y1)
  expect_equivalent(arrays[[1]][[1]][11, ], z1)
  expect_equivalent(arrays[[1]][[1]][12, ], z2)
  expect_equivalent(arrays[[1]][[1]][13, ], z2)
  expect_equivalent(arrays[[1]][[1]][14, ], z1)
  expect_equivalent(arrays[[1]][[1]][15, ], z1)
  
  expect_equivalent(arrays[[1]][[2]][1, , ], mx2)
  expect_equivalent(arrays[[1]][[2]][2, , ], mx2)
  expect_equivalent(arrays[[1]][[2]][3, , ], mx1)
  expect_equivalent(arrays[[1]][[2]][4, , ], mx1)
  expect_equivalent(arrays[[1]][[2]][5, , ], mx1)
  expect_equivalent(arrays[[1]][[2]][6, , ], my1)
  expect_equivalent(arrays[[1]][[2]][7, , ], my2)
  expect_equivalent(arrays[[1]][[2]][8, , ], my2)
  expect_equivalent(arrays[[1]][[2]][9, , ], my1)
  expect_equivalent(arrays[[1]][[2]][10, , ], my1)
  expect_equivalent(arrays[[1]][[2]][11, , ], mz1)
  expect_equivalent(arrays[[1]][[2]][12, , ], mz2)
  expect_equivalent(arrays[[1]][[2]][13, , ], mz2)
  expect_equivalent(arrays[[1]][[2]][14, , ], mz1)
  expect_equivalent(arrays[[1]][[2]][15, , ], mz1)
  
  # 3 targets, target right
  
  gen <- generator_fasta_lm(path_corpus = "fasta_3",
                            batch_size = 5,
                            maxlen = 4,
                            step = 5, 
                            output_format = "target_right",
                            padding = FALSE,
                            target_len = 3)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, , ],
                    matrix(
                      c(1, 0, 0, 0,
                        1, 0, 0, 0, 
                        0, 0, 0, 0,
                        0, 1, 0, 0),
                      byrow = TRUE,  ncol = 4
                    ))
  
  expect_equivalent(arrays[[1]][5, , ],
                    matrix(
                      c(1, 0, 0, 0,
                        0, 1, 0, 0, 
                        0, 0, 1, 0,
                        0, 0, 0, 1),
                      byrow = TRUE,  ncol = 4
                    ))
  
  m1 <- matrix(
    c(0, 1, 0, 0,
      0, 0, 0, 1, 
      0, 0, 0, 0,
      0, 0, 1, 0,
      1, 0, 0, 0),
    byrow = TRUE,  ncol = 4
  )
  m2 <- matrix(
    c(0, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 0, 0, 1,
      0, 1, 0, 0),
    byrow = TRUE,  ncol = 4)
  m3 <- matrix(
    c(0, 0, 1, 0,
      1, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 1,
      0, 0, 1, 0),
    byrow = TRUE,  ncol = 4)
  expect_equivalent(arrays[[2]][ ,1 , ], m1)
  expect_equivalent(arrays[[2]][ ,2 , ], m2)
  expect_equivalent(arrays[[2]][ ,3 , ], m3)
  
  # 3 targets, target middle cnn
  
  gen <- generator_fasta_lm(path_corpus = "fasta_3",
                            batch_size = 5,
                            maxlen = 4,
                            step = 5, 
                            output_format = "target_middle_cnn",
                            padding = FALSE,
                            target_len = 3)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, , ],
                    matrix(
                      c(1, 0, 0, 0,
                        1, 0, 0, 0, 
                        0, 0, 0, 0,
                        0, 0, 1, 0),
                      byrow = TRUE,  ncol = 4
                    ))
  
  expect_equivalent(arrays[[1]][5, , ],
                    matrix(
                      c(1, 0, 0, 0,
                        0, 1, 0, 0, 
                        0, 1, 0, 0,
                        0, 0, 1, 0),
                      byrow = TRUE,  ncol = 4
                    ))
  
  m1 <- matrix(
    c(0, 0, 0, 0,
      0, 0, 1, 0, 
      1, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 1, 0),
    byrow = TRUE,  ncol = 4
  )
  m2 <- matrix(
    c(0, 1, 0, 0,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1),
    byrow = TRUE,  ncol = 4)
  m3 <- matrix(
    c(0, 1, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      0, 0, 1, 0,
      1, 0, 0, 0),
    byrow = TRUE,  ncol = 4)
  expect_equivalent(arrays[[2]][ ,1 , ], m1)
  expect_equivalent(arrays[[2]][ ,2 , ], m2)
  expect_equivalent(arrays[[2]][ ,3 , ], m3)
  
  # 3 targets, target middle lstm
  
  gen <- generator_fasta_lm(path_corpus = "fasta_3",
                            batch_size = 5,
                            maxlen = 4,
                            step = 5, 
                            output_format = "target_middle_lstm",
                            padding = FALSE,
                            target_len = 3)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][[1]][1, , ],
                    matrix(
                      c(1, 0, 0, 0,
                        1, 0, 0, 0),
                      byrow = TRUE,  ncol = 4
                    ))
  
  expect_equivalent(arrays[[1]][[2]][1, , ],
                    matrix(
                      c(0, 0, 1, 0,
                        0, 0, 0, 0),
                      byrow = TRUE,  ncol = 4
                    ))
  
  expect_equivalent(arrays[[1]][[1]][5, , ],
                    matrix(
                      c(1, 0, 0, 0,
                        0, 1, 0, 0),
                      byrow = TRUE,  ncol = 4
                    ))
  
  expect_equivalent(arrays[[1]][[2]][5, , ],
                    matrix(
                      c(0, 0, 1, 0,
                        0, 1, 0, 0),
                      byrow = TRUE,  ncol = 4
                    ))
  
  m1 <- matrix(
    c(0, 0, 0, 0,
      0, 0, 1, 0, 
      1, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 1, 0),
    byrow = TRUE,  ncol = 4
  )
  m2 <- matrix(
    c(0, 1, 0, 0,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1),
    byrow = TRUE,  ncol = 4)
  m3 <- matrix(
    c(0, 1, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      0, 0, 1, 0,
      1, 0, 0, 0),
    byrow = TRUE,  ncol = 4)
  expect_equivalent(arrays[[2]][ ,1 , ], m1)
  expect_equivalent(arrays[[2]][ ,2 , ], m2)
  expect_equivalent(arrays[[2]][ ,3 , ], m3)
  
  # coverage + set learning for label_folder 
  
  directories <- c("coverage_data/x", "coverage_data/y")
  val <- FALSE
  batch_size <- 6
  samples_per_target <- 3
  #new_batch_size <- batch_size/samples_per_target
  path <- directories
  voc_len <- 4
  maxlen <- 7
  reshape_mode <- "time_dist"
  set_learning <- list(reshape_mode = reshape_mode,
                       maxlen = maxlen,
                       samples_per_target = samples_per_target)
  
  gen <- get_generator(path = directories,
                       train_type = "label_folder",
                       val = FALSE,
                       padding = TRUE,
                       format = "fasta",
                       batch_size = batch_size,
                       maxlen = maxlen,
                       vocabulary = c("a", "c", "g", "t"),
                       step = 4,
                       use_coverage = 1,
                       set_learning = set_learning)
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][1, 1, , ], matrix(
    c(7,0,0,0,
      7,0,0,0,
      0,0,0,0,
      0,7,0,0,
      0,7,0,0,
      0,0,0,0,
      0,0,7,0),
    byrow = TRUE,  ncol = 4
  ))
  
  expect_equivalent(arrays[[1]][1, 3, , ], matrix(
    c(11,0,0,0,
      11,0,0,0,
      11,0,0,0,
      11,0,0,0,
      0,0,0,0,
      0,0,0,11,
      0,0,0,11),
    byrow = TRUE,  ncol = 4
  ))
  
  expect_equivalent(arrays[[1]][2, 1, , ], matrix(
    c(0,0,0,0,
      0,0,1,0,
      0,0,1,0,
      0,0,1,0,
      0,0,1,0,
      0,0,0,1,
      0,0,0,1),
    byrow = TRUE,  ncol = 4
  ))
  
  expect_equivalent(arrays[[1]][3, 1, , ], matrix(
    c(0,0,0,0,
      17,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0),
    byrow = TRUE,  ncol = 4
  ))
  
  expect_equivalent(arrays[[1]][3, 2, , ], matrix(
    c(7,0,0,0,
      7,0,0,0,
      0,0,0,0,
      0,7,0,0,
      0,7,0,0,
      0,0,0,0,
      0,0,7,0),
    byrow = TRUE,  ncol = 4
  ))
  
  expect_equivalent(arrays[[1]][4, 1, , ], matrix(
    c(2,0,0,0,
      0,2,0,0,
      0,0,2,0,
      0,0,0,2,
      2,0,0,0,
      2,0,0,0,
      0,2,0,0),
    byrow = TRUE,  ncol = 4
  ))
  
  expect_equivalent(arrays[[1]][5, 3, , ], matrix(
    c(0,0,0,0,
      17,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0),
    byrow = TRUE,  ncol = 4
  ))
  
  expect_equivalent(arrays[[1]][6, 3, , ], matrix(
    c(0,0,0,0,
      0,0,1,0,
      0,0,1,0,
      0,0,1,0,
      0,0,1,0,
      0,0,0,1,
      0,0,0,1),
    byrow = TRUE,  ncol = 4
  ))
  
  expect_equivalent(arrays[[2]], matrix(
    c(1,0,
      1,0,
      1,0,
      0,1,
      0,1,
      0,1),
    byrow = TRUE,  ncol = 2
  ))
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1, , ], matrix(
    c(11,0,0,0,
      11,0,0,0,
      11,0,0,0,
      11,0,0,0,
      0,0,0,0,
      0,0,0,11,
      0,0,0,11),
    byrow = TRUE,  ncol = 4
  ))
  
  expect_equivalent(arrays[[1]][4, 3, , ], matrix(
    c(0,0,0,0,
      17,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0),
    byrow = TRUE,  ncol = 4
  ))
  
  expect_equivalent(arrays[[2]], matrix(
    c(1,0,
      1,0,
      1,0,
      0,1,
      0,1,
      0,1),
    byrow = TRUE,  ncol = 2
  ))
  
  # coverage + set learning for label_folder + normalizing input tensor 
  
  directories <- c("coverage_data/x", "coverage_data/y")
  val <- FALSE
  batch_size <- 6
  samples_per_target <- 3
  #new_batch_size <- batch_size/samples_per_target
  path <- directories
  voc_len <- 4
  maxlen <- 7
  use_coverage <- 17
  reshape_mode <- "time_dist"
  set_learning <- list(reshape_mode = reshape_mode,
                       maxlen = maxlen,
                       samples_per_target = samples_per_target)
  
  gen <- get_generator(path = directories,
                       train_type = "label_folder",
                       val = FALSE,
                       padding = TRUE,
                       format = "fasta",
                       batch_size = batch_size,
                       maxlen = maxlen,
                       vocabulary = c("a", "c", "g", "t"),
                       step = 4,
                       use_coverage = use_coverage,
                       set_learning = set_learning)
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][1, 1, , ], matrix(
    c(7,0,0,0,
      7,0,0,0,
      0,0,0,0,
      0,7,0,0,
      0,7,0,0,
      0,0,0,0,
      0,0,7,0),
    byrow = TRUE,  ncol = 4
  )/use_coverage)
  
  expect_equivalent(arrays[[1]][1, 3, , ], matrix(
    c(11,0,0,0,
      11,0,0,0,
      11,0,0,0,
      11,0,0,0,
      0,0,0,0,
      0,0,0,11,
      0,0,0,11),
    byrow = TRUE,  ncol = 4
  )/use_coverage)
  
  expect_equivalent(arrays[[1]][2, 1, , ], matrix(
    c(0,0,0,0,
      0,0,1,0,
      0,0,1,0,
      0,0,1,0,
      0,0,1,0,
      0,0,0,1,
      0,0,0,1),
    byrow = TRUE,  ncol = 4
  )/use_coverage)
  
  expect_equivalent(arrays[[1]][3, 1, , ], matrix(
    c(0,0,0,0,
      17,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0),
    byrow = TRUE,  ncol = 4
  )/use_coverage)
  
  expect_equivalent(arrays[[1]][3, 2, , ], matrix(
    c(7,0,0,0,
      7,0,0,0,
      0,0,0,0,
      0,7,0,0,
      0,7,0,0,
      0,0,0,0,
      0,0,7,0),
    byrow = TRUE,  ncol = 4
  )/use_coverage)
  
  expect_equivalent(arrays[[1]][4, 1, , ], matrix(
    c(2,0,0,0,
      0,2,0,0,
      0,0,2,0,
      0,0,0,2,
      2,0,0,0,
      2,0,0,0,
      0,2,0,0),
    byrow = TRUE,  ncol = 4
  )/use_coverage)
  
  expect_equivalent(arrays[[1]][5, 3, , ], matrix(
    c(0,0,0,0,
      17,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0),
    byrow = TRUE,  ncol = 4
  )/use_coverage)
  
  expect_equivalent(arrays[[1]][6, 3, , ], matrix(
    c(0,0,0,0,
      0,0,1,0,
      0,0,1,0,
      0,0,1,0,
      0,0,1,0,
      0,0,0,1,
      0,0,0,1),
    byrow = TRUE,  ncol = 4
  )/use_coverage)
  
  expect_equivalent(arrays[[2]], matrix(
    c(1,0,
      1,0,
      1,0,
      0,1,
      0,1,
      0,1),
    byrow = TRUE,  ncol = 2
  ))
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1, , ], matrix(
    c(11,0,0,0,
      11,0,0,0,
      11,0,0,0,
      11,0,0,0,
      0,0,0,0,
      0,0,0,11,
      0,0,0,11),
    byrow = TRUE,  ncol = 4
  )/use_coverage)
  
  expect_equivalent(arrays[[1]][4, 3, , ], matrix(
    c(0,0,0,0,
      17,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0),
    byrow = TRUE,  ncol = 4
  )/use_coverage)
  
  expect_equivalent(arrays[[2]], matrix(
    c(1,0,
      1,0,
      1,0,
      0,1,
      0,1,
      0,1),
    byrow = TRUE,  ncol = 2
  ))
  
  # rds label generator
  
  gen <- generator_rds(rds_folder = "rds", batch_size = 1)
  l_x <- list()
  l_y <- list()
  for (i in 1:40) {
    z <- gen()
    l_x[[i]] <- z[[1]][1,1,1]
    l_y[[i]] <- which.max(z[[2]])
  }
  expect_equivalent(sort(unlist(l_x)), rep(1:20, each=2)) 
  expect_equivalent(sort(unlist(l_y)), rep(1:20, each=2)) 
  
  gen <- generator_rds(rds_folder = "rds", batch_size = 10)
  l_x <- list()
  l_y <- list()
  for (i in 1:4) {
    z <- gen()
    l_x[[i]] <- z[[1]][,1,1]
    l_y[[i]] <- apply(z[[2]], 1, which.max)
  }
  expect_equivalent(sort(unlist(l_x)), rep(1:20, each = 2)) 
  expect_equivalent(sort(unlist(l_y)), rep(1:20, each=2)) 
  
  #  rds lm generator
  
  target_len <- 3
  batch_size <- 1
  gen <- generator_rds(rds_folder = "rds_lm", batch_size = batch_size, target_len = target_len)
  
  for (one_iter in 1:3) {
    first_input <- 1 + (100*(0:4)) 
    for (i in 1:5) {
      z <- gen()
      expect_equivalent(dim(z[[1]]), c(batch_size, 7 - target_len, 4))
      l_x <- z[[1]][1,1,1]
      first_input <- setdiff(first_input, l_x)
      l_y <- NULL
      for (j in 1:target_len) {
        l_y[[j]] <- z[[2]][[j]][1,1]
      }
      expect_equivalent(l_y, l_x + 3 + (1:target_len))
    }
    expect_equivalent(length(first_input), 0)
  }
  
  batch_size <- 5
  gen <- generator_rds(rds_folder = "rds_lm", batch_size = batch_size, target_len = target_len)
  for (one_iter in 1:3) {
    first_input <- 1 + (100*(0:4)) 
    z <- gen()
    expect_equivalent(dim(z[[1]]), c(batch_size, 7 - target_len, 4))
    l_x <- z[[1]][ , 1, 1]
    first_input <- setdiff(first_input, l_x)
    l_y <- NULL
    for (j in 1:target_len) {
      l_y[[j]] <- z[[2]][[j]][,1]
    }
    expect_equivalent(sort(l_y[[1]]), 5 + (100*(0:4)))
    expect_equivalent(sort(l_y[[2]]), 6 + (100*(0:4)))
    expect_equivalent(sort(l_y[[3]]), 7 + (100*(0:4)))
    expect_equivalent(length(first_input), 0)
  }
  
  # n-gram rds
  
  n_gram <- 3
  gen <- generator_rds(rds_folder = "n_gram_rds",
                       batch_size = 1, 
                       target_len = 6,
                       n_gram = n_gram,
                       n_gram_stride = n_gram)
  
  arrays <- gen()
  y <- arrays[[2]]
  y_1_n_gram <- apply(y[[1]], 1, which.max)
  y_2_n_gram <- apply(y[[2]], 1, which.max)
  
  int_seq <- c(1,2,0)
  expect_equivalent(y_1_n_gram[1], 1 + sum(4^((n_gram-1):0) * (int_seq))) # cga
  int_seq <- c(0,0,1)
  expect_equivalent(y_2_n_gram[1], 1 + sum(4^((n_gram-1):0) * (int_seq))) # aac
  
  # set learning concat with coverage encoding
  
  directories <- c("coverage_data/x", "coverage_data/y")
  val <- FALSE
  batch_size <- 8
  samples_per_target <- 3
  #new_batch_size <- batch_size/samples_per_target
  path <- directories
  voc_len <- 4
  maxlen <- 6
  use_coverage <- 17
  reshape_mode <- "concat"
  set_learning <- list(reshape_mode = reshape_mode,
                       maxlen = maxlen,
                       buffer_len = NULL,
                       samples_per_target = samples_per_target)
  buffer_size <- 0
  concat_maxlen <- (maxlen * samples_per_target) + (buffer_size * (samples_per_target - 1))
  
  gen <- get_generator(path = directories,
                       train_type = "label_folder",
                       val = FALSE,
                       padding = TRUE,
                       format = "fasta",
                       batch_size = batch_size,
                       maxlen = maxlen,
                       vocabulary = c("a", "c", "g", "t"),
                       step = maxlen,
                       use_coverage = use_coverage,
                       set_learning = set_learning)
  
  m <- matrix(
    c(0,0,0,0,
      0,0,1/17,0,
      0,0,1/17,0,
      0,0,1/17,0,
      0,0,1/17,0,
      0,0,0,1/17,
      13/17,0,0,0,
      0,13/17,0,0,
      0,0,13/17,0,
      0,0,0,13/17,
      13/17,0,0,0,
      0,13/17,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      1,0,0,0),
    byrow = TRUE,  ncol = 4
  )
  
  m2 <- matrix(
    c(0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      1,0,0,0,
      2/17,0,0,0,
      0,2/17,0,0,
      0,0,2/17,0,
      0,0,0,2/17,
      2/17,0,0,0,
      2/17,0,0,0,
      0,0,3/17,0,
      0,0,3/17,0,
      0,0,3/17,0,
      0,0,3/17,0,
      0,0,0,3/17,
      0,0,0,3/17),
    byrow = TRUE,  ncol = 4
  )
  
  y <- matrix(c(1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1), ncol = 2, byrow = TRUE)
  
  arrays <- gen()
  expect_true(all(arrays[[1]][1,,] == arrays[[1]][3,,]))
  expect_true(all(arrays[[1]][2,,] == arrays[[1]][4,,]))
  expect_equivalent(arrays[[2]], y)
  
  
  expect_equivalent(arrays[[1]][4, , ], m)
  expect_equivalent(arrays[[1]][8, , ], m2)
  
  arrays <- gen()
  expect_true(all(arrays[[1]][1,,] == arrays[[1]][3,,]))
  expect_true(all(arrays[[1]][2,,] == arrays[[1]][4,,]))
  expect_equivalent(arrays[[2]], y)
  expect_equivalent(arrays[[1]][4, , ], m)
  
  arrays <- gen()
  expect_true(all(arrays[[1]][1,,] == arrays[[1]][3,,]))
  expect_true(all(arrays[[1]][2,,] == arrays[[1]][4,,]))
  expect_equivalent(arrays[[2]], y)
  expect_equivalent(arrays[[1]][4, , ], m)
  expect_equivalent(arrays[[1]][5, , ], m2)
  
  # rds generator with multi inputs/outputs
  
  x1 <- array(0, dim = c(9,5,4))
  x2 <- array(0, dim = c(9,5,3))
  y1 <- array(0, dim = c(9,2))
  y2 <- array(0, dim = c(9,6))
  
  for (i in 1:dim(x1)[1]) {
    x1[i,,] <- i
    y1[i, ] <- i
    x2[i,,] <- i + 10
    y2[i, ] <- i + 10
  }
  
  index_1 <- 1:5
  index_2 <- 6:9
  x_list_1 <- list(x1[index_1, , ], x2[index_1, , ])
  x_list_2 <- list(x1[index_2, , ], x2[index_2, , ])
  y_list_1 <- list(y1[index_1, ], y2[index_1, ])
  y_list_2 <- list(y1[index_2, ], y2[index_2, ])
  z1 <- list(x = x_list_1, y = y_list_1)
  z2 <- list(x = x_list_2, y = y_list_2)
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  saveRDS(z1, paste0(temp_dir, "/file_1.rds"))
  saveRDS(z2, paste0(temp_dir, "/file_2.rds"))
  
  gen <- generator_rds(rds_folder = temp_dir,
                       batch_size = 10, path_file_log = NULL,
                       max_samples = NULL,
                       proportion_per_seq = NULL,
                       target_len = NULL,
                       seed = 1,
                       reverse_complement = FALSE,
                       sample_by_file_size = FALSE,
                       n_gram = NULL, n_gram_stride = 1,
                       reverse_complement_encoding = FALSE,
                       add_noise = NULL)
  
  for (k in 1:5) {
    z <- gen()
    x1 <- z[[1]][[1]] %>% as.array()
    x2 <- z[[1]][[2]] %>% as.array()
    y1 <- z[[2]][[1]] %>% as.array()
    y2 <- z[[2]][[2]] %>% as.array()
    
    
    for (i in 1:dim(x1)[1]) {
      expect_equivalent(min(x1[i,,]), max(y1[i,]))
      expect_equivalent(min(x1[i,,]) + 10, max(x2[i,,]))
      expect_equivalent(max(x2[i,,]), min(y2[i,]))
      expect_equivalent(max(y1[i,]) + 10, min(y2[i,]))
    }
  }
  
  # integer encoding label header #
  
  testpath <- file.path("fasta_2")
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = 5, maxlen = 3, step = 2, vocabulary = c("a", "c", "g", "t"),
                                          reverse_complement = FALSE, vocabulary_label = c("w", "x", "y"), return_int = TRUE)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1], 1) # A  
  expect_equivalent(arrays[[1]][1, 2], 1) # A
  expect_equivalent(arrays[[1]][1, 3], 2) # C
  expect_equivalent(arrays[[2]][1, ], c(1, 0, 0)) # W 
  
  expect_equivalent(arrays[[1]][5, 1], 1) # A  
  expect_equivalent(arrays[[1]][5, 2], 1) # A
  expect_equivalent(arrays[[1]][5, 3], 4) # T
  expect_equivalent(arrays[[2]][5, ], c(0, 1, 0)) # W 
  
  
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = 5, maxlen = 8, step = 2, vocabulary = c("a", "c", "g", "t"),
                                          reverse_complement = FALSE, vocabulary_label = c("w", "x", "y"), return_int = TRUE)
  
  arrays <- gen()
  expect_equivalent(arrays[[2]][1, ], c(1, 0, 0))  
  expect_equivalent(arrays[[2]][2, ], c(0, 1, 0))  
  expect_equivalent(arrays[[2]][3, ], c(0, 0, 1))  
  expect_equivalent(arrays[[2]][4, ], c(0, 1, 0))  
  expect_equivalent(arrays[[2]][5, ], c(0, 1, 0))  
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][5, 1], 3) 
  expect_equivalent(arrays[[1]][5, 2], 3) 
  expect_equivalent(arrays[[1]][5, 3], 3) 
  expect_equivalent(arrays[[1]][5, 4], 3) 
  expect_equivalent(arrays[[1]][5, 5], 4) 
  expect_equivalent(arrays[[1]][5, 6], 4) 
  expect_equivalent(arrays[[1]][5, 7], 4) 
  expect_equivalent(arrays[[1]][5, 8], 4) 
  expect_equivalent(arrays[[2]][1, ], c(0, 0, 1))  
  expect_equivalent(arrays[[2]][2, ], c(0, 0, 1))  
  expect_equivalent(arrays[[2]][3, ], c(1, 0, 0))  
  expect_equivalent(arrays[[2]][4, ], c(0, 1, 0))  
  expect_equivalent(arrays[[2]][5, ], c(0, 0, 1))  
  
  
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = 8, maxlen = 7, step = 2, vocabulary = c("a", "c", "g", "t"),
                                          reverse_complement = FALSE, vocabulary_label = c("w", "x", "y"), return_int = TRUE)
  
  arrays <- gen()
  
  # go through a/b.fasta once discard samples with target z
  expect_equivalent(arrays[[1]][8, 1], 1) # A  
  expect_equivalent(arrays[[1]][8, 2], 1) # A
  expect_equivalent(arrays[[1]][8, 3], 2) # C
  expect_equivalent(arrays[[2]][8, ], c(1, 0, 0)) # W 
  
  # label folder with integer encoding
  
  directories <- c("label_folder/x", "label_folder/y", "label_folder/z")
  gen <- get_generator(path = directories,
                       train_type = "label_folder",
                       val = FALSE,
                       padding = TRUE,
                       format = "fasta",
                       batch_size = 6,
                       maxlen = 2,
                       return_int = TRUE,
                       vocabulary = c("a", "c", "g", "t"),
                       step = 2)
  
  arrays <- gen()
  expect_equivalent(arrays[[1]][1, 1], 1)
  expect_equivalent(arrays[[1]][1, 2], 2) 
  expect_equivalent(arrays[[1]][2, 1], 1)
  expect_equivalent(arrays[[1]][2, 2], 2) 
  expect_equivalent(arrays[[1]][3, 1], 3)
  expect_equivalent(arrays[[1]][3, 2], 2) 
  expect_equivalent(arrays[[1]][4, 1], 3)
  expect_equivalent(arrays[[1]][4, 2], 2) 
  expect_equivalent(arrays[[1]][5, 1], 4)
  expect_equivalent(arrays[[1]][5, 2], 4) 
  expect_equivalent(arrays[[1]][6, 1], 4)
  expect_equivalent(arrays[[1]][6, 2], 4) 
  
  expect_equivalent(arrays[[2]][1,  ], c(1, 0, 0)) 
  expect_equivalent(arrays[[2]][2,  ], c(1, 0, 0))
  expect_equivalent(arrays[[2]][3,  ], c(0, 1, 0)) 
  expect_equivalent(arrays[[2]][4,  ], c(0, 1, 0))
  expect_equivalent(arrays[[2]][5,  ], c(0, 0, 1)) 
  expect_equivalent(arrays[[2]][6,  ], c(0, 0, 1))
  
  
  # test skipping file 
  for (i in 1:2) {
    arrays <- gen()
  }
  
  expect_equivalent(arrays[[1]][1, 1], 1)
  expect_equivalent(arrays[[1]][1, 2], 2) 
  expect_equivalent(arrays[[1]][2, 1], 1)
  expect_equivalent(arrays[[1]][2, 2], 3) 
  expect_equivalent(arrays[[1]][3, 1], 2)
  expect_equivalent(arrays[[1]][3, 2], 3) 
  expect_equivalent(arrays[[1]][4, 1], 2)
  expect_equivalent(arrays[[1]][4, 2], 3) 
  expect_equivalent(arrays[[1]][5, 1], 1)
  expect_equivalent(arrays[[1]][5, 2], 1) 
  expect_equivalent(arrays[[1]][6, 1], 1)
  expect_equivalent(arrays[[1]][6, 2], 1) 
  
  expect_equivalent(arrays[[2]][1,  ], c(1, 0, 0)) 
  expect_equivalent(arrays[[2]][2,  ], c(1, 0, 0))
  expect_equivalent(arrays[[2]][3,  ], c(0, 1, 0)) 
  expect_equivalent(arrays[[2]][4,  ], c(0, 1, 0))
  expect_equivalent(arrays[[2]][5,  ], c(0, 0, 1)) 
  expect_equivalent(arrays[[2]][6,  ], c(0, 0, 1))
  
  
  # n-gram integer encoding, label folder #
  
  directories <- c("label_folder/x", "label_folder/y", "label_folder/z")
  gen <- get_generator(path = directories,
                       train_type = "label_folder",
                       batch_size = 6,
                       maxlen = 12,
                       padding = TRUE,
                       n_gram = 3,
                       n_gram_stride = 2, 
                       return_int = TRUE,
                       vocabulary = c("a", "c", "g", "t"),
                       step = 2)
  
  arrays <- gen()
  x <- arrays[[1]]
  y <- arrays[[2]]
  expect_equivalent(dim(x), c(6, 5))
  expect_equivalent(x[1, 1], 0) # padding
  expect_equivalent(x[1, 2], 5) # ACA
  expect_equivalent(unique(x[5, 1:4]), 0) # padding
  expect_equivalent(x[5, 5], 64) # TTT = 4^3
  
  # n-gram one-hot encoding, label folder #
  
  directories <- c("label_folder/x", "label_folder/y", "label_folder/z")
  gen <- get_generator(path = directories,
                       train_type = "label_folder",
                       batch_size = 6,
                       maxlen = 12,
                       padding = TRUE,
                       n_gram = 3,
                       n_gram_stride = 2, 
                       return_int = FALSE,
                       vocabulary = c("a", "c", "g", "t"),
                       step = 2)
  
  arrays <- gen()
  x <- arrays[[1]]
  y <- arrays[[2]]
  expect_equivalent(dim(x), c(6, 5, 64))
  expect_equivalent(unique(x[1, 1, ]), 0) # padding
  expect_equivalent(which.max(x[1, 2, ]), 5) # ACA
  expect_equivalent(unique(as.vector(x[5, 1:4, ])), 0) # padding
  expect_equivalent(which.max(x[5, 5, ]), 64) # TTT = 4^3
  
  ##### masked lm #####
  
  testpath <- file.path("a.fastq")
  masked_lm <- list(mask_rate = 0.25, random_rate = 0.25, identity_rate = 0.25, include_sw = TRUE)
  gen <- get_generator(path = testpath,
                       train_type = "masked_lm",
                       masked_lm = masked_lm,
                       batch_size = 1,
                       maxlen = 200,
                       format = "fastq",
                       padding = TRUE,
                       return_int = TRUE)
  
  z <- gen()
  x <- z[[1]]
  y <- z[[2]]
  sw <- z[[3]]
  
  expect_equivalent(x[1,1:12], rep(0, 12)) # padding
  expect_equivalent(sw[1,1:12], rep(0, 12)) # no sample weights in padding region
  sw_pos <- which(sw[1,] == 1)
  random_pos <- which(x[1,] %in% c(2,3,4))
  masked_pos <- which(x[1,] == 5)
  # masked and random positions must have sw 1
  expect_contains(sw_pos, random_pos)
  expect_contains(sw_pos, masked_pos)
  
  ###
  
  testpath <- file.path("fasta_2/b.fasta")
  masked_lm <- list(mask_rate = 0.25, random_rate = 0.25, identity_rate = 0.25, include_sw = TRUE)
  gen <- get_generator(path = testpath,
                       train_type = "masked_lm",
                       shuffle_input = FALSE,
                       masked_lm = masked_lm,
                       batch_size = 3,
                       maxlen = 10,
                       padding = TRUE,
                       return_int = TRUE)
  
  z <- gen()
  x <- z[[1]]
  y <- z[[2]]
  sw <- z[[3]]
  
  expect_equivalent(sum(x[,1:2]), 0) # padding
  expect_equivalent(sum(sw[,1:2]), 0) # no sample weights in padding region
  for (i in 1:3) {
    sw_pos <- which(sw[i,] == 1)
    masked_pos <- which(x[i,] == 5)
    expect_contains(sw_pos, masked_pos)   # masked positions must have sw 1
  }  
  
  #### test reshape #### 
  
  directories <- c("fasta_2", "fasta_3")
  fx <- function(x) {return(x)}
  reshape_xy <- list(x = fx)
  expect_error(gen <- get_generator(path = directories,
                                    reshape_xy = reshape_xy,
                                    train_type = "label_folder",
                                    batch_size = 4,
                                    maxlen = 3))
  
  
  directories <- c("fasta_2", "fasta_3")
  fx <- function(x = NULL, y = NULL) {
    return(x + 1)
  }
  fy <- function(x = NULL, y = NULL) {
    return(x)
  }
  reshape_xy <- list(x = fx, y = fy)
  gen <- get_generator(path = directories,
                       reshape_xy = reshape_xy,
                       val = FALSE,
                       train_type = "label_folder",
                       format = "fasta",
                       batch_size = 4,
                       maxlen = 3,
                       vocabulary = c("a", "c", "g", "t"),
                       reverse_complement = FALSE, 
                       ambiguous_nuc = "zero",
                       step = 2)
  
  arrays <- gen()
  arrays[[1]][1,,]
  y <- arrays[[2]]
  
  expect_equivalent(arrays[[1]][1, 1,  ], c(1, 0, 0, 0) + 1)
  expect_equivalent(arrays[[1]][1, 2,  ], c(1, 0, 0, 0) + 1)
  expect_equivalent(arrays[[1]][1, 3,  ], c(0, 1, 0, 0) + 1)
  expect_equivalent(arrays[[2]][1, 1,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][1, 2,  ], c(1, 0, 0, 0))
  expect_equivalent(arrays[[2]][1, 3,  ], c(0, 1, 0, 0))
  
  expect_equivalent(arrays[[1]][4, 1,  ], rep(0, 4) + 1)
  expect_equivalent(arrays[[1]][4, 2,  ], c(0, 1, 0, 0) + 1)
  expect_equivalent(arrays[[1]][4, 3,  ], c(0, 1, 0, 0) + 1)
  expect_equivalent(arrays[[2]][4, 1,  ], rep(0, 4))
  expect_equivalent(arrays[[2]][4, 2,  ], c(0, 1, 0, 0))
  expect_equivalent(arrays[[2]][4, 3,  ], c(0, 1, 0, 0))
  
  
  testpath <- file.path("fasta_2")
  label_from_csv <- "output_label.csv"
  fx <- function(x = NULL, y = NULL) {
    return(y + 3)
  }
  fy <- function(x = NULL, y = NULL) {
    return(x + 2)
  }
  reshape_xy <- list(x = fx, y = fy)
  gen <- generator_fasta_label_header_csv(path_corpus = testpath, batch_size = 5,
                                          reshape_xy = reshape_xy,
                                          maxlen = 10, step = 10,
                                          vocabulary = c("a", "c", "g", "t", "Z"),
                                          reverse_complement = FALSE, 
                                          vocabulary_label = c("w", "x", "y"),
                                          shuffle_file_order = FALSE,
                                          seed = 1234,
                                          shuffle_input = FALSE,
                                          padding = TRUE,
                                          concat_seq = "ZZ",
                                          target_from_csv = label_from_csv)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[2]][1, 8, ], c(0, 0, 0, 1, 0) + 2) 
  expect_equivalent(arrays[[2]][1, 9, ], c(0, 0, 0, 0, 1) + 2) 
  expect_equivalent(arrays[[2]][1, 10, ], c(0, 0, 0, 0, 1) + 2) 
  
  expect_equivalent(arrays[[2]][4, 3, ], c(1, 0, 0, 0, 0) + 2) 
  expect_equivalent(arrays[[2]][4, 4, ], c(1, 0, 0, 0, 0) + 2) 
  
  expect_equivalent(arrays[[1]][1, ], 1:4 + 3)
  expect_equivalent(arrays[[1]][2, ], 1:4 + 3)
  expect_equivalent(arrays[[1]][3, ], 1:4 + 3)
  expect_equivalent(arrays[[1]][4, ], 11:14 + 3)
  expect_equivalent(arrays[[1]][5, ], 11:14 + 3)
  
  arrays <- gen()
  
  expect_equivalent(arrays[[2]][1, 8, ], c(1, 0, 0, 0, 0) + 2) 
  expect_equivalent(arrays[[2]][2, 3, ], c(0, 1, 0, 0, 0) + 2)
  
  expect_equivalent(arrays[[1]][1, ], 11:14 + 3)
  expect_equivalent(arrays[[1]][5, ], 11:14 + 3)
  
  
  
})
