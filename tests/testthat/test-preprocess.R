context("preprocess")

test_that("Correct vocabulary extraction", {
  
  expect_equal(getVocabulary("ABC"), c("a", "b", "c"))
  expect_equal(getVocabulary("CBA"), c("a", "b", "c"))
  expect_equal(getVocabulary("AAA"), c("a"))
  expect_equal(getVocabulary("012"), c("0", "1", "2"))
  expect_equal(getVocabulary("()/"), c("(", ")", "/"))
  expect_equal(getVocabulary(" a "), c(" ", "a"))
  expect_equal(getVocabulary("abc"), c("a","b","c"))
  expect_equal(getVocabulary(123),c("1","2","3"))
  expect_equal(getVocabulary(c("A","B","C",1)),c("\n","1","a","b","c"))
  expect_equal(getVocabulary("\n"),c("\n"))
  
  expect_error(getVocabulary(""))
  expect_error(getVocabulary())
  
  expect_is(getVocabulary("abc"),"character")
  
  expect_message(getVocabulary("abc", verbose = T))
  expect_silent(getVocabulary("abc"))
})

test_that("Generating semi-redundant chunks", {
  
  expect_is(preprocessSemiRedundant(char = "abcd", maxlen = 2),"list")
  expect_is(preprocessSemiRedundant(char = "abcd", maxlen = 2)$X,"array")
  expect_is(preprocessSemiRedundant(char = "abcd", maxlen = 2)$Y,"matrix")
  
  expect_equivalent(lengths(preprocessSemiRedundant(char = "abcd", maxlen = 2))[1], 16)
  expect_equivalent(lengths(preprocessSemiRedundant(char = "abcd", maxlen = 2, vocabulary = c("a","b","c","d")))[1], 16)
  expect_equivalent(lengths(preprocessSemiRedundant(char = "abcd", maxlen = 2))[2], 8)
  expect_equivalent(preprocessSemiRedundant(char = "abcd", maxlen = 2)$Y, matrix(c(0,0,1,0,0,0,0,1), byrow = TRUE, nrow = 2))           
  expect_equivalent(length(preprocessSemiRedundant(char="abcd", maxlen = 2)),2)
  
  expect_error(preprocessSemiRedundant(char = "abcd", maxlen = ""))
  expect_error(preprocessSemiRedundant(char = "abcd", maxlen = 0))
  expect_error(preprocessSemiRedundant(char = "abcd", vocabulary = ""))
  expect_error(preprocessSemiRedundant(char = "abcd", vocabulary = 0))
  
  expect_message(preprocessSemiRedundant(char = "abcd", maxlen = 2, verbose = T))
  expect_silent(preprocessSemiRedundant(char = "abcd", maxlen = 2))
  
  expect_type(preprocessSemiRedundant(char = "abcd", maxlen = 2)$X, "double")
  expect_type(preprocessSemiRedundant(char = "abcd", maxlen = 2)$Y, "double")
})

test_that("Generating semi-redundant chunks from Fasta files", {
  
  file <- file.path("fasta/a.fasta")
  
  expect_is(preprocessFasta(file),"list")
  expect_is(preprocessFasta(file)$X,"array")
  expect_is(preprocessFasta(file)$Y,"matrix")
  
  expect_equivalent(lengths(preprocessFasta(file))[2], 39660)
  expect_equivalent(length(preprocessFasta(file)),2)
  
  expect_error(preprocessFasta())
  expect_error(preprocessFasta(""))
  
  expect_silent(preprocessFasta(file))
  
  expect_type(preprocessFasta(file)$X, "double")
  expect_type(preprocessFasta(file)$Y, "double")
})

# TODO: test shuffling
test_that("Checking the generator for the Fasta files", {
  
  expect_error(sequenceToArray())
  expect_error(sequenceToArray(""))
  
  testpath <- file.path("fasta_2")
  vocabulary = c("a", "c", "g", "t")
  batch.size <- 5
  maxlen <- 3
  gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE, vocabulary = vocabulary)
  
  arrays <- gen()
  
  expect_equivalent(dim(arrays[[1]])[1], batch.size)
  expect_equivalent(dim(arrays[[1]])[2], maxlen)
  expect_equivalent(dim(arrays[[1]])[3], length(vocabulary))
  expect_equivalent(dim(arrays[[2]])[1], batch.size)
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
  gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE, vocabulary = vocabulary)
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
  gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE, vocabulary = vocabulary)
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
  gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = 4, maxlen = 3, step = 2, showWarnings = FALSE)
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
  gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = 5, maxlen = 3, step = 2, showWarnings = FALSE, vocabulary = c("c", "g", "t"))
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1, ], c(0, 0, 0)) # a 
  expect_equivalent(arrays[[1]][1, 2, ], c(0, 0, 0)) # a
  expect_equivalent(arrays[[1]][1, 3, ], c(1, 0, 0)) # c
  expect_equivalent(arrays[[2]][1, ], c(1, 0, 0)) # c
  
  ####
  # test padding
  gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = 1, maxlen = 10, step = 4, showWarnings = FALSE,
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
  gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = 2, maxlen = 12, step = 1, showWarnings = FALSE,
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
  expect_error(fastaFileGenerator())
  expect_error(fastaFileGenerator(""))
  
  expect_is(fastaFileGenerator(testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE), "function")
  expect_is(gen(), "list")
  expect_is(gen()[[1]], "array")
  expect_is(gen()[[2]], "matrix")
  
  expect_message(fastaFileGenerator(testpath, batch.size = batch.size, maxlen = maxlen,
                                    showWarnings = FALSE, verbose = T))
  expect_silent(fastaFileGenerator(testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE))
  
  # no T in vocabulary
  expect_warning(fastaFileGenerator(testpath, batch.size = batch.size, maxlen = maxlen, 
                                    vocabulary = c("a", "c", "g"), showWarnings = TRUE))
  
  
  expect_type(gen()[[1]], "double")
  expect_type(gen()[[2]], "double")
  
  ############# Test label generator (header) #############
  testpath <- file.path("fasta_2")
  
  gen <- fastaLabelGenerator(corpus.dir = testpath, batch.size = 5, maxlen = 3, step = 2, showWarnings = FALSE, vocabulary = c("a", "c", "g", "t"),
                             reverseComplement = FALSE, labelVocabulary = c("w", "x", "y"))
  
  arrays <- gen()
  
  expect_equivalent(arrays[[1]][1, 1, ], c(1, 0, 0, 0)) # A  
  expect_equivalent(arrays[[1]][1, 2, ], c(1, 0, 0, 0)) # A
  expect_equivalent(arrays[[1]][1, 3, ], c(0, 1, 0, 0)) # C
  expect_equivalent(arrays[[2]][1, ], c(1, 0, 0)) # W 
  
  expect_equivalent(arrays[[1]][5, 1, ], c(1, 0, 0, 0)) # A  
  expect_equivalent(arrays[[1]][5, 2, ], c(1, 0, 0, 0)) # A
  expect_equivalent(arrays[[1]][5, 3, ], c(0, 0, 0, 1)) # T
  expect_equivalent(arrays[[2]][5, ], c(0, 1, 0)) # W 
  
  
  gen <- fastaLabelGenerator(corpus.dir = testpath, batch.size = 5, maxlen = 8, step = 2, showWarnings = FALSE, vocabulary = c("a", "c", "g", "t"),
                             reverseComplement = FALSE, labelVocabulary = c("w", "x", "y"))
  
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
  
  
  gen <- fastaLabelGenerator(corpus.dir = testpath, batch.size = 8, maxlen = 7, step = 2, showWarnings = FALSE, vocabulary = c("a", "c", "g", "t"),
                             reverseComplement = FALSE, labelVocabulary = c("w", "x", "y"))
  
  arrays <- gen()
  
  # go through a/b.fasta once discard samples with target z
  expect_equivalent(arrays[[1]][8, 1, ], c(1, 0, 0, 0)) # A  
  expect_equivalent(arrays[[1]][8, 2, ], c(1, 0, 0, 0)) # A
  expect_equivalent(arrays[[1]][8, 3, ], c(0, 1, 0, 0)) # C
  expect_equivalent(arrays[[2]][8, ], c(1, 0, 0)) # W 
})

############# Test label generator (folder) #############
directories <- c("label_folder/x", "label_folder/y", "label_folder/z")
val <- FALSE
initializeGenerators(directories = directories,
                     format = "fasta",
                     batch.size = 6,
                     maxlen = 2,
                     vocabulary = c("a", "c", "g", "t"),
                     step = 2)

gen <- labelByFolderGeneratorWrapper(val = val, path = directories)
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
batch.size <- 6
maxlen <- 3
step <- 2
gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE, 
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

gen <- fastaLabelGenerator(corpus.dir = testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE, 
                           vocabulary = vocabulary, ambiguous_nuc = "discard", step = step, reverseComplements = FALSE,
                           labelVocabulary = c("X", "Y"))
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
initializeGenerators(directories = directories,
                     format = "fasta",
                     batch.size = 6,
                     maxlen = 3,
                     ambiguous_nuc = "discard",
                     vocabulary = c("a", "c", "g", "t"),
                     reverseComplements = FALSE, 
                     step = 2)

gen <- labelByFolderGeneratorWrapper(val = FALSE, path = directories)
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
batch.size <- 4
maxlen <- 3
step <- 2
gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE, 
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

gen <- fastaLabelGenerator(corpus.dir = testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE, 
                           vocabulary = vocabulary, ambiguous_nuc = "equal", step = step, reverseComplements = FALSE,
                           labelVocabulary = c("X", "Y"))
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
initializeGenerators(directories = directories,
                     format = "fasta",
                     batch.size = 4,
                     maxlen = 3,
                     vocabulary = c("a", "c", "g", "t"),
                     reverseComplements = FALSE, 
                     ambiguous_nuc = "equal",
                     step = 2)

gen <- labelByFolderGeneratorWrapper(val = FALSE, path = directories)
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
batch.size <- 4
maxlen <- 3
step <- 2
gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE, 
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
batch.size <- 4
maxlen <- 3
step <- 20
gen <- fastaFileGenerator(corpus.dir = testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE, 
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
batch.size <- 4
maxlen <- 3
step <- 2
gen <- fastaLabelGenerator(corpus.dir = testpath, batch.size = batch.size, maxlen = maxlen, showWarnings = FALSE, 
                           vocabulary = vocabulary, ambiguous_nuc = "empirical", step = step, reverseComplements = FALSE,
                           labelVocabulary = c("X", "Y"))
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
initializeGenerators(directories = directories,
                     format = "fasta",
                     batch.size = 4,
                     maxlen = 3,
                     vocabulary = c("a", "c", "g", "t"),
                     reverseComplements = FALSE, 
                     ambiguous_nuc = "empirical",
                     step = 2)

gen <- labelByFolderGeneratorWrapper(val = FALSE, path = directories)
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

gen <- fastaFileGenerator(corpus.dir = "/home/rmreches/deepG/tests/testthat/fasta_3",
                          batch.size = 3,
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

gen <- fastaLabelGenerator(corpus.dir = "/home/rmreches/deepG/tests/testthat/fasta_3",
                           batch.size = 3,
                           maxlen = 15,
                           step = 1,
                           labelVocabulary = c("X", "Y"),
                           reverseComplements = FALSE,
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
initializeGenerators(directories = directories,
                     format = "fasta",
                     batch.size = 6,
                     maxlen = 15,
                     ambiguous_nuc = "equal",
                     vocabulary = c("a", "c", "g", "t"),
                     reverseComplements = FALSE, 
                     step = 1)

equal_vector <- rep(0.25, 4)
gen <- labelByFolderGeneratorWrapper(val = FALSE, path = directories)
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

gen <- fastaFileGenerator(corpus.dir = "fasta_3",
                          batch.size = 8,
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

gen <- fastaLabelGenerator(corpus.dir = "fasta_3",
                           batch.size = 8,
                           maxlen = 12,
                           max_iter = 10000,
                           step = 50, 
                           ambiguous_nuc = "empirical",
                           reverseComplements = FALSE,
                           labelVocabulary = c("X", "Y")
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
initializeGenerators(directories = directories,
                     format = "fasta",
                     batch.size = 20,
                     maxlen = 12,
                     ambiguous_nuc = "empirical",
                     vocabulary = c("a", "c", "g", "t"),
                     reverseComplements = FALSE, 
                     step = 1)

nuc_dist_1 <- 1/18*c(8, 2, 3, 5)
nuc_dist_2 <- 1/17*c(3, 2, 6, 6)
gen <- labelByFolderGeneratorWrapper(val = FALSE, path = directories)
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

gen <- fastaFileGenerator(corpus.dir = "fastq",
                          format = "fastq",
                          batch.size = 10,
                          maxlen = 3,
                          max_iter = 10000,
                          vocabulary = c("a", "c", "g", "t"),
                          verbose = FALSE,
                          randomFiles = FALSE,
                          step = 2, 
                          showWarnings = FALSE,
                          seed = 1234,
                          shuffleFastaEntries = FALSE,
                          numberOfFiles = NULL,
                          fileLog = NULL,
                          reverseComplements = FALSE,
                          output_format = "target_right",
                          ambiguous_nuc = "zeros",
                          use_quality_score = TRUE,    
                          proportion_per_file = NULL,
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

gen <- labelByFolderGenerator(corpus.dir = "fastq",
                              format = "fastq",
                              batch.size = 10,
                              maxlen = 3,
                              max_iter = 10000,
                              vocabulary = c("a", "c", "g", "t"),
                              verbose = FALSE,
                              randomFiles = FALSE,
                              step = 2, 
                              showWarnings = FALSE,
                              seed = 1234,
                              shuffleFastaEntries = FALSE,
                              numberOfFiles = NULL,
                              fileLog = NULL,
                              reverseComplements = FALSE,
                              ambiguous_nuc = "zeros",
                              use_quality_score = TRUE,    
                              proportion_per_file = NULL,
                              numTargets = 2,
                              onesColumn = 1,
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

gen <- labelByFolderGenerator(corpus.dir = "read_data",
                              format = "fastq",
                              batch.size = 5,
                              maxlen = 12,
                              max_iter = 10000,
                              vocabulary = c("a", "c", "g", "t"),
                              verbose = FALSE,
                              randomFiles = FALSE,
                              step = 2, 
                              showWarnings = FALSE,
                              seed = 1234,
                              shuffleFastaEntries = FALSE,
                              numberOfFiles = NULL,
                              fileLog = NULL,
                              read_data = TRUE,
                              reverseComplements = FALSE,
                              ambiguous_nuc = "zeros",
                              use_quality_score = TRUE,    
                              proportion_per_file = NULL,
                              numTargets = 2,
                              onesColumn = 1,
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

gen <- fastaFileGenerator(corpus.dir = "fasta_3",
                          format = "fasta",
                          batch.size = 10,
                          maxlen = 5,
                          vocabulary = c("a", "c", "g", "t"),
                          randomFiles = FALSE,
                          step = 4, 
                          showWarnings = FALSE,
                          shuffleFastaEntries = FALSE,
                          reverseComplements = FALSE,
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
initializeGenerators(directories = dir,
                     format = "fasta",
                     batch.size = 15,
                     maxlen = 4,
                     step = 2, 
                     val = FALSE,
                     padding = FALSE,
                     added_label_path = "label.csv",
                     add_input_as_seq = FALSE)

gen <- labelByFolderGeneratorWrapper(val = FALSE, path = dir) 
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

gen <- labelByFolderGeneratorWrapper(val = FALSE, path = dir) 
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

initializeGenerators(directories = c("read_data_2/label_a", "read_data_2/label_b"),
                     format = "fastq",
                     batch.size = 4,
                     maxlen = 12,
                     max_iter = 10000,
                     vocabulary = c("a", "c", "g", "t"),
                     verbose = FALSE,
                     randomFiles = FALSE,
                     step = 1, 
                     showWarnings = FALSE,
                     seed = 1234,
                     shuffleFastaEntries = FALSE,
                     numberOfFiles = NULL,
                     fileLog = NULL,
                     reverseComplements = FALSE, 
                     val = FALSE,
                     ambiguous_nuc = "zero",
                     proportion_per_file = NULL,
                     target_middle = FALSE,
                     read_data = TRUE,
                     use_quality_score = TRUE,
                     padding = FALSE,
                     added_label_path = NULL,
                     skip_amb_nuc = NULL)

gen <- labelByFolderGeneratorWrapper(val = FALSE, path = c("/read_data_2/label_a", "/read_data_2/label_b"))
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
gen <- fastaLabelGenerator(corpus.dir = testpath, batch.size = 5,
                           maxlen = 10, step = 10, showWarnings = FALSE,
                           vocabulary = c("a", "c", "g", "t", "Z"),
                           reverseComplement = FALSE, 
                           labelVocabulary = c("w", "x", "y"),
                           format = "fasta",
                           max_iter = 10000,
                           verbose = FALSE,
                           randomFiles = FALSE,
                           seed = 1234,
                           shuffleFastaEntries = FALSE,
                           numberOfFiles = NULL,
                           fileLog = NULL,
                           ambiguous_nuc = "zero",
                           proportion_per_file = NULL,
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

gen <- fastaFileGenerator(corpus.dir = "fasta_3",
                          format = "fasta",
                          batch.size = 10,
                          maxlen = 5,
                          vocabulary = c("a", "c", "g", "t"),
                          randomFiles = FALSE,
                          step = 4, 
                          showWarnings = FALSE,
                          shuffleFastaEntries = FALSE,
                          reverseComplements = FALSE,
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
initializeGenerators(directories = dir,
                     format = "fasta",
                     batch.size = 15,
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

gen <- labelByFolderGeneratorWrapper(val = FALSE, path = dir) 
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

gen <- fastaFileGenerator(corpus.dir = "fasta_3",
                          batch.size = 5,
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
expect_equivalent(arrays[[2]][[1]], m1)
expect_equivalent(arrays[[2]][[2]], m2)
expect_equivalent(arrays[[2]][[3]], m3)

# 3 targets, target middle cnn

gen <- fastaFileGenerator(corpus.dir = "fasta_3",
                          batch.size = 5,
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
expect_equivalent(arrays[[2]][[1]], m1)
expect_equivalent(arrays[[2]][[2]], m2)
expect_equivalent(arrays[[2]][[3]], m3)

# 3 targets, target middle lstm

gen <- fastaFileGenerator(corpus.dir = "fasta_3",
                          batch.size = 5,
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
expect_equivalent(arrays[[2]][[1]], m1)
expect_equivalent(arrays[[2]][[2]], m2)
expect_equivalent(arrays[[2]][[3]], m3)

# coverage + set learning for label_folder 

directories <- c("coverage_data/x", "coverage_data/y")
val <- FALSE
batch.size <- 18
samples_per_target <- 3
new_batch_size <- batch.size/samples_per_target
path <- directories
voc_len <- 4
maxlen <- 7
reshape_mode <- "time_dist"
set_learning <- list(reshape_mode = reshape_mode,
                     maxlen = maxlen,
                     samples_per_target = samples_per_target)

initializeGenerators(directories = directories,
                     format = "fasta",
                     batch.size = batch.size,
                     maxlen = maxlen,
                     vocabulary = c("a", "c", "g", "t"),
                     step = 4,
                     use_coverage = 1,
                     set_learning = set_learning)

gen <- labelByFolderGeneratorWrapper(val = val, new_batch_size = new_batch_size,
                                     samples_per_target = samples_per_target,
                                     batch.size = batch.size,
                                     path = path, voc_len = voc_len, maxlen = maxlen,
                                     reshape_mode = reshape_mode)
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
batch.size <- 18
samples_per_target <- 3
new_batch_size <- batch.size/samples_per_target
path <- directories
voc_len <- 4
maxlen <- 7
use_coverage <- 17
reshape_mode <- "time_dist"
set_learning <- list(reshape_mode = reshape_mode,
                     maxlen = maxlen,
                     samples_per_target = samples_per_target)

initializeGenerators(directories = directories,
                     format = "fasta",
                     batch.size = batch.size,
                     maxlen = maxlen,
                     vocabulary = c("a", "c", "g", "t"),
                     step = 4,
                     use_coverage = use_coverage,
                     set_learning = set_learning)

gen <- labelByFolderGeneratorWrapper(val = val, new_batch_size = new_batch_size,
                                     samples_per_target = samples_per_target,
                                     batch.size = batch.size,
                                     path = path, voc_len = voc_len, maxlen = maxlen,
                                     reshape_mode = reshape_mode)
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

gen <- gen_rds(rds_folder = "rds", batch_size = 1)
l_x <- list()
l_y <- list()
for (i in 1:40) {
  z <- gen()
  l_x[[i]] <- z[[1]][1,1,1]
  l_y[[i]] <- which.max(z[[2]])
}
expect_equivalent(sort(unlist(l_x)), rep(1:20, each=2)) 
expect_equivalent(sort(unlist(l_y)), rep(1:20, each=2)) 

gen <- gen_rds(rds_folder = "rds", batch_size = 10)
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
gen <- gen_rds(rds_folder = "rds_lm", batch_size = batch_size, target_len = target_len)

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
gen <- gen_rds(rds_folder = "rds_lm", batch_size = batch_size, target_len = target_len)
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

# n-gram fasta 

test_path <- "fasta_2"
n_gram <- 3
gen <- fastaFileGenerator(corpus.dir = test_path,
                          batch.size = 2,
                          maxlen = 2,
                          step = 1, 
                          padding = FALSE,
                          target_len = 6,
                          n_gram = n_gram,
                          n_gram_stride = n_gram)

arrays <- gen()
y <- arrays[[2]]
y_1_n_gram <- apply(y[[1]], 1, which.max)
y_2_n_gram <- apply(y[[2]], 1, which.max)

int_seq <- c(1,1,2)
expect_equivalent(y_1_n_gram[1], 1 + sum(4^((n_gram-1):0) * (int_seq))) # ccg
int_seq <- c(2,3,3)
expect_equivalent(y_2_n_gram[1], 1 + sum(4^((n_gram-1):0) * (int_seq))) # gtt
int_seq <- c(0,0,3)
expect_equivalent(y_1_n_gram[2], 1 + sum(4^((n_gram-1):0) * (int_seq))) # aat
int_seq <- c(3,3,3)
expect_equivalent(y_2_n_gram[2], 1 + sum(4^((n_gram-1):0) * (int_seq))) # ttt

# n-gram rds

n_gram <- 3
gen <- gen_rds(rds_folder = "n_gram.rds",
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

