# ------------------------------------------------------------------------------
# Objective:
#      Generate 'n' random 'k'-mers by sampling the 20 proteogenic
#      amino acids with replacement
#      ("A","R","N","D","C","Q","E","G","H","I",
#       "L","K","M","F","P","S","T","W","Y","V")
# ------------------------------------------------------------------------------



# Getting started
# ------------------------------------------------------------------------------
# Clear workspace and Load libraries
rm(list=ls())
library('tidyverse')

# Define the 20 proteogenic amino acids
aa = c("A","R","N","D","C","Q","E","G","H","I",
       "L","K","M","F","P","S","T","W","Y","V")

# Set tibble to hold results
results = tibble(sampler = character(), run_no = numeric(), secs = numeric())

# Timing settings
size = 10
n    = 1e5
k    = 9



# Make samplers
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#      _____                       _             __ 
#     / ____|                     | |           /_ |
#    | (___   __ _ _ __ ___  _ __ | | ___ _ __   | |
#     \___ \ / _` | '_ ` _ \| '_ \| |/ _ \ '__|  | |
#     ____) | (_| | | | | | | |_) | |  __/ |     | |
#    |_____/ \__,_|_| |_| |_| .__/|_|\___|_|     |_|
#                           | |                     
#                           |_|                     
# ------------------------------------------------------------------------------

# Sampler_1: The 'math' way
sampler_1 = function(n, k, chars){

  kmers = vector(mode = "character", length = n)
  for( i in 1:n ){
    kmer = vector(mode = "character", length = k)
    for( j in 1:k ){
      kmer[j] = sample(x = chars, size = 1)
    }
    kmers[i] = paste(kmer, collapse = '')
  }
  return(kmers)

}

# Run sampler_1
for( i in 1:size ){
  time_this_run = system.time( sampler_1(n = n, k = k, chars = aa))
  results = results %>%
    bind_rows(tibble(sampler = 'sampler_1',
                     run_no  = i,
                     secs    = unname(time_this_run[3])))
}



# ------------------------------------------------------------------------------
#      _____                       _             ___  
#     / ____|                     | |           |__ \ 
#    | (___   __ _ _ __ ___  _ __ | | ___ _ __     ) |
#     \___ \ / _` | '_ ` _ \| '_ \| |/ _ \ '__|   / / 
#     ____) | (_| | | | | | | |_) | |  __/ |     / /_ 
#    |_____/ \__,_|_| |_| |_| .__/|_|\___|_|    |____|
#                           | |                       
#                           |_|                      
# ------------------------------------------------------------------------------

# But we all know, not to do for loops in R - Do applys!
sampler_2 = function(n, k, chars){

  kmers_mat = sapply(1:n, function(i){
    kmer = sapply(1:k, function(j){
      return( sample(x = chars, size = 1) )
    })
    return(kmer)
  })
  kmers_vec = apply(kmers_mat, 2, function(x){ paste(x, collapse = '') })
  return(kmers_vec)

}

# Run sampler_2
for( i in 1:size ){
  time_this_run = system.time( sampler_2(n = n, k = k, chars = aa))
  results = results %>%
    bind_rows(tibble(sampler = 'sampler_2',
                     run_no  = i,
                     secs    = unname(time_this_run[3])))
}



# ------------------------------------------------------------------------------
#     _____                       _             ____  
#    / ____|                     | |           |___ \ 
#   | (___   __ _ _ __ ___  _ __ | | ___ _ __    __) |
#    \___ \ / _` | '_ ` _ \| '_ \| |/ _ \ '__|  |__ < 
#    ____) | (_| | | | | | | |_) | |  __/ |     ___) |
#   |_____/ \__,_|_| |_| |_| .__/|_|\___|_|    |____/ 
#                           | |
#                           |_|
# ------------------------------------------------------------------------------

# Actually, we should probably do something completely different
sampler_3 = function(n, k, chars){

  # First we generate one long vector with all the samples residues
  smpl_chars = sample(chars, size = n*k, replace = TRUE)
  
  # Then we collapse into one long string
  smpl_string = paste(smpl_chars, collapse = "")
  
  # Now we generate indices corresponding to extracting every 9 characters
  to_index   = seq(k, k*n, by = k)
  from_index = to_index - (k - 1)
  
  # and extract kmers using sub string splitting
  kmers = str_sub(smpl_string, start = from_index, end = to_index)
  
  # Done!
  return(kmers)
}

# Run sampler_3
for( i in 1:size ){
  time_this_run = system.time( sampler_3(n = n, k = k, chars = aa))
  results = results %>%
    bind_rows(tibble(sampler = 'sampler_3',
                     run_no  = i,
                     secs    = unname(time_this_run[3])))
}



# Look at the results
# ------------------------------------------------------------------------------
# Let's visualise it!

# Density plot
results %>%
  ggplot(aes(x = secs, fill = sampler)) +
  geom_density(adjust = 0.5, alpha = 0.5) +
  theme_bw()

# Violin plot
results %>%
  ggplot(aes(x = sampler, y = secs, fill = sampler)) +
  geom_violin(adjust = 0.5, alpha = 0.5, scale = 'width') +
  theme_bw()

# By which factor does each sampler compare?
results %>%
  group_by(sampler) %>%
  summarise(secs_mean = mean(secs)) %>%
  mutate(index_to_min = secs_mean/min(secs_mean))
