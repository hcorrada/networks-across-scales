library(tidyverse)

adj <- read_delim("index.adj", skip=29, col_names=FALSE, delim=' ') %>%
    select(adj=1) %>%
    filter(str_detect(adj, "[-_.0-9]", negate=TRUE)) %>%
    mutate(first_unigram=str_sub(adj,1,1)) %>%
    mutate(first_digram=str_sub(adj,1,2))

 nouns <- read_delim("index.noun", skip=29, col_names=FALSE, delim =' ') %>%
     select(noun=1) %>%
     filter(str_detect(noun,"[-'_.0-9]", negate=TRUE)) %>%
     filter(str_length(noun) > 4) %>%
     mutate(first_unigram=str_sub(noun,1,1)) %>%
     mutate(first_digram=str_sub(noun,1,2)) %>%
     sample_n(100)


nouns %>%
    left_join(adj, by="first_digram") %>%
    group_by(noun) %>%
    sample_n(1) %>%
    ungroup() %>%
    mutate(name = str_c(adj, " ", noun)) %>%
    select(name) %>%
    sample_n(20)
