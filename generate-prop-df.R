data("aaindex")

pub_years <- sapply(aaindex, function(i) i[["D"]]) %>% 
  strsplit(split = "., ") %>% 
  sapply(last) %>% 
  gsub(pattern = ")", replacement = "", x = .) %>% 
  as.numeric 

prop_df <- data.frame(ID = sapply(aaindex, function(i) i[["H"]]),
                      desc = sapply(aaindex, function(i) i[["D"]]),
                      year = pub_years,
                      stringsAsFactors = FALSE)

