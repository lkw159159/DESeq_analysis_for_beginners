install.packages("wordcloud")
install.packages('data.table')
library(wordcloud)


# Read data
Up_res=data.table::fread('UP.txt') %>% data.frame
Down_res=data.table::fread('DN.txt') %>% data.frame

# Check the data
head(Up_res)

Up_res$Category %>% table
Down_res$Category %>% table

##etc

# Select the significantly enriched Terms

## Check the significance
Up_res$PValue<0.05
Down_res$PValue<0.05

## select only significant terms
Up_res_sig=Up_res[Up_res$PValue<0.05,]
Down_res_sig=Down_res[Down_res$PValue<0.05,]


# make them into several word

## check the data
Up_res_sig$Term

## divide by "~"
strsplit(Up_res_sig$Term,"~") 

## select only last value
strsplit(Up_res_sig$Term,"~") %>% sapply(function(x) x[2])

## divide by each word
Up_words=strsplit(Up_res_sig$Term,"~") %>% sapply(function(x) x[2]) %>% 
  strsplit(" ") %>% unlist

## same as Down
Down_words=strsplit(Down_res_sig$Term,"~") %>% sapply(function(x) x[2]) %>% 
  strsplit(" ") %>% unlist

#check the data
Up_words
Down_words

## exclude not-wanted words
Up_words %>% gsub('of',NA,.) %>% gsub('to',NA,.) %>% gsub('in',NA,.)
Down_words %>% gsub('of',NA,.) %>% gsub('to',NA,.) %>% gsub('in',NA,.)

## exclude the NA value
Up_words=Up_words %>% gsub('of',NA,.) %>% gsub('to',NA,.) %>% gsub('in',NA,.) %>% na.omit
Down_words=Down_words %>% gsub('of',NA,.) %>% gsub('to',NA,.) %>% gsub('in',NA,.) %>% na.omit



# make a wordcloud
Up_tab=table(Up_words)
Down_tab=table(Down_words)

Up_tab



x11()

wordcloud(names(Up_tab),freq=Up_tab,min.freq=1,colors='brown',random.order = F)

wordcloud(names(Down_tab),freq=Down_tab,min.freq=1,colors='blue',random.order = F)

wordcloud(names(Down_tab),freq=Down_tab,min.freq=3,colors='blue',random.order = F)

