options("install.lock"=FALSE)
install.packages("wordcloud")
library(wordcloud)
library(openxlsx)
sig=read.xlsx('GOenrichment table simple30.xlsx',2)
table(sig$module)

##make them a word
brown2=unlist(strsplit(sig$term.name[61:90]," "))
coral2=unlist(strsplit(sig$term.name[121:150]," "))
firebrick4=unlist(strsplit(sig$term.name[1:30]," "))
lightsteelblue1=unlist(strsplit(sig$term.name[31:60]," "))
orangered3=unlist(strsplit(sig$term.name[91:120]," "))

brown2=gsub("of",NA,brown2)
coral2=gsub("of",NA,coral2)
firebrick4=gsub("of",NA,firebrick4)
lightsteelblue1=gsub("of",NA,lightsteelblue1)
orangered3=gsub("of",NA,orangered3)

brown2=gsub("to",NA,brown2)
coral2=gsub("to",NA,coral2)
firebrick4=gsub("to",NA,firebrick4)
lightsteelblue1=gsub("to",NA,lightsteelblue1)
orangered3=gsub("to",NA,orangered3)

brown2=gsub("in",NA,brown2)
coral2=gsub("in",NA,coral2)
firebrick4=gsub("in",NA,firebrick4)
lightsteelblue1=gsub("in",NA,lightsteelblue1)
orangered3=gsub("in",NA,orangered3)


A_tab=table(brown2)
B_tab=table(coral2)
C_tab=table(firebrick4)
D_tab=table(lightsteelblue1)
E_tab=table(orangered3)
pal=brewer.pal(8,"Accent")
x11()
wordcloud(names(A_tab),freq=A_tab,min.freq=1,colors='Brown',random.order = F)
wordcloud(names(B_tab),freq=B_tab,min.freq=1,colors='Coral',random.order = F)
wordcloud(names(C_tab),freq=C_tab,min.freq=1,colors='Firebrick',random.order = F)
wordcloud(names(D_tab),freq=D_tab,min.freq=1,colors='Lightsteelblue',random.order = F)
wordcloud(names(E_tab),freq=E_tab,min.freq=1,colors='Orangered',random.order = F)

?wordcloud
