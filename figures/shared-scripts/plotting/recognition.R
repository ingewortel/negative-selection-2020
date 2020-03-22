source("../shared-scripts/plotting/mytheme.R")
library(scales)

argv <- commandArgs( trailingOnly = TRUE )
alph_size <- as.numeric( argv[1] ) # the alphabet size; 20 for AAs, 27 (letters + _ ) for latin.
filesuffix <- argv[2]

#Build a matrix of all 64 possible combinations of 0/1 for the 6 positions (mismatch/match):
m <- 6
matches <- c(0,1)

for(p in 2:m){
  n <- nrow(as.matrix(matches))
  matches <- rbind(cbind(as.matrix(matches),rep(0,n)),cbind(as.matrix(matches),rep(1,n)))
}

#Function to determine max length of 1-stretches from bin vector
stretch.length <- function(bin.vec){
  #find stretches
  s <- rle(bin.vec)
  
  #check if there are any matches, then subset the stretches of '1'
  if(sum(bin.vec)>0){
    s1 <- s$lengths[s$values==1]
    
    #return max
    max(s1)
  } else {
    #return 0
    0
  }
  
  
}

#Now take into account the number of possibilities for a mismatch. 
#Positions with a zero (mismatch) can have 19 amino acids, positions with a 1 (match) only 1.
possibilities <- matches
possibilities[possibilities==0] <- alph_size - 1

#Make vectors: kvec with the max stretch for each of the 64 patterns, num the number of sequences that fit the pattern.
#Make a dataframe of those two.
kvec <- apply(matches,1,stretch.length)
num <- apply(possibilities,1,prod)
data <- data.frame(k=kvec,num=num)

for(p in 1:m){
  pvec <- matches[,p]==1
  data <- cbind(data,pvec)
  names(data)[p+2]<-paste0("p",p)
}
total.binders <- function(k,data){
  sum(data$num[data$k>=k])
}
pos.specificity <- function(data,k,pos){
  100*sum(data$num[data$k>=k & data[,2+pos]==TRUE])/total.binders(k,data)
}
pos.entropy <- function(data,k,pos){
  nseqs <- total.binders(k,data)
  en <- ((1/log(2))*(alph_size - 1 )/(2*nseqs))
  
  matches <- sum(data$num[data$k>=k & data[,2+pos]==TRUE])
  mismatches <- sum(data$num[data$k>=k & data[,2+pos]==FALSE])
  
  f.mismatch <- (mismatches/(alph_size - 1 ))/nseqs
  f.match <- (matches/nseqs)
  if(f.mismatch==0){
    Hi <- -f.match*log2(f.match)
  } else {
    Hi <- -((alph_size - 1 )*(f.mismatch*log2(f.mismatch))+(f.match*log2(f.match)))
  }

  Ri <- log2(alph_size) -(Hi)
  Ri
}



pos.data <- data.frame()
for(k in 1:6){
  for(p in 1:m){
    pos.data <- rbind(pos.data, 
                      data.frame(k=paste0("k = ",k),p=p,spec=pos.entropy(data,k,p)))
  }
}

pos.data$k <- factor(pos.data$k,levels=sort(levels(pos.data$k),decreasing=TRUE))

plot <- ggplot(pos.data, aes(x=p,
                             y=spec,
                             colour=k,
                             group=k))+
  geom_point(size=1.2,show.legend=FALSE)+
  geom_line(show.legend=FALSE)+
  labs(x=paste0("Position in ",m,"-mer"),
       y="Information (bits)")+
  scale_colour_manual(values=kcolors)+
  scale_x_continuous(limits=c(1,6),breaks=seq(1,6),labels=seq(1,6))+
  #scale_y_continuous(limits=c(0,100))+
  mytheme + theme(
    legend.position="right",
    axis.title.y=element_text(hjust=0.5)
  )
ggsave(paste0("plots/Fig-specificity",filesuffix,".pdf"),height=4.1,width=3.7,units="cm")


knums <- sapply(seq(1,6), total.binders, data)

num.data <- data.frame(k=paste0("t = ",(seq(1,6))),k2=as.character(seq(1,6)),n=knums)

plotcolors <- c("1"="gray90","2"="gray90","3"="gray90", "4"="lightskyblue2",
"5"="gray90", "6"="gray90") 

plot2 <- ggplot(num.data,aes(x=k,
                             y=n,
                             fill=k,
                             colour=k2))+
  geom_bar(stat="identity",show.legend=FALSE,size=0,width=0.70, aes(fill=k2))+
  labs(x="k",
       y="# matching peptides") +
  scale_fill_manual(values=plotcolors)+
  scale_colour_manual(values=plotcolors)+
  scale_y_log10(expand=c(0,0),breaks=c(1,1e3,1e6),labels= trans_format("log10", math_format(10^.x))) +
  coord_flip()+
  mytheme + theme(
    axis.line=element_blank(),
    axis.ticks=element_blank(),
    axis.title=element_blank(),
    axis.text.x=element_blank()
  )
ggsave(paste0("plots/Fig-crossreactivity",filesuffix,".pdf"),height=4.6,width=4,units="cm")

print( pos.data )
print( num.data )
