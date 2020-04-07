library( ggplot2 )
source( "../shared-scripts/plotting/mytheme.R" )

argv <- commandArgs( trailingOnly = TRUE )

concfile <- argv[1]
controlfile <- argv[2]
id <- argv[3]
outfile <- argv[4]

d <- read.table( concfile, header=TRUE )

if( id == "lang"){
  
  lookup_classes <- c(xh="Xhosa",
                      ta="Tagalog",
                      hi="Hiligaynon",
                      la="Latin",
                      pd="Plautdietsch",
                      me="M. English",
                      en="English 2",
                      en2="English 2",
                      self="English")
  
  d <- d[d$class != "en" & d$class != "en2" & d$class != "en_t2", ]
  levs <- c("English","M. English","Latin","Plautdietsch","Tagalog","Hiligaynon", "Xhosa")

  
} else if (id=="peps"){
  
  lookup_classes <- c(hiv="HIV",
                      hepb="Hep. B",
                      hepc="Hep. C",
                      zika="Zika",
                      mal="Malaria",
                      lis="Listeria",
                      hcmv="HCMV",
                      ebola="Ebola",
                      self="Self",
                      vac="Vaccinia")
  
  levs <- c("Self", "Ebola", "HCMV", "Hep. B","Hep. C","HIV", "Listeria", "Malaria", "Vaccinia", "Zika")
  
} else {
  stop("Unknown string id. Please choose either 'lang' or 'peps'.")
}


# add control
dcontrol <- read.table( controlfile, header = TRUE )
d <- rbind( d, data.frame( class = "self", conc = mean(dcontrol$concordance), sd=NA, sem=NA))

# lookup language
d$lang <- lookup_classes[as.character(d$class)]
d$lang <- factor( d$lang, levels = levs)

# manual coloring
n.classes <- length(d$class)
fill.colors <- setNames( rep("black", n.classes), d$class )
fill.colors["self"] <- pathogencolors["self"]           # pathogencolors defined in ../scripts/plotting/mytheme.R
if( id == "lang" ){
#   fill.colors["self"] <- pathogencolors["self"]
#   fill.colors["xh"] <- pathogencolors["xh"]
#   fill.colors["la"]<- pathogencolors["la"]
} else {
  fill.colors["self"] <- pathogencolors["self"]
  fill.colors["hiv"] <- pathogencolors["hiv"]
}

p <- ggplot( d, aes(x = lang,
                    y = conc,
                    fill = class)) +
  geom_col(  position="dodge" ) +
  scale_y_continuous( limits=c(0,1), expand=c(0,0), breaks=seq(0.5,1,by=0.1)) +
  coord_cartesian( ylim = c(0.48,1)) +
  labs( x = "", y = "concordance") +
  mytheme +
  scale_fill_manual( values=fill.colors ) +
  guides( fill = FALSE, colour=FALSE ) +
  theme(
    axis.title.x= element_blank(),
    axis.line.x=element_blank(),
    axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)
  ) 


# p2 <- ggplot( d, aes(x = lang,
#                     y = conc,
#                     fill = class)) +
#   geom_bar( stat = "identity" ) +
#   scale_y_continuous( limits=c(0,1), expand=c(0,0)) +
#   labs( x = "", y = "class concordance") +
#   mytheme +
#   scale_fill_manual( values=fill.colors ) +
#   guides( fill = FALSE ) +
#   theme(
#     axis.title.y= element_blank(),
#     axis.line.y=element_blank()
#   ) + coord_flip()
  

pwidth <- nrow(d)*0.4
p <- set_panel_size(p, width=unit(pwidth,"cm"), height=unit(2.5,"cm")) 
ggsave(outfile, plot=p, width=pwidth + 1.5, height=4.9, units="cm")
