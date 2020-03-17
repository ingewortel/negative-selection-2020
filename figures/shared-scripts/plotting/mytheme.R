library(ggplot2)
library(grid)

set_panel_size <- function(p=NULL, g=ggplotGrob(p), width=unit(3, "cm"), height=unit(3, "cm")){
  panel_index_w<- g$layout$l[g$layout$name=="panel"]
  panel_index_h<- g$layout$t[g$layout$name=="panel"]
  g$widths[[panel_index_w]] <- width
  g$heights[[panel_index_h]] <- height
  class(g) <- c("fixed", class(g), "ggplot")
  g
}

#General plotting theme
mytheme <-  theme_classic() +
  theme(
    line=element_line(size=2),
    text=element_text(size=9),
    axis.text=element_text(size=9),
    legend.position="top",
    legend.title=element_blank(),
    axis.line.x=element_line(size=0.25),
    axis.line.y=element_line(size=0.25),
    axis.ticks=element_line(size=0.25),
    plot.margin = unit(c(0.3,0.3,0.3,0.3),"cm")#trbl
  )

#Theme for matrix plots:
matrixtheme <-  theme_classic() +
  theme(
    legend.title=element_blank(),
    axis.line.x=element_blank(),
    axis.line.y=element_blank(),
    panel.border=element_blank(),
    axis.ticks=element_blank(),
    legend.position="bottom"
  )

#Function to generate color gradient for plotting
plot.colors <- function(catnum,ctrlcat=NULL){
  library(RColorBrewer)
  if(is.null(ctrlcat)){
    if(catnum<3){
      colorvec <- brewer.pal(5,"Reds")
      colorvec <- colorvec[4:(catnum+3)]
    } else {
      colorvec <- brewer.pal(catnum,"Reds")
    }
  } else {
    if((catnum-1)<3){
      colorvec <- brewer.pal(5,"Reds")
      colorvec <- colorvec[4:(catnum+2)]
    } else {
      colorvec <- brewer.pal(catnum-1,"Reds")
    }
    
    colorvec <- append(colorvec,"#666666",after=(ctrlcat-1))
  }
  colorvec
}
plot.colorsblue <- function(catnum,ctrlcat=NULL){
  library(RColorBrewer)
  if(is.null(ctrlcat)){
    if(catnum<3){
      colorvec <- brewer.pal(5,"Blues")
      colorvec <- colorvec[4:(catnum+3)]
    } else {
      colorvec <- brewer.pal(catnum,"Blues")
    }
  } else {
    if((catnum-1)<3){
      colorvec <- brewer.pal(5,"Blues")
      colorvec <- colorvec[4:(catnum+2)]
    } else {
      colorvec <- brewer.pal(catnum-1,"Blues")
    }
    
    colorvec <- append(colorvec,"#666666",after=(ctrlcat-1))
  }
  colorvec
}


#predefined colorschemes
kcolors <- c("black","darkorange",plot.colors(1),"blue2","forestgreen","gray70")
    names(kcolors) <- paste0("k = ",seq(1,6))
epcolors <- c(plot.colors(1),"gray70",plot.colors(1))
    names(epcolors) <- c("hiv","self","for")
epcolors.t <- alpha(epcolors,0.15)
    names(epcolors.t) <- names(epcolors)
mylinetypes <- c( 	R = 1,
			B1 = 2,
			O = 2,
			B5 = 3 )
shufcolors <- c(plot.colors(1),plot.colors(1),plot.colors(1),"gray70","blue2",plot.colors(1),plot.colors(1),"darkorange",plot.colors(1),"blue")
    names(shufcolors) <- c("Random","Biased","none","both","self","Biased1","Biased2","Biased3","Biased4","Biased5")
pathogencolors <- c(hiv="dodgerblue3",
                    hcmv="blue",
                    hepb="darkorange",
                    hepc="forestgreen",
                    zika="purple",
                    ebola="black",
		    mal="deeppink",
		    lis="deepskyblue",
		    vac="deepskyblue4",
		    xh="dodgerblue3",
		    ta="deeppink",
		    hi="darkorange",
		    la="cornflowerblue",
		    pd="purple",
		    en2="black",
		    en="black",
		    moby="black",
		    me="forestgreen",
		    self="gray60",
		    self2="black",
		    en_t="gray60")
