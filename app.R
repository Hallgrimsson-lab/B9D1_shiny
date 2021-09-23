library(shiny)
library(Morpho)
library(ggplot2)
library(rgl)
library(Rvcg)


# load("/mnt/Hallgrimsson/Users/Rebecca_Green/Michelle Leong/B9d1 embryos/Landmarked files/B9d1_morphs.Rdata")
# #points, ply, cva.som
# 
# ply <- file2mesh("~/shiny/B9D1_morphs/ml_b9d1e115_low_26_06_dec.ply")
# 
# B9d1.landmark.coordinates.May2018 <- read.csv("/mnt/Hallgrimsson/Users/Rebecca_Green/Michelle Leong/B9d1 embryos/Landmarked files/B9d1 landmark coordinates May2018.csv", row.names=1)

# save(cva.som, cvadf, data, LMS, ply, points, diet, diet_genotype, genotype, tail_somite, file = "~/shiny/B9D1_morphs/B9D1_data.Rdata")
load("/srv/shiny-server/B9D1_morphs/B9D1_data.Rdata")


#define new morpho plotting method until he releases it on cran

render <- function(x,...) UseMethod("render")

#' @rdname render
#' @method render meshDist
#' @export
render.meshDist <- function(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,uprange=NULL,tol=NULL,tolcol=NULL,rampcolors=NULL,NAcol=NULL,displace=FALSE,shade=TRUE,sign=NULL,add=FALSE,scaleramp=NULL,...) {
  clost <- x$clost
  dists <- x$dists
  distsOrig <- dists
  colorall <- x$cols
  colramp <- x$colramp
  params <- x$params
  distqual <- x$distqual
  if (!is.null(tolcol))
    tolcol <- colorRampPalette(tolcol)(1)
  if (!add) {
    if (rgl.cur() !=0)
      rgl.clear()
  }
  if (!is.null(from) || !is.null(to) || !is.null(uprange) ||  !is.null(tol)  ||  !is.null(sign) || !is.null(steps) || !is.null(rampcolors) || !is.null(NAcol) || !is.null(tolcol) || !is.null(scaleramp)) {
    neg=FALSE
    colMesh <- x$colMesh
    if(is.null(steps))
      steps <- x$params$steps
    if (is.null(rampcolors))
      rampcolors <- x$params$rampcolors
    if (is.null(NAcol))
      NAcol <- x$params$NAcol
    if (is.null(tolcol))
      tolcol <- x$params$tolcol
    if (is.null(tol))
      tol <- x$params$tol
    if(is.null(sign))
      sign <- x$params$sign
    if (!sign) {
      distsOrig <- dists
      dists <- abs(dists)
    }
    if(is.null(ceiling))
      ceiling <- x$params$ceiling
    if(is.null(uprange))
      uprange <- x$params$uprange
    
    if (is.null(from)) {
      mindist <- min(dists)
      if (sign && mindist < 0 ) {
        from <- quantile(dists,probs=(1-uprange)) 
        neg <- TRUE            
      } else {
        from <- 0
      }             
    }
    if (is.null(scaleramp))
      scaleramp <- x$params$scaleramp
    
    if (from < 0)
      neg <- TRUE
    if (is.null(to))
      to <- quantile(dists,probs=uprange)    
    if(ceiling)
      to <- ceiling(to)
    
    to <- to+1e-10
    #ramp <- blue2green2red(maxseq*2)
    ramp <- colorRampPalette(rampcolors)(steps-1)
    colseq <- seq(from=from,to=to,length.out=steps)
    coldif <- colseq[2]-colseq[1]
    if (neg && sign) {
      
      negseq <- length(which(colseq<0))
      poseq <- steps-negseq
      maxseq <- max(c(negseq,poseq))
      if (scaleramp) {
        ramp <- colorRampPalette(rampcolors)(maxseq*2)
        ramp <- ramp[c(maxseq-negseq+1):(maxseq+poseq)]
        
      }
      else
        ramp <- colorRampPalette(rampcolors)(steps-1)
      distqual <- ceiling(((dists+abs(from))/coldif)+1e-14)
      #distqual[which(distqual < 1)] <- steps+10
    } else if (from > 0) {
      distqual <- ceiling(((dists-from)/coldif)+1e-14)
    } else {
      distqual <- ceiling((dists/coldif)+1e-14)
    }
    distqual[which(distqual < 1)] <- steps+10
    colorall <- ramp[distqual]
    if (!is.null(tol)) {
      if ( length(tol) < 2 ) {
        if (sign) {
          tol <- c(-tol,tol)
        } else {
          tol <- c(0,tol)
        }
      }
      good <- which(abs(dists) < tol[2])
      colorall[good] <- tolcol
    }
    colfun <- function(x){x <- colorall[x];return(x)}
    colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
    colMesh$material$color[is.na(colMesh$material$color)] <- NAcol
    #colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
    colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
  } else {
    if (is.null(tol))
      tol <- x$params$tol
    colramp <- x$colramp
    colMesh <- x$colMesh
  }
  if (is.null(tolcol))
    tolcol <- x$params$tolcol
  
  if (shade)
    shade3d(vcgUpdateNormals(colMesh),specular="black",meshColor="legacy",...)
  if (displace) {
    dismesh <- colMesh
    vl <- dim(colMesh$vb)[2]
    dismesh$vb <- cbind(colMesh$vb,rbind(clost,1))
    dismesh$it <- rbind(1:vl,1:vl,(1:vl)+vl)
    dismesh$material$color <- rbind(colorall,colorall,colorall)
    wire3d(dismesh,lit=FALSE,meshColor="legacy")
  }
  diffo <- ((colramp[[2]][2]-colramp[[2]][1])/2)
  image(colramp[[1]],colramp[[2]][-1]-diffo,t(colramp[[3]][1,-1])-diffo,col=colramp[[4]],useRaster=TRUE,ylab="Distance in mm",xlab="",xaxt="n")
  if (!is.null(tol)) {
    if (sum(abs(tol)) != 0)
      image(colramp[[1]],c(tol[1],tol[2]),matrix(c(tol[1],tol[2]),1,1),col=tolcol,useRaster=TRUE,add=TRUE)
  }
  params <- list(steps=steps,from=from,to=to,uprange=uprange,ceiling=ceiling,sign=sign,tol=tol,rampcolors=rampcolors,NAcol=NAcol,tolcol=tolcol)
  out <- list(colMesh=colMesh,dists=distsOrig,cols=colorall,colramp=colramp,params=params,distqual=distqual,clost=clost)
  
  class(out) <- "meshDist"
  invisible(out)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("B9D1 embryo morphs"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      c("Double click the plot to morph the embryo"),
      tags$br(),
      tags$br(),
      numericInput("xaxis", label = "Dimension on X", value = 1, min = 1, max = 3),
      numericInput("yaxis", label = "Dimension on Y", value = 2, min = 1, max = 3), 
      sliderInput("transparency", label = "Mesh transparency", min = 0, max = 1, step = NULL, value = .3),
      width = 2,
      tags$head(
        tags$style("body {background-color: #1a1a1a; }
                   .well {background-color: #404040;}
                   .control-label {color: white}
                   .irs-min, .irs-max {color: white}
                   .checkbox {color: white}
                   .container-fluid {color: white}
                   pre {
                   width: 250px; 
                   margin: auto;
                   text-align: center;
                   }
                   ")
      )
    ),
    
    mainPanel(
      fluidRow(
        splitLayout(
          plotOutput("scoreplot", dblclick = "plot_click"),
          rglwidgetOutput("embryo_warp")
        ),
        verbatimTextOutput("CVscores")
      )
    , width = 10)
  )
)

server <- function(input, output){
  
  # data <- B9d1.landmark.coordinates.May2018
  # LMS <- data[,-(1:3)]
  # genotype <- data[,1]
  # diet <- data[,3]
  # tail_somite <- data[,2]
  # diet_genotype <- as.factor(paste(diet,genotype, sep="_"))
  # cvadf <- as.data.frame(cva.som$CVscores[,1:3])
  
  output$scoreplot <- renderPlot({
    
    y <- ggplot(data=cvadf, aes(x=cvadf[,input$xaxis], y=cvadf[,input$yaxis]))
    
    y  + geom_point(aes(fill = diet_genotype), colour="black", pch=21,size=I(5), alpha=I(.8)) +
      labs(x= paste0("CV", input$xaxis), y= paste0("CV", input$yaxis)) + 
      guides(fill =guide_legend(title="Series")) +
      theme_dark() + 
      theme(plot.background=element_rect(fill = "#1a1a1a"), 
            legend.background = element_rect(fill = "black", color = NA),
            panel.background = element_rect(fill = '#404040'),
            legend.key = element_rect(color = "#404040", fill = "#404040"),
            legend.title = element_text(color = "white"),
            legend.text = element_text(color = "white"))
    
  })
  
  output$embryo_warp <- renderRglwidget({
    if (names(dev.cur()) != "null device") dev.off()
    pdf(NULL)
    
    #input score
    tmp_scorex <- if(is.null(input$plot_click$x)){0} else{input$plot_click$x}
    tmp_scorey <-  if(is.null(input$plot_click$y)){0} else{input$plot_click$y}
    
    #multiply loadings
    shape.mean <-  tps3d(ply, refmat = points, tarmat = cva.som$Grandm)
    shape.estimate <- tmp_scorex*matrix(cva.som$CVvis[,input$xaxis],nrow(cva.som$Grandm),ncol(cva.som$Grandm)) + tmp_scorey*matrix(cva.som$CVvis[,input$yaxis],nrow(cva.som$Grandm),ncol(cva.som$Grandm))  + cva.som$Grandm
    
    #diagonal view
    #par3d(userMatrix = matrix(c(.22,-.21,.95,0,-.97,-.03,.219,0,-.01,-.97,-.2,0,0,0,0,1),ncol =4,nrow = 4))
    #front view
    par3d(userMatrix = matrix(c(.69,-.2,.7,0,-.72,-.15,.67,0,-.037,-.97,-.23,0,0,0,0,1),ncol =4,nrow = 4))
    par3d(zoom = .7)
    bg3d(color = "#1a1a1a")
    #tps warp to default mesh
    shape.transformed <- tps3d(shape.mean, refmat = cva.som$Grandm, tarmat = shape.estimate)
    # shape.warp <-  plot3d(tps3d(shape.mean, refmat = cva.som$Grandm, tarmat = shape.estimate), col = "grey", alpha = .3, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "")
    mD <- meshDist(shape.mean, shape.transformed, plot = F, scaleramp = F, displace = T, alpha = input$transparency)
    a <- render(mD, displace = T, alpha = input$transparency)
    spheres3d(shape.estimate, radius = .0025, color = "red")
    # if(is.null(input$plot_click$x)){} else{for(i in 1:38) arrow3d(shape.estimate[i,], cva.som$Grandm[i,], type = "lines", col = "red")}
    
    rglwidget()
    
    
  })
  
  output$CVscores <- renderText({
    paste0("CV x score = ",  if(is.null(input$plot_click$x)){0} else{round(input$plot_click$x, digits = 1)}, "\nCV y score = ",  if(is.null(input$plot_click$y)){0} else{round(input$plot_click$y, digits = 1)})
  })
  
}

shinyApp(ui = ui, server = server)

