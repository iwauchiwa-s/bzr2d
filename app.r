library(shiny)
library(plotly)
library(shinyWidgets)

# Define UI
ui <- shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("2-D Belousov-Zhabotinsky Reaction Simulator"),
  
  # Sidebar with sliders
  sidebarPanel(
    
    tags$h5("Caution! Selecting large numbers for grid and steps will take a long integration time. It will take about a few minutes for 100x100 grids and 100 steps."),
    
    actionButton("sim", "Simulate"),
    
    tags$h4("Parameter Settings"),
    
    sliderInput("ngrid", label = h5("Grid number:"), min = 10, 
                max = 100, value = 20),
    
    sliderInput("nstep", label = h5("Integration steps:"), min = 10, 
                max = 100, value = 30),
    
    sliderInput("rlxp", label = h5("Diffusion parameter:"), min = 0, 
                max = 1, value = 1),
    
    sliderInput("racp", label = h5("Chemial reaction parameter:"), min = 0, 
                max = 1, value = 1),
    
    br(),
    
  ),
  
  # Show plots
  mainPanel(
  	tags$h3("A"),
    plotlyOutput("graphA"),
  	tags$h3("B"),
    plotlyOutput("graphB"),
  	tags$h3("C"),
    plotlyOutput("graphC"),
  	tags$h3("A-B-C phase space"),
    plotlyOutput("graphD"),
    tags$h3("A, B, C at final step"),
    plotlyOutput("graphSa"),
    plotlyOutput("graphSb"),
    plotlyOutput("graphSc")
  )
  
))

server <- function(input, output, session) {
  
  mySim <- eventReactive(input$sim, {
    
    # define data length
    n <- input$ngrid
    ni <- n
    nj <- n
    nn <- ni * nj
    
    # relaxation for diffusion (0 to 1) 0:no diffusion
    rlx <- input$rlxp
    
    # chemical constants
    ca <- input$racp
    cb <- input$racp
    cc <- input$racp
    
    # integration cycle number
    nc <- input$nstep
    
    # define 2-D data and initial random values
    a <- array(runif(nn,0,1), dim = c(ni,nj))
    b <- array(runif(nn,0,1), dim = c(ni,nj))
    c <- array(runif(nn,0,1), dim = c(ni,nj))
    
    # working matrix
    aw <- array(0, dim = c(ni,nj))
    bw <- array(0, dim = c(ni,nj))
    cw <- array(0, dim = c(ni,nj))
    
    # prepare data frame for matrix a,b,c
    t <- c(1:nc+1)
    df <- 
    p <- c(1:ni)
    q <- c(1:nj)
    df <- expand.grid(p=p,q=q)
    
    # add data to data frame
    da <- a
    db <- b
    dc <- c
    dim(da) <- c(ni*nj,1)
    dim(db) <- c(ni*nj,1)
    dim(dc) <- c(ni*nj,1)
    
    df$sa <- da
    df$sb <- db
    df$sc <- dc
    df$t <- 1
    adf <- df
    
    
    for (istep in 1:nc) { # integration steps
      
      for (i in 1:ni) {
        for (j in 1:nj) {
            # diffusion (local 3x3 grid box averaging)
            a1 <- 0
            b1 <- 0
            c1 <- 0
            # 3x3 grid summation
            for (x in (i-1):(i+1)) {
              for (y in (j-1):(j+1)) {
                  a1 <- a1 + a[ (( x - 1 ) %% ni ) + 1, (( y - 1 ) %% nj ) + 1]
                  b1 <- b1 + b[ (( x - 1 ) %% ni ) + 1, (( y - 1 ) %% nj ) + 1]
                  c1 <- c1 + c[ (( x - 1 ) %% ni ) + 1, (( y - 1 ) %% nj ) + 1]
              }
            }
            # averaging 
            a1 <- a1 / 9
            b1 <- b1 / 9
            c1 <- c1 / 9
            
            # relaxation effect for diffusion process
            a2 <- a[i,j] + ( a1 - a[i,j] ) * rlx
            b2 <- b[i,j] + ( b1 - b[i,j] ) * rlx
            c2 <- c[i,j] + ( c1 - c[i,j] ) * rlx
            
            # reactions
            a3 <- a2 + a2 * ( ca * b2 - cc * c2 )
            b3 <- b2 + b2 * ( cb * c2 - ca * a2 )
            c3 <- c2 + c2 * ( cc * a2 - cb * b2 )
            
            # constrain within 0-1
            if ( a3 < 0 ) {
              a3 <- 0
            }
            if ( a3 > 1 ) {
              a3 <- 1
            }
            if ( b3 < 0 ) {
              b3 <- 0
            }
            if ( b3 > 1 ) {
              b3 <- 1
            }
            if ( c3 < 0 ) {
              c3 <- 0
            }
            if ( c3 > 1 ) {
              c3 <- 1
            }
            
            # store change values at grid [i,j] into working matrix
            aw[i,j] <- a3
            bw[i,j] <- b3
            cw[i,j] <- c3
        }
      }
      # change matrix data after each step calculation
      a <- aw
      b <- bw
      c <- cw
      
    # add data to data frame
    da <- a
    db <- b
    dc <- c
    dim(da) <- c(ni*nj,1)
    dim(db) <- c(ni*nj,1)
    dim(dc) <- c(ni*nj,1)
    
    df$sa <- da
    df$sb <- db
    df$sc <- dc
    df$t <- istep + 1
    adf <- rbind(adf, df)

    } ## sim cycle nc
    ## end of simulation
    
    return(list(adf,a,b,c))
  })
  
  # make graph
  output$graphA <- renderPlotly({
    fig <- plot_ly(mySim()[[1]], x = ~p, y = ~q, z = ~sa, frame = ~t,
                   marker = list(color = ~sa, colorscale = "Rainbow", cmin = 0, cmax = 1,size = 2, symbol = "circle", showscale = TRUE))
    fig <- fig %>% add_markers()
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'X'),
                                       yaxis = list(title = 'Y'),
                                       zaxis = list(title = 'Conc.',range = list(0.0, 1.0))))
  }) ## renderPlotly

  output$graphB <- renderPlotly({
    fig <- plot_ly(mySim()[[1]], x = ~p, y = ~q, z = ~sb, frame = ~t,
                   marker = list(color = ~sb, colorscale = "Rainbow", cmin = 0, cmax = 1,size = 2, symbol = "circle", showscale = TRUE))
    fig <- fig %>% add_markers()
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'X'),
                                       yaxis = list(title = 'Y'),
                                       zaxis = list(title = 'Conc.',range = list(0.0, 1.0))))
  }) ## renderPlotly

  output$graphC <- renderPlotly({
    fig <- plot_ly(mySim()[[1]], x = ~p, y = ~q, z = ~sc, frame = ~t,
                   marker = list(color = ~sc, colorscale = "Rainbow", cmin = 0, cmax = 1,size = 2, symbol = "circle", showscale = TRUE))
    fig <- fig %>% add_markers()
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'X'),
                                       yaxis = list(title = 'Y'),
                                       zaxis = list(title = 'Conc.',range = list(0.0, 1.0))))
  }) ## renderPlotly
  
  output$graphD <- renderPlotly({
    fig <- plot_ly(mySim()[[1]], x = ~sa, y = ~sb, z = ~sc, frame = ~t,
                   marker = list(color = ~sa, colorscale = "Rainbow", cmin = 0, cmax = 1,size = 2, symbol = "circle", showscale = TRUE))
    fig <- fig %>% add_markers()
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'A',range = list(0.0, 1.0)),
                                       yaxis = list(title = 'B',range = list(0.0, 1.0)),
                                       zaxis = list(title = 'C',range = list(0.0, 1.0))))
  }) ## renderPlotly
  
  output$graphSa <- renderPlotly({
	fig <- plot_ly(showscale = TRUE)
	fig <- fig %>% add_surface(z = mySim()[[2]], opacity = 0.8, colorscale = "Rainbow", cmin = 0, cmax = 1)
    fig <- fig %>% layout(scene = list(zaxis = list(title = 'A conc.',range = list(0.0, 1.0))))
  }) ## renderPlotly

  output$graphSb <- renderPlotly({
	fig <- plot_ly(showscale = TRUE)
	fig <- fig %>% add_surface(z = mySim()[[3]], opacity = 0.8, colorscale = "Rainbow", cmin = 0, cmax = 1)
    fig <- fig %>% layout(scene = list(zaxis = list(title = 'B conc.',range = list(0.0, 1.0))))
  }) ## renderPlotly

  output$graphSc <- renderPlotly({
	fig <- plot_ly(showscale = TRUE)
	fig <- fig %>% add_surface(z = mySim()[[4]], opacity = 0.8, colorscale = "Rainbow", cmin = 0, cmax = 1)
    fig <- fig %>% layout(scene = list(zaxis = list(title = 'C conc.',range = list(0.0, 1.0))))
  }) ## renderPlotly


    
} ## server

shinyApp(ui, server)
