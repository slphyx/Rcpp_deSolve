#
# saralamba@gmail.com
#

library(shiny)
library(deSolve)
#if you got this error message '"unlock_solver" not resolved from current namespace (deSolve)'
#use these two lines instead of 'library(deSolve)'

#library.dynam.unload("deSolve", libpath=paste(.libPaths()[1], "//deSolve", sep="")) 
#library.dynam("deSolve", package="deSolve", lib.loc=.libPaths()[1]) 

library(Rcpp)

sourceCpp("SIR.cpp")

SolveODE <- function(parms){
  initial <- c(S = 0.9, I = 0.1, R = 0)
  times <- seq(0,300,by = 0.1)
  
  out <- ode(y = initial, times = times, func = SIRmodel, parms = parms)
  return(out)
}

ui <- fluidPage(
   
   titlePanel("SIR model"),
   
   sidebarLayout(
      sidebarPanel(
         sliderInput("birth",
                     "birth rate",
                     min = 0,
                     max = 1,
                     value = 0.001, step = 0.0001),
         
         sliderInput("death",
                     "death rate",
                     min = 0,
                     max = 1,
                     value = 0.001, step = 0.0001),
         
         sliderInput("recovery",
                     "recovery rate",
                     min = 0,
                     max = 1,
                     value = 0.005, step = 0.0001),
         
         sliderInput("beta",
                     "beta",
                     min = 0,
                     max = 10,
                     value = 0.1, step = 0.01)
         
      ),
      
      mainPanel(
         plotOutput("graphs")
      )
   )
)

server <- function(input, output) {

   parms <- reactive(c(birth = input$birth,
                       death = input$death,
                       recovery = input$recovery,
                       beta = input$beta))

   
   odeout <- reactive(SolveODE(parms()))
   
   plotgraphs <- function(){
     out <- odeout()
     plot(out)
   }
  
   output$txt <- renderText(head(odeout()))   
   output$graphs <- renderPlot({plotgraphs()})
}


shinyApp(ui = ui, server = server)

