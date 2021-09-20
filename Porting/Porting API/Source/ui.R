
# Define UI for application that draws a histogram
fluidPage(
  
  shinyjs::useShinyjs(),
  # Application title
  # titlePanel("SENTINEL-1/2/3 DOWNLOADER"),
  tags$head(
    # tags$title(title),  
    tags$script(src="https://maps.googleapis.com/maps/api/js?key=&libraries=places"), 
    tags$script(src="js/myfuncs.js"),  
    tags$script(src="global/js/autocomplete_street.js"), 
    tags$link(rel="stylesheet", type="text/css", href="css/extra.css?v=2")
  ),
  
  sidebarLayout(
    sidebarPanel( width = 6, style="background-color:#ffffff;", 
                  tabsetPanel(id="mytabset",
                    tabPanel("Scihub search", value="searchtab",
                             hr(),
                             fluidRow(
                               column(width=6,  textInput("user", label =NULL,placeholder = "SciHub username")),
                               column(width=6,  passwordInput("password",label =NULL, placeholder = "password"))
                             ),
                             fluidRow(
                               column(width=12, 
                                      searchInput(inputId = "searchBar", label="Search address:",  
                                                  width="100%", 
                                                  placeholder = "Street address ....",  
                                                  btnSearch = icon("search"), btnReset = icon("remove") )   
                                      
                               )
                             ),
                             fluidRow (
                               column(width=5, selectInput("platform", "Platform:", 
                                                           choices = list("", "Sentinel-1", "Sentinel-2"))
                               ),
                               
                               column(width=7, dateRangeInput("daterange",
                                                              "Date range:", start  = "2018-01-01") 
                               )
                             ),
                             
                             fluidRow (
                               column( width=3, selectizeInput("month",
                                                               "Month:", choices=c( months.list ), 
                                                               multiple = T, options = list(placeholder="select month(s)")) ),
                               
                               column( width=3, selectizeInput("weekday",
                                                               "Weekday(s):", choices=format(ISOdate(2000, 1, 3:9), "%A"), 
                                                               multiple = T, options = list(placeholder="select weekday(s)")) ),
                               
                               column( width=2, selectizeInput("day",
                                                               "Day(s):", choices=c('', as.character(1:31)), 
                                                               multiple = T, options = list(placeholder="select day(s)")) ),
                               column(width=4, prettyCheckboxGroup("orbit.dir", "Orbit direction", 
                                                                   choices = orbit.direction, 
                                                                   selected = orbit.direction,
                                                                   inline=T) )
                             ),
                             fluidRow( 
                               conditionalPanel(
                                 condition = "input.platform == 'Sentinel-2'",
                                 column(width=12,
                                     sliderInput("cloud",
                                               "Max Cloud Cover (%):", 10, max = 100, step=1, min = 0  ) 
                               ) ),
                               conditionalPanel(
                                 condition = "input.platform == 'Sentinel-1'",
                                 column(width=4,
                                       selectInput("polarisationmode", multiple=T,
                                           "Polarisation:", choices=polarization.mode  ) 
                                 ),
                                 column(width=2,
                                        selectInput("sensoroperationalmode", multiple=T,
                                                    "Mode:", choices=operationalmode  ) 
                                 ),
                                 column(width=6,
                                        selectInput("producttype", multiple=T,
                                                    "Product type:", choices=product.type  ) 
                                 )
                               
                               ) 
                             ), 
                             actionButton("search", "Search", icon = icon("search"))    
                    ),
                    tabPanel("Results", value="tabletab",
                             
                             downloadLink("downloadData", "Download script to wget selected images")   ,
                             DT::dataTableOutput("mytable")
                    ) 
                  )      
                  
    ),
    
    # Show a plot of the generated distribution
    mainPanel(width = 6,  
              leafletOutput("mymap") )
  )
)