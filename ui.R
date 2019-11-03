library(shiny)

loadingBar <- tags$div(class="progress progress-striped active",
                       tags$div(class="bar", style="width: 100%;"))
# Code for loading message
loadingMsg <- tags$div(class="modal", tabindex="-1", role="dialog", 
                       "aria-labelledby"="myModalLabel", "aria-hidden"="true",
                       tags$div(class="modal-header",
                                tags$h3(id="myModalHeader", "Loading...")),
                       tags$div(class="modal-footer",
                                loadingBar))
# The conditional panel to show when shiny is busy
loadingPanel <- conditionalPanel(paste("input.goButton > 0 &&", 
                                       "$('html').hasClass('shiny-busy')"),
                                 loadingMsg)

shinyUI(navbarPage(
  id="index",
  title=HTML('<div style="position: relative ; right: 10px ;  bottom: 5px; line-height:20px ; color:black ; font-weight:bold">SpadeR<img src="SpadeR.jpg" width = "53px" height = "53px"></img></div>'),
  windowTitle = "SpadeR", 
  tabPanel(HTML("<div style = 'height:40px;line-height:40px';>Species</div>"),h2("Species"),value="Species"),
  
  tabPanel(HTML("<div style = 'height:40px;line-height:40px';>Diversity Profile Estimation</div>"),h2("Diversity Profile Estimation"),value="Diversity"),
  
  tabPanel(HTML("<div style = 'height:40px;line-height:40px';>Shared Species</div>"),h2("Shared Species"),value="Shared Species"),
  
  #tabPanel(title="Interpolation and Extrapolation",h2("Interpolation and Extrapolation"),value="IE"),
  
  tabPanel(HTML("<div style = 'height:40px;line-height:40px';>Two-Community Measures</div>"),h2("Two-Community Measures"),value="Two-Community Measures"),
  
  tabPanel(HTML("<div style = 'height:40px;line-height:40px';>Multiple-Community Measures</div>"),h2("Multiple-Community Measures"),value="Multiple-Community Measures"),
  
  tabPanel(HTML("<div style = 'height:40px;line-height:40px';>Genetics Measures</div>"), h2("Genetics Measures"), value = "Genetics Measures"),#######2015.09.26-(S.W.Wei)
  
  fluidPage(
    tags$head(tags$link(rel = "icon", type = "image/x-icon", href = "http://s6d4.turboimg.net/sp/0a1023688902103219468ecfb3cbe3e5/_.png")),
    fluidRow(
      column(3,
      # tags$em("//",style="color:white;"),
      # tags$img(src="http://s6d4.turboimg.net/sp/0a1023688902103219468ecfb3cbe3e5/_.png" ,alt="W3Schools.com",width = "65px", height = "65px"),
      # # tags$figure(src = "http://s6d4.turboimg.net/sp/0a1023688902103219468ecfb3cbe3e5/_.png", width = "65px", height = "65px"),
      # tags$em("////////////////////////",style="color:white;"), 
      tags$head(      
        tags$style(type='text/css', ".span4 { max-width: 300px; }"),
        tags$style(type="text/css", "label.radio { display: inline-block; }", ".radio input[type=\"radio\"] { float: none; }")
      ),
      h4("Data Setting"),
      wellPanel(
        #selectInput("type","Select data type:",choices=c("Abundance data"="abun","Frequency Data"="abun_infr","Incidence data"="inci")),
        conditionalPanel("input.index=='Species' | input.index=='Diversity'",
                         selectInput("type2","Select data type:",choices=c("Abundance data"="abun_speci","Abundance-frequency counts"="abun_infr","Incidence-frequency data"="inci_speci","Incidence-frequency counts"="inci_count","Incidence-raw data"="inci_raw")),
                         conditionalPanel("input.index=='Species'",
                         numericInput("cut", "Cut off point (For ACE and ICE):", min = 1,  max = 100, value = 10)
                         )),
        conditionalPanel("input.index!='Species' & input.index!='Genetics Measures' & input.index!='Diversity' ",
                         selectInput("type","Select data type:",choices=c("Abundance data"="abun","Incidence-frequency data"="inci","Incidence-raw data"="inci_raw"))),
        conditionalPanel("input.index=='Genetics Measures'",
                         selectInput("typeG","Select data type:",choices=c("Hypothetical Data"="abun1","Human Allele Data"="abun2"))),
        radioButtons("source","Choose one:",choices=c("Demo data"="demo","Upload data"="upload")),
        conditionalPanel("input.source=='upload'",
                         fileInput("files","Choose File"),
                         conditionalPanel(" (input.index=='Two-Community Measures' |input.index=='Multiple-Community Measures'|input.index=='Shared Species') & (input.type=='inci_raw') ",
                                          h5("Number of sampling units in each community:", style = "font-weight: bold;"),
                                          tags$textarea(id="keyin_t", rows = 1, cols=40))),
        conditionalPanel("input.index=='IE' |input.index=='Multiple-Community Measures'|input.index=='Genetics Measures'",
                         radioButtons("order","Choose one order to compute pairwise similarity:",
                                     choices=c("q=0"="0","q=1"="1","q=2"="2"),
                                     selected = 0)),
        conditionalPanel(" (input.index=='IE' |input.index=='Multiple-Community Measures') & (input.order=='1' | input.order=='2')  ",
                         selectInput("method","Select method:",
                                      choices=c("Comparing species relative abundances"="relative","Comparing species absolute abundances"="absolute"))),
        ################################################################2015.09.23-(H.W.Hsu)
        conditionalPanel(condition="input.index == 'Diversity'", 
                         h5("Order q of Hill numbers:"),
                         tags$textarea(id="orderq", rows = 2, cols=40,"0.00 0.25 0.5 0.75 1.00 1.25 1.50 1.75 2.00 2.25 2.50 2.75 3.00")
                         
        )
        ################################################################2015.09.23
      ),
      h4("General Setting"),
      wellPanel(
        numericInput(inputId="nboot", label="Number of bootstraps", value=100, min=1, max=1000, step=1),
        conditionalPanel(condition="input.index == 'IE'",
                         MynumericInput(inputId="endpt", label="Endpoint Setting", value=NULL, placeholder="Double sample size in default"),
                         numericInput(inputId="knot", label="Number of knots", value=40)                       
        )
        #conditionalPanel(condition="input.index =='Shared Species'",
        #                 numericInput(inputId="conf", label="Confidence level", value=0.95, min=0, max=1, step=0.01)
        #)
      ),
      actionButton("goButton",span("Run!",style="font-size: 30px"),icon("rocket","fa-3x"))
    ),
    column(9,
      tabsetPanel(
        tabPanel("Estimation",icon=icon("list-alt"),
                 h3("Estimation"),
                 loadingPanel,
                 verbatimTextOutput("estimation"),
                 downloadLink("dlest", "Download as txt file")
        ),
        
        tabPanel("Visualization",icon=icon("picture-o"),
                 h3("Visualization"),
                 conditionalPanel(" (input.index=='Multiple-Community Measures') | (input.index=='Two-Community Measures')  ",
                                  selectInput("method2",h5("Select method:"), choices=c("Comparing species relative abundances"="relative","Comparing species absolute abundances"="absolute"))),
                 loadingPanel,
                 plotOutput("visualization"),
                 downloadButton('dlvis','Download as PNG File'),
                 plotOutput("visualization2"),
                 #loadingPanel,
                 conditionalPanel(" (input.index=='Multiple-Community Measures') | (input.index=='Two-Community Measures')  ",
                 
                 downloadButton('dlvis2','Download as PNG File')
                 )
                 
                 #downloadLink("dlvis", "Download as pdf file")
        ),
        
        
        tabPanel("Data Summary",icon=icon("file-text-o"),
                 h3("Data Summary"),
                 loadingPanel,
                 h4("Data (First ten species frenquency)"),
                 tableOutput("table"),
                 downloadLink("dltab", "Download as txt file"),
                 uiOutput('rawdata')
        ),
        tabPanel("Introduction",icon=icon("question-circle"),
                 includeMarkdown("src/about.md")
        ),
        tabPanel("User Guide",icon=icon("question-circle"),
                 tags$a("User Guide Link",href="http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/SpadeR_UserGuide.pdf")
        )
      )
    )
    ))
))