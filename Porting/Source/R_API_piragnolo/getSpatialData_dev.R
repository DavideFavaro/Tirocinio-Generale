library(httr)
library(xml2)
source("getSpatialData_dev_internal.R")



cophub_api<-function (x, p, user, pw) 
{
  if (x == "auto") {
    if (p == "Sentinel-1" | p == "Sentinel-2") {
      x <- "operational"
    }
    else {
      x <- "pre-ops"
    }
  }
  if (x == "operational") {
    x <- getOption("gSD.api")$dhus
  }
  if (x == "pre-ops") {
    x <- getOption("gSD.api")$s3
    user <- "s3guest"
    pw <- "s3guest"
  }
  return(c(user, pw, x))
}




set_archive<-function (dir_archive) 
{
  if (!is.character(dir_archive)) {
    out(paste0("Argument 'dir_archive' needs to be of type 'character'."), 
        type = 3)
  }
  if (!dir.exists(dir_archive)) 
    out("The defined directory does not exist.", type = 3)
  options(gSD.archive = dir_archive)
  options(gSD.archive_set = TRUE)
}



login_CopHub<-function (username, password = NULL) 
{
  if (is.null(password)) {
    password <- getPass()
  }
  char_args <- list(username = username, password = password)
  for (i in 1:length(char_args)) {
    if (!is.character(char_args[[i]])) {
      out(paste0("Argument '", names(char_args[i]), "' needs to be of type 'character'."), 
          type = 3)
    }
  }
  options(gSD.cophub_user = username)
  options(gSD.cophub_pass = password)
  options(gSD.cophub_set = TRUE)
}


getSentinel_query2<-function (time_range, platform, aoi = NULL, username = NULL, 
          password = NULL, hub = "auto", verbose = TRUE, ingestiondate=FALSE, extra=NULL) 
{
  if (is.TRUE(getOption("gSD.cophub_set"))) {
    if (is.null(username)) 
      username <- getOption("gSD.cophub_user")
    if (is.null(password)) 
      password <- getOption("gSD.cophub_pass")
  }
  if (!is.character(username)) 
    out("Argument 'username' needs to be of type 'character'. You can use 'login_CopHub()' to define your login credentials globally.", 
        type = 3)
  if (!is.null(password)) {
    password = password
  }
  else {
    password = getPass()
  }
  if (inherits(verbose, "logical")) 
    options(gSD.verbose = verbose)
  # if (is.null(aoi)) {
  #   if (is.TRUE(getOption("gSD.aoi_set"))) {
  #     aoi <- getOption("gSD.aoi")
  #   }
  #   else {
  #     out("Argument 'aoi' is undefined and no session AOI could be obtained. Define aoi or use set_aoi() to define a session AOI.", 
  #         type = 3)
  #   }
  # }
  #aoi <- make_aoi(aoi, type = "matrix")
  char_args <- list(time_range = time_range, platform = platform)
  for (i in 1:length(char_args)) if (!is.character(char_args[[i]])) 
    out(paste0("Argument '", names(char_args[i]), "' needs to be of type 'character'."), 
        type = 3)
  if (length(time_range) != 2) {
    out("Argument 'time_range' must contain two elements (start and stop time).", 
        type = 3)
  }
  cop.url <- function(ext.xy, url.root, platform, time.range, 
                      row.start) {
   if(!ingestiondate) qs <- list(url.root = paste0(url.root, "/"), search = c("search?start=", 
                                                            "&rows=100&q=("), 
               and = "%20AND%20", 
               aoi.poly = c("footprint:%22Intersects(POLYGON((", ")))%22"), 
               platformname = "platformname:", 
               time = list(`[` = "beginposition:%5b", to = "%20TO%20", `]` = "%5d"))
   else  qs <- list(url.root = paste0(url.root, "/"), search = c("search?start=", 
                                                                 "&rows=100&q=("), 
                    and = "%20AND%20", 
                    aoi.poly = c("footprint:%22Intersects(POLYGON((", ")))%22"), 
                    platformname = "platformname:", 
                    time = list(`[` = "ingestiondate:%5b", to = "%20TO%20", `]` = "%5d"))
   # aoi.str <- paste0(apply(ext.xy, MARGIN = 1, function(x) paste0(x, 
   #                                                                collapse = "%20")), collapse = ",")
    time.range <- sapply(time.range, function(x) {
      tt<-as.POSIXct(time.range)
     # paste0(x, "T00:00:00.000Z")
      format(tt, "%Y-%m-%dT%H:%M:%S.000Z" )
      }, USE.NAMES = F)
 
     
  if(is.null(extra)) urlstring<-paste0(qs$url.root, qs$search[1], toString(row.start),
                      qs$search[2], qs$platformname, platform,  qs$and, qs$time$`[`, 
                      time.range[1], qs$time$to, time.range[2], qs$time$`]`, 
                      ")")
     
  else  urlstring<-paste0(qs$url.root, qs$search[1], toString(row.start),
                      qs$search[2], qs$platformname, platform,  qs$and, extra,  qs$and, qs$time$`[`, 
                      time.range[1], qs$time$to, time.range[2], qs$time$`]`, 
                      ")")
    print(urlstring)
    print(as.character(Sys.time()))
    return(urlstring)
  } 
  cred <- cophub_api(hub, platform, username, password)
  row.start <- -100
  re.query <- T
  give.return <- T
  query.list <- list()
  while (is.TRUE(re.query)) {
    row.start <- row.start + 100
    query <- gSD.get(cop.url(ext.xy = aoi, url.root = cred[3], 
                             platform = platform, time.range = time_range, row.start = row.start), 
                     cred[1], cred[2])
    query.xml <- suppressMessages(xml_contents(content(query)))
    query.list <- c(query.list, lapply(query.xml[grep("entry", 
                                                      query.xml)], function(x) xml_contents(x)))
    if (length(query.list) == 0 & row.start == 0) {
      out("No results could be obtained for this request.", 
          msg = T)
      re.query <- F
      give.return <- F
    }
    if (length(query.list) != row.start + 100) 
      re.query <- F
  }
  if (is.TRUE(give.return)) {
    field.tag <- c("title", "link href", "link rel=\"alternative\"", 
                   "link rel=\"icon\"", "summary", "name")
    field.names <- c("title", "url", "url.alt", "url.icon", 
                     "summary")
    query.cont <- lapply(query.list, function(x, field = field.tag) unlist(lapply(field, 
                                                                                  function(f, y = x) grep(f, y, value = T))))
    query.names <- lapply(query.cont, function(x, field.n = field.names) c(field.n, 
                                                                           sapply(x[(length(field.n) + 1):length(x)], function(y) strsplit(y, 
                                                                                                                                           "\"")[[1]][2], USE.NAMES = F)))
    query.fields <- lapply(query.cont, function(x) sapply(x, 
                                                          function(y) strsplit(strsplit(y, ">")[[1]][2], "<")[[1]][1], 
                                                          USE.NAMES = F))
    query.fields <- lapply(1:length(query.fields), function(i, 
                                                            qf = query.fields, qc = query.cont, qn = query.names) {
      x <- qf[[i]]
      x[which(is.na(x) == T)] <- sapply(qc[[i]][which(is.na(x) == 
                                                        T)], function(y) strsplit(strsplit(y, "href=\"")[[1]][2], 
                                                                                  "\"/>")[[1]][1], USE.NAMES = F)
      names(x) <- qn[[i]]
      return(x)
    })
    return.names <- unique(unlist(query.names))
    return.df <- as.data.frame(stats::setNames(replicate(length(return.names), 
                                                         numeric(0), simplify = F), return.names))
    return.df <- do.call(rbind.data.frame, lapply(query.fields, 
                                                  function(x, rn = return.names, rdf = return.df) {
                                                    rdf[1, match(names(x), rn)] <- sapply(x, as.character)
                                                    return(rdf)
                                                  }))
  }
  if (is.TRUE(give.return)) 
    return(return.df)
}
