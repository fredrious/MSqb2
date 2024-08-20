
## Read and save 'config_template' function to sourceable R script.
## Config file will be saved in Script folder and will be sourced by msqb_workflow().
makeConfigFile <- function(conf.file, whichConf) {
  tempfl <- conf.file # file.path(Scripts.path, conf.file)
  # sink(file = paste0(tempfl, ".R"))
  sink(file = tempfl)
  if (whichConf == "workflow") print(config_template)
  if (whichConf == "viz") print(config_viz_template)
  sink()

  x <- readLines(tempfl)
  
  il <- max(grep("## start config", x)) # was: ## start config
  nl <- max(grep("## end config", x))
  x <- readLines(tempfl, n = nl - 1)
  x[1:il] <- ""
  sink(file = tempfl)
  cat(x, sep = "\n")
  sink()
}
