#' Package Startup Message
#'
#' This function displays a custom startup message when the package is loaded.
#' It is automatically called when the package is attached via `library()` or `require()`.
#'
#' @importFrom glue glue
#' @importFrom utils packageVersion
#' @noRd
.onAttach <- function(libname, pkgname) {

logo1 <- "
 _______  _______         __
|   |   ||     __|.-----.|  |--.
|       ||__     ||  _  ||  _  |
|__|_|__||_______||__   ||_____|
                     |__|
"

logo2 <- "
{__       {__  {__ __          {__
{_ {__   {___{__    {__        {__
{__ {__ { {__ {__        {__   {__
{__  {__  {__   {__    {_  {__ {__ {__
{__   {_  {__      {__ {_  {__ {__   {__
{__       {__{__    {__ {__{__ {__   {__
{__       {__  {__ __      {__ {__ {__
                           {___
"

logo3 <- "______/\\/\\______/\\/\\____/\\/\\/\\/\\/\\____________/\\/\\________
_____/\\/\\/\\__/\\/\\/\\__/\\/\\__________/\\/\\/\\/\\__/\\/\\________
____/\\/\\/\\/\\/\\/\\/\\____/\\/\\/\\/\\__/\\/\\__/\\/\\__/\\/\\/\\/\\____
___/\\/\\__/\\__/\\/\\__________/\\/\\__/\\/\\/\\/\\__/\\/\\__/\\/\\__
__/\\/\\______/\\/\\__/\\/\\/\\/\\/\\________/\\/\\__/\\/\\/\\/\\____
___________________________________/\\/\\______________
"
  if (!interactive()) return()

  v = packageVersion("MSqb2")
  d = read.dcf(system.file("DESCRIPTION", package = "MSqb2"),
               fields = c("Packaged", "Built", "Revision"))

  if (is.na(d[1L])) {
    if (is.na(d[2L])) {
      return() # Packaged, Built fields does not exists
    } else {
      d = unlist(strsplit(d[2L], split = "; "))[3L]
      g = if (is.na(d[3L])) "" else paste0(" (",d[3L],")")
    }
  } else {
    d = d[1L]
    g = ""
  }

  packageStartupMessage(glue::glue("\n\nloading ...\n{logo3}MSqb2 v{v} built {d}.\n\n"), appendLF = FALSE)
}
