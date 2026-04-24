# Data loading functions

library(readxl)

load_bodyfat <- function() {
  readxl::read_excel("data/raw/bodyfat.xlsx")
}
