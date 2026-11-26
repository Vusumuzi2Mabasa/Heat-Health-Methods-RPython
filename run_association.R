# Rscript --vanilla run_association.R [--analysis-only] [--render]
args<-commandArgs(trailingOnly=TRUE);file_arg<-grep('^--file=',commandArgs(),value=TRUE)
root<-if(length(file_arg))dirname(normalizePath(sub('^--file=','',file_arg[1]),winslash='/')) else normalizePath('.',winslash='/')
options(repos=c(CRAN='https://cloud.r-project.org'))
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
suppressPackageStartupMessages({library(jsonlite);library(sandwich);library(car);library(lpSolve)})
source(file.path(root,'R/01_preprocess.R'),encoding='UTF-8')
source(file.path(root,'R/02_associations.R'),encoding='UTF-8')
source(file.path(root,'R/validate_association.R'),encoding='UTF-8')
if(!'--analysis-only' %in% args){
 source(file.path(root,'R/03_tables_figures.R'),encoding='UTF-8')
 source(file.path(root,'R/04_manuscript.R'),encoding='UTF-8')
 if('--render' %in% args){
  Sys.setenv(R_ENVIRON_USER=if(.Platform$OS.type=='windows')'NUL' else '/dev/null',R_PROFILE_USER=if(.Platform$OS.type=='windows')'NUL' else '/dev/null')
  quarto<-Sys.getenv('QUARTO_PATH',unset=Sys.which('quarto'))
  if(quarto=='')stop('Set QUARTO_PATH to the Quarto executable')
  code<-system2(quarto,c('render',shQuote(file.path(root,'report/Heat_Health_Manuscript.qmd')),'--to','docx'))
  if(code!=0)stop('Quarto rendering failed')
 }
}
cat('R association workflow finished\n')
