library(rTRIPLEXSIF)
library(testthat)

test_that("basic model run", {
  path <- system.file(Inputvariable, package = "rTRIPLEXSIF")
  path <- system.file(Inputpara, package = "rTRIPLEXSIF")

  examdata<-data.frame()
  for(i in 1:12){
    subdata<-subset(Inputvariable,Inputvariable$Month==i)
    mondata<-subset(subdata,subdata$Day<=2)
    examdata<-rbind(examdata,mondata)
  }

  out <- TRIPLEX_SIF(
    Input_variable=examdata,Input_parameter=Inputpara)

  #expect_type(out,'list')
})
