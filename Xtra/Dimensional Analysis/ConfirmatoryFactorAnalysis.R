# http://www.statpower.net/Content/312/Handout/Confirmatory%20Factor%20Analysis%20with%20R.pdf
# http://methodsconsultants.com/tutorial/9/Confirmatory-Factor-Analysis-Using-the-SEM-Package-in-R
# http://pareonline.net/getvn.asp?v=18&n=4

usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("sem")

visr.applyParameters()