devtools::install_github("kgoldfeld/simstudy")

library(brms)
library(projpred)
library(simstudy)

set.seed(87261)


gen.school <- defData(varname = "s0", dist = "normal", formula = 0, variance = 3,
                      id = "idSchool")
gen.school <- defData(gen.school, varname = "nClasses", dist = "noZeroPoisson", formula = 3)

dtSchool <- genData(8, gen.school)
dtSchool <- trtAssign(dtSchool, n = 2)

gen.class <- defDataAdd(varname = "c0", dist = "normal", formula = 0, variance = 2)
gen.class <- defDataAdd(gen.class, varname = "nStudents", dist = "noZeroPoisson",
                        formula = 20)

dtClass <- genCluster(dtSchool, "idSchool", numIndsVar = "nClasses", level1ID = "idClass")
dtClass <- addColumns(gen.class, dtClass)

gen.student <- defDataAdd(varname = "Male", dist = "binary",
                          formula = 0.5)
gen.student <- defDataAdd(gen.student, varname = "age", dist = "uniform",
                          formula = "9.5; 10.5")
gen.student <- defDataAdd(gen.student, varname = "test", dist = "normal",
                          formula = "50 - 5 * Male + s0 + c0 + 8 * trtGrp", 
                          variance = 5)

dtStudent <- genCluster(dtClass, cLevelVar = "idClass", numIndsVar = "nStudents",
                        level1ID = "idChild")
dtStudent <- addColumns(gen.student, dtStudent)
dtStudent$Male <- as.factor(dtStudent$Male)
dtStudent$idClass <- as.factor(dtStudent$idClass)
dtStudent$idSchool <- as.factor(dtStudent$idSchool)
dtStudent

# write the table
setwd("~/Desktop/projpred-workflow")
csv_name <- paste0("./data/test_results.csv")
ff <- file(csv_name, open="w")
write.table(dtStudent, file = ff, sep = ",", row.names = FALSE)
close(ff)


# fit the BRMS model
fit <- brm(
  formula = "test ~ C(Male) + C(idSchool) + C(idClass) + age", 
  family = "gaussian", 
  data = dtStudent
)
fit

# perform projection predictive inference
vs <- varsel(
  fit, 
  method = "forward"
)
plot(vs)
vs$solution_terms

