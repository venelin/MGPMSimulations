library(PCMBase)
library(TestPCMFit)
TestPCMFit:::GeneratePCMModels()

library(fpc)
trueModel1 <- testData_t2$model[[1]]
tree1 <- testData_t2$treeWithRegimes[[1]]

trueMu1 <- PCMMean(tree1, trueModel1)
trueVar1 <- PCMVar(tree1, trueModel1)

load("local-data/ResultsTestData_t2/FinalResult_MixedGaussian_testData_t2_id_1_t2_.RData")

tree1Inferred <- fitMappings$tree
model1Inferred <- PCMFit::RetrieveBestFitAIC(fitMappings)$inferredModel


inferredMu1 <- PCMMean(model1Inferred$, model1Inferred$inferredModel)
inferredVar1 <- PCMMean(tree1Inferred, model1Inferred$inferredModel)

bhattacharyya.dist()
