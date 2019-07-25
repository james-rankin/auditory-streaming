
FullData=read.table("data-all-switches.csv",sep=',',header=T)

colnames(FullData)=c('trialnum','subject','percept','TMod','dur','ontime','offtime','direction','phase','SwTGroup','meandur','kappaIS','kappaSI','thetaIS','thetaSI');
# meandur, kappas and thetas are computed on a subject-by-subject basis for each TMod condition.
# Hence all rows for the same subject and with same TMod value have the same meandur, kappas and thetas.
frameFD=as.data.frame(FullData)

NumSubj = 17
NumCond = 4

DataBySub = data.frame()
colnames(DataBySub)=c('subject','TMod','meandur','kappaIS','kappaSI','thetaIS','thetaSI','SwTGroup');
DataBySub=as.data.frame(DataBySub)

# Columns with relevant subject summary data
SubjSummInd = c(2,4,11,12,13,14,15,10)
  
TMods = c(0,5,10,20)

# Take single entry for each subject for each TMod condition
for (i in 1:NumSubj) {
  for (j in 1:NumCond) {
    SubjSubset = subset(frameFD,subject==paste("s",i,sep = '') & TMod==paste("TMod",TMods[j],sep = '') )
    temp = SubjSubset[1,SubjSummInd]
    DataBySub = rbind(DataBySub,temp)
  }
}

head(DataBySub)
tail(DataBySub)

library(ez)

# ANOVA on the mean durations w.r.t. to TMod conditions
ezANOVA(data=DataBySub,dv=meandur,wid=subject,within=.(TMod),detailed=TRUE,return_aov = TRUE)
# Corresponding pairwise t-test
pairwise.t.test(DataBySub[,3],DataBySub[,2],p.adj = "bonf")

# Anova between TMod conditions on von mises shape parameter (kappa)
ezANOVA(data=DataBySub,dv=kappaIS,wid=subject,within=.(TMod),detailed=TRUE,return_aov = TRUE)
ezANOVA(data=DataBySub,dv=kappaSI,wid=subject,within=.(TMod),detailed=TRUE,return_aov = TRUE)

# Testing significance between TMod cases on von mises shape parameter (kappa)
pairwise.t.test(DataBySub[,4],DataBySub[,2],p.adj = "bonf")
pairwise.t.test(DataBySub[,5],DataBySub[,2],p.adj = "bonf")

