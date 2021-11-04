library(tidyverse)
library(FactoMineR)

for(tmp in commandArgs()){
	tmp = strsplit(tmp, '=')
	if(tmp[[1]][1] == 'projectpath'){ projectpath = tmp[[1]][2] }
	if(tmp[[1]][1] == 'binpath'){ binpath = tmp[[1]][2] }
	if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nIterations'){ nIterations = as.numeric(tmp[[1]][2]) }
}

setwd(paste(projectpath, '/', timeStamp, sep=''))

obs = read.table('ABCstat.txt', h=T)
toReject = sort(unique(c(grep('pearson', colnames(obs)), grep('min', colnames(obs)), grep('max', colnames(obs)), grep('dataset', colnames(obs)))))


# models
model_names = c()
for(i in c('ab', 'AB')){
	for(j in c('abc', 'ABC', 'AbC', 'aBC')){
		model_names = c(model_names, paste(i, j, sep='_'))
	}
}


models = c()
statistics = NULL
for(model in model_names){
	for(i in 0:(nIterations-1)){
		tmp = read.table(paste('modelComp/', model, '_1M_1N_', i, '/ABCstat.txt', sep=''), h=T)[, -toReject]
		statistics = rbind(statistics, tmp)
		models = c(models, rep(model, nrow(tmp)))
		
		tmp = read.table(paste('modelComp/', model, '_1M_2N_', i, '/ABCstat.txt', sep=''), h=T)[, -toReject]
		statistics = rbind(statistics, tmp)
		models = c(models, rep(model, nrow(tmp)))
		
		tmp = read.table(paste('modelComp/', model, '_2M_1N_', i, '/ABCstat.txt', sep=''), h=T)[, -toReject]
		statistics = rbind(statistics, tmp)
		models = c(models, rep(model, nrow(tmp)))
		
		tmp = read.table(paste('modelComp/', model, '_2M_2N_', i, '/ABCstat.txt', sep=''), h=T)[, -toReject]
		statistics = rbind(statistics, tmp)
		models = c(models, rep(model, nrow(tmp)))
	}
}



models = c('observation', models)
statistics = rbind(obs[-toReject], statistics)

dataset = cbind(models, statistics)

res.pca = PCA(dataset[, -1], scale.unit=TRUE, ncp=3, graph=F)

write.table(res.pca$ind$coord, 'modelComp/prior_checking.txt', col.names=T, row.names)

