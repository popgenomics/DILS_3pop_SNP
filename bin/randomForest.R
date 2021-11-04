library(abcrf)

for(tmp in commandArgs()){
	tmp = strsplit(tmp, '=')
	if(tmp[[1]][1] == 'projectpath'){ projectpath = tmp[[1]][2] }
	if(tmp[[1]][1] == 'binpath'){ binpath = tmp[[1]][2] }
	if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nIterations'){ nIterations = as.numeric(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'ncores'){ ncores = as.numeric(tmp[[1]][2]) }
}

setwd(paste(projectpath, '/', timeStamp, sep=''))

obs = read.table('ABCstat.txt', h=T)
toReject = sort(unique(c(grep('min', colnames(obs)), grep('max', colnames(obs)), grep('dataset', colnames(obs)))))

epsilon=0.00001
ntree=1000


### first comparison : ab versus AB
ab = c()
for(i in c('ab')){
	for(j in c('abc', 'ABC', 'AbC', 'aBC')){
		for(M in c('1M', '2M')){
			for(N in c('1N', '2N')){
				ab = c(ab, paste(i,j,M,N,sep='_'))
			}
		}
	}
}

AB = c()
for(i in c('AB')){
	for(j in c('abc', 'ABC', 'AbC', 'aBC')){
		for(M in c('1M', '2M')){
			for(N in c('1N', '2N')){
				AB = c(AB, paste(i,j,M,N,sep='_'))
			}
		}
	}
}

ab_sim = NULL
for(model in ab){
	for(i in 0:(nIterations-1)){
		tmp = read.table(paste('modelComp/', model, '_', i, '/ABCstat.txt', sep=''), h=T)[1:200,]
		ab_sim = rbind(ab_sim, tmp)
	}
}

AB_sim = NULL
for(model in AB){
	for(i in 0:(nIterations-1)){
		tmp = read.table(paste('modelComp/', model, '_', i, '/ABCstat.txt', sep=''), h=T)[1:200,]
		AB_sim = rbind(AB_sim, tmp)
	}
}

simulations = rbind(ab_sim, AB_sim)
toReject_tmp = toReject
for(i in 1:ncol(simulations)){
	sd_tmp = sd(simulations[,i])
	if(sd_tmp<epsilon){
		toReject_tmp = c(toReject_tmp, i)
	} 
}
toReject_tmp = unique(toReject_tmp)
print(colnames(simulations)[toReject_tmp])

modIndex = as.factor(c(rep('ab', nrow(ab_sim)), rep('AB', nrow(AB_sim))))
data_simulations = data.frame(modIndex, simulations[,-toReject_tmp])

# build the forest
model_rf = abcrf(modIndex~., data=data_simulations, ntree=ntree, paral=T, ncores=ncores)

# predict the best model
pred_1 = predict(model_rf, obs[-toReject_tmp], data_simulations, ntree=ntree, paral.predict=T, ncores.predict=ncores)

res = paste(c(as.character(pred_1$allocation), as.character(pred_1$post.prob)), collapse=': ')


### Second comparison : abc versus ABC versus AbC versus aBC
abc = c()
for(i in c('ab', 'AB')){
	for(j in c('abc')){
		for(M in c('1M', '2M')){
			for(N in c('1N', '2N')){
				abc = c(abc, paste(i,j,M,N,sep='_'))
			}
		}
	}
}
ABC = c()
for(i in c('ab', 'AB')){
	for(j in c('ABC')){
		for(M in c('1M', '2M')){
			for(N in c('1N', '2N')){
				ABC = c(ABC, paste(i,j,M,N,sep='_'))
			}
		}
	}
}
AbC = c()
for(i in c('ab', 'AB')){
	for(j in c('AbC')){
		for(M in c('1M', '2M')){
			for(N in c('1N', '2N')){
				AbC = c(AbC, paste(i,j,M,N,sep='_'))
			}
		}
	}
}
aBC = c()
for(i in c('ab', 'AB')){
	for(j in c('aBC')){
		for(M in c('1M', '2M')){
			for(N in c('1N', '2N')){
				aBC = c(aBC, paste(i,j,M,N,sep='_'))
			}
		}
	}
}
abc_sim = NULL
for(model in abc){
	for(i in 0:(nIterations-1)){
		tmp = read.table(paste('modelComp/', model, '_', i, '/ABCstat.txt', sep=''), h=T)
		abc_sim = rbind(abc_sim, tmp)
	}
}
ABC_sim = NULL
for(model in ABC){
	for(i in 0:(nIterations-1)){
		tmp = read.table(paste('modelComp/', model, '_', i, '/ABCstat.txt', sep=''), h=T)
		ABC_sim = rbind(ABC_sim, tmp)
	}
}
AbC_sim = NULL
for(model in AbC){
	for(i in 0:(nIterations-1)){
		tmp = read.table(paste('modelComp/', model, '_', i, '/ABCstat.txt', sep=''), h=T)
		AbC_sim = rbind(AbC_sim, tmp)
	}
}
aBC_sim = NULL
for(model in aBC){
	for(i in 0:(nIterations-1)){
		tmp = read.table(paste('modelComp/', model, '_', i, '/ABCstat.txt', sep=''), h=T)
		aBC_sim = rbind(aBC_sim, tmp)
	}
}

modIndex = as.factor(c(rep('abc', nrow(abc_sim)), rep('ABC', nrow(ABC_sim)), rep('AbC', nrow(AbC_sim)), rep('aBC', nrow(aBC_sim))))
simulations = rbind(abc_sim, ABC_sim, AbC_sim, aBC_sim)

toReject_tmp = toReject
for(i in 1:ncol(simulations)){
	sd_tmp = sd(simulations[,i])
	if(sd_tmp<epsilon){
		toReject_tmp = c(toReject_tmp, i)
	} 
}
toReject_tmp = unique(toReject_tmp)
print(colnames(simulations)[toReject_tmp])

data_simulations = data.frame(modIndex, simulations[, -toReject_tmp])

# build the forest
model_rf = abcrf(modIndex~., data=data_simulations, ntree=ntree, paral=T, ncores=ncores)

# predict the best model
pred_2 = predict(model_rf, obs[-toReject_tmp], data_simulations, ntree=ntree, paral.predict=T, ncores.predict=ncores)

### Third comparison : 1M_1N versu 2M_2N versu 1M_2N versus 2M_1N
model = paste(as.character(pred_1$allocation), as.character(pred_2$allocation), sep='_')

M1_N1 = M2_N1 = M1_N2 = M2_N2 = NULL

for(i in 0:(nIterations-1)){
	# 1M_1N
	tmp = read.table(paste('modelComp/', model, '_1M_1N_', i, '/ABCstat.txt', sep=''), h=T)
	M1_N1 = rbind(M1_N1, tmp)

	# 2M_1N
	tmp = read.table(paste('modelComp/', model, '_2M_1N_', i, '/ABCstat.txt', sep=''), h=T)
	M2_N1 = rbind(M2_N1, tmp)

	# 1M_2N
	tmp = read.table(paste('modelComp/', model, '_1M_2N_', i, '/ABCstat.txt', sep=''), h=T)
	M1_N2 = rbind(M1_N2, tmp)

	# 2M_2N
	tmp = read.table(paste('modelComp/', model, '_2M_2N_', i, '/ABCstat.txt', sep=''), h=T)
	M2_N2 = rbind(M2_N2, tmp)
}

modIndex = as.factor(c(rep('1M_1N', nrow(M1_N1)), rep('2M_1N', nrow(M2_N1)), rep('1M_2N', nrow(M1_N2)), rep('2M_2N', nrow(M2_N2))))
simulations = rbind(M1_N1, M2_N1, M1_N2, M2_N2)

toReject_tmp = toReject
for(i in 1:ncol(simulations)){
	sd_tmp = sd(simulations[,i])
	if(sd_tmp<epsilon){
		toReject_tmp = c(toReject_tmp, i)
	} 
}
toReject_tmp = unique(toReject_tmp)
print(colnames(simulations)[toReject_tmp])

data_simulations = data.frame(modIndex, simulations[, -toReject_tmp])

# build the forest
model_rf = abcrf(modIndex~., data=data_simulations, ntree=ntree, paral=T, ncores=ncores)

# predict the best model
pred_3 = predict(model_rf, obs[-toReject_tmp], data_simulations, ntree=ntree, paral.predict=T, ncores.predict=ncores)


header = c('ab vs AB', 'abc vs ABC vs AbC vs aBC', '1M_1N vs 1M_2N vs 2M_1N vs 2M_2N')
line1 =  c(as.character(pred_1$allocation), as.character(pred_2$allocation), as.character(pred_3$allocation))
line2 =  c(as.character(pred_1$post.prob), as.character(pred_2$post.prob), as.character(pred_3$post.prob))
write(header, 'modelComp/hierarchical_models.txt', append=F, sep='\t', ncolumns=3)
write(line1, 'modelComp/hierarchical_models.txt', append=T, sep='\t', ncolumns=3)
write(line2, 'modelComp/hierarchical_models.txt', append=T, sep='\t', ncolumns=3)

# estimation of parameters of best parameters
best = paste(model, as.character(pred_3$allocation), sep='_')
param = stats = NULL
for(i in 0:(nIterations-1)){
		tmp = read.table(paste('modelComp/', best, '_', i, '/ABCstat.txt', sep=''), h=T)
		stats = rbind(stats, tmp)

		tmp = read.table(paste('modelComp/', best, '_', i, '/priorfile.txt', sep=''), h=T)
		param = rbind(param, tmp)
}


toReject_tmp = toReject
for(i in 1:ncol(stats)){
	sd_tmp = sd(stats[,i])
	if(sd_tmp<epsilon){
		toReject_tmp = c(toReject_tmp, i)
	} 
}
toReject_tmp = unique(toReject_tmp)
print(colnames(stats)[toReject_tmp])

nParams = ncol(param)

res_params = c()
res_expectation = c()
res_med = c()
res_variance = c()
res_varianceCDF = c()
res_quantile_0025 = c()
res_quantile_0975 = c()
res_NMAEmean = c()

for(i in 1:nParams){
	parameter = param[,i]
	data = data.frame(parameter, stats[, -toReject_tmp])
	mod = regAbcrf(parameter~., data, ntree=ntree, paral=T, ncores=ncores)
	estimate = predict(mod, obs[-toReject_tmp], data, paral=T, ncores=ncores)

	res_params = c(res_params, colnames(param)[i])
	res_expectation = c(res_expectation, estimate$expectation)
	res_med = c(res_med, estimate$med)
	res_variance = c(res_variance, estimate$variance)
	res_varianceCDF = c(res_varianceCDF, estimate$variance.cdf)
	res_quantile_0025 = c(res_quantile_0025, estimate$quantiles[1])
	res_quantile_0975 = c(res_quantile_0975, estimate$quantiles[2])
	res_NMAEmean = c(res_NMAEmean, estimate$post.NMAE.mean)
}


res_RF = data.frame(parameters = res_params, expectation = res_expectation, median = res_med, variance = res_variance, variance_CDF = res_varianceCDF, quantile_0025 = res_quantile_0025, quantile_0975 = res_quantile_0975, NMAE_mean = res_NMAEmean)
write.table(res_RF, 'modelComp/posterior_RF.txt', col.names=T, row.names=F, sep='\t', quote=F)

write(best, 'modelComp/best_model.txt')

