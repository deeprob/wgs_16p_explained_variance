#!/usr/bin/R


library(xlsx)         # to load excel sheets
library(glue)         # to format strings


#----------------------------
# load in data
#----------------------------



pheno_cols = c('Full-scale IQ','ABCL_CBCL external','ABCL_CBCL internal','SRS raw','BSI','BMI Z-score','Head circumference Z-score')
binary_cols  = c('Sex')
numeric_cols = c('Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'Genes DEL', 'DELs LOEUF<0.35', 'Genes DUP', 'DUPs LOEUF<0.35', 'SCZ_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'autism_PRS')
other_numeric = c('Missense CADD>25', 'LOF CADD25_NA', 'Splice CADD25', 'Genes DEL', 'Genes DUP')


cohorts = c('16p-deletion', '16p-duplication')

variant_cols = c('Missense CADD>25', 'LOF CADD25_NA', 'LOF CADD LOEUF<0.35', 'Splice CADD25', 'Splice CADD25 LOEUF<0.35', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'Genes DEL', 'Genes DUP', 'DELs LOEUF<0.35', 'DUPs LOEUF<0.35', 'SCZ_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'autism_PRS')


linear_regression = function(df, formula, outfilename) {
	# fit the linear model
	mod = lm(formula, data=df)
	
	# get beta coefficients and confidennce intervals
	summ = confint(mod)
	summ = cbind(summ, coef(mod))
	
	pvalues = coef(summary(mod))[,'Pr(>|t|)']
	summ = cbind(summ, pvalues)
	
	# format output table
	summ = data.frame(summ)
	colnames(summ) = c('2.5% C.I.', '97.5% C.I.', 'Log Odds Ratio', 'P-value')
	summ[,'Variable'] = rownames(summ)
	cols = c('Variable', 'Log Odds Ratio', '2.5% C.I.', '97.5% C.I.', 'P-value')
	summ = summ[,cols]
	
	summ[,'Test'] = 'Linear regression'
	summ[,'metric'] = 'Beta coefficient'
	# get number of observations
	num_observations = nobs(mod)
	summ['Num_samples'] = num_observations
	
	# save table	
	write.table(summ, outfilename, sep='\t', row.names=F)
}



for (cohort in cohorts) {
	filename = glue('../1_variant_preperation/SVIP.csv')
	df = read.csv(filename, check.names=F)
	df = df[df[,'Family type'] == cohort,]
	colnames(df) = gsub('/', '_', colnames(df))
	rownames(df) = df$Sample	
	
	df['All_rare_del_var'] = df[,'Missense CADD>25'] + df[,'LOF CADD25_NA'] +
	df[,'Splice CADD25'] + df[, 'Genes DEL'] + df[, 'Genes DUP']
	
	# convert binary columns from char to numeric
	for (col in binary_cols){
		df[, col] <- as.factor(df[, col])
		df[, col] <- as.numeric(df[, col]) - 1
	}
	
	# scale the numerical variables
	for (col in c(numeric_cols, pheno_cols, other_numeric, c('All_rare_del_var'))){
		df[, col] = as.numeric(df[, col])
		df[, col] = scale(df[, col])
	}

	
	print(cohort)
	
	print('MODEL 1')
	# for each phenotypic column, fit a linear model
	for (col in pheno_cols) {
		print(col)

		formula = glue('`{col}` ~ Sex + SCZ_PRS + All_rare_del_var')
		outfilename = glue('statistics/{cohort}_{col}_model1.tsv')
		linear_regression(df, formula, outfilename)

	}
	
	print('MODEL 2')
	for (col in pheno_cols) {
		print(col)

		formula = glue('`{col}` ~ Sex + SCZ_PRS + Rare_Deleterious_SNVs + `Genes DEL` + `Genes DUP`')
		outfilename = glue('statistics/{cohort}_{col}_model2.tsv')
		linear_regression(df, formula, outfilename)

	}

	print('MODEL 3')
	for (col in pheno_cols) {
		print(col)

		formula = glue('`{col}` ~ Sex + SCZ_PRS + Rare_Deleterious_SNVs_LOEUF + `DELs LOEUF<0.35` + `DUPs LOEUF<0.35`')
		outfilename = glue('statistics/{cohort}_{col}_model3.tsv')
		linear_regression(df, formula, outfilename)

	}
}







	
	
	
	
	
