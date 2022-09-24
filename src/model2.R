#!/usr/bin/R


library(xlsx)         # to load excel sheets
library(glue)         # to format strings


#----------------------------
# load in data
#----------------------------

dropbox_location = '~/Dropbox/'

filename = glue('{dropbox_location}/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx')
df = read.xlsx(filename, sheetIndex=1, check.names=FALSE)

# rename rows
rownames(df) = df$Sample

# FILTER 1
# only keep Probands
df = df[df$Relationship == 'P',]


#----------------------------
# Prepare phenotypes
#----------------------------

pheno_cols = c('Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial')

# for each phenotype split the cohort in two
numerical2binary = function(a, split_value=0) {
	if (is.na(a)) {
		return(NA)
	} else if (a<=split_value) {
		return(0)
	} else if (a>split_value) {
		return(1)
	} 
}

col = 'Child_ID_DD'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=2))

col = 'Child_behav'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=1))

col = 'Child_psych'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=0))

col = 'Child_nervous_system'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=0))

col = 'Child_congenital'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=1))

col = 'Child_craniofacial'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=2))

col = 'Family_history_max'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=1))



#----------------------------
# Prepare predictor variables
#----------------------------


# select columns
binary_cols  = c('Sex')
numeric_cols = c('Rare_Deleterious_SNVs', 'genes_del', 'genes_dup', 'STRs_exonic','SCZ_PRS')



# convert binary columns from char to numeric
for (col in binary_cols){
	df[, col] <- as.factor(df[, col])
	df[, col] <- as.numeric(df[, col]) - 1
}

# convert numeric columns from char to numeric
for (col in numeric_cols){
	df[, col] <- as.numeric(df[, col])
}

# scale the numerical variables
for (col in numeric_cols){
	df[, col] = scale(df[, col])
}

#----------------------------
# Logistic regression
#----------------------------

# print(dim(df))
# for each phenotypic column, fit a linear model
for (col in pheno_cols) {
	print(col)
	formula = glue('{col} ~ Sex + SCZ_PRS + Rare_Deleterious_SNVs + genes_del + genes_dup + STRs_exonic')

	# fit the linear model (prints warning message if exists)
	mod = glm(formula, data=df, family=binomial, x=T)
	
	# print(summary(mod))
	
	# get beta coefficients and confidennce intervals
	summ = confint(mod)
	summ = cbind(summ, coef(mod))
	
	pvalues = coef(summary(mod))[,'Pr(>|z|)']
	summ = cbind(summ, pvalues)
	
	# format output table
	summ = data.frame(summ)
	colnames(summ) = c('2.5% C.I.', '97.5% C.I.', 'Log Odds Ratio', 'P-value')
	summ[,'Variable'] = rownames(summ)
	cols = c('Variable', 'Log Odds Ratio', '2.5% C.I.', '97.5% C.I.', 'P-value')
	summ = summ[,cols]
	
	summ[,'Test'] = 'Logistic regression'
	summ[,'metric'] = 'Log odds ratio'
	summ[,'model'] = 'model2'
	num_samples = dim(mod$x)[1]
	summ[,'Num_samples'] = num_samples
	print(num_samples)
	
	# save table
	outfilename = glue('statistics/model2_{col}.tsv')
	write.table(summ, outfilename, sep='\t', row.names=F)
	
	
}







