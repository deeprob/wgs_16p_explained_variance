#!/usr/bin/R


library(openxlsx)         # to load excel sheets
library(glue)         # to format strings
library(DescTools)


#----------------------------
# load in data
#----------------------------

filename = glue('../data/16p12_cohort_summary_v17.xlsx')
df = read.xlsx(filename, sheet=1, check.names=FALSE)

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




# # print phenotype value counts
# for (col in pheno_cols) {
# 	print(col)
# 	print(table(df[,col]))
# }

#----------------------------
# Prepare predictor variables
#----------------------------


# variant_groups = ['Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF']


# select columns
binary_cols  = c('Sex')
numeric_cols = c('Rare_Deleterious_SNVs_LOEUF', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic_LOEUF035','SCZ_PRS')


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


# for each phenotypic column and all genetic variables together, fit a linear model
for (col in pheno_cols) {
	formula = glue('{col} ~ Sex + SCZ_PRS + Rare_Deleterious_SNVs_LOEUF + dels_loeuf + dups_loeuf + STRs_exonic_LOEUF035')
	print(formula)
	# fit the linear model
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
	summ[,'model'] = 'model3'
	num_samples = dim(mod$x)[1]
	summ[,'Num_samples'] = num_samples
	# add Nagelkerke's R2
	summ[, "R2"] = PseudoR2(mod, which="Nagelkerke")
	# print(summ)
	# save table
	outfilename = glue('../data/statistics/model3_{col}.tsv')
	write.table(summ, outfilename, sep='\t', row.names=F)
	
	
}


# for each phenotypic column and each genetic variable individually, fit a linear model
for (col in pheno_cols) {
	for (gv in numeric_cols) {

		formula = glue('{col} ~ Sex + {gv}')
		print(formula)

		# fit the linear model
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
		summ[,'model'] = 'model3'
		num_samples = dim(mod$x)[1]
		summ[,'Num_samples'] = num_samples
		# add Nagelkerke's R2
		summ[, "R2"] = PseudoR2(mod, which="Nagelkerke")
		# print(summ)
		# save table
		outfilename = glue('../data/statistics/model3_{col}_{gv}.tsv')
		write.table(summ, outfilename, sep='\t', row.names=F)
	}

}





