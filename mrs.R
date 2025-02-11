library(data.table)
library(matrixStats)
gsub2 <- \(x,pattern,replacement) gsub(pattern,replacement,x) 

ms <- fread('1_NGT_QC_final_240801.csv')
ws <- fread('ukb_nightingale_biomarker_disease_association_atlas.csv')[icd10=='E11']

ws[, biomarker_name := biomarker_name |>
  gsub2(' %','_pct') |>
  gsub2('/','_by_') |>
  gsub2('-| ','_') |>
  gsub2('_particle_size','_size') |>

  gsub2('Alanine','Ala') |>
  gsub2('Glutamine','Gln') |>
  gsub2('Glycine','Gly') |>
  gsub2('Glycoprotein_acetyls','GlycA') |>
  gsub2('Histidine','His') |>
  gsub2('Isoleucine','Ile') |>
  gsub2('Leucine','Leu') |>
  gsub2('Phenylalanine','Phe') |>
  gsub2('Phosphatidylcholines','Phosphatidylc') |>
  gsub2('Phosphoglycerides','Phosphoglyc') |>
  gsub2('3_Hydroxybutyrate','bOHbutyrate') |>
  gsub2('Tyrosine','Tyr') |>
  gsub2('Total_cholines','Cholines') |>
  gsub2('Total_fatty_acids','Total_FA') |>
  gsub2('Total_triglycerides','Total_TG') |>
  gsub2('Valine','Val')
]

setdiff(names(ms), ws$biomarker_name)
setdiff(ws$biomarker_name, names(ms))
# Metabolite names harmonized!

mtx <- as.matrix(ms[,!c('Sample id','Biobank Subject ID')])
mtx[is.na(mtx)] <- (colMins(mtx,na.rm=T)/2)[col(mtx)[is.na(mtx)]] # Half-min imputation

dir.create('hists')
for(i in 1:ncol(mtx)) {
  nm <- colnames(mtx)[i]
  png(file.path('hists',nm))
  hist(mtx[,i],main=nm)
  dev.off()
}

inv_norm <- \(v) qnorm( (rank(v)-0.5)/length(v) )
for(i in 1:ncol(mtx)) mtx[,i] <- inv_norm(mtx[,i])


ws <- ws[
 ][, .(weight_vec = list(   estimate   ),
         name_vec = list(biomarker_name)),
     by=c('age_group','sex','endpoint_type','analysis_method')
 ][, set_nm := paste(age_group,sex,endpoint_type,analysis_method,sep='-')
]

mtx <- mtx[,ws$name_vec[[1]]] # Rearrange
all(sapply(ws$name_vec, identical, colnames(mtx)))
# Metabolites in same order, ready for matrix * vec multiplication!

# This is the important part right here: |----↓--|
ws[, score_vec := lapply(weight_vec,\(v) mtx %*% v |> scale())
 ][, score_vec_no_gluc := lapply(weight_vec, \(v) mtx[,-which(colnames(mtx)=='Glucose')] %*% v[-which(name_vec[[1]]=='Glucose')] |> scale())
]

to_write <-
  as.data.table(ws$score_vec) |>
  #as.data.table(ws$score_vec_no_gluc) |>
  setnames(ws$set_nm)
to_write[
  ][, sample_id          := ms$`Sample id`
  ][, biobank_subject_id := ms$`Biobank Subject ID`
] |> setcolorder(c('sample_id','biobank_subject_id'))

#fwrite(to_write,'mrs_no_gluc.tsv',sep='\t')
fwrite(to_write,'mrs.tsv',sep='\t')



# Sanity check to make sure I know how mtx mult works!
#mtx[1:10,1:3]*tmp$weight_vec[[1]][1:3]
#tmp$weight_vec[[1]][1:3]
#mtx[1:6,1:3]
#tmp$weight_vec[[1]][1:3]%*%mtx[1:6,1:3]
#mtx[1:6,1:3] %*% tmp$weight_vec[[1]][1:3]
