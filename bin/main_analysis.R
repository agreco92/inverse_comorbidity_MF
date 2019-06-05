# reproducible example
#### libraries ####
#### input ####
working_directory = './'
setwd(working_directory)
datasets_folder = './data/datasets_IC2/' 
# this folder must contain separate folders with the metagenes and metasample matrices to analyze.
# The metagenes and metasample matrices must be in the format: dataset name + m_samples_tag/m_genes_tag
# For example, the dataset 'pippo' will have to follow the following path: ./data/datasets_IC2/pippo/pippo_A.xls (analogously for S)
m_samples_tag = '_A.xls' # identifier for meta-samples (samples x ic) matrix inside single datasets folders
m_genes_tag = '_S.xls'  # identifier for meta-genes matrix (genes x ic) inside single datasets folders
case_tag = 'Case' # identifier of case samples for orienting components
diseases = c('Alzheimer','Lung') # strings that identify the conditions in the network

# network parameters
disease_pointing_orientation = TRUE # if TRUE, orients the components maximising the projection of case metasamples onto metagenes [ref paper]

#### disease pointing metagenes ####
flipping_matrix = matrix(1, nrow = length(dir(datasets_folder)), ncol = 100, 
                        dimnames = list(dir(datasets_folder),paste0('IC',seq(1,100))))

if(disease_pointing_orientation){
  for(folder in dir(datasets_folder)){ # loop over datasets to orient components
    # print(folder)
    m_samples_file = grep(m_samples_tag, dir(paste0(datasets_folder,folder),full.names = T),value = T, fixed = T)
    m_genes_file = grep(m_genes_tag, dir(paste0(datasets_folder,folder),full.names = T),value = T, fixed = T)

    m_samples = read.table(m_samples_file, header = TRUE,row.names = 1)
    m_genes = read.table(m_genes_file, header = TRUE,row.names = 1)
    
    cases_cols = rownames(m_samples)[grep(case_tag,x = rownames(m_samples))] # columns containing case tag
    cases = m_samples[cases_cols,] # subsetting full matrix
    # flip ICs in which the average projection of the cases is negative
    for (col in colnames(m_genes)){
      if (mean(cases[,col]) < 0){
        flipping_matrix[folder,col] = -1
      }
    }
  }
}
#### correlation network ####
d1l = list()
d2l = list()
weight = list()
for (n1 in dir(datasets_folder)[-length(dir(datasets_folder))]){
  #n1 = dir(datasets_folder)[1]
  print(n1)
  name1 = paste0(datasets_folder,n1,'/',n1,m_genes_tag)
  d1 = read.table(name1,header = TRUE,row.names = 1)
  dis1 = grep(pattern = "Alzheimer|Lung",x = unlist(strsplit(n1,split="-|_")),value = TRUE)
  for (n2 in dir(datasets_folder)[(match(x = n1,dir(datasets_folder))+1):length(dir(datasets_folder))]){ #cycle on remaining folders
    #n2 = dir(datasets_folder)[2]
    print(c(n1,n2))
    name2 <- paste0(datasets_folder,n2,'/',n2,m_genes_tag)
    d2 = read.table(name2,header = TRUE,row.names = 1)
    dis2 = grep(pattern = "Alzheimer|Lung",x = unlist(strsplit(n2,split="-|_")),value = TRUE)
    subd1 = t(t(d1)*flipping_matrix[n1,1:ncol(d1)])
    subd2 = t(t(d2)*flipping_matrix[n2,1:ncol(d2)])
    int = intersect(rownames(d1),rownames(d2)) # common genes
    subd1 = subd1[int,]
    subd2 = subd2[int,]
    
    corr_mat = cor(subd1,subd2, method = 'pearson')
    for (i in rownames(corr_mat)){  #loop on rows (d1)
      minrow = min(corr_mat[i,])
      maxrow = max(corr_mat[i,])
      for (j in colnames(corr_mat)){  #loop on columns (d2)
        if(abs(corr_mat[i,j])>0.05){
          mincol = min(corr_mat[,j])
          maxcol = max(corr_mat[,j])
          if ((((mincol==minrow) && (minrow<0))||((maxcol==maxrow) && (maxrow>0)))){ # RBH condition
            IC1 = paste0(n1,'_',i)
            IC2 = paste0(n2,'_',j)
            d1l = c(d1l,IC1)
            d2l = c(d2l,IC2)
            weight = c(weight,corr_mat[i,j])
          }
        }
      }
    }
  }
}

d1l = unlist(d1l)
d2l = unlist(d2l)
weight = unlist(weight)
matrice = matrix(NA,nrow = length(d1l),ncol = 3)
matrice[,1]= unlist(d1l)
matrice[,2]= unlist(d2l)
matrice[,3]= unlist(weight)
tabla = as.table(matrice)

#### prune links #####
disease_table = tabla # this table will be used to prune the links

for(disease in diseases){
  disease_table = gsub(paste0('.*',disease,'.*'), disease ,disease_table)
}

# keep negative links between metagenes from different conditions
mask = (disease_table[,1]!= disease_table[,2]) & disease_table[,3]<0
tmp_tabla = tabla[mask,] # filtered edge list
#### community detection ####
nodes = unique(c(tmp_tabla[,1], tmp_tabla[,2]))
adj_matrix = matrix(0, nrow = length(nodes), ncol = length(nodes))
colnames(adj_matrix) = rownames(adj_matrix) = nodes

for(row in 1:nrow(tmp_tabla)){ #fill the adjacency matrix
  adj_matrix[tmp_tabla[row,1],tmp_tabla[row,2]] =1
  adj_matrix[tmp_tabla[row,2],tmp_tabla[row,1]] =1
}

# community analysis
library(MCL)
labels = MCL::mcl(adj_matrix,addLoops = F,allow1 = T)$Cluster

#### average metagenes ####
library(tidyverse)
sigmas = 3 # how many standard deviations from the mean should a gene be to enter the list?
community_table = tibble(metagene = nodes, cluster = labels)
for (cluster in unique(community_table$cluster)){
  metagenes = unlist(community_table[community_table$cluster == cluster,1])
  metagenes = lapply(metagenes, function(mg){
    num = as.integer(gsub('.*IC','',mg)) # metagene number
    dset.string = gsub('_IC.*','',mg)
    dset = read.table(paste0('./data/datasets_IC2/', dset.string, '/', dset.string,'_S.xls'), header = TRUE, row.names = 1)
    full_mgene = tibble(gene = rownames(dset), mg = dset[,num]*flipping_matrix[dset.string,num]) #isolate full metagene and orient accordingly to flipping_matrix
    colnames(full_mgene) = c('gene', mg)
    return(full_mgene)
  })
  mat = Reduce(inner_join, metagenes)
  # flip the Lung components (substitute with appropriate tag
  mat[,grepl(diseases[2], colnames(mat))] = - mat[,grepl(diseases[2], colnames(mat))]
  # rownames(mgenes) = rownames(dset)
  average.metagene = rowMeans(mat %>% select(-gene))
  # enrichment thing
  names(average.metagene) = mat$gene
  top_mergene = average.metagene[average.metagene> mean(average.metagene) + sigmas *sqrt(var(average.metagene))]
  bot_mergene = average.metagene[average.metagene< mean(average.metagene) - sigmas *sqrt(var(average.metagene))]
  top_mergene = names(top_mergene[order(top_mergene,decreasing=TRUE)])
  bot_mergene = names(bot_mergene[order(bot_mergene,decreasing=TRUE)])
  # write results
  write.table(top_mergene,file = paste('./data/top_genes/',cluster,'_',sigmas,'_sigmas','_up.csv',sep = ''),col.names =FALSE,quote = FALSE,sep = '\t',row.names=FALSE)
  write.table(bot_mergene,file = paste('./data/top_genes/',cluster,'_',sigmas,'_sigmas','_dn.csv',sep = ''),col.names = FALSE,quote = FALSE,sep = '\t',row.names=FALSE)
}
