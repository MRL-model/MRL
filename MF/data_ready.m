%%%%%%%%%%%%getdata%%%%%%
load('expression_matrix.mat');%mRNA expression profile
X1 = expression_Matrix;
load('target.mat');%
C1 = Tar_Matrix;
%%%%%%%input%%%%%%
RCA_batch(X1,C1,'median',0,1,lnc_list,gene_list,sample);
