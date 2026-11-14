function RCA_m(X1, C1, zMode, nmfFlag, keepOrilncFlag,lnc_list,gene_list,sample)
 X = X1;
 A = C1;
 Z = zMode;
 RCA_batch(X,A,Z,0,1,lnc_list,gene_list,sample);

exit;
