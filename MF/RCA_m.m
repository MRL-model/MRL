function RCA_m(X1, C1, zMode, nmfFlag, keepOrilncFlag,lnc_label,gene_label,sample_label)
 X = X1;
 A = C1;
 Z = zMode;
 RCA_batch(X,A,Z,0,1,lnc_label,gene_label,sample_label);

exit;%%%自动退出
