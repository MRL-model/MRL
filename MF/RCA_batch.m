function RCA_batch(X1, C1, zMode, nmfFlag, keepOrilncFlag,lnc_label,gene_label,sample_label)
% this program create output1.dat, not sort y matrix  (Bin Yan)
% (1) Open C-T binding data file
tokens = C1;
lncSym = lnc_label;
% Get data
subDataX = X1;

[N, M] = size(subDataX);  % N: number of genes. 
                        % M: number of experiments.

% Get gene symbol of array data

geneSymID = gene_label;

% Get sample label of array data
sampleID = sample_label;

% construct connectionComp sparse matrix from conserved motif information
numC = length(lncSym);
connectionC = C1';

% get Comp subdataset from subdataset subDataX
sublncDataX = [];
if strcmpi(zMode, 'pls') == 1
    xmean = mean(connectionC);
    xstd = std(connectionC);
    ymean = mean(subDataX);
    ystd = std(subDataX);
    connectionlnccentered = (connectionC - xmean(ones(N,1),:))./xstd(ones(N,1),:);%%%./点除,如果a、b是矩阵，a./b就是a、b中对应的每个元素相除，得到一个新的矩阵
    subDataXcentered = (subDataX - ymean(ones(N,1),:))./ystd(ones(N,1),:);
    [Xloadings,Yloadings,Xscores,Yscores,Weights] = simpls(connectionlnccentered,subDataXcentered,numC);
    sublncDataX = Xloadings*Yloadings';
else
    for iiii = 1:numC
        inx_targetgene = find(connectionC(:,iiii)==1);
        if isempty(inx_targetgene),
            oneRowlncDataX = subDataX(1,:); %randomly assign
        elseif length(inx_targetgene) == 1,
            oneRowlncDataX = subDataX(inx_targetgene,:);
        else
            switch lower(zMode)
                case 'median'
                    oneRowlncDataX = median(subDataX(inx_targetgene,:));
                 case 'kmean'
                    [kmeansInx,Centriod] = kmeans(subDataX(inx_targetgene,:),2, 'distance','correlation','replicates', 10);
                    if std(Centriod(1,:)) > std(Centriod(2,:))
                        oneRowlncDataX = Centriod(1,:);
                    else
                        oneRowlncDataX = Centriod(2,:);
                    end
                case 'pca'
                    [U,S,V] = svd(subDataX(inx_targetgene,:));
                    X2 = U*S;
                    oneRowlncDataX = X2(1,:);
                case 'pls'
                    [U,S,V] = svd(subDataX(inx_targetgene,:));
                    X2 = U*S;
                    oneRowlncDataX = X2(1,:);
            end
        end
        sublncDataX = [sublncDataX; oneRowlncDataX];%%%%构建矩阵Z
    end
end
%%%%%% output Z file

str_file_name1111 = [strtok('expression','\.') '_' strtok('lncTargerList','\.') '_' zMode '_' 'ZMat_' num2str(numC) '.txt'];
fid1111 = fopen(str_file_name1111, 'wt');
fprintf(fid1111,'%s\t','lnc_NAME');
for ii = 1:M-1,
    fprintf(fid1111,'%s\t',sampleID{ii,1});
end
fprintf(fid1111,'%s\n',sampleID{M,1}); %finish firsr row

for i = 1:numC,
    fprintf(fid1111,'%s\t',lncSym{i});
    for j = 1:M-1,
        fprintf(fid1111,'%f\t',sublncDataX(i,j));
    end
    fprintf(fid1111,'%f\n',sublncDataX(i,M)); %finish one inter row        
end
fclose(fid1111);
%%%%%%%%

% rank standard derivation of Component profile across conditions and select top
% lnc with higher std
stdsublncDataX = std(sublncDataX, 0, 2);
[sortedStdsublncDataX,sortedInxsublncDataX] = sort(stdsublncDataX,'descend');
%%if numC < 100,
%%   numToplncs = numC;
%%else 
%%    numToplncs = 100;
%%end
numToplncs = numC;
 sortedsublncYSYM = lncSym(sortedInxsublncDataX(1:numToplncs));
 sortedsublncDataX = sublncDataX(sortedInxsublncDataX(1:numToplncs),:);
 sortedconnectionlnc = connectionC(:,sortedInxsublncDataX(1:numToplncs));
%sortedsublncYSYM = lncSym(1:numToplncs);
%sortedsublncDataX = sublncDataX(1:numToplncs,:);
%sortedconnectionlnc = connectionC(:,1:numToplncs);

%%%% (3) Using SCA to infer regulatory network 

 Z0 = sortedsublncDataX;
 L = numToplncs;
 sW = 0.65; %0.967; %0.96; %define sparseness
 iteration = 2000;
 fname = ['results_SCA_' num2str(L) '.mat']; 
 %[Y,Z,objhist] = nmfscfixedH( subDataX, Z, L, sW, fname, iteration);
 [Y,Z,objhist] = nmfscfixedHInitW( subDataX, sortedconnectionlnc, Z0, L, sW, nmfFlag, keepOrilncFlag, fname, iteration);
 
 
 %plot(objhist);
 
% construct clusters based on the Y  matrix nonzero elements

[gene_inx,tf_inx_repeat,val_weight] = find(Y);

MCell = cell(L,1); % MCell is a cell array in which each element contains a subcluster
clusterInx = [];
clusterInxStart = 0;
for ii = 1:L
   temp_gene_inx_array = [];
   for jj = 1:length(gene_inx)
       if tf_inx_repeat(jj) == ii
           temp_gene_inx_array = [temp_gene_inx_array gene_inx(jj)];
           clusterInx = [clusterInx;clusterInxStart];
       end
   end
   if length(temp_gene_inx_array) ~= 0
       clusterInxStart = clusterInxStart + 1;
   end
   MCell{ii} = temp_gene_inx_array;
   clear temp_gene_inx_array;
end


geneSYM = geneSymID(gene_inx);


% output cluster file
str_file_name3 = [strtok('expression','\.') '_' strtok('lnctargetlist','\.') '_' zMode '_' 'NNTopology_' num2str(L) '.txt']; 
fid3 = fopen(str_file_name3, 'wt');
fprintf(fid3,'%s\t%s\t%s\t%s\t%s\n','Source','Target','IntType','TrueNNVeri','Weight');
for i = 1:length(geneSYM),
    if Y(gene_inx(i),tf_inx_repeat(i))>0
        fprintf(fid3,'%s\t%s\t%s\t%d\t%f\n',sortedsublncYSYM{tf_inx_repeat(i)},geneSYM{i},'Activate', sortedconnectionlnc(gene_inx(i),tf_inx_repeat(i)), Y(gene_inx(i),tf_inx_repeat(i)));
    elseif Y(gene_inx(i),tf_inx_repeat(i))< 0
        fprintf(fid3,'%s\t%s\t%s\t%d\t%f\n',sortedsublncYSYM{tf_inx_repeat(i)},geneSYM{i},'Inhabit', sortedconnectionlnc(gene_inx(i),tf_inx_repeat(i)), Y(gene_inx(i),tf_inx_repeat(i)));
    else
        fprintf(fid3,'%s\t%s\t%s\t%d\t%f\n',sortedsublncYSYM{tf_inx_repeat(i)},geneSYM{i},'Noeffect', sortedconnectionlnc(gene_inx(i),tf_inx_repeat(i)), Y(gene_inx(i),tf_inx_repeat(i)));
    end
end
fclose(fid3);

% Load prior total lnc in human
%[humanlncRowID, humanlncGeneIDCell]=lnc_label;
%humanlncGeneID = humanlncGeneIDCell(:,1);
%[overlaplncListSYM, inxlncListInsubYSYM, inxInhumanlnc]= intersect(geneSymID, humanlncGeneID);
[overlaplncListSYM, inxlncListInsubYSYM, inxInhumanlnc]= intersect(geneSymID, lnc_label);
% output node attribute file
meansubDataX = mean(subDataX, 2);
typelncinsubDataX = zeros(size(geneSymID));
typelncinsubDataX(inxlncListInsubYSYM) = 1;
str_file_name1 = [strtok('expression','\.') '_' strtok('lnctargetlist','\.') '_' zMode '_' 'NodeAtt_' num2str(L) '_3.txt']; 
fid1 = fopen(str_file_name1, 'wt');
fprintf(fid1,'%s\t%s\t%s\n','Symbol','TYPE','AveExVal');
for i = 1:length(geneSymID)
    fprintf(fid1,'%s\t%d\t%d\n',geneSymID{i},typelncinsubDataX(i),meansubDataX(i));
end
fclose(fid1);

% output Y file
%str_file_name11 = [strtok('lnctargetlist','\.') '.dat'];
str_file_name11 = [strtok('lnctargetlist','\.') '.txt'];
fid11 = fopen(str_file_name11, 'wt');
fprintf(fid11,'%s\t','NAME');
for ii = 1:L-1,
    fprintf(fid11,'%s\t',sortedsublncYSYM{ii});
end
fprintf(fid11,'%s\n',sortedsublncYSYM{L}); %finish firsr row

for i = 1:N,
    fprintf(fid11,'%s\t',geneSymID{i});
    for j = 1:L-1,
        fprintf(fid11,'%f\t',Y(i,j));
    end
    fprintf(fid11,'%f\n',Y(i,L)); %finish one inter row        
end
fclose(fid11);
exit;