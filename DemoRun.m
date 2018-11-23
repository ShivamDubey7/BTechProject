clear;
##load('./SyntheticData/data2.txt');
##load('./DATA/Demodata.mat');
##load('./DATA/Slashdot.mat');
##load('./DATA/Epinion.mat');
##load('./DATA/Wikielect.mat');
load('./DATA/Wordnet.mat');

## Run KOCG_main.m
global alpha beta;
B = -A;
A(A < 0) = 0;
B(B < 0) = 0;

alpha = 0.9;
beta = 0.5;
k_start = 10;
k_iter = 10;
k_limit = 11;


log_file = 'graph_demorun.txt';
version = 0; %0 for old and 1 for new
if version == 0
  log_file = strcat('LogFile_main_prev_',log_file);
else
  log_file = strcat('LogFile_main2_',log_file);
end
fileID = fopen(log_file,'a');
fprintf(fileID, '%f\t%f\t%f\n',NV,alpha,beta);
fclose(fileID);
NG = k_start;
iter = k_iter;
while NG < k_limit
  if version == 0
    [X_enumKOCG_cell, time,finalAns1] = KOCG_main_prev(A, B, NG, alpha, beta);
  else
    [X_enumKOCG_cell, time,finalAns1] = KOCG_main2(A, B, NG, alpha, beta);
  end
  X = zeros(NV,NG);
  for row=1:size(X_enumKOCG_cell,1)
    for col=1:size(X_enumKOCG_cell,2)
      [x,y,z] = find(X_enumKOCG_cell(row,col){1,1});
      for k=1:length(x)
        X(x(k),y(k))=z(k);
      end
    end
  end
  U = ones(NG,NG) - eye(NG);
  fx_val = trace(X'*A*X) + alpha*trace(X'*B*X*U) - beta*trace(X'*X*U);
  fileID = fopen(log_file,'a');
  fprintf(fileID,'%d\t%f\t%f\n',NG,fx_val,time); 
  fclose(fileID);
  NG = NG + k_iter;
end


log_file = 'graph_demorun.txt';
version = 1; %0 for old and 1 for new
if version == 0
  log_file = strcat('LogFile_main_prev_',log_file);
else
  log_file = strcat('LogFile_main2_',log_file);
end
fileID = fopen(log_file,'a');
fprintf(fileID, '%f\t%f\t%f\n',NV,alpha,beta);
fclose(fileID);
NG = k_start;
iter = k_iter;
while NG < k_limit
  if version == 0
    [X_enumKOCG_cell, time,finalAns1] = KOCG_main_prev(A, B, NG, alpha, beta);
  else
    [X_enumKOCG_cell, time,finalAns1] = KOCG_main2(A, B, NG, alpha, beta);
  end
  X = zeros(NV,NG);
  for row=1:size(X_enumKOCG_cell,1)
    for col=1:size(X_enumKOCG_cell,2)
      [x,y,z] = find(X_enumKOCG_cell(row,col){1,1});
      for k=1:length(x)
        X(x(k),y(k))=z(k);
      end
    end
  end
  U = ones(NG,NG) - eye(NG);
  fx_val = trace(X'*A*X) + alpha*trace(X'*B*X*U) - beta*trace(X'*X*U);
  fileID = fopen(log_file,'a');
  fprintf(fileID,'%d\t%f\t%f\n',NG,fx_val,time); 
  fclose(fileID);
  NG = NG + k_iter;
end




##temp=sum(X,1);
##temp
##for i=1:size(X,2)
##  fprintf(fileID,'%d\t',temp(i));
##end




##clear;
##NV = 20;
##NG = 2;
##B = [];
##sparse_thres = 0.99;
##[A, GT_X] = genSignedGraphOptEqSize(NV, NG, sparse_thres);
##save data3.txt NV NG sparse_thres A B;







