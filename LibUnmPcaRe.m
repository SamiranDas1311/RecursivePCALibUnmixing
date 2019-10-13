function [LibEmIdx,RecErr]=LibUnmPcaRe(Data,Lib,NoEm)
[NoPix,NoBand]=size(Data);
[NoLibEle,NoBand]=size(Lib);
%% Calculate Nuclear Norm of the original Library Abundance Matrix
[XPca,RecErr]=HSIPca(Data,NoEm-1);
DataCov=cov(Data);
for i=1:NoLibEle
%[~,~,~,DataPca,Re1(i,1)]=EigDecRank1Mod(Data,DataCov,Lib(i,:),NoEm-1);
%[N,L]=size(Data);
%Y=[Data;Lib(i,:)];
[XPca1,Re1(i,1)]=HSIPca(Data,NoEm-1);
%[XPca2,RE2(i,1)]=HSIPca(Y,NoEm);
EvalDif(i,1)=abs(RecErr-Re1(i,1));
i=i+1;
end
RecErr=EvalDif;
[Er,idx]=sort(EvalDif,'ascend');
LibEmIdx=idx(1:NoEm,1);
end