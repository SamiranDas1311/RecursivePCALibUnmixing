function [LibEmIdx,EvalDif]=LibUnmPcaDif(Data,Lib,NoEm)
[NoPix,NoBand]=size(Data);
[NoLibEle,NoBand]=size(Lib);
%% Calculate Nuclear Norm of the original Library Abundance Matrix
[XPca,Re]=HSIPca(Data,NoEm-1);
for i=1:NoLibEle
Y=[Data;Lib(i,:)];
[XPca1,Re1(i,1)]=HSIPca(Y,NoEm-1);
%[XPca2,RE2(i,1)]=HSIPca(Y,NoEm);
EvalDif(i,1)=Re-Re1(i,1);
i=i+1;
end
[Er,idx]=sort(EvalDif,'descend');
LibEmIdx=idx(1:NoEm,1);
end