function [Loc,Year]=sparseSquenceInformation(sqArray,delimeter,ordLoc,ordYear)
for t=1:1:length(sqArray);
B{t,1}=strsplit(sqArray{t,1},delimeter);
% ID(t)=cellstr(B{t,1}{1,ordID});
Loc(t)=cellstr(B{t,1}{1,ordLoc});
y(t)=cellstr(B{t,1}{1,ordYear});
Year(t) = str2num(y{1,t});
end
end