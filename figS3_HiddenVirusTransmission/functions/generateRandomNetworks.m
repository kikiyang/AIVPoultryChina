% �������ڽӹ�ϵ�����������е�·��������󣬰�����ԭʼ����ʦ���㷨���м���õ�����Ծ�͵�·��
% when no examination on maximum diversity, no mean evolutionary rate input [on conditions of internal genes] 
function route_corr2gd=generateRandomNetworks(networkNum,geolocs,geodis,genedis,ordLocInSqName,ordYearInSqName,varargin)
narginchk(6,7)

route_corr2gd = zeros(networkNum,1);
if(length(geolocs)==length(geodis))
    n=length(geolocs);
end
% genedis = readtable(genedisfilename,'ReadRowNames',true,'ReadVariableNames',false);
% range = strcat('A1:',transNumber2ExcelColumn(width(genedis),width(genedis)+1));
% genedis = readtable(genedisfilename,'Range',range,'ReadRowNames',true,'ReadVariableNames',false);
GeneticDis=table2array(genedis);
[loc,year]=sparseSquenceInformation(genedis.Properties.RowNames,'/',ordLocInSqName,ordYearInSqName);
sqTimeInterval = zeros(width(genedis),width(genedis));
for r=1:1:width(genedis)
    for c=1:1:width(genedis)
       sqTimeInterval(r,c) = abs(year(r)-year(c));
    end
end

divergence = GeneticDis ./ (sqTimeInterval + 1);
Div = divergence(:);

% TODO #1 ȷ���Ƿ���Ҫ���diversity [No examination for internal genes]
% ind: divergence > maximum possible div (in one year)������
if (isempty(varargin))
    ind = isnan(Div);
else
    TimInt = sqTimeInterval(:);
    ind = find((TimInt ==0 & Div>varargin{1})| isnan(Div));
end

GeneDis_v = GeneticDis(:);
GeneDis_v(ind) = [];

for i=1:1:networkNum
    % flag = 1;
    % if edge not exist, randperm again
    % while (flag~=0)
    idx=randperm(n);
    t=1:1:n-1;
%     idx2=idx(t+1);
%     idx2(n)=idx(1);
    geodisnet = geodis((idx(t)),(idx(t+1)));
    NetworkDis = diag(geodisnet);
    EdgesLen=sum(NetworkDis);
    %     flag = nnz(NetworkDis==0);
    % end
    
    NetDis=diag(0,n-1);
    t= 2:(n+1):(n^2-n);
    k=1:1:n-1;
    NetDis(t)=NetworkDis(k);
    for k=3:n
        for t=k:1:n
            NetDis(t,t-k+1)=NetDis(t,t-1)+NetDis(t-1,t-k+1);
        end
    end
    %�Գƾ���
    NetDistance=NetDis+NetDis';
    %���Գƾ���
    % NetDistance_in=NetDis-NetDis'+triu(repmat(EdgesLen,n,n));
    % NetDistance_in(1:n+1:end)=diag(NetDistance);
    
    t=1:1:n;
    GeoLocationNames(t) = geolocs(idx(t));
    NetDistance(t*n-n+t) = EdgesLen;
    %��loc���ֵ��Ӧ��GeoLocationNames��orderȥ
    loc_idx = arrayfun(@(x)find(strcmp(GeoLocationNames,x),1),loc);
    
    sqRouteDis = zeros(width(genedis),width(genedis));
    % sqTimeInterval = zeros(width(genedis),width(genedis));
    for r=1:1:width(genedis)
        for c=1:1:width(genedis)
            %        sqRouteDis(r,c) = NetDistance_in(loc_idx(r),loc_idx(c));
            sqRouteDis(r,c) = NetDistance(loc_idx(r),loc_idx(c));
            %        sqTimeInterval(r,c) = abs(year(r)-year(c));          
        end
    end
    % ����cycle -> R2 ���ϵ����� ��ΪNetworkDistance���и��������ͬ����
    sqCycleDis = sqTimeInterval * EdgesLen;
    NetworkDistance = sqRouteDis + sqCycleDis;
    % ����cycle,�㣺ʱ����������ܵ�distance
%     NetworkDistance = sqRouteDis .* (sqTimeInterval+1);
    
    NetworkDis_v = NetworkDistance(:);
    NetworkDis_v(ind) = [];
    R = corrcoef(GeneDis_v,NetworkDis_v);
%     [R,P] = corrcoef(GeneDis_v,NetworkDis_v);
    route_corr2gd(i) = R(1,2);
%     correlation(i,2) = P(2,1);
    t =1:1:n;
    route_corr2gd(i,t+1) = idx(t);
%     correlation(i)=corr(GeneticDis(:),NetworkDistance(:),'rows','all');
end
end