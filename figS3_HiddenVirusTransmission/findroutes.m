function [geneFlow,wnl,adjlist] = findroutes(N,distance,adjacentM,geolocs,loc_ord,gdtable,ordLocInSqName,ordYearInSqName)
disM = distance .* adjacentM;
disM(disM==0)=Inf;
[fwdisM,fwpath]=all_shortest_paths(disM,struct('full2sparse',1,'algname','floyd_warshall'));
fwpath_ul = generate_underlying_paths(N,fwpath);
dictlist =fwpath_ul(loc_ord,loc_ord);
geolocs =geolocs(loc_ord);
loc_dis =fwdisM(loc_ord,loc_ord);
route_corr2gd = generateRandomNetworks(10000,geolocs,loc_dis,gdtable,ordLocInSqName,ordYearInSqName);
[wnl,freqm,adjlist] = orignodelist2adjm(route_corr2gd,dictlist);
geneFlow = freqm ./ sum(sum(freqm));