function [cdffreq,cdfval] = ecdf(in)

cdffreq = (1:length(in))./length(in);
cdfval = sort(in);

