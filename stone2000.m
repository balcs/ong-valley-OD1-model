function out = stone2000(lat,P, Fsp);% Calculates the geographic scaling factor for cosmogenic-nuclide prodction as % a function of site latitude and atmospheric pressure, according to:%% Stone, J., 2000, Air Pressure and Cosmogenic Isotope Production. JGR 105:B10, % p. 23753. %% Syntax: scalingfactor = stone2000(latitude,pressure,fsp)%% Units: % latitude in decimal degrees% pressure in hPa% fsp is the fraction (between 0 and 1) of production at sea level % and high latitude due to spallation (as opposed to muons). % This argument is optional and defaults to 0.978, which is the value % used by Stone (2000) for Be-10. The corresponding value for Al-26% is 0.974. Note that using 0.844 for Be-10 and 0.826 for Al-26 will % closely reproduce the Lal, 1991 scaling factors as long as the standard% atmosphere is used to convert sample elevation to atmospheric pressure.% Also note that this function will yield the scaling factor for spallation% only when fsp=1, and that for muons only when fsp=0.  %% Elevation can be converted to pressure with the functions% stdatm.m (general use) and antatm.m (Antarctica). % % Vector argments are OK. All arguments must be the same size. %% Written by Greg Balco -- UW Cosmogenic Nuclide Lab% balcs@u.washington.edu% First version, Feb. 2001% checked March, 2006% Part of the CRONUS-Earth online calculators: %      http://hess.ess.washington.edu/math%% Copyright 2001-2007, University of Washington% All rights reserved% Developed in part with funding from the National Science Foundation.%% This program is free software; you can redistribute it and/or modify% it under the terms of the GNU General Public License, version 2,% as published by the Free Software Foundation (www.fsf.org).% check for obvious errorsif ~isempty(find(abs(lat) > 90));	error('Latitudes below 90, please');end;if length(lat) ~= length(P);	error('Vectors the same size, please');end;% default Fspif nargin == 2;	Fsp = 0.978;end;% Spallogenic production at index latitudes;% enter constants from Table 1a = [31.8518 34.3699 40.3153 42.0983 56.7733 69.0720 71.8733];b = [250.3193 258.4759 308.9894 512.6857 649.1343 832.4566 863.1927];c = [-0.083393 -0.089807 -0.106248 -0.120551 -0.160859 -0.199252 -0.207069];d = [7.4260e-5 7.9457e-5 9.4508e-5 1.1752e-4 1.5463e-4 1.9391e-4 2.0127e-4];e = [-2.2397e-8 -2.3697e-8 -2.8234e-8 -3.8809e-8 -5.0330e-8 -6.3653e-8 -6.6043e-8];ilats = [0 10 20 30 40 50 60];% calculate index latitudes at given P'slat0 = a(1) + (b(1) .* exp(P./(-150))) + (c(1).*P) + (d(1).*(P.^2)) + (e(1).*(P.^3));lat10 = a(2) + (b(2) .* exp(P./(-150))) + (c(2).*P) + (d(2).*(P.^2)) + (e(2).*(P.^3));lat20 = a(3) + (b(3) .* exp(P./(-150))) + (c(3).*P) + (d(3).*(P.^2)) + (e(3).*(P.^3));lat30 = a(4) + (b(4) .* exp(P./(-150))) + (c(4).*P) + (d(4).*(P.^2)) + (e(4).*(P.^3));lat40 = a(5) + (b(5) .* exp(P./(-150))) + (c(5).*P) + (d(5).*(P.^2)) + (e(5).*(P.^3));lat50 = a(6) + (b(6) .* exp(P./(-150))) + (c(6).*P) + (d(6).*(P.^2)) + (e(6).*(P.^3));lat60 = a(7) + (b(7) .* exp(P./(-150))) + (c(7).*P) + (d(7).*(P.^2)) + (e(7).*(P.^3));% initialize outputcorrection = zeros(size(P));% northernize southern-hemisphere inputslat = abs(lat);% set high lats to 60;lat(find(lat > 60)) = (zeros(size(find(lat > 60))) + 60);% loop b =1;while b <= length(lat);			%interpolate for actual elevation:	S(b) = interp1(ilats,[lat0(b) lat10(b) lat20(b) lat30(b) lat40(b) lat50(b) lat60(b)], lat(b));		% continue loop		b = b+1;	end;% Production by muons% constantsmk = [0.587 0.600 0.678 0.833 0.933 1.000 1.000];% index latitudes at given P'sml0 = mk(1) .* exp((1013.25 - P)./242);ml10 = mk(2) .* exp((1013.25 - P)./242);ml20 = mk(3) .* exp((1013.25 - P)./242);ml30 = mk(4) .* exp((1013.25 - P)./242);ml40 = mk(5) .* exp((1013.25 - P)./242);ml50 = mk(6) .* exp((1013.25 - P)./242);ml60 = mk(7) .* exp((1013.25 - P)./242);% loop b =1;while b <= length(lat);			%interpolate for actual elevation:	M(b) = interp1(ilats,[ml0(b) ml10(b) ml20(b) ml30(b) ml40(b) ml50(b) ml60(b)], lat(b));		% continue loop		b = b+1;	end;% Combine spallogenic and muogenic production; returnFm = 1 - Fsp;out_1 = ((S .* Fsp) + (M .* Fm));% make vectors horizontalif size(out_1,1) > size(out_1,2);	out = out_1';else;	out = out_1;end;