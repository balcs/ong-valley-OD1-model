function out = P_of_z(z,P,nuclide)

% Returns production rate (atoms/g/yr) at depth z (g/cm2)
% 
% Syntax: out = P_of_z(z,P,nuclide)
% z is the depth you want. Vector should work.
% P is a structure with P.zg (depth) and P.P (production rates) for
% interpolation.
% Nuclide is 10,21,or 26; this tells it which row in P to use. 
% 
% Accepts vector z and returns corresponding size P
%
% Greg Balco
% June 2019

if any(z > P.zg(end))
    error('P_of_z.m: z outside available range');
end

% Figure out which row of P to use
nindex = find(nuclide == [10 21 26]);
if isempty(nindex); error('P_of_z.m: nuclide??'); end

% Choose appropriate row and interpolate
out = interp1(P.zg,P.P(nindex,:),z);


