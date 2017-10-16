function scale = spect_correction(NDB)
%
% Calculates correction factor for fat peaks under the water peak
% such that FFcorr = scale * FF.
%

% fat chemical shifts in ppm
d = [5.29 5.19 4.2 2.75 2.2 2.02 1.6 1.3 0.9];

% Mark's heuristic fordR2las
CL = 16.8+0.25*NDB;
NDDB = 0.093*NDB^2;

% Gavin's formulas (no. protons per molecule)
awater = 2;
a(1) = NDB*2;
a(2) = 1;
a(3) = 4;
a(4) = NDDB*2;
a(5) = 6;
a(6) = (NDB-NDDB)*4;
a(7) = 6;
a(8) = (CL-4)*6-NDB*8+NDDB*2;
a(9) = 9;

% peaks 1,2,3 under water peak
fat_under_water = a(1)+a(2)+a(3);

scale = 1+fat_under_water/sum(a);

end
