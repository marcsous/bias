function [A psi] = fat_basis(te,NDB,H2O,Tesla,precession,units)
% [A psi] = fat_basis(te,NDB,H2O,Tesla,precession,units)
%
% Function that produces the fat-water matrix. Units can be adjusted
% to reflect the relative masses of water and fat. Default is proton
% suitable for PDFF calculation. For backward compatibility, NDB=-1
% uses the old fat spectrum.
%
% Inputs:
%  te (echo times in seconds)
%  NDB (number of double bonds or -1 for old spectrum)
%  H2O (water freq in ppm)
%  Tesla (field strength in Tesla)
%  precession (+1 clockwise or -1 anticlockwise)
%  units (default 'proton' or can be 'mass')
%
% Ouputs:
%  A = matrix of basis vectors [water fat]
%  psi = best fit to exp(-2*pi*psi*te): R2*=Re(psi) B0=Im(psi) in Hz
%
% Reference:
% Bydder M, Girard O, Hamilton G. Magn Reson Imaging. 2011;29:1041

% argument checks
if ~exist('units','var') || isempty(units)
    units = 'proton';
end
if ~exist('NDB','var') || isempty(NDB)
    NDB = -1;
    disp(['Warning: fat_basis.m assuming NDB=' num2str(NDB)])
    if ~isequal(units,'proton'); error('NDB=-1 requires units=proton'); end
end
if ~exist('H2O','var') || isempty(H2O)
    H2O = 4.7;
    disp(['Warning: fat_basis.m assuming H2O=' num2str(H2O)])
end
if ~exist('Tesla','var') || isempty(Tesla)
    Tesla = 3;
    disp(['Warning: fat_basis.m assuming Tesla=' num2str(Tesla)])
end
if ~exist('precession','var') || isempty(precession)
    precession = -1; % GE default is -1
    disp(['Warning: fat_basis.m assuming precession=' num2str(precession)])
end

% backward compatibility
if NDB == -1
    % old spectrum values
    d = [1.300 2.100 0.900 5.300 4.200 2.750];
    a = [0.700 0.120 0.088 0.047 0.039 0.006];
    awater = 1;
else
    % fat peak chemical shifts in ppm
    d = [0.90 1.30 1.60 2.02 2.20 2.75 4.20 5.19 5.29];
    
    % heuristic formulas (chain length and double-double bonds)
    CL = 16.8 + 0.25 * NDB;
    NDDB = 0.093 * NDB^2;
    
    % formulas for no. protons per molecule
    a(9) = NDB*2;
    a(8) = 1;
    a(7) = 4;
    a(6) = NDDB*2;
    a(5) = 6;
    a(4) = (NDB-NDDB)*4;
    a(3) = 6;
    a(2) = (CL-4)*6-NDB*8+NDDB*2;
    a(1) = 9;
    awater = 2;
end

% units above are in no. molecules, i.e. 1 unit of water = 2 protons
% and 1 unit of fat = (2+6*CL-2*NDB) protons
if isequal(units,'proton')
    % convert to proton units: 1 unit = 1 proton (used in PDFF)
    awater = awater/sum(awater);
    a = a/sum(a);
else
    % convert to mass units, i.e. 18 units of water = 2 protons and
    % (134+42*CL-2*NDB) units of fat = (2+6*CL-2*NDB) protons
    awater = awater/18;
    a = a/(134+42*CL-2*NDB);
end

% time evolution matrix (te in ms)
nte = numel(te);
te = reshape(te,nte,1);
water = repmat(awater,nte,1);
fat = zeros(nte,1);
larmor = 42.576*Tesla; % larmor freq in MHz
for j = 1:numel(d)
    freq = larmor*(d(j)-H2O); % Hz relative to water
    fat = fat + a(j)*exp(precession*2*pi*i*freq*te);
end
A = [water fat];

% best fit to exp(-2*pi*psi*te): R2*=Re(psi) B0=Im(psi) in Hz
if nargout>1
    psi = [30 precession*150*Tesla]; % initial guesses
    psi = fminsearch(@(psi)myfun(psi,te,fat),psi);
    psi = complex(psi(1),psi(2));
    %cplot(te,fat,'o');
    %te = linspace(min(te),max(te),10*numel(te));
    %hold on; cplot(te,exp(-2*pi*psi*te)); hold off;
end

% fitting function
function resnorm = myfun(psi,te,fat)
psi = complex(psi(1),psi(2)); % R2* and B0 in Hz
f = exp(-2*pi*psi*te); % complex decay
v = (f'*fat)/(f'*f); % varpro scaling
resnorm = norm(v*f-fat); % residual norm
