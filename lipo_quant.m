function [F R r rsq] = lipo_quant(te,data,Tesla,init,H2O,NDB,NOISE,C,BETA)
%[F R r rsq] = lipo_quant(te,data,Tesla,init,H2O,NDB,NOISE,C,BETA)
%
% Estimates fat fraction and R2* on magnitude data (single pixel only).
%
% Inputs:
%  te is a vector of echo times (in seconds)
%  data is a vector of signal (same size as te)
%  Tesla is a the field strength (in Tesla)
%
% Inputs (optional):
%  init (scalar) is an initial estimate of FF
%  NDB (2-vector) is the no. double bonds [NDB(1)+NDB(2)*FF]
%  H2O (scalar) is the water ppm
%  NOISE (scalar) is quadrature noise bias
%  C (2-vector) is a delta R2* between water and fat [C(1)+C(2)*FF]
%  BETA (scalar) is Gaussian decay [exp(-BETA*te.^2)]
%
% Outputs:
%  F is the fitted fat fraction (range 0-1)
%  R is the fitted R2* (unit seconds)
%  r is the residual for all pixels
%  rsq is the r-squared
%
%% set up parameters
if ~exist('H2O','var') || isempty(H2O)
    H2O = 4.7; % water frequency in ppm
end
if ~exist('NDB','var') || isempty(NDB)
    NDB = -1; % no. double bonds (-1 is a flag to use old spectrum)
end
if ~exist('init','var') || isempty(init)
    init = 0.05; % initial estimate of FF
end
if ~exist('NOISE','var') || isempty(NOISE)
    NOISE = 0; % quadrature offset for noise bias
end
if ~exist('C','var') || isempty(C)
    C = [0 0]; % offset for R2* decay
end
if ~exist('BETA','var') || isempty(BETA)
    BETA = 0; % non-exp decay term
end

%% check compatiblity

% must be double precision
te = double(te);
Tesla = double(Tesla);
H2O = double(H2O);
NDB = double(NDB);
init = double(init);
NOISE = double(NOISE);
C = double(C);
BETA = double(BETA);

% catch size errors
nte = numel(te);
te = reshape(te,nte,1);
data = reshape(data,nte,1);
Tesla = reshape(Tesla,1,1);
H2O = reshape(H2O,1,1);
if numel(NDB)==1 || NDB(1)==-1
    NDB = [NDB(1) 0]; % allow scalar NDB for backward compatiblity
end
NDB = reshape(NDB,2,1);
init = reshape(init,1,1);
NOISE = reshape(NOISE,1,1);
C = reshape(C,2,1);
BETA = reshape(BETA,1,1);

% parameter ranges
if max(te)>1
    error('Argument ''te'' should be in seconds.');
end
if init<0 || init>1
    warning('Argument ''init'' should be from 0 to 1.');
end

%% setup fitting with 3 parameters: M=scale F=ff R=exp decay

% initial estimates
M = 1.25*max(data); % signal at te=0
F = init; % fat fraction (range 0-1)
R = double(25); % unit s^-1

tol = 1e-10; % fitting tolerance
residual = @(x)myfun(x,te,Tesla,H2O,NDB,NOISE,C,BETA)-data; % function to minimize
opts = optimset('display','off','algorithm','levenberg-marquardt','tolfun',tol,'tolx',tol);

%% perform fitting

[x,resnorm2,r] = lsqnonlin(residual,[M F R],[],[],opts);

% returned variables
F = x(2);
R(1) = x(3); %+C(2)*F; % R2* of water
%R(2) = x(3)+C(1); % R2* of fat

% r-squared
rsq = 1 - resnorm2./sum((data-repmat(mean(data),nte,1)).^2);

% display single pixel fit
if nargout==0
    plot(1e3*te,data,'o','color',[0 0.5 0])
    dte = min(diff(te))/3;
    te = linspace(min(te)-dte,max(te)+dte,1000)';
    signal = myfun(x,te,Tesla,H2O,NDB,NOISE,C,BETA);
    hold on; plot(1e3*te,signal,'blue'); hold off
    txt = sprintf('FF=%.3f R2*=%.1f R^2=%.3f',F,R,rsq);
    text(0.16,0.94,txt,'FontName','FixedWidth','Units','Normalized')
    ylabel('signal (a.u.)','FontName','FixedWidth');
    xlabel('TE (ms)','FontName','FixedWidth');
end


%% signal function, 3 parameters: M=scale F=ff R=exp decay
function signal = myfun(x,te,Tesla,H2O,NDB,NOISE,C,BETA)

M = x(1); % scale term
F = x(2); % fat fraction
R = x(3); % R2* decay

% apply NDB linear term
NDB = NDB(1) + NDB(2)*F;

% time evolution matrix
A = fat_basis(te,NDB,H2O,Tesla,-1);

% apply Gaussian decay
if BETA
    gaussian_decay = @(t)exp(-BETA*t.^2); % gaussian
    %sinc_decay = @(t)sinc(sqrt(6*BETA)*t/pi); % sinc
    A = bsxfun(@times,gaussian_decay(te),A);
end

% fat and water signals
f = A(:,2) * M * F;
w = A(:,1) * M * [1-F];

% R2* of fat and water (such that R2w-R2f = C1+C2*F)
Rf = R - C(1);
Rw = R + C(2)*F;

% apply R2* decay
f = f .* exp(-Rf*te);
w = w .* exp(-Rw*te);

% magnitude
signal = hypot(f+w,NOISE);
