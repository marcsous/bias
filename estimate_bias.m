function estimate_bias()
%
% Program to estimate bias parameters in PDFF quantification.
% It reads in 2 data sets with 6 or 16 echos and spectroscopy.
%
% Ref. Sources of Systematic Error in Proton Density Fat Fraction
% Quantification in the Liver Evaluated From Magnitude Images With
% Different Numbers of Echoes. NMR Biomed 2017
%
%% load 2 data sets

% structure containing subject data
subject = []; 

for k = 1:2

    if k==1; load data16echo.mat; end
    if k==2; load data6echo.mat; end

    for n = 1:size(data,1)
        subject(end+1).fa  = fa(n);     % flip angle (deg)
        subject(end).tr    = tr(n);     % repetition time (ms)
        subject(end).te    = te(n,:);   % echo times (s)
        subject(end).Tesla = Tesla;     % field strength (T)
        subject(end).data  = data(n,:); % magnitude signal
        subject(end).noise = noise(n);  % background noise
        subject(end).spect = spect(n);  % uncorrected spect PDFF
    end
   
end

% exclude subject with very high R2*
subject(95) = [];


%% clean up workspace

close all
clearvars -except subject

N = numel(subject);
M = numel([subject(:).data]);

disp('SUBJECTS')
fprintf(' no. subjects = %i\n',N);
fprintf(' no. data points = %i\n',M);
fprintf(' tr = %.1f-%.1f ms\n',min([subject(:).tr]),max([subject(:).tr]))
fprintf(' te = %.2f-%.2f ms\n',1e3*min([subject(:).te]),1e3*max([subject(:).te]))
fprintf(' fa = %.1f-%.1f deg\n',min([subject(:).fa]),max([subject(:).fa]))
fprintf(' Tesla = %.1f-%.1f T\n',min([subject(:).Tesla]),max([subject(:).Tesla]))


%% normalize signal by background noise

for n = 1:N
   
    % background noise is 2x higher at edge of FOV than at center
    noise_fudge_factor = 2;
    subject(n).noise = subject(n).noise / noise_fudge_factor;

    % divide by background noise to normalize variance between subjects
    subject(n).data = subject(n).data / subject(n).noise;

end


%% parameter fitting - see myfunc() and lipo_quant() for details
disp('PARAMETER FITTING')

tol = 1e-10; % fitting tolerance
x0 = [1000 0 0 4.7 1 1.92 0]; % initial estimates [Gaussian C0 C1 H20 NOISE NDB0 NDB1]
opts = optimset('display','off','tolfun',tol,'tolx',tol,'algorithm','levenberg-marquardt','findifftype','central','findiffrelstep',sqrt(tol));

% step through discrete values of NDB1
NDB1 = 1:-1:-6;

for k = 1:numel(NDB1)

    % use eval string to apply constraint
    constraint = ['x(7)=' num2str(NDB1(k))];
    
    % nonlinear fitting (returns resnorm^2)
    [x resnorm2] = lsqnonlin(@(x)myfunc(x,subject,constraint),x0,[],[],opts);
    
    % overwrite whatever lsqnonlin left in the constrained variable(s)
    evalc(constraint);
    
    % save fitting results
    xsave(k,:) = [x sqrt(resnorm2)];
    
    % display table
    fprintf('\n constraint: %s\n',constraint);
    disp(array2table(xsave,'VariableNames',{'BETA','C0','C1','H2O','NOISE','NDB0','NDB1','rnorm'}));
    
end

% RESULT FROM PREVIOUS FITTING (GAUSSIAN)
%1005.592  -47.09085   218.1014   4.687844   0.7071270   2.27142   1   5.385011
%1007.019  -45.11461   209.6835   4.687813   0.7040253   2.50437   0   5.367468
%1008.279  -43.21209   201.4985   4.687797   0.7005217   2.73707  -1   5.353403
%1009.379  -41.38465   193.5591   4.687793   0.6966518   2.96959  -2   5.342758
%1010.326  -39.63384   185.8779   4.687800   0.6924447   3.20197  -3   5.335472
%1011.124  -37.96132   178.4671   4.687817   0.6879236   3.43423  -4   5.331481
%1011.775  -36.36882   171.3382   4.687843   0.6831065   3.66636  -5   5.330717
%1012.281  -34.85810   164.5020   4.687877   0.6780070   3.89837  -6   5.333110

% fit again to get 95% confidence intervals
constraint = 'x(7)=-1';
[x,resnorm2,~,~,~,~,J] = lsqnonlin(@(x)myfunc(x,subject,constraint),x0,[],[],opts);

v = M-3*N-numel(x)-1; % degrees of freedom
cov = pinv(full(J'*J))*resnorm2/v;
ci95 = sqrt(max(diag(cov),0))*1.96;

% clean up results
x(7)=-1; ci95(7)=0;

% display
fprintf('\n')
disp([' x(1): BETA  = ' num2str(x(1),'%f') ' ± ' num2str(ci95(1),'%f')])
disp([' x(2): C0    = ' num2str(x(2),'%f') ' ± ' num2str(ci95(2),'%f')])
disp([' x(3): C1    = ' num2str(x(3),'%f') ' ± ' num2str(ci95(3),'%f')])
disp([' x(4): H2O   = ' num2str(x(4),'%f') ' ± ' num2str(ci95(4),'%f')])
disp([' x(5): NOISE = ' num2str(x(5),'%f') ' ± ' num2str(ci95(5),'%f')])
disp([' x(6): NDB0  = ' num2str(x(6),'%f') ' ± ' num2str(ci95(6),'%f')])
disp([' x(7): NDB1  = ' num2str(x(7),'%f') ' ± ' num2str(ci95(7),'%f')])
fprintf(' resnorm = %.9f\n',sqrt(resnorm2));


%% add effects one at a time

% fit using up to 16 echos
for e = 1:16
    
    parfor n = 1:N
        
        % current subject
        te = subject(n).te;
        data = subject(n).data;
        Tesla = subject(n).Tesla;
        init = subject(n).spect;
        
        if e < 3 || e > numel(te)
            
            % mark as missing data
            F1(n,e) = NaN; R1(n,e) = NaN;
            F2(n,e) = NaN; R2(n,e) = NaN;
            F3(n,e) = NaN; R3(n,e) = NaN;
            F4(n,e) = NaN; R4(n,e) = NaN;
            F5(n,e) = NaN; R5(n,e) = NaN;
            F6(n,e) = NaN; R6(n,e) = NaN;
            
        else

            % starting values (original model)
            H2O = 4.7;
            NDB = -1; % -1 is a flag to use the old Hamilton spectrum
            NOISE = 0;
            C = [0 0];
            BETA = 0;
            
            % original
            [F1(n,e) R1(n,e) tmp1 tmp2] = lipo_quant(te(1:e),data(1:e),Tesla,init,H2O,NDB,NOISE,C,BETA);
            if e==numel(te); r1{n} = tmp1; rsq1(n) = tmp2; end % residual and r-squared with max. echos

            % Gaussian
            BETA = x(1);
            [F2(n,e) R2(n,e) tmp1 tmp2] = lipo_quant(te(1:e),data(1:e),Tesla,init,H2O,NDB,NOISE,C,BETA);
            if e==numel(te); r2{n} = tmp1; rsq2(n) = tmp2; end % residual and r-squared with max. echos
            
            % NDB
            NDB = [x(6) x(7)];
            [F3(n,e) R3(n,e) tmp1 tmp2] = lipo_quant(te(1:e),data(1:e),Tesla,init,H2O,NDB,NOISE,C,BETA);
            if e==numel(te); r3{n} = tmp1; rsq3(n) = tmp2; end % residual and r-squared with max. echos
                        
            % R2* offset
            C = [x(2) x(3)];
            [F4(n,e) R4(n,e) tmp1 tmp2] = lipo_quant(te(1:e),data(1:e),Tesla,init,H2O,NDB,NOISE,C,BETA);
            if e==numel(te); r4{n} = tmp1; rsq4(n) = tmp2; end % residual and r-squared with max. echos
            
            % H2O
            H2O = x(4);
            [F5(n,e) R5(n,e) tmp1 tmp2] = lipo_quant(te(1:e),data(1:e),Tesla,init,H2O,NDB,NOISE,C,BETA);
            if e==numel(te); r5{n} = tmp1; rsq5(n) = tmp2; end % residual and r-squared with max. echos
            
            % noise bias
            NOISE = x(5);
            [F6(n,e) R6(n,e) tmp1 tmp2] = lipo_quant(te(1:e),data(1:e),Tesla,init,H2O,NDB,NOISE,C,BETA);
            if e==numel(te); r6{n} = tmp1; rsq6(n) = tmp2; end % residual and r-squared with max. echos
            
        end
        
    end
    
end


%% resnorm - fitted using all echos
disp('RESIDUAL NORMS')

% residual norms after each correction (using all echos)
r1 = cell2mat(r1(:));
r2 = cell2mat(r2(:));
r3 = cell2mat(r3(:));
r4 = cell2mat(r4(:));
r5 = cell2mat(r5(:));
r6 = cell2mat(r6(:));
resnorm = [norm(r1) norm(r2) norm(r3) norm(r4) norm(r5) norm(r6)];

fprintf(' original   = %f\n',resnorm(1));
fprintf(' optimized  = %f\n',resnorm(end));
fprintf(' difference = %f (%.1f%%)\n',resnorm(1)-resnorm(end),100*(resnorm(1)-resnorm(end))/resnorm(1));

disp([' norm   = ' num2str(resnorm,'%.4f ')])
disp([' diff   = ' num2str([0 -diff(resnorm)],'%.4f ')])
disp([' diff(%)=  ' num2str(100*[0 diff(resnorm)./(resnorm(end)-resnorm(1))],' %.1f%% ')])


%% r-squared
disp('QUALITY OF FIT')

fprintf(' original R^2  = %f\n', mean(rsq1));
fprintf(' optimized R^2 = %f\n', mean(rsq6));
fprintf(' improvement in R^2: %i out of %i\n',sum(rsq6>rsq1),N)


%% spect ndb correction
disp('SPECT NDB CORRECTION')

% ndb correction for fat peaks under water
old_correction = 1.087; % Gavin's spectrum
new_correction = spect_correction(x(6)); % just NDB0, ignore NDB1 

% undo old correction and apply new correction to spect
for n = 1:N
    spect_ndb_corr(n) = subject(n).spect * new_correction / old_correction;
end

disp([' Old spectroscopy fat under water peak: ' num2str(100*(old_correction-1)) '%'])
disp([' New spectroscopy prior knowledge: ndb = ' num2str(x(6))])
disp([' New spectroscopy fat under water peak: ' num2str(100*(new_correction-1)) '%'])


%% T1 correction
disp('T1 CORRECTION')

% assume T1s in liver
T1f = 260; % fat
T1w = 790; % water

% spect T1 correction
for n = 1:N

    fa = 90 * (pi/180); % spect fa is always 90
    TR = 3500; % spect tr is always 3500
    
    % theoretical fat and water signals
    ftrue = -0.1:0.01:0.5;
    wtrue = 1-ftrue;
    f = ftrue * (1-exp(-TR/T1f))*sin(fa)./(1-exp(-TR/T1f)*cos(fa));
    w = wtrue * (1-exp(-TR/T1w))*sin(fa)./(1-exp(-TR/T1w)*cos(fa));
    
    % theoretical fat fractions
    Fmeas = f./(f+w);
    Ftrue = ftrue./(ftrue+wtrue);

    % T1-corrected spect fat fractions
    spect_T1_corr(n) = interp1(Fmeas,Ftrue,subject(n).spect);
    spect_ndb_T1_corr(n) = interp1(Fmeas,Ftrue,spect_ndb_corr(n));

    % absolute error in PDFF
    E(n) = spect_ndb_corr(n) - spect_ndb_T1_corr(n);

end
fprintf(' Max. absolute PDFF error due to T1 (spect) = %f%%\n',100*norm(E,Inf));

% imaging T1 correction
for n = 1:N
    
    fa = subject(n).fa * (pi/180);
    TR = subject(n).tr;
    
    % theoretical fat and water signals
    ftrue = -0.1:0.01:0.5;
    wtrue = 1-ftrue;
    f = ftrue * (1-exp(-TR/T1f))*sin(fa)./(1-exp(-TR/T1f)*cos(fa));
    w = wtrue * (1-exp(-TR/T1w))*sin(fa)./(1-exp(-TR/T1w)*cos(fa));
    
    % theoretical fat fractions
    Fmeas = f./(f+w);
    Ftrue = ftrue./(ftrue+wtrue);

    % amplification factor (A=Fmeas/Ftrue) evaluated at spect_ndb_T1_corr
    A = interp1(Ftrue,Fmeas,spect_ndb_T1_corr(n)) / spect_ndb_T1_corr(n);
    
    % T1-corrected imaging fat fractions
    F1(n,:) = F1(n,:) / A;
    F2(n,:) = F2(n,:) / A;
    F3(n,:) = F3(n,:) / A;
    F4(n,:) = F4(n,:) / A;
    F5(n,:) = F5(n,:) / A;
    F6(n,:) = F6(n,:) / A;
    
    % absolute error in PDFF
    E(n) = (A-1) * spect_ndb_T1_corr(n);

end
fprintf(' Max. absolute PDFF error due to T1 (imaging) = %f%%\n',100*norm(E,Inf));


%% correlations with no. echos
disp('CORRELATION PDFF vs NECHO')

% system matrix
X = ones(N,1)*(1:16);
X = reshape(X,[],1);
X = [ones(size(X)) X];

% regress for each model (treats NaN as missing data)
for m = 1:6
    if m==1; Y = reshape(F1,[],1); end % original
    if m==2; Y = reshape(F2,[],1); end % Gaussian
    if m==3; Y = reshape(F3,[],1); end % NDB
    if m==4; Y = reshape(F4,[],1); end % R2*
    if m==5; Y = reshape(F5,[],1); end % H2O
    if m==6; Y = reshape(F6,[],1); end % noise bias
     [B,~,~,~,stats] = regress(Y,X);
    fprintf(' PDFF: intercept=%f slope=%f r=%f p=%f\n',B(1),B(2),sqrt(stats(1)),stats(3));
end

disp('CORRELATION R2* vs NECHO')
for m = 1:6
    if m==1; Y = reshape(R1,[],1); end % original
    if m==2; Y = reshape(R2,[],1); end % Gaussian
    if m==3; Y = reshape(R3,[],1); end % NDB
    if m==4; Y = reshape(R4,[],1); end % R2*
    if m==5; Y = reshape(R5,[],1); end % H2O
    if m==6; Y = reshape(R6,[],1); end % noise bias
    [B,~,~,~,stats] = regress(Y,X); 
    fprintf(' R2*: intercept=%f slope=%f r=%f p=%f\n',B(1),B(2),sqrt(stats(1)),stats(3));
end


%% display plots

% trends vs no. echoes (show 16 echo only)
figure('Position',[19 598 934 380]);

% get mask for the 16 echo data
k = false(N,1);
for n = 1:N
    k(n) = numel(subject(n).te)==16;
end

subplot(1,2,1)
h = plot(100*[mean(F1(k,:));mean(F2(k,:));mean(F3(k,:));mean(F4(k,:));mean(F5(k,:));mean(F6(k,:))]');
xlabel('No. Echoes'); ylabel('Mean PDFF (%)');
axis([2 17 9.55 12.3]);
legend({'Original','Gaussian','No. double bonds','R2* offset','Water freq.','Noise bias'},'location','southwest');
title('PDFF vs No. Echoes')

subplot(1,2,2)
g = plot([mean(R1(k,:));mean(R2(k,:));mean(R3(k,:));mean(R4(k,:));mean(R5(k,:));mean(R6(k,:))]');
xlabel('No. Echoes'); ylabel('Mean R2* (s^{-1})');
title('R2* vs No. Echoes')
axis([2 17 12 63]);

h(1).Marker = 'o'; g(1).Marker = 'o';
h(2).Marker = '.'; g(2).Marker = '.';
h(3).Marker = '*'; g(3).Marker = '*';
h(4).Marker = '+'; g(4).Marker = '+';
h(5).Marker = 'diamond'; g(5).Marker = 'diamond';
h(6).Marker = 'square'; g(6).Marker = 'square';

% example fits (show 16 echo only)
figure('Position',[10 12 727 995]);

examples = [74 22 28 61];

subplot(numel(examples),2,1)
for j = 1:numel(examples)
    
    te = subject(examples(j)).te;
    data = subject(examples(j)).data;
    Tesla = subject(examples(j)).Tesla;
    subplot(numel(examples),2,1+2*(j-1))
    init = subject(examples(j)).spect;
    
    % original values
    H2O = 4.7;
    NDB = -1; % -1 is a flag to use the old spectrum
    NOISE = 0;
    C = [0 0];
    BETA = 0;
    lipo_quant(te,data,Tesla,init,H2O,NDB,NOISE,C,BETA);
    
    axis tight; ax1 = axis;
    subplot(numel(examples),2,2+2*(j-1))
    
    % optimized values
    H2O = x(4);
    NDB = [x(6) x(7)];
    NOISE = x(5);
    C = [x(2) x(3)];
    BETA = x(1);
    lipo_quant(te,data,Tesla,init,H2O,NDB,NOISE,C,BETA); 
    
    axis tight; ax2 = axis;
    subplot(numel(examples),2,1+2*(j-1))
    xlim([0 20]); ylim([min(ax1(3),ax2(3)),max(ax1(4),ax2(4))]);
    subplot(numel(examples),2,2+2*(j-1))
    xlim([0 20]); ylim([min(ax1(3),ax2(3)),max(ax1(4),ax2(4))]);
    
end
subplot(numel(examples),2,1); title('Original');
subplot(numel(examples),2,2); title('Optimized');


%% compare PDFF with spect
disp('BLAND ALTMAN')

% use PDFF and R2* from max no. echos
for n = 1:N
    nte = numel(subject(n).te);
    F1_all_echos(n) = F1(n,nte);
    F6_all_echos(n) = F6(n,nte);
    R1_all_echos(n) = R1(n,nte);
    R6_all_echos(n) = R6(n,nte);
end

% original PDFF must be compared with original spect values (same ndb)
figure('Position',[19 598 934 380]);

subplot(1,2,1);
bland_altman(100*spect_T1_corr,100*F1_all_echos,'.');
title('Original')
xlabel('(MRS + MRI) / 2 (%)')
ylabel('MRS - MRI (%)')
axis([0 40 -5 6])

subplot(1,2,2);
bland_altman(100*spect_ndb_T1_corr,100*F6_all_echos,'.');
title('Optimized')
xlabel('(MRS + MRI) / 2 (%)')
ylabel('MRS - MRI (%)')
axis([0 40 -5 6])

%% compare R2* with PDFF
disp('R2* VERSUS PDFF')

% original R2 must be compared with original spect values (same ndb)
figure('Position',[19 598 934 380]);

subplot(1,2,1);
fit_linear(100*spect_T1_corr,R1_all_echos,'.');
title('Original'); xlabel('MRS PDFF (%)'); ylabel('R2* (s^{-1})')
axis([-1 41 -50 150])

subplot(1,2,2);
fit_linear(100*spect_ndb_T1_corr,R6_all_echos,'.');
title('Optimized'); xlabel('MRS PDFF (%)'); ylabel('R2* (s^{-1})')
axis([-1 41 -50 150])



keyboard

end

%% residual generating function
function r = myfunc(x,subject,constraint)

% apply constraint, e.g. 'x(7)=0'
if exist('constraint','var')
    evalc(constraint);
end

% bias parameters
BETA  = x(1);
C     = [x(2) x(3)];
H2O   = x(4);
NOISE = x(5);
NDB   = [x(6) x(7)];

% fit subject specific parameters (scale, PDFF, R2*)
parfor n = 1:numel(subject)
    te = subject(n).te;
    data = subject(n).data;
    Tesla = subject(n).Tesla;
    init = subject(n).spect; % initial estimate for PDFF
    [~,~,r{n}] = lipo_quant(te,data,Tesla,init,H2O,NDB,NOISE,C,BETA);
end

% total residual in vector form
r = cell2mat(r(:));

end
