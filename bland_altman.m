function [ave err] = bland_altman(x,y,marker)
%[ave err] = bland_altman(x,y,marker)
%
% Creates Bland-Altman plots.
% Treats NaN as missing data.

% ignore NaN
x = reshape(x,[],1);
y = reshape(y,[],1);
ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);

% marker
if ~exist('marker','var')
    marker = 'o';
end

% difference and mean
d = (x-y);
m = (x+y)/2;

plot(m,d,marker);

ave = mean(d);
err = std(d)*1.96;

line([floor(min(m)) ceil(max(m))],[ave ave],'color','r')
line([floor(min(m)) ceil(max(m))],[ave+err ave+err],'color','r','linestyle',':','linewidth',1)
line([floor(min(m)) ceil(max(m))],[ave-err ave-err],'color','r','linestyle',':','linewidth',1)

title([num2str(ave) ' ± ' num2str(err) ' (95% CI)'])
xlabel('(arg1+arg2) / 2')
ylabel('arg1-arg2')

%text(0.05,0.9,str,'Units','Normalized','FontName','FixedWidth')

axis tight;
ax = axis;
offsetx = -(ax(2)-ax(1))/8;
offsety = err/10;
offsetn = offsetx/4; % extra space for negative sign

text(double(max(m)+offsetx+offsetn*(ave<0)),double(ave+offsety),num2str(ave,'%.2f'),'FontName','FixedWidth');
text(double(max(m)+offsetx+offsetn*(ave+err<0)),double(ave+err+offsety),num2str(ave+err,'%.2f'),'FontName','FixedWidth');
text(double(max(m)+offsetx+offsetn*(ave-err<0)),double(ave-err+offsety),num2str(ave-err,'%.2f'),'FontName','FixedWidth');

if nargout==0
    fprintf(' bias=%f ci95=%f\n',ave,err);
    clear;
end
