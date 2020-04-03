function [r2,rv,Chi] = mf_stats
%----------------------------------------------------------------
% function mf_stats
%    Purpose : calculate statistics of fit from current params
%   Calls to : mf_figs, fitfunction
%
% M. Zinkin 30.11.84
%----------------------------------------------------------------
tic;
mf_msg('Calculating stats...');

%-------------- Extract data from figure -----------------------
h=findobj('tag','mf_DataWindow');
if isempty(h)
  disp('MFIT error:No data window open');
  r2 = 0; rv = 0; Chi = [ 0 0 ];
  return
end
data=get(h,'userdata');
index=data(:,4);
index = find(index ~= 0);
x=data(index,1);
y=data(index,2);
err=data(index,3);
i = find(err);
wt = 0*x;
wt(i)=1./err(i);

%-------------- Extract fit function name and dir ------------------
fitfun=get(findobj('tag','mf_FitFuncFile'),'string');

%-------------- Calculate fit -----------------------------------
if (nargout == 0)
  disp([ '* File: ' get(findobj('Tag','mf_DataFile'),'string') ]);
  fprintf(1,'  Mean=%g; Integral=%g; Std. Dev=%g\n', mean(y), trapz(x,y), std(y));
  [ym,xm]=min(y); xm=x(xm);
  fprintf(1,'  Min=%g at x=%g\n', ym, xm);
  [ym,xm]=max(y); xm=x(xm);
  fprintf(1,'  Max=%g at x=%g\n', ym, xm);
end

[p,dp,fixed]=mf_rpars; 
if isempty(p)
	return
end
hfit=findobj('Tag','mf_fitline');
if ~isempty(hfit)
	mf_msg('Getting fit results');
	yfit= get(hfit,'Ydata');
	xfit= get(hfit,'Xdata');
	if iscell(yfit)
		yfit = yfit{1};
		xfit = xfit{1};
	end
%	ycalc=interp1(xfit,yfit,x);
	ycalc=reshape(mf_interp(xfit,yfit,x),size(x));
	

else
	mf_msg('Evaluating fit function...');
	[ycalc]=feval(fitfun,x,p);
end
t=toc;

index = find(~isnan(x) & ~isnan(y) & ~isnan(err)& ~isnan(ycalc));
x=x(index);
y=y(index);
err=err(index);
ycalc=ycalc(index);
wt = wt(index);
%-------- Work out Chi squared ------------------------
% NB number of free pars isn't right: should be sum(1-dp) not length(p)
v=length(y)-length(find(~fixed));		% # of degrees of freedom (#sel points - #free pars)
ChiSq = sum( ((y-ycalc)./err).^2 );
Q=gammainc(0.5*v,0.5*ChiSq);
Chi=[ChiSq/v Q];

%---------- Display -----------------------------------
h(1)=findobj('tag','ChiSq');
h(2)=findobj('tag','QChiSq');

set(h(1),'string',num2str(Chi(1),6));
set(h(2),'string',num2str(Chi(2),2));

r=corrcoef(y.*wt,ycalc.*wt);
r2=r.*r';
r2=r2(1,2);
rv=sum(((y-ycalc).*wt).^2/length(y));
if (nargout == 0)
  disp(['* Fit with model: ' get(findobj('tag','mf_FitFuncName'),'string') '(' get(findobj('tag','mf_FitFuncFile'),'string') ')' ]);
  fprintf(1,'  Chi Squared: %g', Chi(1));
	fprintf(1,'  Correlation Coefficient for fit R^2: %g', r2)
end

mf_msg(['Done (' num2str(t) 's)']);



