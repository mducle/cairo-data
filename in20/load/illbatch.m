function [x, y, err, xlab, ylab,monitor,optpars]=illbatch(filespec)
%function [x, y, err, xlab, ylab,monitor,optpars]=illbatch(filespec)
%
% MVIEW/MFIT load routine for ILL triple axis data files. 
% This routine either accepts a simple filename input 
% argument, or alternatively a compound file
% specification of the form:
%
%      filename,{options,...}
%
% the 'setpar' option tries to set Rescal parameters and mfit 'Temperature',
% 'Energy' and 'H/K/Lzero' parameters.
%
% Valid options are =
%      X=string     name for X, use '-' for auto setting, '#' for point number
%      Y=string     name for Y
%      M=string     name for monitor, can also be : none, use '-' for auto setting
%      E=string     name for Y error, can also be : sqrt(y) or none
%      N=number     normalisation value, can be an expression
%      S=number     scan number
%      P=string     name for Polarisation Analysis. 'PAL' is default.
%      F=number     flip number in PAL column
%      F1=number (0 or 1)  first flipper state
%      F2=number (0 or 1)  second flipper state
%      gui          ask user with GUI
%      setpar       set parameters if possible (in rescal or mfit)
%      silent       set no output mode
%      
% example :
%      batch('foo.txt,X=2:10,M=-,E=sqrt(y),Y=Counts,N=norm(y),gui')
% optpars = [ ki kf ti mn PALval ];

% MZ 11.10.95 and DFM 2.11.95, EF 4.09.97 (update), ARW 21.10.98(polarisation)
%
% this file is composed of those different parts :
% * filespec analysis, retrieval of options
% * load data and datastr (plus possible options)
% * display file info for user                 <- file format dependent
% * extract column guess or choice              <- file format dependent
% * possibly call GUI                          <- programer's choice
% * last column check
% * extract x,y,error,monitor,labels
% * options : automatic set params,...         <- file format dependent
%
% for GUI, you MIGHT use mf_coldg (x,y,error,monitor,scan,norm selector) at your choice.


%===== Parse filespec, open data file, and get column names ============================

%---- Set default column names and flags

scan='';
xname='-';
yname='CNTS';
mname='none';
ename='-';
normf=1; 
PALname='PAL';
%F1name='F1';
%F2name='F2';
F1name='f1';
F2name='f2';
F1val = [];
F2val = [];
PALval = [];

setpars = 0;  % automatic set of rescal pars
gui = 0;    % no gui choice

x=[]; y= []; err=[]; ylab=''; xlab=''; monitor = ''; Fstate = []; optpars = [];

%----- Parse filespec --------------------------------------

[fspec filespec]=strtok(filespec,',');
while ~isempty(filespec)
   [s filespec]=strtok(filespec,',');
   s = s(find(~isspace(s)));
   fspec=str2mat(fspec,s);
end
[nargs,nchars]=size(fspec);

%----- Update scan parameters from filespec---------------------------

i=strmatch('X=',fspec);
if ~isempty(i)
  xname=deblank(fspec(i(end),3:nchars)); end

i=strmatch('Y=',fspec);
if ~isempty(i)
  yname=deblank(fspec(i(end),3:nchars)); end

i=strmatch('M=',fspec);
if ~isempty(i)
  mname=deblank(fspec(i(end),3:nchars)); end

i=strmatch('N=',fspec);
if ~isempty(i)
  normf=deblank(fspec(i(end),3:nchars)); end

i=strmatch('E=',fspec);
if ~isempty(i)
  ename=deblank(fspec(i(end),3:nchars)); end


i=strmatch('S=',fspec);
if ~isempty(i)
  scan=deblank(fspec(i(end),3:nchars)); end

i=strmatch('P=',fspec);
if ~isempty(i)
  PALname=deblank(fspec(i(end),3:nchars)); end

i=strmatch('F=',fspec);
if ~isempty(i)
  PALval=sscanf(fspec(i(end),3:nchars),'%d'); end

i=strmatch('F1=',fspec);
if ~isempty(i)
  F1val=sscanf(fspec(i(end),4:nchars),'%d'); end

i=strmatch('F2=',fspec);
if ~isempty(i)
  F2val=sscanf(fspec(i(end),4:nchars),'%d'); end

filename=deblank(fspec(1,:));

if strmatch('gui',fspec)
  gui = 1;
end

if strmatch('setpar',fspec)
  setpars=1;
end

if strmatch('silent',fspec)
  silent = 1;
else
  silent = 0;
end

if strmatch('file_list',fspec)
  file_index = fspec;
else
  file_index = [];
end

%===== Get data from file ============================

if ~silent, disp([ 'Loading ILL data : ' filename ]); end
[data,datastr,pres,pscan,stepvec,stepstr,com,header,dathead,flip]=illdata(filename);
if isempty(data) return; end    % load error
if isempty(datastr)
  [nchar,ncolumns] = size(data);
  datastr=[];
  for i=1:ncolumns
    datastr=str2mat(datastr,sprintf('col%i',i));
  end
  datastr(1,:)=[];
end

PALcol=eval(PALname,'strmatch(lower(PALname),lower(datastr),''exact'')');

%===== Establish polarisation states

if flip & isempty(PALcol)
  pal = length(data(:,1))/length(unique(data(:,1)));
  PALcol = repmat(1:pal, [1 max(data(:,1))]);
  PALcol = PALcol';
elseif ~isempty(PALcol)
  pal=max(data(:,PALcol));
  PALcol = data(:,PALcol);
end

if length(datastr) ~= size(data,2)
  disp(['Warning: ' filename 'the column names do not match the data column number']);
  foundPAL = 0;
  for i=1:size(data,2)
    if all(data(:,i) == PALcol)
      foundPAL = i;
    end
  end
  if foundPAL
    if foundPAL > 1
      datastr = str2mat(datastr(1:(foundPAL-1),:), 'PAL', datastr(foundPAL:end,:));
    else 
      datastr = str2mat('PAL', datastr(1:end,:));
    end
  end
end

F1col=eval(lower(F1name),'strmatch(lower(F1name),lower(datastr),''exact'')');
F2col=eval(lower(F2name),'strmatch(lower(F2name),lower(datastr),''exact'')');

[ncolumns,nchar] = size(datastr);
ndatastr = '';
for i=1:ncolumns          % set a header line
  ndatastr = [ ndatastr ' ' datastr(i,:) ];
end

% ===== print informations ============================

if ~silent
  disp(com);
  if (any(pscan(1:4) - pscan(5:8))) & ~silent
    disp('HKLE scan :');
    fprintf(1,'%g\t',pscan(1:4));
    fprintf(1,'\n\t\t...\n');
    fprintf(1,'%g\t',pscan(5:8));
    fprintf(1,'\n');
  end
end

t=0;
for i=1:ncolumns
  if (findstr(datastr(i,:),'TT'))
    t=mean(data(:,i)); % there is a TT column
    if ~silent
        fprintf(1,'T = %.2g K (+/- %.2g K)\n',t,std(data(:,i)));
    end
    break;
  end
end
if (~t)  
  ttpos = findstr(header,'TT'); % looks for TT in header
  if ~isempty(ttpos)
    tt = findstr(header(ttpos(1):length(header)),'=');
    t=sscanf(header((ttpos(1)+tt):length(header)),'%f');
    if ~silent
        fprintf(1,'T = %.2g K\n',t);
    end
  end
end
kfix = 0;
if (pres(10) == 2)
  if ~silent
      fprintf(1,'KF');
  end
  ki = 0; kf = pres(9); kfix = 2;
else 
  if ~silent
      fprintf(1,'KI');
  end
  kf = 0; ki = pres(9); kfix = 1;
end
if ~silent
  fprintf(1,' = %g A-1 (fixed) \n',pres(9));
end
for i=1:ncolumns
  if ki == 0 & (findstr(datastr(i,:),'KI'))
    ki=mean(data(:,i)); % there is a KI column
    if ~silent
        fprintf(1,'KI = %.2g +/- %.2g A-1\n',ki,std(data(:,i)));
    end
    break;
  end
  if kf == 0 & (findstr(datastr(i,:),'KF'))
    kf=mean(data(:,i)); % there is a KF column
    if ~silent
        fprintf(1,'KF = %.2g +/- %.2g A-1\n',kf,std(data(:,i)));
    end
    break;
  end
end
if (ki == 0)  
  kipos = findstr(header,'KI'); % looks for KI in header
  if ~isempty(kipos)
    ki = findstr(header(round(kipos(1)):length(header)),'=');
    ki=sscanf(header(round(kipos(1)+ki):length(header)),'%f');
    if ~silent
        fprintf(1,'KI = %.2g A-1\n',t);
    end
  end
end
if (kf == 0)  
  kfpos = findstr(header,'KF'); % looks for KI in header
  if ~isempty(kfpos)
    kf = findstr(header(round(kfpos(1)):length(header)),'=');
    kf=sscanf(header(round(kfpos(1)+ki):length(header)),'%f');
    if ~silent
        fprintf(1,'KF = %.2g A-1\n',t);
    end
  end
end
[n,c] = size(data);

if ~silent
  fprintf(1,'Data is ( %ix%i ).\n',n,c);
  fprintf(1,'DATA_:\n%s\n',ndatastr);
  nc = min(10,ncolumns);
  for i=1:min(2,n)
    fprintf(1,'%g ',data(i,1:nc));
    fprintf(1,'\n');
  end
  fprintf(1,'    ...\n');
  for i=max(1,n-1):n
    fprintf(1,'%g ',data(i,1:nc));
    fprintf(1,'\n');
  end
end
if ~isempty(stepvec)
  steps = find(stepvec ~= 0);
else
  steps = 1;
end
[n,c] = size(stepstr);
if ~silent & ~isempty(stepvec)
  if any(stepvec), fprintf(1,'STEP vector is :\n'); end
  for i = 1:n
    if stepvec(i), fprintf(1,'\t%i - %s = %g\n',i,stepstr(i,1:c),stepvec(i)); end
  end
end
monflag = -1; mn = 0; ti = 0;
i = findstr(' mn=',lower(header));
if ~isempty(i)
  i=i(1)+4; %length of item to search
  mon = sscanf(header(i:length(header)),'%f');
else
  i=findstr(' mn ',lower(com));  
  if ~isempty(i)
    i=i(1)+4; %length of item to search
    mon = sscanf(com(i:length(com)),'%f');
  end
end
if ~isempty(i)
  monflag = 1; mn = mon;
  if ~silent
      fprintf(1,'Constant monitor MN= %.2g\n',mon);
  end
end
i = findstr(' ti=',lower(header));
if ~isempty(i)
  i=i(1)+4; %length of item to search
  mon = sscanf(header(i:length(header)),'%f');
else
  i=findstr(' ti ',lower(com));  
  if ~isempty(i)
    i=i(1)+4; %length of item to search
    mon = sscanf(com(i:length(com)),'%f');
  end
end
if ~isempty(i)
  monflag = 0; ti = mon;
  if ~silent
      fprintf(1,'Constant time TI= %.2g\n',mon);
  end
end

if flip & ~silent
  fprintf(1,'Polarisation analysis was used for this scan\n');
end

%===== extract x,y,mon columns guesses ============================

if (strcmp(xname,'-'))      % auto X name set 
  if (length(steps) == 1)    % ok, only one variable in scan...
    xname = deblank(stepstr(steps,1:c));
    if lower(xname(1)) == 'd', xname = xname(2:end); end
    i = findstr(xname,'=');
    if ~isempty(i) & i>1
      xname = xname(1:(i-1));
    end
    if ~silent
        fprintf(1,'Scan variable is : %s\n',xname);
    end
  elseif ~isempty(stepvec)
    [j,i] = max(stepvec(steps));
    xname = stepstr(i,2:c);
  end
end

n = ncolumns;    
xcol=eval(xname,'strmatch(xname,datastr,''exact'')');
if isempty(xcol) xcol = ncolumns+1; end

ycol=eval(yname,'strmatch(yname,datastr,''exact'')');

mcol=eval(mname,'strmatch(mname,datastr,''exact'')');
if isempty(mcol) mcol = ncolumns+1; end

if (strcmp(ename,'-'))
  ename='sqrt(y)';
end

ecol=eval(ename,'strmatch(ename,datastr,''exact'')');
if isempty(ecol)
  if strcmp(ename,'sqrt(y)'), ecol = ncolumns+1; end
end

% polarisation

if ~isempty(PALcol)
  if isempty(F1col) & isempty(F2col)
    Fstate = 1:pal;
  else
  for i = 1:pal
    pstate = find(PALcol == i);
    if ~isempty(F1col), fl1 = data(pstate(1),F1col); else fl1 = 0; end
    if ~isempty(F2col), fl2 = data(pstate(1),F2col); else fl2 = 0; end
    Fstate = [Fstate; fl1 fl2 ];  % list of polarisation states used during scan
  end;
  end
else Fstate = [];
end

if length([ F1val F2val ]) == 2  % a set of flippers is asked in filespec
  i = find((Fstate(:,1) == F1val) & (Fstate(:,2) == F2val));
  if ~isempty(i), PALval = i(1); end
end

%===== GUI option call  ============================

if gui
  [p,name,ext,ver] = fileparts(filename);
  my_title = [ name,ext,ver ' (ILL TAS)' ];
  [xcol, ycol, ecol, mcol, xlab, ylab,scantmp,normf,keep,palo]=mf_coldg(datastr, ncolumns,xcol,ycol,ecol,mcol,scan,normf,Fstate, my_title, file_index);
  if isempty(xcol) return; end   % cancel button
  if isempty(scantmp) scantmp='1'; end
  if (~strcmp(scantmp,scan)) % need to reach a new part of file 
    scan = scantmp;
    if (xcol<=ncolumns) 
      xcol = datastr(xcol,:); 
    else
      xcol='#';
    end
    if isempty(PALval) & ~isempty(palo), PALval = palo; end
    if (ycol<=ncolumns) ycol = datastr(ycol,:); 
    else ycol = num2str(ycol); end
    if (mcol<=ncolumns) mcol = datastr(mcol,:); 
    else mcol = 'none'; end
    if (ecol<=ncolumns) ecol = datastr(ecol,:); 
    else
      if (ecol == ncolumns+1) ecol = 'sqrt(y)';
      else ecol = 'none'; end
    end
    if (PALcol<=ncolumns) PALcol = datastr(PALcol,:);
    else PALcol = 'none'; end
    filespec = sprintf('%s,X=%s,Y=%s,M=%s,E=%s,S=%s,N=%s',filename,xcol,ycol,mcol,ecol,scan,normf);
    if ~isempty(PALcol) & (PALcol<=ncolumns)
      filespec = sprintf('%s,P=%s',filespec,datastr(PALcol,:));
    end
    if ~isempty(PALval)
      filespec = sprintf('%s,F=%i',filespec,PALval);
    end
    if setpars
      filespec = [ filespec ',setpar' ];
    end
    if keep
      filespec = [ filespec ',gui' ];
    end
    if silent | ~keep
      filespec = [ filespec ',silent' ];
    end
    [x, y, err, xlab, ylab,monitor,optpars]=illbatch(filespec);
    return
  end
end

%===== last check for columns =============

if isempty(ecol) ecol=ncolumns+2; end
if isempty(mcol) mcol=ncolumns+1; end

%===== Check that there's not multiple columns of x,y,m =====
if length(xcol)>1 xcol=xcol(1); end
if length(ycol)>1 ycol=ycol(1); end
if length(ecol)>1 ecol=ecol(1); end
if length(mcol)>1 mcol=mcol(1); end
   
if (isempty(xcol)) % ----- first column that varies
  xcol=find(std(data));
  xcol=xcol(1);
end

if (isempty(ycol)) %----- Find first column that varies for y next to x
  ycol=find(std(data(:,(xcol+1):ncolumns)));
  ycol=ycol(1);
end


%===== extract x,y values ============================

%---- Find correct polarisation state if applicable
if flip 
  if (length(PALcol) == size(data,1)) & ~isempty(PALval)
    data = data(find(round(PALcol) == round(PALval)),:);
  elseif ~isempty(PALcol) & ~isempty(PALval)
    data = data( find(round(data(:,PALcol))==round(PALval)) ,:);
  elseif ~silent & size(Fstate,1) > 1
    disp('Polarization available, but not selected correctly !')
  end
end;

if isempty(data)
  disp('No data could be extracted !')
  disp([ 'Extraction command is ' filespec ]);
  return
end
%---- Work out which columns to extract, and extract --------

% first group columns in case x/ycol is a group of columns

if (xcol <= ncolumns)
  x=data(:,xcol);
  [n,c]=size(x);
  x=reshape(x',n*c,1);
  ncx=n*c;
else
  ncx = -1; % causes x = 1:length(y)
end

y=data(:,ycol);
[n,c]=size(y);
y=reshape(y',n*c,1);

if (ncx ~= n*c)
  x=1:(n*c);
  x=x';
end


if mcol <= ncolumns
  monitor=data(:,mcol);
else
  monitor=ones(length(x),1);
end
[n,c]=size(monitor);
monitor=reshape(monitor',n*c,1);

[monzeros]=find(monitor==0);  %----- Test to see if selected monitor is zero
if ~isempty(monzeros);
  disp(' ')
  disp(['  Warning: ' filename ' selected monitor has some zeros'])
  disp(' ')
  monitor=ones(length(x),1);
end

if (ecol == ncolumns+1)
  err=sqrt(abs(y));      % ...root(y) errors requested        
elseif (ecol >= ncolumns+2)
  err=1e-3*max(abs(y));    % no error bars - so equal weights
else
  err=data(:,ecol);      % specified errors
end
[n,c]=size(err);
err=reshape(err',n*c,1);
err = err./monitor;

%----- Normalise counts

normfv=eval(num2str(normf));

y=normfv*y ./monitor; 
err=normfv*err;



%----- Create labels

   if (xcol <= ncolumns)
     xlab=datastr(xcol(1),:); 
   else
  xlab='Point number';
   end
   if ~isempty(t) & (t ~= 0)
  xlab=sprintf('%s T=%.2f K',xlab,t);
   end
   if (monflag == 0)
  xlab=sprintf('%s ti=%.1f',xlab,mon);
  ymon=[' (per ' num2str(mean(monitor)) ' TI)' ];
   elseif (monflag == 1)
  xlab=sprintf('%s mn=%.1f',xlab,mon);
  ymon=[' (per ' num2str(mean(monitor)) ' MN)' ];
   else
  ymon = '';
   end

   if (mean(normfv) ~= 1) & (strcmp(normf,'1') ~= 1)
  ymon=sprintf('%s*%.1f', ymon, mean(normfv));
   end
   if kfix == 1
  xlab=sprintf('%s ki=%.3f',xlab,pres(9));
   elseif kfix == 2
  xlab=sprintf('%s kf=%.3f',xlab,pres(9));
   end
   
   ylab2=ymon;
   ylab1=fliplr(deblank(fliplr(deblank(datastr(ycol(1),:)))));
   ylab=[ylab1 ylab2];
   if ~isempty(PALval) & ~isempty(PALcol)
  if (size(Fstate,1) >= PALval) & (size(Fstate,2) == 2)
    ylab3 = sprintf(' F=%d:%d',Fstate(PALval,1),Fstate(PALval,2));
  else
    ylab3 = sprintf(' P=%d',PALval);
  end
  ylab = [ ylab ylab3 ];
   end

% check for xdata limits and POSQE agreement

  xmin = min(x);

  if isempty(find(pscan(1:4) == xmin)) & xcol(1) <= ncolumns
    % looks for X variable in stepstr(HKLE)
  [n,c] = size(stepstr);
  i = strmatch(deblank(datastr(xcol(1),:)),str2mat('QH','QK','QL','EN'));
  if ~isempty(i) & ~silent
    pscan(i) = xmin;
    fprintf(1,'Adjust %s limit in PSCAN\n',datastr(xcol(1),:))
  end
  end

  optpars = [ ki kf ti mn PALval ];


%===== Set params ============================

if setpars

% pres : pfile format '.par' with 42 values. --------------

if length(pres) ~=42 & ~silent
  disp([' Warning: ' filename ' Incorrect number of rescal parameters'])
else
  hrc_paras = findobj('Tag','hrc_paras');
  if (~isempty(hrc_paras))
    hrc_paras=get(hrc_paras,'Userdata');
    toupdate=[ 9 10 19 20 21 31:42 ];
    for i=toupdate
      if ~isnan(pres(i))
        set(hrc_paras(i),'String',num2str(pres(i)));
      end
    end
    if (~isempty(hrc_paras)) & ~silent
      fprintf(1,'Transfert main rescal parameters .\n');
    end
  end
end

% pscan : tfile format '.scn' with 8 values. --------------

if length(pscan) ~=8 & ~silent
  disp([' Warning: '  filename ' Incorrect number of scan parameters'])
end 
htrix_scan_paras = findobj('Tag','htrix_scan');
if (~isempty(htrix_scan_paras))
  htrix_scan_paras = get(htrix_scan_paras,'Userdata');
  for i=1:length(htrix_scan_paras)
    set(htrix_scan_paras(i),'String',num2str(pscan(i)));
  end
  if (~isempty(htrix_scan_paras)) & ~silent
    fprintf(1,'Transfert scan parameters.\n');
  end
end

% transfert to mfit params
hmf_pars = findobj('Tag','mf_ParWindow');
h=[];
if ~isempty(hmf_pars)
  h=get(hmf_pars,'userdata');
end
[npars d]=size(h);
funfile = findobj('Tag','mf_FitFuncFile');
if ~isempty(funfile) 
  funfile = get(funfile,'string');
else 
  funfile='';
end


if ~silent,
    fprintf(1,'Transfert mfit parameters for %s fit function.\n',funfile);
end
bragg = floor(pscan(1:3)+0.5);
dscan = (pscan(5:8) + pscan(1:4))/2;
dscan(1:3) = dscan(1:3) - bragg;
dsn = norm(dscan(1:3));
xlab = sprintf('%s q=%.2f',xlab,dsn);
fprintf(1,'|q| = %.3f  ',dsn);
if (dsn == 0) 
  if ~silent, fprintf(1,'(Zone center scan) '); end
  V0 = [ 0 1 0];
  dsn = 1;
else
  fprintf(1,'-- ');
  V0 = dscan/dsn;
end
typ = abs(sum(bragg.*dscan(1:3))/norm(bragg)/dsn);
if ~silent
  if (typ <= 0.01)
    fprintf(1,'Pure transverse or center scan\n');
  elseif (typ >= 0.99)
    fprintf(1,'Pure longitudinal scan\n');
  else
    fprintf(1,'Mixed T/L scan (%.1f %% L,angle %.1f deg)\n',typ*100,asin(typ)*180/pi);
  end
end
p=1:npars;
for i=1:npars
  c=str2num(get(h(i,1),'string'));
  if (isempty(c))
    p(i) = 0;
  else
    p(i) = c;
  end
  pnames = lower(get(h(i,3),'string'));
  if     findstr(pnames,'backg')
    p(i) = min(y);
  elseif findstr(pnames,'temp') & ~isempty(t)
    p(i) = t;
  elseif findstr(pnames,'energy')
    p(i) = dscan(4);
  elseif findstr(pnames,'hzero')
    p(i) = bragg(1);
  elseif findstr(pnames,'kzero')
    p(i) = bragg(2);
  elseif findstr(pnames,'lzero')
    p(i) = bragg(3);
  elseif findstr(pnames,'ph:xo')
    p(i) = dscan(1);
  elseif findstr(pnames,'ph:yo')
    p(i) = dscan(2);
  elseif findstr(pnames,'ph:zo')
    p(i) = dscan(3);
  elseif findstr(pnames,'dc:xdir')
    p(i) = V0(1);
  elseif findstr(pnames,'dc:ydir')
    p(i) = V0(2);
  elseif findstr(pnames,'dc:zdir')
    p(i) = V0(3);
  elseif findstr(pnames,'mon flag') & (monflag >= 0)
    p(i) = monflag;
  end
  
end
% make modifs

for i=1:npars
  set(h(i,1),'string',num2str(p(i),6));
end

end % if setpars

