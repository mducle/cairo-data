function [x,y,mn,setindex,emn,ey] = flexmerge(runnums,datafolder,norm_mn,temperature_threshold)
% Function to look through a set of datafiles, determine which have the same experimental parameters and sums them into sets
% Will sort by kfix, temperature and magnetic field. Will only parse HKLE scans.
% 
% Syntax:  [x,y,mn,setindex] = flexmerge(runnums,datafolder,norm_mn);
%
% Inputs:  runnums    - An array of run numbers
%          datafolder - A string with the path to the data
%          norm_mn    - _OPTIONAL_ if given this is the monitor value to normalise the data sets to. 
%
% Outputs: x          - A cell array of the independent variable
%          y          - A cell array of the dependent variable (the total counts, or normalised counts)
%          mn         - A cell array of the total monitor counts corresponding to each x,y, or the error if norm_mn is given
%          setindex   - A cell array with:
%                       {independent_var,   % cell   - The dependent variable name (e.g. QH, QK, QL or EN)
%                           scan_numbers,   % cell   - The numbers of the scans with the same parameters
%                                   Qfix,   % cell   - The fixed wavevector/energy transfer. (e.g. '(1 0 L), EN=0meV')
%                               kfix_var,   % cell   - Which of the incident or scattered wavevector is fixed
%                               kfix_val,   % vector - The value of the fixed wavevector in 1/A
%                            temperature,   % vector - The mean temperature during the scan
%                         magnetic_field}   % vector - The mean magnetic field during the scan
%
% If you only specify a single, or no output variable [e.g. x=flexmerge() or just flexmerge()], then the function will output 
% a list of the Q's, Counts, temperature and field to the screen.
%
% Note that if you have not outputed MAG, then the function assumes you measured at zero field!
% If you do not give a value for norm_mn, then the function will sum the counts and monitor
% Also note that this function requires the ILL Matlab library. You can download this from: 
%    http://www.ill.eu/de/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/
%
% Some examples:
%
% >> flexmerge(start:stop,'DataFolder');                     % Will output a list of summed scans with Q/E, Counts, Monitor.
%
% >> [x,y,e,sets] = flexmerge(start:stop,'Data',mon); 
% >> for iset=1:length(x); figure;                           % Will plot a graph of all summed datasets.
%    errorbar(x{iset},y{iset},m{iset},'.');          
%    xlabel(sets{1}{iset}); title(sets{3}{iset}); end

% Duc Le - Fri Aug 14 00:00:54 CEST 2009 - duc.le@helmholtz-berlin.de

% TODO: Save a set of dummy datafiles for use in plotting dispersion...
% TODO: Option to save output text to a word file...

if exist(datafolder)~=7; error('datafolder does not exist'); end

% ---------------------------------------------------------------------------------------------------------------------------------- %
% Some initial variables - Change these to suit your needs!
% ---------------------------------------------------------------------------------------------------------------------------------- %
monitor = 'M1';                   % The name of the monitor to use (either M1 or M2, etc.)
if ~exist('temperature_threshold','var')
temperature_threshold = 20;       % Percent of the temperature to categorise scans - e.g. if =5, then scans at 2K and 1.96K will be 
end                               %    considered to be at the same temperature.
step_threshold = 10;              % The percent of the step size below which data points will be considered to be repeated. Eg. if =11
                                  %    then points at EN=0.1 and EN=0.11 will be considered to be the same.
max_time = 600;                   % Longest counting time per point to consider - points with counting time longer than this (e.g.
                                  %    when reactor goes off, shutter goes down accidentally, are ignored.

% ---------------------------------------------------------------------------------------------------------------------------------- %
% Don't change below this...
% ---------------------------------------------------------------------------------------------------------------------------------- %
nothkle = zeros(size(runnums));   % Flags to show which scans are not HKLE scans
% ---------------------------------------------------------------------------------------------------------------------------------- %
% Determines the operating system
if ispc
   dirsep = '\';  % Windows
else
   dirsep = '/';  % Mac or Unix
end
% ---------------------------------------------------------------------------------------------------------------------------------- %

% ---------------------------------------------------------------------------------------------------------------------------------- %
% Loops through all the runs, determining their parameters, for subsequent processing
% ---------------------------------------------------------------------------------------------------------------------------------- %
for irun=1:length(runnums)
   [data,datastr,pres,pscan,stepvec,stepstr,com,header,head]=illdata([datafolder dirsep sprintf('eiger2013n%06.f.scn',runnums(irun))]);
   pscan = round(pscan*100)/100;
   % Checks that this is an HKLE scan
   if isempty(findstr(upper(com),'QH'))
      nothkle(irun) = 1; 
      continue;   % Skips the rest of the loop...
   end
   data(:,5)=round(data(:,5)*100)/100;
   % Determines whether ki or kf is fixed
   ifx = strfind(header,'FX='); fx=str2num(header((ifx+3):(ifx+5)));
   if(fx==1) 
      kfix_var{irun} = 'ki Fixed'; 
   elseif(fx==2) 
      kfix_var{irun} = 'kf Fixed'; 
   else 
      kfix_var{irun} = 'Unknown kfix!';
   end
   % Determines the value of kfix
   tstr = header((ifx+5):(ifx+100)); kfix_val(irun) = str2num(tstr((strfind(tstr,'KFIX=')+5):(strfind(tstr,'PARAM')-1)));
   % Determines which variable out of QH,QK,QL,EN is varied
   istepstr = findstr(header,'STEPS:'); 
   icolonstr = findstr(header,':'); ii=1;
   for iii = 1:length(istepstr)
      inextline = icolonstr(find(icolonstr>istepstr(iii)));
      stepline = header(istepstr(iii):inextline(2));
      iequals = findstr(stepline,'=');
      for iiii = 1:length(iequals)
         mystepvec(ii) = str2num(stepline((iequals(iiii)+1):(iequals(iiii)+8))); 
         mystepstr{ii} = stepline((iequals(iiii)-2):(iequals(iiii)-1));
         ii=ii+1; 
      end
   end
   stepvec = mystepvec;
   Qvar = find(mystepvec); %mystepstr = cellstr(stepstr);
   if length(Qvar)>1
      mystepvec = stepvec./stepvec(Qvar(1));
      for iQvar=2:length(Qvar)
         if mystepvec(iQvar)~=1
            mystepstr{iQvar} = [num2str(mystepvec(iQvar)) mystepstr{Qvar(1)}]; 
         else
            mystepstr{iQvar} = mystepstr{Qvar(1)};
         end
      end
   elseif isempty(Qvar)  % Zero step size
      Qvar = 1; mystepvec = stepvec;
   end
   independent_var{irun} = mystepstr{Qvar(length(Qvar))}; 
   % Determines the fixed wavevector transfer, and fixed energy transfer (if any)
   Qinit = pscan(1:4);
   xvar = '(';
   for iQ=1:3
      if isempty(find(Qvar==iQ))
         xvar = [xvar num2str(Qinit(iQ)) ' ']; 
      else 
         xvar = [xvar mystepstr{iQ} ' ']; 
      end
   end
   xvar = [xvar ')']; 
   if isempty(find(Qvar==4))
      xvar = [xvar ', EN=' num2str(Qinit(4))];
   end
   Qfix{irun} = xvar;
   % Determines the mean temperature
   mydatastr = cellstr(datastr);
   [tmpvar,ttColInd] = intersect(mydatastr,'TT'); 
   if ~isempty(ttColInd)
      tt = data(:,ttColInd); tt(tt>9e3)=[]; temperature(irun) = mean(tt);
      % Checks to see if the TT sample thermometer is offline - TUNDER - as
      %    it will report an unphysically low (and constant) temperature
      if sum(diff(tt))==0 && (temperature(irun)>1 && temperature(irun)<1.3)
         ttColInd = findstr(header,', RT='); tstr = header((ttColInd+10):(ttColInd+50)); ifx = findstr(tstr,',');
         if isempty(ifx); temperature(irun)=300; else; temperature(irun) = str2num(tstr(1:ifx(1))); end
      end
   else  % TT was not outputted, so we try the header...
      ttColInd = findstr(header,'PARAM: TT='); tstr = header((ttColInd+10):(ttColInd+50)); ifx = findstr(tstr,',');
      if isempty(ifx)       % DTI was not turn on - assume we are at room temperature
         temperature(irun) = 300;
      else
         temperature(irun) = str2num(tstr(1:ifx(1)));
      end
   end
   % Determines the mean magnetic field
   [tmpvar,mgColInd] = intersect(mydatastr,'MAG');
   if ~isempty(mgColInd)
      magnetic_field(irun) = mean(data(:,mgColInd)); 
   else
      magnetic_field(irun) = 0; 
   end
   % Stores the independent variable values, counts and monitor
   [tmpvar,tColInd] = intersect(mydatastr,'TIME'); tColInd=data(:,tColInd)<max_time;
   [tmpvar,xColInd] = intersect(mydatastr,mystepstr{Qvar(length(Qvar))});
   xval{irun} = data(tColInd,xColInd);
   [tmpvar,yColInd] = intersect(mydatastr,'CNTS');
   yval{irun} = data(tColInd,yColInd);
   [tmpvar,mColInd] = intersect(mydatastr,monitor);
   mval{irun} = data(tColInd,mColInd);
   scan_numbers(irun) = runnums(irun);
end

% ---------------------------------------------------------------------------------------------------------------------------------- %
% Begins parsing the list of scans
% ---------------------------------------------------------------------------------------------------------------------------------- %
hkleInd = ~nothkle;                      % We only want the HKLE scans, not scans of any particular angle...
kfixvar = {'ki Fixed', 'kf Fixed'};      % To allow looping

% Determines the number of unique kfix values, temperature and magnetic field
ukfix_val = unique(round(kfix_val*100)/100); ukfix_val(find(ukfix_val==0))=[];
for iii=1:length(kfix_val); ukfix_val = round(ukfix_val*100)/100; end;    % Round to nearest centimal(?) place
umagnetic_field = unique(round(magnetic_field*100)/100);
utemperature = unique(temperature); utemperature(find(utemperature==0))=[];
% For the temperature, if the value is within the temperature threshold, we count it as a single temperature
for iii=1:2; itemp = 2;   % Do two passes!
while(itemp<=length(utemperature))
   if abs(utemperature(itemp)-utemperature(itemp-1)) < (utemperature(itemp-1)*(temperature_threshold/100))
      utemperature(itemp)=[];
   else
      itemp = itemp + 1;
   end
end; end; utemperature = unique(utemperature); 
utemperature = unique(round(utemperature*1000)/1000);
% Determine all the combinations of kfix, temperature and field
numphyschar = 1;
for ikfx = 1:length(ukfix_val); for itemp = 1:length(utemperature); for imag = 1:length(umagnetic_field);
   physchar(numphyschar,:) = [ukfix_val(ikfx) utemperature(itemp) umagnetic_field(imag)]; numphyschar=numphyschar+1;
end; end; end
% Finds those scans with the same independent variable
iscans = 1; mystepstr = {'QH' 'QK' 'QL' 'EN'};
for itype=1:4
   scanInd = find(hkleInd).*ismember(independent_var(hkleInd),mystepstr{itype}); % 1=QH, 2=QK, 3=QL, 4=EN
   scanInd(find(scanInd==0))=[]; %scantype{itype} = scan_numbers(scanInd);
   % Of these split into those with ki fixed or kf fixed
   if ~isempty(scanInd)
      for ikfx=1:2
         scnkfxInd = scanInd.*ismember(kfix_var(scanInd),kfixvar{ikfx}); scnkfxInd(find(scnkfxInd==0))=[];
         % Split into those with particular values of kfix, temperature and fields
         if ~isempty(scnkfxInd)
            for iii = 1:(numphyschar-1)
               scnkxIndFound = find( (abs(kfix_val(scnkfxInd)-physchar(iii,1))<1e-2) ... 
                           .*(abs(temperature(scnkfxInd)-physchar(iii,2))<(physchar(iii,2)*(temperature_threshold/100))) ...
                           .*(abs(magnetic_field(scnkfxInd)-physchar(iii,3))<1e-3) );
               sameInd = scnkfxInd(scnkxIndFound); sameInd(find(sameInd==0))=[];
               if ~isempty(sameInd)
                  scnkfxInd(scnkxIndFound) = [];  % So that we don't count scans twice or more!
                  usame = unique(Qfix(sameInd));  % Finds the unique fixed Q-vectors (i.e. repeated scans)
                  for iSame = 1:length(usame)
                     indices = sameInd.*ismember(Qfix(sameInd),usame{iSame}); indices(find(indices==0))=[];
                     s_independent_var{iscans} = independent_var{indices(1)};
                     s_scan_numbers{iscans} = scan_numbers(indices);
                     s_Qfix{iscans} = Qfix{indices(1)};
                     s_kfix_var{iscans} = kfix_var{indices(1)};
                     s_kfix_val(iscans) = physchar(iii,1);
                     s_temperature(iscans) = physchar(iii,2);
                     s_magnetic_field(iscans) = physchar(iii,3);
                     s_indices{iscans} = indices;
                     iscans = iscans + 1;
                  end
               end
            end
         end
      end
   end
end
% Finished sorting, sets the output setindex variables
setindex = {s_independent_var, s_scan_numbers, s_Qfix, s_kfix_var, s_kfix_val, s_temperature, s_magnetic_field};

% List of scans for debug purposes
for i=1:length(setindex{1})
   disp(sprintf('%s scan at %s,  \t%s at %5.3fA^-1,  \tT=%5.3fK  \tB=%5.3fT',setindex{1}{i},setindex{3}{i},setindex{4}{i},setindex{5}(i),setindex{6}(i),setindex{7}(i)));
end

% ---------------------------------------------------------------------------------------------------------------------------------- %
% Now loop throught each set of repeated scan, summing counts and monitor values
% ---------------------------------------------------------------------------------------------------------------------------------- %
for iset = 1:length(s_scan_numbers)
   x=[]; y=[]; mn=[];
   for irun = 1:length(s_scan_numbers{iset})
      step = mean(diff(xval{s_indices{iset}(irun)}))*step_threshold/100; % The step length in this scan
      if ((length(xval{s_indices{iset}(irun)}))==1)
         step=xval{s_indices{iset}(irun)}(1)*step_threshold/100; 
      end
      for ix = 1:length(xval{s_indices{iset}(irun)})                     % Loop through each value of x (e.g. EN)
         xInd = find(abs(x-xval{s_indices{iset}(irun)}(ix))<abs(step));
         if isempty(xInd)                                                % Does not match a value yet -> append to end
            x = [x xval{s_indices{iset}(irun)}(ix)]; 
            y = [y yval{s_indices{iset}(irun)}(ix)]; 
            mn = [mn mval{s_indices{iset}(irun)}(ix)];
         else                                                            % Does match, sum y (CNTS) and MN
            y(xInd) = y(xInd) + yval{s_indices{iset}(irun)}(ix);
            mn(xInd) = mn(xInd) + mval{s_indices{iset}(irun)}(ix);
         end
      end
   end
   [s_x{iset},ix] = sort(x); s_y{iset} = y(ix); s_mn{iset} = mn(ix);     % Sorts the data into ascending order in x
end

% ---------------------------------------------------------------------------------------------------------------------------------- %
% Determine if the user wants error values or monitor values (i.e. if they have given a monitor value to normalise to).
% ---------------------------------------------------------------------------------------------------------------------------------- %
if nargin>2 
   emn = s_mn; ey = s_y;
   for iset = 1:length(s_scan_numbers)
      s_e{iset} = (norm_mn./s_mn{iset}) .* sqrt(s_y{iset});
      s_y{iset} = (norm_mn./s_mn{iset}) .* s_y{iset};
   end
   x = s_x; y = s_y; mn = s_e;
else
   x = s_x; y = s_y; mn = s_mn;
end

% ---------------------------------------------------------------------------------------------------------------------------------- %
% Now determine if the user wants the data or just a list of points
% ---------------------------------------------------------------------------------------------------------------------------------- %
%setindex = {s_independent_var, s_scan_numbers, s_Qfix, s_kfix_var, s_kfix_val, s_temperature, s_magnetic_field};
if nargout<1
   for iset = 1:length(s_scan_numbers)
      if(nargin>2)
         countstr = ['Counts/MN=' num2str(norm_mn)]; mnstr = 'Error';
      else
         countstr = 'Total Counts'; mnstr = 'Total Monitor';
      end
      num_scn = length(s_scan_numbers{iset}); 
      if num_scn==1
         scanstr = num2str(s_scan_numbers{iset}); 
      else
         scanstr = sprintf('%i, ',s_scan_numbers{iset}(1:(length(s_scan_numbers{iset}))-1));
         scanstr = [scanstr num2str(s_scan_numbers{iset}(length(s_scan_numbers{iset})))];
      end
      disp('|------------------------------------------------------------|');
      disp(['| ' s_kfix_var{iset} ' = ' num2str(s_kfix_val(iset)) ',  ' s_Qfix{iset} ',  ' ...
            'T = ' num2str(s_temperature(iset)) 'K, B = ' num2str(s_magnetic_field(iset)) 'T.']);
      disp(['| Scans:  ' scanstr]);
      disp('|------------------------------------------------------------|');
      disp([s_independent_var{iset} sprintf('\t') countstr sprintf('\t') mnstr]);
      for ix = 1:length(x{iset})
         disp([num2str(x{iset}(ix)) sprintf('\t\t') num2str(y{iset}(ix)) sprintf('\t\t') num2str(mn{iset}(ix))]);
      end
   end
   x = [];
elseif nargout==1
   x = {x y mn setindex};
end
