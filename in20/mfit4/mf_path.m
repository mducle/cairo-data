function h=mf_path(direc)
%
% MFIT function mf_path
%  Add directory to matlab path, if not there already
% M. Zinkin 5.12.94
%

if ~exist('direc')
	direc = '';
end

%---- Work out what path separator is on this machine ---------------
pathsep = ':';
c = computer;
if strcmp(c(1:2),'PC') |  strcmp(c(1:2),'MA')
	pathsep = ';';
elseif strcmp(c(1:2),'VA')
	if strcmp(c(1:7),'VAX_VMS')
		pathsep = ',';
	end
end

%---- Get old mfit path ---------------------------------------------
h=findobj('tag','mf_path');
if isempty(h)
	h=uicontrol('tag','mf_path','visible','off','string','');
end
mfpath=get(h,'string');

%----- Delete old mfit path -----------------------------------------
mpath=path;
p=findstr(mpath,mfpath);
if ~isempty(p)
	mpath(p:(p+length(mfpath)-1))=[];
end

%----- Add new mfit path ---------------------------------------------
mfpath = '';
if ~isempty(direc) & isempty(findstr(direc,mpath))
	mfpath = [ mfpath pathsep direc ];
end	
loadroutdir =get(findobj('Tag','mf_LoadRoutineDir'),'string');
direc = loadroutdir;
if ~isempty(direc) & isempty(findstr(direc,mpath))
	mfpath = [ mfpath pathsep direc ];
end
fitroutdir =get(findobj('Tag','mf_FitRoutineDir'),'string');
direc = fitroutdir;
if ~isempty(direc) & isempty(findstr(direc,mpath))
	mfpath = [ mfpath pathsep direc ];
end
fitfuncdir =get(findobj('Tag','mf_FitFuncDir'),'string');
direc = fitfuncdir;
if ~isempty(direc) & isempty(findstr(direc,mpath))
	mfpath = [ mfpath pathsep direc ];
end

if isempty(mfpath)
	i = [];
else
	i=find(mfpath==' ');
end
if ~isempty(i)
	mfpath(i)=[];
end
%mf_path=lower(mf_path);
set(h,'string',mfpath);
eval('addpath([mfpath ]);','');
h = mfpath;

