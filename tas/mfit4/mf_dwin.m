function hmf_data=mf_dwin(xlab, ylab)
%
% MFIT function  hmf_data=mf_dwin
%		Create new data window
% 		MZ 29.11.94
%
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

%========== Create data window ===================

if ~exist('xlab')
 xlab='';
end
if ~exist('ylab')
 ylab='';
end

if (hmf_data~=0)
	figure(hmf_data);
%	cla;
else
	pos=eval([ '[' get(findobj('tag','mf_FigurePosition'),'string') ' ' get(findobj('tag','mf_FigureSize'),'string') ']' ] );
	if length(pos) < 4, pos = [ 200 20 400 350 ]; end
	pos=pos(:);
	bgcolor =get(findobj('tag','mf_FigureBgColor'),'string');
	if (isempty(bgcolor) | strcmp(bgcolor, 'gray')), bgcolor = [0.8 0.8 0.8]; end
	hmf_data = figure('Position',pos,...
						'tag','mf_DataWindow',...
						'color',bgcolor,...
						'MenuBar','none',...
						'Name','MFIT:Data',...
						'visible','off',...
						'HandleVisibility','on',...
						'WindowButtonDownFcn',...
						'mf_btprs',...
						'NextPlot','Add',...
				      'Interruptible','on');

	if ~isempty(bgcolor) 
	end
	     
   uicontrol('Tag','mf_zoom_list','visible','off');      % Zoom list
   uicontrol('Tag','mf_rbbox_action','visible','off');   % Action of rbbox

%--------- Zoom menu ------------------------------
	hzoom=uimenu(hmf_data,'Label','&Zoom');
	uimenu(hzoom,...
			'Label','Zoom in',...
			'Tag','hmf_mvzoom',...
			'Checked','on',...
			'Accelerator','z',...
			'Callback','mf_rbbox(''zoomin'')');
	uimenu(hzoom,...
			'Label','Zoom out',...
			'Callback','mf_rbbox(''zoomout'')');
	uimenu(hzoom,...
			'Label','Reset Zoom',...
			'Accelerator','r',...
			'Callback','mf_rbbox(''zoomreset'')');


%---------Select/Deselect data menu ----------------
	hsel=uimenu(hmf_data,'Label','&Data');
	uimenu(hsel,...
			'Label','Select Data',...
			'Tag','hmf_mvsel',...
			'Checked','off',...
			'accelerator','s',...
			'Callback','mf_rbbox(''select'')');
	uimenu(hsel,...
			'Label','Deselect Data',...
			'Tag','hmf_mvusel',...
			'Checked','off',...
			'accelerator','d',...
			'Callback','mf_rbbox(''deselect'')');
	uimenu(hsel,...
			'Label','Select All',...
			'Callback','mf_rbbox(''selectall'')');
	uimenu(hsel,...
			'Label','Deselect All',...
			'Callback','mf_rbbox(''deselectall'')');

		uimenu(hsel,...
			'Label','Change X Axis...',...
			'accelerator','a',...
			'Callback','mf_chgx(''mfit+show'');');
	if exist('frombase') & exist('mf_exprdg')
		uimenu(hsel,...
			'Label','Transform...',...
			'accelerator','f',...
			'Callback','[mf_x, mf_y, mf_err] = frombase(''mfit''); tomfit(mf_x,mf_y,mf_err); clear mf_x mf_y mf_err');
	end
	uimenu(hsel,...
			'Label','Free all params',...
			'separator','on',...
			'Callback','mf_batch(''free all'');');
	uimenu(hsel,...
			'Label','Fix all params',...
			'Callback','mf_batch(''fix all'');');

    if isnumeric(hmf_data); hmf_data_n = num2str(hmf_data); else hmf_data_n = num2str(hmf_data.Number); end;
	uimenu(hsel,...
			'Label','Close Data',...
			'accelerator','w',...
			'Callback',[ 'close(' hmf_data_n ');' ]);


%---------- Text menu-----------------------------
	htxt=uimenu(hmf_data,'Label','&Text');
	uimenu(htxt,...
	      		'Label','Add Text with Pars',...
			'accelerator','t',...
			'CallBack','mf_text(0,''mf_text_add'');');

	uimenu(htxt,...
	      		'Label','Add Text',...
			'CallBack','mf_opt(''addtext'');');

	uimenu(htxt,...
	      		'Label','Remove Text',...
			'CallBack','mf_text([],''mf_text_deleteall''); mf_gdata(''noload'');');

%---------- Axes menu -------------------------------
   haxes=uimenu(hmf_data,'Label','&Axes');

   uimenu(haxes,...
        			'Label','Error bars',...
       			'Tag','mf_ebarsonoff',...
      			'Checked','on',...
         			'Callback','mf_opt(''ebars'')');

	uimenu(haxes,...
			'Label','Grid',...
         			'Tag','mf_gridonoff',...
         			'Checked','off',...
			'Callback','mf_opt(''grid'')');
	uimenu(haxes,...
			'Label','Log y axis',...
			'Tag','mf_logliny',...
			'Checked','off',...
			'Callback','mf_opt(''logliny'')');

	uimenu(haxes,...
			'Label','Log x axis',...
			'Tag','mf_loglinx',...
			'Checked','off',...
			'Callback','mf_opt(''loglinx'')');

	uimenu(haxes,...
			'Label','Axes limits',...
			'Callback','mf_xylim');
	uimenu(haxes,...
			'Label','Hide Fit',...
			'Callback','set(findobj(''Tag'',''mf_fitline''),''Visible'',''off'');');

%---------- Print menu-----------------------------
	hprt=uimenu(hmf_data,'Label','&Print');
	uimenu(hprt,'Label','Print figure',...
			'accelerator','p',...
			 'Callback','mf_print');
	uimenu(hprt,'Label','Configure...',...
			 'Callback','printdlg');
	uimenu(hprt,'Label','Save figure',...
			 'Callback','[f,p] = uiputfile(''*.m'',''Save Figure As...'');if f,eval([''print -dmfile '' p f ]); disp(''Figure saved''); end');


end

%=========== Setup axes ===========================

axescolor =get(findobj('tag','mf_AxesColor'),'string');
if isempty(axescolor), axescolor = 'black'; end
axesfont =get(findobj('tag','mf_AxesFont'),'string');
if isempty(axesfont), axesfont = 'Helvetica'; end
axesfontsize =str2num(get(findobj('tag','mf_AxesFontSize'),'string'));
if isempty(axesfontsize), axesfontsize = '12'; end

hax=get(hmf_data,'CurrentAxes');
hgrid=findobj('Tag','mf_gridonoff');
if ~isempty(hgrid)
	hgrid = get(hgrid,'Checked');
else
	hgrid = 'off';
end

set(hax,...
	'box','on',...
	'drawmode','fast',...
	'position',[ 0.175 0.175 0.8 0.7],...
	'userdata','',...
	'xgrid',hgrid,...
	'ygrid',hgrid,...
	'xcolor',axescolor,...
	'ycolor',axescolor,...
	'Fontname',axesfont,...
	'Fontsize',axesfontsize);


%----- Axis labels and title ----------------------------------------

labelfont=get(findobj('tag','mf_LabelFont'),'string');
if isempty(labelfont), labelfont = 'Helvetica'; end
labelfontsize=str2num(get(findobj('tag','mf_LabelFontSize'),'string'));
if isempty(labelfontsize), labelfontsize = '12'; end

datadir=get(findobj('Tag','mf_DataDir'),'string');
datafile=get(findobj('Tag','mf_DataFile'),'string');

titl=[datafile]; % titl = [datadir datafile] can be ok for Unix or VMS

if isempty(xlab)
	xlab = get(findobj('tag','mf_text_xlabel'),'String');
	if isempty(xlab)
		xlab = get(get(findobj(hmf_data,'type','axes'),'Xlabel'),'string');		if isempty(xlab)
			xlab = 'x';
		end
	end
end
if isempty(ylab)
	ylab = get(findobj('tag','mf_text_ylabel'),'String');
	if isempty(ylab)
		ylab = get(get(findobj(hmf_data,'type','axes'),'Ylabel'),'string');
		if isempty(ylab)
			ylab = 'x';
		end
	end
end

xlabl=xlab;
ylabl=ylab;

set(get(hax,'Title'),...               % Make title
   'String',titl,...
   'Tag','mf_text_title',...
   'Fontname',labelfont,...
   'Fontsize',labelfontsize,...
   'handlevisibility','on',...
	'color',axescolor,...
   'userdata',get(hax,'Title'));
set(get(hax,'Xlabel'),...                % Make xlabel
   'String',xlabl,...
   'Tag','mf_text_xlabel',...
   'Fontname',labelfont,...
   'Fontsize',labelfontsize,...
   'handlevisibility','on',...
	'color',axescolor,...
   'userdata',get(hax,'Xlabel'));
set(get(hax,'Ylabel'),...                % Make ylabel
   'String',ylabl,...
   'Tag','mf_text_ylabel',...
   'Fontname',labelfont,...
   'Fontsize',labelfontsize,...
   'handlevisibility','on',...
	'color',axescolor,...
   'userdata',get(hax,'Ylabel'));

set(findobj('Tag','mf_zoom_list'),'userdata','');           % Clear zoom list
set(findobj('Tag','mf_rbbox_action'),'userdata','zoomin');  % Set default action

set(hmf_data,'visible','on');


