classdef Plot
% Plot contains customized routines to plot.
    
    properties
    end
    
    methods (Static)
        function lineBands(ks,es,nbins,marker,labels)
            kpts = size(ks,1);
            if nargin<4
                plot(1:kpts,es)
            else
                plot(1:kpts,es,marker)
            end
            xlim([1,kpts])

            xticks = 1+[0,cumsum(nbins)];
            xticklabels = cell(length(nbins)+1,1);
            hold on
            for ii = 0:length(nbins)
                xii = xticks(ii+1);
                plot(xii*[1,1],ylim,'--k')
                xticklabels{ii+1} = ['[',num2str(ks(xii,:),'%5.2f'),']'];
            end
            hold off

            ax = gca;
            ax.XTick = xticks;
            if nargin<5
                ax.XTickLabel = xticklabels;
                ax.XTickLabelRotation = -60;
            else
                ax.XTickLabel = labels;
            end
        end
        
        function spectrum(es,varargin)
        % varargin contains ks
        % Examples (suppose system is 3D):
        % e.g.1: if size(es) = [nes, kxpts, kypts, kzpts]
        %        then Plot.spectrum(es, kxs, kys, kzs)
        %        or   Plot.spectrum(es, {kxs, kys, kzs})
        % e.g.2: if size(es) = [nes, kxpts]
        %        then Plot.spectrum(es, kxs)
        %        or   Plot.spectrum(es, kxs, [], []) % to tell it's 3D
        %        or   Plot.spectrum(es, {kxs, [], []})
        % e.g.3: if size(es) = [nes, kzpts]
        %        then Plot.spectrum(es, [], [], kzs)
        %        or   Plot.spectrum(es, {[], [], kzs})
        % TODO: add parameter control (es contains iterations over pars)
            if nargin<2
                figure, plot(es,'o')
                return
            elseif length(varargin)==1
                if iscell(varargin{1})
                    ks = reshape(varargin{1},[],1);
                else
                    ks = varargin;
                end
            else
                ks = reshape(varargin,[],1);
            end
            dks = cellfun(@length,ks);
            dks = reshape(dks+(dks==0),1,[]);
            es = reshape(es,[size(es,1),dks]);
            es = permute(es,[1,fliplr(2:length(dks)+1)]); % a la plot_sp
            plot_sp(ks,es)
        end
        
        function wavefuncs(wfs,es,ns,varargin)
        % ns is size of the lattice
        % varargin contains ks
        % Examples (suppose system is 3D):
        % e.g.1: if size(wfs) = [ndof*ny*nz, nes, kxpts]
        %           and size(es) = [nes, kxpts]
        %        then Plot.wavefuncs(wfs,es,[0,ny,nz],kxs,[],[])
        %        or   Plot.wavefuncs(wfs,es,[0,ny,nz],{kxs,[],[]})
        % e.g.2: if size(wfs) = [ndof*ny, nes, kxpts, kzpts]
        %           and size(es) = [nes, kxpts, kzpts]
        %        then Plot.wavefuncs(wfs,es,[0,ny,0],kxs,[],kzs)
        %        or   Plot.wavefuncs(wfs,es,[0,ny,0],{kxs,[],kzs})
        % TODO: add parameter control; customized operators
            if nargin<4
                ks = cell(length(ns),1);
            elseif length(varargin)==1
                if iscell(varargin{1})
                    ks = reshape(varargin{1},[],1);
                else
                    ks = varargin;
                end
            else
                ks = reshape(varargin,[],1);
            end
            dks = cellfun(@length,ks);
            dks = reshape(dks+(dks==0),1,[]);
            nes = size(es,1);
            es = reshape(es,[nes,dks]);
            es = permute(es,[1,fliplr(2:length(dks)+1)]); % a la plot_wfs
            ns = reshape(ns+(ns==0),1,[]);
            ndof = numel(wfs)/(numel(es)*prod(ns));
            wfs = reshape(wfs,[ndof,ns,nes,dks]);
            nn = length([ndof,ns,nes]);
            wfs = permute(wfs,[1:nn,nn+fliplr(1:length(dks))]); % a la plot_wfs
            plot_wfs(ks,es,wfs)
        end
        
        function manipulate(varargin)
        % manipulate plotting function with varying parameters
        % eg 1: manipulate('Figure_Name',@(ip)plot(xs,ys(:,ip)),p)
        % eg 2: manipulate(@(pind)plot(xs,ys(:,pind(1),pind(2))),p1,p2)
        % eg 3: manipulate({@(ip)plot(xs,ys(:,ip)),{'YLim'},{[-0.2,0.2]}},p)

            if isa(varargin{1},'function_handle') || iscell(varargin{1})
                fname = 'Manipulate';
                fplot = varargin{1};
                poffset = 1;
            elseif isa(varargin{1},'char') && ...
                   (isa(varargin{2},'function_handle') || iscell(varargin{2}))
                fname = varargin{1};
                fplot = varargin{2};
                poffset = 2;
            else
                error('Invalid input!')
            end

            if ~iscell(fplot)
                f = fplot;
            else
                f = fplot{1};
                apn = fplot{2}; % axis property name (cell array)
                apv = fplot{3}; % axis property value (cell array)
            end

            p = varargin(poffset+1:end); % parameters
            npar = length(p);
            nval = cellfun(@length,p);
            pind = ones(npar,1);

            sht = 15; % slider and text height
            sph = 5; % space height
            aht = 520; % axis height
            fht = aht+sph+npar*(sht+sph);

            %  Create and then hide the GUI as it is being constructed.
            hf = figure('Visible','off','Position',[360,500,560,fht],...
                'NumberTitle', 'off');

            %  Construct the components.
            hctls(npar,3) = uicontrol('Visible','off'); % handles for sliders and texts
            for ipar = 1:npar
                cbm = fht-ipar*(sht+sph);
                hctls(ipar,1) = uicontrol('Style','text',...
                                'BackgroundColor',get(hf,'Color'),...
                                'String',inputname(poffset+ipar),...
                                'Position',[20,cbm,50,sht]);
                if nval(ipar)>1
                hctls(ipar,2) = uicontrol('Style','slider',...
                                'Max',nval(ipar),'Min',1,'Value',pind(ipar),...
                                'SliderStep',[1/nval(ipar)+1e-6 5/nval(ipar)+1e-6],...
                                'Position',[80,cbm,420,sht],...
                                'Callback',{@slider_Callback,ipar});
                else
                hctls(ipar,2) = uicontrol('Style','text',...
                                'BackgroundColor',get(hf,'Color'),...
                                'String','',...
                                'Position',[80,cbm,420,sht]);
                end
                hctls(ipar,3) = uicontrol('Style','text',...
                                'BackgroundColor',get(hf,'Color'),...
                                'String',num2str(p{ipar}(pind(ipar))),...
                                'Position',[505,cbm,50,sht]);
            end
            ha = axes('Units','Pixels','OuterPosition',[0,0,560,aht]);

            % Initialize the GUI.
            % Change units to normalized so components resize automatically.
            set([hf,ha,reshape(hctls,1,[])],...
                'Units','normalized');

            %Create a plot in the axes.
            update_plot

            % Assign the GUI a name to appear in the window title.
            set(hf,'Name',fname)
            % Move the GUI to the center of the screen.
            movegui(hf,'center')
            % Make the GUI visible.
            set(hf,'Visible','on');

            % Subfuntions implementing the callbacks
            % ------------------------------------------

            function slider_Callback(source,eventdata,ipar) 
                val = round(get(source,'Value'));
                set(source,'Value',val);
                pind(ipar) = val;
                set(hctls(ipar,3),'String',num2str(p{ipar}(val)));
                update_plot
            end

            % Subfuntions implementing other functions
            % ------------------------------------------

            function update_plot
                set(0,'CurrentFigure',hf)
                f(pind);
                axis tight
                if iscell(fplot)
                    set(gca,apn,apv)
                end
            end
        end
        
        function manipulate2(varargin)
        % manipulate plotting function with varying parameters
        % separated control and plot windows
        % eg 1: manipulate('Figure_Name',@(ip)plot(xs,ys(:,ip)),p)
        % eg 2: manipulate(@(pind)plot(xs,ys(:,pind(1),pind(2))),p1,p2)
        % eg 3: manipulate({@(ip)plot(xs,ys(:,ip)),{'YLim'},{[-0.2,0.2]}},p)

            if isa(varargin{1},'function_handle') || iscell(varargin{1})
                fname = 'Manipulate';
                fplot = varargin{1};
                poffset = 1;
            elseif isa(varargin{1},'char') && ...
                   (isa(varargin{2},'function_handle') || iscell(varargin{2}))
                fname = varargin{1};
                fplot = varargin{2};
                poffset = 2;
            else
                error('Invalid input!')
            end

            if ~iscell(fplot)
                f = fplot;
            else
                f = fplot{1};
                apn = fplot{2}; % axis property name (cell array)
                apv = fplot{3}; % axis property value (cell array)
            end

            p = varargin(poffset+1:end); % parameters
            npar = length(p);
            nval = cellfun(@length,p);
            pind = ones(npar,1);

            sht = 15; % slider and text height
            sph = 5; % space height
            aht = 520; % axis height
            fht = sph+npar*(sht+sph);

            %  Create and then hide the GUI as it is being constructed.
            hf1 = figure('Visible','off','Position',[360,500,560,fht],...
                'NumberTitle', 'off');
            hf2 = figure('Visible','off','Position',[360,500,560,aht],...
                'NumberTitle', 'off','CloseRequestFcn',@f2_closereq);

            %  Construct the components.
            hctls(npar,3) = uicontrol(hf1,'Visible','off'); % handles for sliders and texts
            for ipar = 1:npar
                cbm = fht-ipar*(sht+sph);
                hctls(ipar,1) = uicontrol(hf1,'Style','text',...
                                'BackgroundColor',get(hf1,'Color'),...
                                'String',inputname(poffset+ipar),...
                                'Position',[20,cbm,50,sht]);
                if nval(ipar)>1
                hctls(ipar,2) = uicontrol(hf1,'Style','slider',...
                                'Max',nval(ipar),'Min',1,'Value',pind(ipar),...
                                'SliderStep',[1/nval(ipar)+1e-6 5/nval(ipar)+1e-6],...
                                'Position',[80,cbm,420,sht],...
                                'Callback',{@slider_Callback,ipar});
                else
                hctls(ipar,2) = uicontrol(hf1,'Style','text',...
                                'BackgroundColor',get(hf1,'Color'),...
                                'String','',...
                                'Position',[80,cbm,420,sht]);
                end
                hctls(ipar,3) = uicontrol(hf1,'Style','text',...
                                'BackgroundColor',get(hf1,'Color'),...
                                'String',num2str(p{ipar}(pind(ipar))),...
                                'Position',[505,cbm,50,sht]);
            end
            ha = axes(hf2,'Units','Pixels','OuterPosition',[0,0,560,aht]);

            % Initialize the GUI.
            % Change units to normalized so components resize automatically.
            set([hf1,reshape(hctls,1,[])],'Units','normalized');
            set([hf2,ha],'Units','normalized');

            %Create a plot in the axes.
            update_plot

            % Assign the GUI a name to appear in the window title.
            set(hf1,'Name',['Control: ',fname])
            set(hf2,'Name',['Figure: ',fname])
            % Move the GUI to the center of the screen.
            movegui(hf1,'west')
            movegui(hf2,'east')
            % Make the GUI visible.
            set(hf1,'Visible','on');
            set(hf2,'Visible','on');

            % Subfuntions implementing the callbacks
            % ------------------------------------------

            function slider_Callback(source,eventdata,ipar) 
                val = round(get(source,'Value'));
                set(source,'Value',val);
                pind(ipar) = val;
                set(hctls(ipar,3),'String',num2str(p{ipar}(val)));
                update_plot
            end

            function f2_closereq(src,callbackdata) 
                delete(hf2)
                delete(hf1)
            end

            % Subfuntions implementing other functions
            % ------------------------------------------

            function update_plot
                set(0,'CurrentFigure',hf2)
                f(pind);
                axis tight
                if iscell(fplot)
                    set(gca,apn,apv)
                end
            end

        end
        
        function surfall(X,Y,Zs,N)
            cla
            hold on
            if ndims(Zs)==3
                for ii = 1:size(Zs,1)
                    surf(X,Y,squeeze(Zs(ii,:,:)));
                end
            elseif ismatrix(Zs) && length(X)==prod(N)
                for ii = 1:size(Zs,2)
                    surf(reshape(X,N),reshape(Y,N),reshape(Zs(:,ii),N).');
                end
            end
            hold off
            axis equal
            axis tight
        end
    end
    
end

function plot_sp(ks,es)
% Plot spectrum
% ks is a column cell array of length (number of system dimensions)
% es(nes,[kz_pts,ky_pts,]kx_pts)
% inside [] is dimension dependent

nds = size(ks,1); % number of system dimensions
nes = size(es,1); % number of energy eigenstates

kds = [];
kks = {};
for id = 1:nds
    if length(ks{id})>1
        kds = [kds id]; %#ok<AGROW>
        kks = cat(1,kks,ks(id));
    end
end
nkds = length(kds);
ds = 2+nds-kds; % nontrivial dimensions in es

klabels = {'kx'; 'ky'; 'kz'};
axis_choices = klabels(kds);
switch nkds
    case 0
        plot_choices = {'dots'};
        plot_type = 'dots';
        xa_sel = {};
        ya_sel = {};
    case 1
        plot_choices = {'dots','line'};
        plot_type = 'line';
        xa_sel = axis_choices(1);
        ya_sel = {};
    otherwise
        plot_choices = {'dots','line','surf','pcolor'};
        plot_type = 'surf';
        xa_sel = axis_choices(1);
        ya_sel = subsref(setdiff(axis_choices,xa_sel),substruct('()',{1}));
end
fds = setdiff(axis_choices,cat(1,xa_sel,ya_sel));

xs = [];
ys = [];
iixa = [];
iiya = [];
iifd1 = [];
iifd2 = [];
indxa = {};
indya = {};
indfd1 = {};
indfd2 = {};

%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[360,500,640,450],...
    'NumberTitle', 'off');

%  Construct the components.
th = 320; % top height
htext0 = uicontrol('Style','text',...
                   'ForegroundColor','b',...
                   'String','Plot Type',...
                   'Position',[510,th,100,15]);
htextfn = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','func.',...
                'Position',[480,th-30,30,15]);
hpopupfn = uicontrol('Style','popupmenu',...
                'String',plot_choices,...
                'Value',find(strcmp(plot_type,plot_choices)),...
                'Position',[510,th-30,100,20],...
                'Callback',{@popupfn_Callback});
htextxa = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','xaxis',...
                'Position',[480,th-60,30,15]);
hpopupxa = uicontrol('Style','popupmenu',...
                'String',popup_str(axis_choices),...
                'Position',[510,th-60,100,20],...
                'Callback',{@popupax_Callback,1});
htextya = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','yaxis',...
                'Position',[480,th-90,30,15]);
hpopupya = uicontrol('Style','popupmenu',...
                'String',popup_str(setdiff(axis_choices,xa_sel)),...
                'Position',[510,th-90,100,20],...
                'Callback',{@popupax_Callback,2});
update_axis_sel(0);

htext1 = uicontrol('Style','text',...
                   'ForegroundColor','b',...
                   'String','Select Slice',...
                   'Position',[510,th-120,100,15]);
htextfd1 = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','',...
                'Position',[480,th-150,30,15]); 
hsliderfd1 = uicontrol('Style','slider',...
                'Max',10,'Min',1,'Value',1,...
                'SliderStep',[0.1 0.3],...
                'Position',[510,th-150,100,20],...
                'Callback',{@sliderfd_Callback,1}); 
htextfdv1 = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','',...
                'Position',[510,th-165,100,15]);
htextfd2 = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','',...
                'Position',[480,th-195,30,15]); 
hsliderfd2 = uicontrol('Style','slider',...
                'Max',10,'Min',1,'Value',1,...
                'SliderStep',[0.1 0.3],...
                'Position',[510,th-195,100,20],...
                'Callback',{@sliderfd_Callback,2}); 
htextfdv2 = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','',...
                'Position',[510,th-210,100,15]);
update_sliderfds;

ha = axes('Units','Pixels','Position',[60,60,400,350]); 

% Initialize the GUI.
% Change units to normalized so components resize automatically.
set([f,ha,...
    htext0,htextfn,hpopupfn,htextxa,hpopupxa,htextya,hpopupya,...
    htext1,htextfd1,hsliderfd1,htextfdv1,...
    htextfd2,hsliderfd2,htextfdv2],...
    'Units','normalized');

%Create a plot in the axes.
update_plot;

% Assign the GUI a name to appear in the window title.
set(f,'Name','Spectrum')
% Move the GUI to the center of the screen.
movegui(f,'center')
% Make the GUI visible.
set(f,'Visible','on');

% Subfuntions implementing the callbacks
% ------------------------------------------

    function popupfn_Callback(source,eventdata)  %#ok<*INUSD>
     % Determine the selected data set.
     str = get(source,'String');
     val = get(source,'Value');
     plot_type = str{val};
     update_axis_sel(0);
     update_sliderfds;
     update_plot;
    end

    function popupax_Callback(source,eventdata,ia)  %#ok<*INUSL>
     update_axis_sel(ia);
     update_sliderfds;
     update_plot;
    end

    function sliderfd_Callback(source,eventdata,ifd) 
     val = round(get(source,'Value'));
     set(source,'Value',val);
     switch ifd
         case 1
             fd_data = kks{iifd1};
             set(htextfdv1,'String',num2str(fd_data(val)));
             indfd1 = {val};
         case 2
             fd_data = kks{iifd2};
             set(htextfdv2,'String',num2str(fd_data(val)));
             indfd2 = {val};
     end
     update_plot;
    end

% Subfuntions implementing other functions
% ------------------------------------------

    function update_axis_sel(ia)
        if ia<2
            if sum(strcmp(plot_type,{'dots','line','surf','pcolor'}))
                set([htextxa,hpopupxa],'Visible','on');
                str = get(hpopupxa,'String');
                val = get(hpopupxa,'Value');
                xa_sel = str(val);
                iixa = find(strcmp(xa_sel,axis_choices));
                xs = kks{iixa};
                indxa = {':'};
                set(hpopupya,'String',...
                    popup_str(setdiff(axis_choices,xa_sel)));
            else
                set([htextxa,hpopupxa],'Visible','off');
                xa_sel = {};
                xs = [];
                indxa = {};
                iixa = [];
            end
        end
        if sum(strcmp(plot_type,{'surf','pcolor'}))
            set([htextya,hpopupya],'Visible','on');
            str = get(hpopupya,'String');
            val = get(hpopupya,'Value');
            ya_sel = str(val);
            iiya = find(strcmp(ya_sel,axis_choices));
            ys = kks{iiya};
            indya = {':'};
        else
            set([htextya,hpopupya],'Visible','off');
            ya_sel = {};
            ys = [];
            indya = {};
            iiya = [];
        end
    end

    function update_sliderfds
        fds = setdiff(axis_choices,cat(1,xa_sel,ya_sel));
        nfd = length(fds);
        if nfd>0
            iifd1 = find(strcmp(fds{1},axis_choices));
            fd_data = kks{iifd1};
            set(htextfd1,'String',fds{1});
            set(hsliderfd1,'Max',length(fd_data)+0.0001,'Value',1,....
                'SliderStep',[1/length(fd_data) 2/length(fd_data)]);
            set(htextfdv1,'String',num2str(fd_data(1)));
            set([htext1,htextfd1,hsliderfd1,htextfdv1],'Visible','on');
            indfd1 = {1};
        else
            set([htext1,htextfd1,hsliderfd1,htextfdv1],'Visible','off');
            indfd1 = {};
            iifd1 = [];
        end
        if nfd>1
            iifd2 = find(strcmp(fds{2},axis_choices));
            fd_data = kks{iifd2};
            set(htextfd2,'String',fds{2});
            set(hsliderfd2,'Max',length(fd_data)+0.0001,'Value',1,....
                'SliderStep',[1/length(fd_data) 2/length(fd_data)]);
            set(htextfdv2,'String',num2str(fd_data(1)));
            set([htextfd2,hsliderfd2,htextfdv2],'Visible','on');
            indfd2 = {1};
        else
            set([htextfd2,hsliderfd2,htextfdv2],'Visible','off');
            indfd2 = {};
            iifd2 = [];
        end
    end

    function update_plot
        ind = [iiya iixa iifd1 iifd2];
        dsa = [1 ds(ind)]; % all nontrivial dimensions in es
        sp = permute(es,[dsa, setdiff(1:(1+nds),dsa)]);
        indsp = cat(2,':',indya,indxa,indfd1,indfd2);
        sp = subsref(sp,substruct('()',indsp));
        set(0,'CurrentFigure',f)
        cla(ha)
        switch plot_type
            case 'dots'
                if nkds<1
                    hold on
                    for ine = 1:nes
                        scatter(ha,0,sp(ine),'ok');
                    end
                    hold off
                    set(ha,'XTick',[]);
                    set(ha,'XTickLabel',{});
                else
                    plot(xs,sp,'ok');
                    xlabel(xa_sel);
                end
                ylabel(ha,'E');
            case 'line'
%                 plot(xs,sp,'o')
                plot(xs,sp,'k')
                xlabel(xa_sel)
                ylabel('E');
                axis tight
            case 'surf'
                hold on
                for ine = 1:nes
                    surf(ha,xs,ys,squeeze(sp(ine,:,:)));
                end
                hold off
                xlabel(xa_sel)
                ylabel(ya_sel)
                zlabel('E');
                axis tight
                shading interp
            case 'pcolor'
                hold on
                for ine = 1:nes
                    pcolor(ha,xs,ys,squeeze(sp(ine,:,:)));
                end
                hold off
                xlabel(xa_sel)
                ylabel(ya_sel)
                axis equal
                axis tight
                shading interp
        end
    end

end

function plot_wfs(ks,es,wfs)
% Plot wavefunctions
% ks is a column cell array of length (number of system dimensions)
% es(nes,[kz_pts,ky_pts,]kx_pts)
% wfs(dof,Lx,[Ly,Lz,]nes,[kz_pts,ky_pts,]kx_pts)
% inside [] is dimension dependent

nds = size(ks,1); % number of system dimensions
nes = size(wfs,nds+2); % number of energy eigenstates
dof = size(wfs,1); % number of components
% if (2+nds*2)~=length(size(wfs))
%     errordlg('Dimension inconsistency!');
%     return
% end

% reshape es to the same form as wfs
ens = reshape(es,[ones(1,nds+1) size(es)]);

rds = [];
kds = [];
rs = cell(nds,1);
for id = 1:nds
    rs{id} = (1:size(wfs,id+1))';
    if ~isempty(ks{id})
        kds = [kds id]; %#ok<AGROW>
    else
        rds = [rds id]; %#ok<AGROW>
    end
end
ds = [1+rds, 2+nds*2-kds+1]; % nontrivial dimensions in wfs
rks = cat(1,rs(rds),ks(kds));
rlabels = {'x'; 'y'; 'z'};
klabels = {'kx'; 'ky'; 'kz'};
axis_choices = cat(1,rlabels(rds),klabels(kds));

iixa = 1;
xs = rks{iixa};
xa_sel = axis_choices(iixa);
iiya = [];
ys = [];
ya_sel = {};

fds = setdiff(axis_choices,cat(1,xa_sel,ya_sel));
iifd1 = []; iifd2 = [];

indxa = {':'};
indya = {};
inden = {1};
indfd1 = {};
indfd2 = {};
part = 'abssq';
comp = 'sum';
scal = 'linear';

if nds==1
    plot_choices = {'line'};
else
    plot_choices = {'line','surf','pcolor'};
end
plot_type = 'line';

%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[360,500,640,450],...
    'NumberTitle', 'off');

%  Construct the components.
th = 400; % top height
htext0 = uicontrol('Style','text',...
                   'ForegroundColor','b',...
                   'String','Plot Type',...
                   'Position',[510,th,100,15]);
htextfn = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','func.',...
                'Position',[480,th-30,30,15]);
hpopupfn = uicontrol('Style','popupmenu',...
                'String',plot_choices,...
                'Value',find(strcmp(plot_type,plot_choices)),...
                'Position',[510,th-30,100,20],...
                'Callback',{@popupfn_Callback});
htextxa = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','xaxis',...
                'Position',[480,th-60,30,15]);
hpopupxa = uicontrol('Style','popupmenu',...
                'String',axis_choices,...
                'Value',iixa,...
                'Position',[510,th-60,100,20],...
                'Callback',{@popupax_Callback,1});
htextya = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','yaxis',...
                'Position',[480,th-90,30,15]);
hpopupya = uicontrol('Style','popupmenu',...
                'String',popup_str(setdiff(axis_choices,xa_sel)),...
                'Position',[510,th-90,100,20],...
                'Callback',{@popupax_Callback,2});
update_yaxis_sel;

htext1 = uicontrol('Style','text',...
                   'ForegroundColor','b',...
                   'String','Select Eigenstate',...
                   'Position',[510,th-120,100,15]);
htexten = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','#E',...
                'Position',[480,th-150,30,15]);
hpopupen = uicontrol('Style','popupmenu',...
                'String',num2str((1:nes)'),...
                'Value',1,...
                'Position',[510,th-150,100,20],...
                'Callback',{@popupen_Callback});
htextfd1 = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','',...
                'Position',[480,th-180,30,15]); 
hsliderfd1 = uicontrol('Style','slider',...
                'Max',10,'Min',1,'Value',1,...
                'SliderStep',[0.1 0.3],...
                'Position',[510,th-180,100,20],...
                'Callback',{@sliderfd_Callback,1}); 
htextfdv1 = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','',...
                'Position',[510,th-195,100,15]);
htextfd2 = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','',...
                'Position',[480,th-225,30,15]); 
hsliderfd2 = uicontrol('Style','slider',...
                'Max',10,'Min',1,'Value',1,...
                'SliderStep',[0.1 0.3],...
                'Position',[510,th-225,100,20],...
                'Callback',{@sliderfd_Callback,2}); 
htextfdv2 = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','',...
                'Position',[510,th-240,100,15]);
update_sliderfds;

htext2 = uicontrol('Style','text',...
                   'ForegroundColor','b',...
                   'String','Data Type',...
                   'Position',[510,th-270,100,15]);
htextpar = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','part',...
                'Position',[480,th-300,30,15]);
hpopuppar = uicontrol('Style','popupmenu',...
                'String',{'abssq','real','imag'},...
                'Position',[510,th-300,100,20],...
                'Callback',{@popuppar_Callback});
htextcom = uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','comp.',...
                'TooltipString','component',...
                'Position',[480,th-330,30,15]);
hpopupcom = uicontrol('Style','popupmenu',...
                'String',cat(1,'sum',cellstr(num2str((1:dof)'))),...
                'Position',[510,th-330,100,20],...
                'Callback',{@popupcom_Callback});
htextsca= uicontrol('Style','text',...
                'BackgroundColor',get(f,'Color'),...
                'String','scale',...
                'Position',[480,th-360,30,15]);
hpopupsca = uicontrol('Style','popupmenu',...
                'String',{'linear','log10'},...
                'Position',[510,th-360,100,20],...
                'Callback',{@popupsca_Callback});

ha = axes('Units','Pixels','Position',[60,60,400,350]); 

% Initialize the GUI.
% Change units to normalized so components resize automatically.
set([f,ha,...
    htext0,htextfn,hpopupfn,htextxa,hpopupxa,htextya,hpopupya,...
    htext1,htexten,hpopupen,htextfd1,hsliderfd1,htextfdv1,...
    htextfd2,hsliderfd2,htextfdv2,...
    htext2,htextpar,hpopuppar,htextcom,hpopupcom,htextsca,hpopupsca],...
    'Units','normalized');

%Create a plot in the axes.
update_plot;

% Assign the GUI a name to appear in the window title.
set(f,'Name','Wavefunctions')
% Move the GUI to the center of the screen.
movegui(f,'center')
% Make the GUI visible.
set(f,'Visible','on');

% Subfuntions implementing the callbacks
% ------------------------------------------

    function popupfn_Callback(source,eventdata)  %#ok<*INUSD>
     str = get(source,'String');
     val = get(source,'Value');
     plot_type = str{val};
     update_yaxis_sel;
     update_sliderfds;
     update_plot;
    end

    function popupax_Callback(source,eventdata,ia)  %#ok<*INUSL>
     str = get(source,'String');
     val = get(source,'Value');
     switch ia
         case 1
             xa_sel = str(val);
             iixa = find(strcmp(xa_sel,axis_choices));
             xs = rks{iixa};
             set(hpopupya,'String',popup_str(setdiff(axis_choices,xa_sel)));
             update_yaxis_sel;
         case 2
             ya_sel = str(val);
             iiya = find(strcmp(ya_sel,axis_choices));
             ys = rks{iiya};
     end
     update_sliderfds;
     update_plot;
    end

    function popupen_Callback(source,eventdata)
     inden = {get(source,'Value')};
     update_plot;
    end

    function sliderfd_Callback(source,eventdata,ifd)
     val = round(get(source,'Value'));
     set(source,'Value',val);
     switch ifd
         case 1
             fd_data = rks{iifd1};
             set(htextfdv1,'String',num2str(fd_data(val)));
             indfd1 = {val};
         case 2
             fd_data = rks{iifd2};
             set(htextfdv2,'String',num2str(fd_data(val)));
             indfd2 = {val};
     end
     update_plot;
    end

    function popuppar_Callback(source,eventdata)
     str = get(source,'String');
     val = get(source,'Value');
     part = str{val};
     update_plot;
    end

    function popupcom_Callback(source,eventdata)
     str = get(source,'String');
     val = get(source,'Value');
     comp = str{val};
     update_plot;
    end

    function popupsca_Callback(source,eventdata)
     str = get(source,'String');
     val = get(source,'Value');
     scal = str{val};
     update_plot;
    end

% Subfuntions implementing other functions
% ------------------------------------------

    function update_yaxis_sel
        if sum(strcmp(plot_type,{'surf','pcolor'}))
            set([htextya,hpopupya],'Visible','on');
            str = get(hpopupya,'String');
            val = get(hpopupya,'Value');
            ya_sel = str(val);
            iiya = find(strcmp(ya_sel,axis_choices));
            ys = rks{iiya};
            indya = {':'};
        else
            set([htextya,hpopupya],'Visible','off');
            ya_sel = {};
            ys = [];
            indya = {};
            iiya = [];
        end
    end

    function update_sliderfds
        fds = setdiff(axis_choices,cat(1,xa_sel,ya_sel));
        nfd = length(fds);
        if nfd>0
            iifd1 = find(strcmp(fds{1},axis_choices));
            fd_data = rks{iifd1};
            set(htextfd1,'String',fds{1});
            set(hsliderfd1,'Max',length(fd_data)+0.0001,'Value',1,....
                'SliderStep',[1/length(fd_data) 2/length(fd_data)]);
            set(htextfdv1,'String',num2str(fd_data(1)));
            set([htextfd1,hsliderfd1,htextfdv1],'Visible','on');
            indfd1 = {1};
        else
            set([htextfd1,hsliderfd1,htextfdv1],'Visible','off');
            indfd1 = {};
            iifd1 = [];
        end
        if nfd>1
            iifd2 = find(strcmp(fds{2},axis_choices));
            fd_data = rks{iifd2};
            set(htextfd2,'String',fds{2});
            set(hsliderfd2,'Max',length(fd_data)+0.0001,'Value',1,....
                'SliderStep',[1/length(fd_data) 2/length(fd_data)]);
            set(htextfdv2,'String',num2str(fd_data(1)));
            set([htextfd2,hsliderfd2,htextfdv2],'Visible','on');
            indfd2 = {1};
        else
            set([htextfd2,hsliderfd2,htextfdv2],'Visible','off');
            indfd2 = {};
            iifd2 = [];
        end
    end

    function update_plot
        ind = [iiya iixa iifd1 iifd2];
        dsa = [1 ds(ind) nds+2]; % all nontrivial dimensions in wfs
        wf = permute(wfs,[dsa, setdiff(1:(2+nds*2),dsa)]);
        en = permute(ens,[dsa, setdiff(1:(2+nds*2),dsa)]);
        if ~isempty(iifd1)
            if ds(iifd1)<nds+2 % is space dimension
                indfd1e = {':'}; % ens does not have space dimension
            else
                indfd1e = indfd1;
            end
        else
            indfd1e = {};
        end
        if ~isempty(iifd2)
            if ds(iifd2)<nds+2 % is space dimension
                indfd2e = {':'}; % ens does not have space dimension
            else
                indfd2e = indfd2;
            end
        else
            indfd2e = {};
        end
        indwf = cat(2,':',indya,indxa,indfd1,indfd2,inden);
        inde = cat(2,':',indya,indxa,indfd1e,indfd2e,inden);
        wf = subsref(wf,substruct('()',indwf));
        en = subsref(en,substruct('()',inde));
        switch part
            case 'abssq'
                wf = abs(wf).^2;
            case 'real'
                wf = real(wf);
            case 'imag'
                wf = imag(wf);
        end
        switch comp
            case 'sum'
                wf = squeeze(sum(wf,1));
            otherwise
                indwf = cat(2,str2num(comp),indya,indxa); %#ok<ST2NM>
                wf = squeeze(subsref(wf,substruct('()',indwf)));
        end
        switch scal
            case 'linear'
            case 'log10'
                wf = log10(abs(wf));
        end
        set(0,'CurrentFigure',f)
        switch plot_type
            case 'line'
                plot(xs,wf,'x-')
                xlabel(xa_sel)
                axis tight
            case 'surf'
                if length(ys)==1
                    warndlg('ys is singular! Change to line plot.')
                    set(0,'CurrentFigure',f)
                    plot(xs,wf,'x-')
                    xlabel(xa_sel)
                    axis tight
                elseif length(xs)==1
                    warndlg('xs is singular! Change to line plot.')
                    set(0,'CurrentFigure',f)
                    plot(ys,wf,'x-')
                    xlabel(ya_sel)
                    axis tight
                else
                    surf(xs,ys,wf)
                    xlabel(xa_sel)
                    ylabel(ya_sel)
                    shading interp
                end
            case 'pcolor'
                if length(ys)==1
                    warndlg('ys is singular! Change to line plot.')
                    set(0,'CurrentFigure',f)
                    plot(xs,wf,'x-')
                    xlabel(xa_sel)
                    axis tight
                elseif length(xs)==1
                    warndlg('xs is singular! Change to line plot.')
                    set(0,'CurrentFigure',f)
                    plot(ys,wf,'x-')
                    xlabel(ya_sel)
                    axis tight
                else
                    pcolor(xs,ys,wf)
                    xlabel(xa_sel)
                    ylabel(ya_sel)
                    shading interp
                end
        end
        if isscalar(squeeze(en))
            title(ha,['E = ',num2str(en)])
        end
    end

end

function str = popup_str(str0)
    if isempty(str0)
        str = {'N/A'};
    else
        str = str0;
    end
end