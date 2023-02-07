classdef bdControlVectorDialog < handle
    %bdControlVectorDialog  Dialog box for editing a vector control widget
    %   This class specialised to work in tandem with the control panel.
    %   It should not be called directly by the user. 
    %
    %AUTHORS
    %  Stewart Heitmann (2017c-d,2018a-b,2019a)

    % Copyright (C) 2016-2019 QIMR Berghofer Medical Research Institute
    % All rights reserved.
    %
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions
    % are met:
    %
    % 1. Redistributions of source code must retain the above copyright
    %    notice, this list of conditions and the following disclaimer.
    % 
    % 2. Redistributions in binary form must reproduce the above copyright
    %    notice, this list of conditions and the following disclaimer in
    %    the documentation and/or other materials provided with the
    %    distribution.
    %
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    % "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    % LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    % FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    % COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
    % INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
    % BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    % LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    % CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    % LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
    % ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    % POSSIBILITY OF SUCH DAMAGE.
    
    properties (Access=private)
        control         % handle to control panel
        dialogfig       % handle to dialog box figure
        datatable       % handle to the table widget
        baraxes         % handle to bar graph axes
        bargraph        % handle to bar graph widget
        histogrm        % handle to histogram widget
        minbox          % handle to minbox
        maxbox          % handle to maxbox
        zerobutton      % handle to zero button
        perbbutton      % handle to perb button
        randbutton      % handle to rand button
        haltbutton      % handle to halt button
        listener1       % handle to listener1
        listener2       % handle to listener2
    end
    
    methods
        % Constructs a bdControlVectorDialog dialog box where 
        % control = handel to the bdControl object
        % xxxdef = 'pardef' or 'vardef' or 'lagdef' (string).
        % xxxindx is an index of the sys.xxxdef array.
        function this = bdControlVectorDialog(control,xxxdef,xxxindx,titlestr)
            % remember the control panel handle
            this.control = control; 

            % extract the relevant fields from control.sys.xxxdef
            xxxname  = control.sys.(xxxdef)(xxxindx).name;
            xxxvalue = control.sys.(xxxdef)(xxxindx).value;
            xxxlim   = control.sys.(xxxdef)(xxxindx).lim;
            
            % construct dialog box (at the current mouse position)
            xypos = get(groot,'PointerLocation'); 
            this.dialogfig = figure('Units','pixels', ...
                'Position',[xypos(1) xypos(2), 400, 330], ...
                'MenuBar','none', ...
                'Name',titlestr, ...
                'NumberTitle','off', ...
                'ToolBar', 'none', ...
                'Resize','off', ...
                'DeleteFcn', @(~,~) delete(this) );

            % data table 
            this.datatable = uitable(this.dialogfig,'Position',[10 10 140, 300], ...
                'Data',reshape(xxxvalue,[],1), ...       % ensure data is column vector
                'ColumnName',{xxxname}, ...
        ...        'ColumnWidth',{75}, ...
                'ColumnEditable',true, ...
                'CellEditCallback', @(src,~) this.CellEditCallback(xxxdef,xxxindx));

            % axes for vector values (bar graph)
            this.baraxes = axes('parent',this.dialogfig, ...
                'Units','pixels', ...
                'Position',[195 250 190 60]); 
            
            % axes for histogram
            ax2 = axes('parent',this.dialogfig, ...
                'Units','pixels', ...
                'Position',[195 135 190 60]);
            
            % bar graph
            n = numel(xxxvalue);
            this.bargraph = bar(xxxvalue, 'parent',this.baraxes);
            xlim(this.baraxes,[0.5 n+0.5]);
            ylim(this.baraxes,xxxlim+[-1e-6,1e-6]);
            xlabel('parent',this.baraxes,'index');
            ylabel('parent',this.baraxes,'value');
            title(this.baraxes,xxxname);
            grid(this.baraxes,'on');
            
            % histogram with 10 bins (11 edges)
            this.histogrm = histogram(xxxvalue, 11, ...
                'parent',ax2, ...
              ...  'BinLimits', xxxlim, ...
              ...  'BinMethod','sturges', ...
                'Normalization','probability');
            xlabel('parent',ax2,'value');
           %ylabel('parent',ax2,'proportion');
            title(ax2,'Histogram');
           
            % 'ZERO' button
            this.zerobutton = uicontrol('Style','pushbutton', ...
                'String','ZERO', ...
                'HorizontalAlignment','center', ...
                'FontUnits','pixels', ...
                'FontSize',12, ...
                'Parent', this.dialogfig, ...
                'Callback', @(~,~) this.ZeroCallback(xxxdef,xxxindx), ...
                'Position',[195 65 60 20], ...
                'ToolTipString','Zero the data');            

            % 'PERB' button
            this.perbbutton = uicontrol('Style','pushbutton', ...
                'String','PERB', ...
                'HorizontalAlignment','center', ...
                'FontUnits','pixels', ...
                'FontSize',12, ...
                'Parent', this.dialogfig, ...
                'Callback', @(~,~) this.PerbCallback(xxxdef,xxxindx), ...
                'Position',[260 65 60 20], ...
                'ToolTipString','Uniform perturbation (5%)');            

            % 'RAND' button
            this.randbutton = uicontrol('Style','pushbutton', ...
                'String','RAND', ...
                'HorizontalAlignment','center', ...
                'FontUnits','pixels', ...
                'FontSize',12, ...
                'Parent', this.dialogfig, ...
                'Callback', @(~,~) this.RandCallback(xxxdef,xxxindx), ...
                'Position',[325 65 60 20], ...
                'ToolTipString','Uniform random data');            

            % min box
            this.minbox = uicontrol('Parent',this.dialogfig, ...
                'Style', 'edit', ...
                'Units','pixels',...
                'Position',[195 40 60 20], ...
                'String',num2str(xxxlim(1),'%0.4g'), ...
                'Value',xxxlim(1), ...
                'HorizontalAlignment','center', ...
                'Visible','on', ...
                'Callback', @(~,~) this.minboxCallback(xxxdef,xxxindx), ...
                'ToolTipString','Lower limit');

            % max box
            this.maxbox = uicontrol('Parent',this.dialogfig, ...
                'Style', 'edit', ...
                'Units','pixels',...
                'Position',[260 40 60 20], ...
                'String',num2str(xxxlim(2),'%0.4g'), ...
                'Value',xxxlim(2), ...
                'HorizontalAlignment','center', ...
                'Visible','on', ...
                'Callback', @(~,~) this.maxboxCallback(xxxdef,xxxindx), ...
                'ToolTipString','Upper limit');

            % HALT button
            this.haltbutton = uicontrol('Style','radio', ...
                'String','HALT', ...
                'Value',this.control.sys.halt, ...
                'HorizontalAlignment','left', ...
                'FontUnits','pixels', ...
                'FontSize',12, ...
                'FontWeight','bold', ...
                'ForegroundColor', 'r', ...
                'Parent', this.dialogfig, ...
                'ToolTipString', 'Halt the solver', ...
                'Callback', @(src,~) this.HaltCallback(src), ...
                'Position',[195 10 60 20]);

            % 'Close' button
            uicontrol('Style','pushbutton', ...
                'String','Close', ...
                'HorizontalAlignment','center', ...
                'FontUnits','pixels', ...
                'FontSize',12, ...
                'Parent', this.dialogfig, ...
                'Callback', @(~,~) this.visible('off'), ...
                'Position',[325 10 60 20], ...
                'ToolTipString','Close the dialog box');
 
            % force a refresh at startup
            this.refreshListener(xxxdef,xxxindx);   

            % listen to the control panel for widget refresh events (incuding those generate by this dialog box)
            this.listener1 = addlistener(control,'refresh',@(~,~) this.refreshListener(xxxdef,xxxindx));   
            this.listener2 = addlistener(control,xxxdef,@(~,~) this.refreshListener(xxxdef,xxxindx));   
        end
        
        % Destructor (called when this object is no longer referenced)
        function delete(this)
            delete(this.listener2);
            delete(this.listener1);
            delete(this.dialogfig);
        end
        
        % Make the dialog box visible/invisible
        function visible(this,flag)
            figure(this.dialogfig);
            this.dialogfig.Visible = flag;
        end
    end
    
    methods (Access=private)
        % TABLE cell edit callback
        function CellEditCallback(this,xxxdef,xxxindx)
            %disp('CellEditCallback'); 
            % get the data from the table
            data = get(this.datatable,'data');
            
            % reshape the data to match the sys.xxxdef.value.
            data = reshape(data,size(this.control.sys.(xxxdef)(xxxindx).value));
            
            % update the control panel.
            this.control.sys.(xxxdef)(xxxindx).value = data;

            % update bar graph
            xxxlim = this.control.sys.(xxxdef)(xxxindx).lim;
            this.bargraph.YData = data;
            this.baraxes.YLim = xxxlim + [-1e-6 1e-6];

            % update histogram
            this.histogrm.Data = data;
            this.histogrm.BinLimitsMode='auto';

            % notify all widgets (except ourself) that sys.xxxdef has changed
            this.listener1.Enabled = 0;
            this.listener2.Enabled = 0;
            notify(this.control,xxxdef);
            this.listener1.Enabled = 1;
            this.listener2.Enabled = 1;
            
            % tell the solver to recompute the solution
            if ~this.control.sys.halt
                notify(this.control,'recompute');
            end
        end
        
        % ZERO button callback
        function ZeroCallback(this,xxxdef,xxxindx)
            % update the control panel.
            valsize = size(this.control.sys.(xxxdef)(xxxindx).value);
            this.control.sys.(xxxdef)(xxxindx).value = zeros(valsize);
            
            % notify all widgets (which includes ourself) that sys.xxxdef has changed
            notify(this.control,xxxdef);
            
            % tell the solver to recompute the solution
            if ~this.control.sys.halt
                notify(this.control,'recompute');
            end
        end

        % RAND button callback
        function RandCallback(this,xxxdef,xxxindx)
            % determine the limits of the random values
            xxxlim = this.control.sys.(xxxdef)(xxxindx).lim;
            lo = xxxlim(1);
            hi = xxxlim(2);
            
            % update the control panel.
            valsize = size(this.control.sys.(xxxdef)(xxxindx).value);
            this.control.sys.(xxxdef)(xxxindx).value = (hi-lo)*rand(valsize) + lo;
            
            % notify all widgets (which includes ourself) that sys.xxxdef has changed
            notify(this.control,xxxdef);
            
            % tell the solver to recompute the solution
            if ~this.control.sys.halt
                notify(this.control,'recompute');
            end
        end

        % PERB button callback
        function PerbCallback(this,xxxdef,xxxindx)
            % determine the limits of the random values
            xxxlim = this.control.sys.(xxxdef)(xxxindx).lim;
            lo = xxxlim(1);
            hi = xxxlim(2);
            
            % update the control panel with a perturned version of the data
            valsize = size(this.control.sys.(xxxdef)(xxxindx).value);
            this.control.sys.(xxxdef)(xxxindx).value =  ...
                this.control.sys.(xxxdef)(xxxindx).value + ...
                0.05*(hi-lo)*(rand(valsize)-0.5);

            % notify all widgets (which includes ourself) that sys.xxxdef has changed
            notify(this.control,xxxdef);
            
            % tell the solver to recompute the solution
            if ~this.control.sys.halt
                notify(this.control,'recompute');
            end
        end
        
        % min box callback function
        function minboxCallback(this,xxxdef,xxxindx)
            % read the minbox string and convert to a number
            str = this.minbox.String;
            minval = str2double(str);
            if isnan(minval)
                hndl = errordlg(['Invalid number: ',str], 'Invalid Number', 'modal');
                uiwait(hndl);
                % restore the minbox string to its previous value
                this.minbox.String = num2str(this.minbox.Value,'%0.4g');                 
            else           
                % adjust the max box if necessary
                maxval = max(this.maxbox.Value, minval);
                
                % update control.sys
                this.control.sys.(xxxdef)(xxxindx).lim = [minval maxval];
                
                % notify all widgets (which includes ourself) that sys.xxxdef has changed
                notify(this.control,xxxdef);

                % notify all display panels to redraw themselves
                notify(this.control,'redraw');
            end
        end        
        
        % max box callback function
        function maxboxCallback(this,xxxdef,xxxindx)
            % read the maxbox string and convert to a number
            str = this.maxbox.String;
            maxval = str2double(str);
            if isnan(maxval)
                hndl = errordlg(['Invalid number: ',str], 'Invalid Number', 'modal');
                uiwait(hndl);
                % restore the minbox string to its previous value
                this.maxbox.String = num2str(this.maxbox.Value,'%0.4g');                 
            else           
                % adjust the min box if necessary
                minval = min(this.minbox.Value, maxval);
                
                % update control.sys
                this.control.sys.(xxxdef)(xxxindx).lim = [minval maxval];
                
                % notify all widgets (which includes ourself) that sys.xxxdef has changed
                notify(this.control,xxxdef);

                % notify all display panels to redraw themselves
                notify(this.control,'redraw');
            end
        end
        
        % HALT button callback
        function HaltCallback(this,haltbutton)
            this.control.sys.halt = haltbutton.Value;   % get the HALT button state
            notify(this.control,'refresh');             % notify all widgets to refresh themselves
            if ~this.control.sys.halt
                notify(this.control,'recompute');       % tell the solver to recompute
            end
        end
        
        % Listener for widget refresh events from the control panel
        function refreshListener(this,xxxdef,xxxindx)
            %disp(['bdControlVectorDialog.refreshListener:' xxxdef])         
            
            % extract the data from control.sys.xxxdef
            xxxvalue = this.control.sys.(xxxdef)(xxxindx).value;
            xxxlim = this.control.sys.(xxxdef)(xxxindx).lim;
          
            % rehsape the value data to a column vector
            data = reshape(xxxvalue,[],1);
            
            % update the data table
            this.datatable.Data = data;
            
            % update bar graph
            this.bargraph.YData = data;
            this.baraxes.YLim = xxxlim + [-1e-6 1e-6];

            % update histogram
            this.histogrm.Data = data;
            this.histogrm.BinLimitsMode='auto';

            % update the min box
            this.minbox.Value = xxxlim(1);
            this.minbox.String = num2str(xxxlim(1),'%0.4g');
           
            % update the max box
            this.maxbox.Value = xxxlim(2);
            this.maxbox.String = num2str(xxxlim(2),'%0.4g');

            % update the HALT button
            this.haltbutton.Value = this.control.sys.halt; 
            
            % special case: if this is a vardef control and the evolve button
            % is ON then disable the edit buttons.
            switch xxxdef
                case 'vardef'
                    if this.control.sys.evolve
                        % disable the edit buttons
                        this.datatable.Enable = 'off';
                        this.zerobutton.Enable = 'off';
                        this.perbbutton.Enable = 'off';
                        this.randbutton.Enable = 'off';
                    else
                        % enable the buttons
                        this.datatable.Enable = 'on';
                        this.zerobutton.Enable = 'on';
                        this.perbbutton.Enable = 'on';
                        this.randbutton.Enable = 'on';
                    end
            end

        end
        
    end
    
end

