function guidlg_time(varargin);

persistent hfig_dlg tstart ancho;

if (nargin == 0)
    opcion = 'inicializar';
else
    opcion = varargin{1};
end

if strcmp(opcion, 'inicializar')
    tstart = clock;
    ancho = 0.90;
    
    set(0, 'Units', 'pixels');
    rect_figure = get(0, 'ScreenSize');
    rect_figure(1) = rect_figure(3)/3.1;
    rect_figure(2) = rect_figure(4)/2.2;
    rect_figure(3) = rect_figure(3)-2*rect_figure(1);
    rect_figure(4) = rect_figure(4)-2*rect_figure(2);
    Option.Background_color = [0.8 0.8 0.8];
    hfig_dlg = figure('Color', Option.Background_color, 'MenuBar', 'none', 'Position', rect_figure);
    set(hfig_dlg, 'Units', 'normalized', 'NumberTitle', 'off', 'Name', '');
    % Create handle struct
    handles = guihandles(hfig_dlg);
    
%% controls for dialog
    uicontrol('Style', 'Text', 'String', 'Please wait...', 'HorizontalAlignment', 'center', ...
        'Units', 'normal', 'Position', [0.05 0.55 0.90 0.25], ...
        'FontSize', 10, 'BackgroundColor', Option.Background_color);
    handles.clockedit = uicontrol('Style', 'edit', ...
        'Units', 'normal', 'Position', [0.05 0.20 ancho 0.25], 'BackgroundColor', [1 1 1]);
    handles.timeedit = uicontrol('Style', 'edit', ...
        'Units', 'normal', 'Position', [0.05 0.20 ancho/1000 0.25], 'BackgroundColor', [1 0 0]);
    guidata(hfig_dlg, handles);
elseif strcmp(opcion, 'finalizar')
    delete(hfig_dlg);
else
    handles = guidata(hfig_dlg);
    
    if (opcion < 0) || (opcion > 1)
        error('Argument must be a real value in the range [0, 1].');
    end
    rect = get(handles.timeedit, 'Position');
    rect(3) = opcion*ancho;
    set(handles.timeedit, 'Position', rect);
    set(hfig_dlg, 'Name', [num2str(etime(clock, tstart)) ' seconds']);
end
return;