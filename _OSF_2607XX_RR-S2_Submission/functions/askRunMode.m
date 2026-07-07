function mode = askRunMode()
%ASKRUNMODE  Fast modal dialog (legacy HG) with bigger fonts + multi-monitor.
%   mode = askRunMode()

fontSize = 25;

% Headless? Fallback to CLI.
if ~usejava('desktop')
    txt  = input('Process all subjects? [y/N]: ','s');
    mode = ternary(strcmpi(strtrim(txt),'y'),'all','new');
    return
end

% Monitor selection
mp = get(0,'MonitorPositions');  % [x y w h] per monitor
%monitorIndex = min(2, size(mp,1));  % prefer second if available
monitorIndex = 1;
mon = mp(monitorIndex,:);

% Dialog geometry
W = max(420, 28*fontSize);    % width scales with font
H = max(180, 12*fontSize);    % height scales with font
px = mon(1) + (mon(3)-W)/2;
py = mon(2) + (mon(4)-H)/2;

% Build modal figure (classic, fast)
bg = 0.97*[1 1 1];

% Build figure invisible, normal first
f = figure('Name','Run mode','NumberTitle','off', ...
    'MenuBar','none','ToolBar','none','Resize','off', ...
    'Color',bg, 'Visible','off', 'Units','pixels',...
    'DefaultUicontrolFontSize',fontSize, ...
    'DefaultUicontrolFontName','Helvetica');

% Position BEFORE making it modal
set(f,'Position',[px py W H]);  % uses your computed [px py W H]
drawnow;                        % let HG apply the position

% Now make it modal and ensure it's on-screen
set(f,'WindowStyle','modal');
movegui(f,'onscreen');          % guard against OS nudges

% Layout numbers (pixels)
pad   = round(0.08*W);
btnH  = max(28, round(2.2*fontSize));
btnW  = round((W - 2*pad - 2*12)/3);
lblY  = H - pad - 2.8*fontSize;
subY  = lblY - 1.8*fontSize;
btnY  = pad;

% Message
uicontrol(f,'Style','text','String','Process all subjects or only new ones?', ...
    'Units','pixels','Position',[pad lblY W-2*pad 2*fontSize], ...
    'BackgroundColor',bg,'HorizontalAlignment','center', ...
    'FontWeight','bold');

% Return value holder
setappdata(f,'mode','new');

% Buttons
x1 = pad;               x2 = x1 + btnW + 12;     x3 = x2 + btnW + 12;
uicontrol(f,'Style','pushbutton','String','ALL', ...
    'Units','pixels','Position',[x1 btnY btnW btnH], ...
    'Callback',@(h,~) choose('all'));
uicontrol(f,'Style','pushbutton','String','NEW', ...
    'Units','pixels','Position',[x2 btnY btnW btnH], ...
    'Callback',@(h,~) choose('new'));
uicontrol(f,'Style','pushbutton','String','Cancel', ...
    'Units','pixels','Position',[x3 btnY btnW btnH], ...
    'Callback',@(h,~) choose('new'));

% Keyboard shortcuts: Enter/Esc default to 'new', 'a' -> all, 'n' -> new
%set(f,'KeyPressFcn',@(~,e) keyChoice(e));

% Show and block
set(f,'Visible','on'); drawnow;
uiwait(f);

% Output
if isvalid(f)
    mode = getappdata(f,'mode');
    delete(f);
else
    mode = 'new';
end

if strcmp(mode, 'all')
    clc
    disp(upper('processing all subjects...'));
else
    clc
    disp(upper('processing only new subjects...'));
end

% Nested helpers
    function choose(m)
        if isvalid(f), setappdata(f,'mode',m); uiresume(f); end
    end
    function keyChoice(e)
        k = lower(e.Key);
        switch k
            case {'return','enter','escape'}
                choose('new');
            case 'a'
                choose('all');
            case 'n'
                choose('new');
        end
    end
end

function out = ternary(cond,a,b)
if cond, out = a; else, out = b; end
end
