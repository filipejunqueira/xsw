function hh = polarerror(varargin)
% Modified ERRORBAR function to work with polarplots
% hh = polarerror(theta,r,theta error,r error,opts)


%ERRORBAR Plot error bars along curve
%   ERRORBAR(Y,E) plots Y and draws a vertical error bar at each element of
%   Y. The error bar is a distance of E(i) above and below the curve so
%   that each bar is symmetric and 2*E(i) long.
%
%   ERRORBAR(X,Y,E) plots Y versus X with symmetric vertical error bars
%   2*E(i) long. X, Y, E must be the same size. When they are vectors, each
%   error bar is a distance of E(i) above and below the point defined by
%   (X(i),Y(i)). When they are matrices, each error bar is a distance of
%   E(i,j) above and below the point defined by (X(i,j),Y(i,j)).
%
%   ERRORBAR(X,Y,NEG,POS) plots X versus Y with vertical error bars
%   NEG(i)+POS(i) long specifying the lower and upper error bars. X and Y
%   must be the same size. NEG and POS must be the same size as Y or empty.
%   When they are vectors, each error bar is a distance of NEG(i) below and
%   POS(i) above the point defined by (X(i),Y(i)). When they are matrices,
%   each error bar is a distance of NEG(i,j) below and POS(i,j) above the
%   point defined by (X(i,j),Y(i,j)). When they are empty the error bar is
%   not drawn.
%
%   ERRORBAR( ___ ,Orientation) specifies the orientation of the error
%   bars. Orientation can be 'horizontal', 'vertical', or 'both'. When the
%   orientation is omitted the default is 'vertical'.
%
%   ERRORBAR(X,Y,YNEG,YPOS,XNEG,XPOS) plots X versus Y with vertical error
%   bars YNEG(i)+YPOS(i) long specifying the lower and upper error bars and
%   horizontal error bars XNEG(i)+XPOS(i) long specifying the left and
%   right error bars. X and Y must be the same size. YNEG, YPOS, XNEG, and
%   XPOS must be the same size as Y or empty. When they are empty the error
%   bar is not drawn.
%
%   ERRORBAR( ___ ,LineSpec) specifies the color, line style, and marker.
%   The color is applied to the data line and error bars. The line style
%   and marker are applied to the data line only.
%
%   ERRORBAR(AX, ___ ) plots into the axes specified by AX instead of the
%   current axes.
%
%   H = ERRORBAR( ___ ) returns handles to the errorbarseries objects
%   created. ERRORBAR creates one object for vector input arguments and one
%   object per column for matrix input arguments.
%
%   Example: Draws symmetric error bars of unit standard deviation.
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      errorbar(x,y,e)

%   L. Shure 5-17-88, 10-1-91 B.A. Jones 4-5-93
%   Copyright 1984-2017 MathWorks, Inc.

% Look for a parent among the input arguments
[~, cax, args] = parseplotapi(varargin{:},'-mfilename',mfilename);

% We require at least one input
narginchk(1, inf);

% Separate Name/Value pairs from data inputs, convert LineSpec to
% Name/Value pairs, and filter out the orientation flag.
[pvpairs,args,nargs,orientation] = parseargs(args);

pvpairs = matlab.graphics.internal.convertStringToCharArgs(pvpairs);
% Check that we have the correct number of data input arguments.
if nargs < 2
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargs > 4 && ~isempty(orientation)
    error(message('MATLAB:errorbar:InvalidUseOfOrientation'));
elseif nargs == 5
    error(message('MATLAB:errorbar:InvalidNumberDataInputs'));
elseif nargs > 6
    error(message('MATLAB:narginchk:tooManyInputs'));
end


% Grab the X data if present.
if nargs >= 3
    % errorbar(x,y,e,...)
    x = args{1};
    args = args(2:end);
else
    % errorbar(y,e)
    x = [];
end

% Grab the Y data
y = checkSingleInput(args{1}, [], 'YData');
sz = size(y);
n = numel(y);

% Now that we have the size of the YData, validate the size of the XData
x = checkSingleInput(x, sz, 'XData');

xpos = args{2};
ypos = args{3};

% Handle vectorized data sources and display names
extrapairs = cell(n,0);


% Determine the Color and LineStyle property names
% If the Color/LineStyle is not specified use the _I property names so that
% the ColorMode or LineStyleMode properties are not toggled.
if ~any(strcmpi('color',pvpairs(1:2:end)))
    c = 'k';
else
    c = pvpairs{find(strcmpi(pvpairs,'color'))+1};
end

if ~any(strcmpi('marker',pvpairs(1:2:end)))
    m = 'x';
else
    m = pvpairs{find(strcmpi(pvpairs,'marker'))+1};
end

if ~any(strcmpi('linestyle',pvpairs(1:2:end)))
    ls = 'none';
else
    ls = pvpairs{find(strcmpi(pvpairs,'linestyle'))+1};
end


removeidx = find(strcmpi(pvpairs,'color'));
pvpairs = pvpairs([1:removeidx-1,removeidx+2:end]);
removeidx = find(strcmpi(pvpairs,'marker'));
pvpairs = pvpairs([1:removeidx-1,removeidx+2:end]);
removeidx = find(strcmpi(pvpairs,'linestyle'));
pvpairs = pvpairs([1:removeidx-1,removeidx+2:end]);

% Create the ErrorBar objects
h = gobjects(1,n);
h2 = gobjects(1,n);
for k = 1:n
    % extract data from vectorizing over columns
    
    
      
    tick = 0.05;

    r = y(k);
    t = x(k);
    
    re = ypos(k);
    te = xpos(k);
    t_tick = tick/2;
    r_tick_m = tick./(r-re)/2;
    r_tick_p = tick./(r+re)/2;
    
    precision = 5;
    
    t_plot = [[0,linspace(0,-r_tick_m,precision+1),linspace(-r_tick_m,r_tick_m,2*precision+1),linspace(r_tick_m,0,precision+1),0]+t,...
        [0,linspace(0,-r_tick_p,precision+1),linspace(-r_tick_p,r_tick_p,2*precision+1),linspace(r_tick_p,0,precision+1),0]+t,...
        t,...
        [linspace(0,-te,precision+1),ones(1,3)*-te]+t,...
        [linspace(-te,te,2*precision+1),ones(1,3)*te]+t];
    
    r_plot = [[r,ones(1,4*precision+3)*r,r]-re,...
        [r,ones(1,4*precision+3)*r,r]+re,...
        r,...
        ones(1,precision+1)*r,[-t_tick,t_tick,0]+r,...
        ones(1,2*precision+1)*r,[-t_tick,t_tick,0]+r];
    
    h(k) = polarplot(t_plot,r_plot,'color',c,'LineStyle','-','marker','none',...
        pvpairs{:},extrapairs{k,:});
    hold on

    

   
    
end

h2 = polarplot(x,y,'color',c,'marker',m,'linestyle',ls,...
        pvpairs{:});


hold off

if nargout>0, hh = [h,h2]; end
disp(pvpairs)
end


%-------------------------------------------------------------------------%
function [pvpairs,args,nargs,orientation] = parseargs(args)
% separate pv-pairs from opening arguments
[args,pvpairs] = parseparams(args);

% Check for LineSpec or Orientation strings
% Allow the orientation flag to occur either before or after the LineSpec
% Allow LineSpec and Orientation to occur at most once each.
validOrientations = {'horizontal','vertical','both'};
orientation = '';
keepArg = true(1,numel(pvpairs));
extraPairs = {};
for a = 1:min(2,numel(pvpairs))
    if matlab.graphics.internal.isCharOrString(pvpairs{a})
        % Check for partial matching of the orientation flag using a
        % minimum of 3 characters.
        tf = strncmpi(pvpairs{a},validOrientations,max(3,numel(pvpairs{a})));
        if isempty(orientation) && any(tf)
            orientation = validOrientations{tf};
            keepArg(a) = false;
        else
            % Check for LineSpec string
            [l,c,m,tmsg]=colstyle(pvpairs{a},'plot');
            if isempty(tmsg) && isempty(extraPairs)
                keepArg(a) = false;
                if ~isempty(l)
                    extraPairs = {'LineStyle',l};
                end
                if ~isempty(c)
                    extraPairs = [{'Color',c},extraPairs]; %#ok<AGROW>
                end
                if ~isempty(m)
                    extraPairs = [{'Marker',m},extraPairs]; %#ok<AGROW>
                end
            else
                break;
            end
        end
    else
        % Not a string, so stop looking.
        break
    end
end

pvpairs = [extraPairs, pvpairs(keepArg)];
nargs = numel(args);

end

%-------------------------------------------------------------------------%
function val = checkSingleInput(val, sz, propName)

if isvector(val)
    val = val(:);
end

if ~isempty(sz) && ~isempty(val) && ~isequal(sz, size(val))
    if strcmp(propName,'XData')
        error(message('MATLAB:errorbar:XDataSizeMismatch'));
    else
        error(message('MATLAB:errorbar:DeltaSizeMismatch', propName));
    end
end

end

%-------------------------------------------------------------------------%
function col = getColumn(val, k)
    if isempty(val)
        col = val;
    else
        col = datachk(val(:,k));
    end
end
