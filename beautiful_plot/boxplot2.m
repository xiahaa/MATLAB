function h = boxplot2(varargin)
%BOXPLOT2 Enhanced boxplot plots
% 
% h = boxplot2(y)
% h = boxplot2(y,x)
% h = boxplot2(..., p1, v1, ...)
%
% I don't like the original boxplot function... it messes with axis
% dimensions, replaces tick labels with text (even if no rotation is
% needed), and in general hijacks the axis too much.  Plus it doesn't
% return the handles to the plotted objects in an easy-to-use way (I can
% never remember which row corresponds to which part of the boxplot).  This
% version creates more light-handed boxplots, assuming that any cosmetic
% changes (color, tick labels, line specs, etc) can be added by the user
% afterwards if necessary.  It also allows one to create clustered
% boxplots, similar to an unstacked bar graph.
%
% Input variables:
%
%   y:              either a ndata x nx array (as in boxplot) or nx x ny x
%                   ndata array where nx indicates the number of
%                   x-clusters, ny the number of boxes per cluster, and
%                   ndata the number of points per boxplot.
%
%   x:              vector of x locations for each box cluster
%
% Optional input variables (passed as parameter/value pairs)
%
%   notch:          'on' or 'off' ['off']
%
%   orientation:    'vertical' or 'horizontal' ['vertical']
%
%   barwidth:       Barwidth value used to position boxes (see bar) [0.8]
%
%   whisker:        whisker length factor (see boxplot) [1.5]
%
%   axes:           axis to plot to [gca]
%
% Output variables:
%
%   h:              structure of handles to boxplot, with the following
%                   fields: 
%                   'box':      box
%                   'ladj':     lower adjacent value
%                   'lwhis':    lower whisker
%                   'med':      median
%                   'out':      outliers
%                   'uadj':     upper adjacent value
%                   'uwhis':    upper whisker    
%

% Copyright 2012 Kelly Kearney       

% Parse input

p = inputParser;
p.addRequired('y', @isnumeric);
p.addOptional('x', [], @isnumeric);
p.addParamValue('notch', 'off', @ischar);
p.addParamValue('orientation', 'vertical', @ischar);
p.addParamValue('axes', gca, @(x) isscalar(x) && ishandle(x) && strcmp(get(x,'type'),'axes'));
p.addParamValue('barwidth', 0.8, @(x) isscalar(x) && x > 0 && x <= 1);
% p.addParamValue('boxwidth', [], @(x) isscalar(x));
p.addParamValue('whisker', 1.5, @(x) isscalar(x));


p.parse(varargin{:});

In = p.Results;
In.notch = validatestring(In.notch, {'on', 'off'});
In.orientation = validatestring(In.orientation, {'vertical', 'horizontal'});

if ndims(In.y) == 2
    In.y = permute(In.y, [2 3 1]);
end
[nx, ny, ndata] = size(In.y);

if isempty(In.x)
    In.x = 1:nx;
end

ybox = reshape(In.y, [], ndata)';

% Use bar graph to get x positions

figtmp = figure('visible', 'off');
try
    hax = axes;
    hb = bar(In.x, In.y(:,:,1), In.barwidth);
    for ib = 1:length(hb)
        if verLessThan('matlab','8.4.0')
            xbar = get(get(hb(ib), 'children'), 'xdata');
            xb(:,ib) = mean(xbar,1);
        else
            xb(:,ib) = hb(ib).XData + hb(ib).XOffset;
        end
    end
    if verLessThan('matlab', '8.4.0')
        boxwidth = diff(minmax(xbar(:,1)));
    else
        if ny > 1
            boxwidth = diff([hb(1:2).XOffset])*In.barwidth;
        else
            mindx = min(diff(In.x));
            boxwidth = mindx .* In.barwidth;
        end
    end
    delete(hb);

    boxplot(ybox, 'positions', xb(:), ...
                  'notch', In.notch, ...
                  'orientation', In.orientation, ...
                  'symbol', '+', ...
                  'widths', boxwidth, ...
                  'whisker', In.whisker);

    h.box   = copyobj(findall(hax, 'tag', 'Box'), In.axes);
    h.ladj  = copyobj(findall(hax, 'tag', 'Lower Adjacent Value'), In.axes);
    h.lwhis = copyobj(findall(hax, 'tag', 'Lower Whisker'), In.axes);
    h.med   = copyobj(findall(hax, 'tag', 'Median'), In.axes);
    h.out   = copyobj(findall(hax, 'tag', 'Outliers'), In.axes);
    h.uadj  = copyobj(findall(hax, 'tag', 'Upper Adjacent Value'), In.axes);
    h.uwhis = copyobj(findall(hax, 'tag', 'Upper Whisker'), In.axes);

    close(figtmp);
catch ME
    close(figtmp);
    rethrow(ME);
end

    h = structfun(@(x) reshape(flipud(x), ny, nx), h, 'uni', 0);

end

function varargout = minmax(a, type, whis)
    %MINMAX Returns minimum and maximum value in the given array
    %
    % [minval maxval] = minmax(a)
    % lims = minmax(a);
    % lims = minmax(a, type);
    % lims = minmax(a, type, w);
    %
    % Computes the minimum and maximum value in entire array (all dimensions).
    %
    % Input variables:
    %
    %   a:      numeric array
    %
    %   type:   'all':          absolute minimum and maximum (default)   
    %           'noout':        discards outliers
    %           'center':       centers on zero
    %           'centernoout':  centers on 0 and eliminates outliers
    %           'expand':       wider that 'all' by a specifed fraction
    %
    %   w:      for no-outlier version, factor defining an outlier.  Point is
    %           considered an outlier if larger than q3+w*(q3-q1) or smaller
    %           than q1-w*(q3-q1) [1.5]  
    %           
    %           for expand version, fraction of actual range to add onto each
    %           end [0.1]
    %
    % Output variables:
    %
    %   minval: minimum value in a
    %
    %   maxval: maximum value in a

    % Copywrite 2005 Kelly Kearney

    if nargin == 1
        type = 'all';
    end

    switch type
        case 'all'
            minval = min(a(:));
            maxval = max(a(:));
        case 'noout'
            if nargin < 3
                whis = 1.5;
            end
            pctiles = prctile(a(:),[25 75]);
            q1 = pctiles(1);
            q3 = pctiles(2);
            vhi = q3+whis*(q3-q1);
            vlo = q1-whis*(q3-q1);
            minval = min(a(a > vlo));
            maxval = max(a(a < vhi));
            if isempty(minval)
                minval = vlo;
            end
            if isempty(maxval)
                maxval = vhi;
            end
        case 'center'
            temp = max(abs([min(a(:)) max(a(:))]));
            minval = -temp;
            maxval = temp;
        case 'centernoout'
            if nargin < 3
                whis = 1.5;
            end
            pctiles = prctile(a(:),[25 75]);
            q1 = pctiles(1);
            q3 = pctiles(2);
            vhi = q3+whis*(q3-q1);
            vlo = q1-whis*(q3-q1);
            temp = max(abs([min(a(a > vlo)) max(a(a < vhi))]));
            minval = -temp;
            maxval = temp;
        case 'expand'
            if nargin < 3
                whis = 0.1;
            end
            minval = min(a(:));
            maxval = max(a(:));
            da = maxval - minval;
            minval = minval - whis*da;
            maxval = maxval + whis*da;
        otherwise
            error('Unrecognized option: %s', type);
    end


    if nargout == 2
        varargout{1} = minval;
        varargout{2} = maxval;
    elseif nargout == 1
        varargout{1} = [minval maxval];
    elseif nargout == 0
        [minval maxval]
    else
        error('Wrong number of output arguments');
    end
end

