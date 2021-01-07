% eegplugin_std_envtopo() - std_envtopo plugin for updated std_envtopo
%
% Usage:
%   >> eegplugin_std_envtopo(fig, try_strings, catch_strings);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   To create a new plugin, simply create a file beginning with "eegplugin_"
%   and place it in your eeglab folder. It will then be automatically
%   detected by eeglab. See also this source code internal comments.
%   For eeglab to return errors and add the function's results to
%   the eeglab history, menu callback must be nested into "try" and
%   a "catch" strings. For more information on how to create eeglab
%   plugins, see http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% See also: pop_std_envtopo, std_envtopo, envtopo, pop_envtopo

%123456789012345678901234567890123456789012345678901234567890123456789012

% Author: Makoto Miyakoshi & Arnaud Delorme, JSPS/SCCN, INC, UCSD
% History:
% 02/05/2014 ver 1.5 by Makoto. Bug fixed; cb_clustedit was preventing hisotry function. 
% 10/09/2013 ver 1.4 by Makoto. Took care of 'args' 'LASTCOM'.
% 05/09/2013 ver 1.3 by Makoto. Gave up 'com'.
% 01/12/2012 ver 1.2 by Makoto. Button renamed with caution.
% 12/23/2011 ver 1.1 by Makoto. Renamed the function.
% 08/01/2011 ver 1.0 by Makoto. Created.

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function vers = eegplugin_std_evntopo( fig, trystrs, catchstrs);

vers = 'std_envtopo 4.01';
if nargin < 3
    error('eegplugin_std_evntopo requires 3 arguments');
end;

% create menu
std = findobj(fig, 'tag', 'study');
uimenu( std, 'label', 'Run Envtopo analysis', 'callback', 'gui_std_envtopo', 'userdata', 'startup:off;study:on');
