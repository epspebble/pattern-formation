function options = abmset(varargin)
%ABMSET Create/alter ABM OPTIONS structure.
%   OPTIONS = ABMSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an options
%   structure OPTIONS in which the named properties have the specified
%   values. Any unspecified properties have default values. It is
%   sufficient to type only the leading characters that uniquely identify
%   the property. Case is ignored for property names.  
%   
%   OPTIONS = ABMSET(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS.
%   
%   OPTIONS = ABMSET(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties
%   overwrite corresponding old properties. 
%   
%   ABMSET with no input arguments displays all property names and their
%   possible values.
%   
%ABMSET PROPERTIES
%
% AGENT PARAMETERS
%
%ntype - Number of type of cells 
%   Integer value, default to 2, meaning melanophore and xanthophore.
%
%domx, domy - number of evenly spaced grid points in x and y direction
%   Default to 2000 and 1000.
%
%domxt, domyt - 
%
%rm, rx - 
%
%gamma_loc - 
%
%gamma_long_inner, gamma_long_outer - Cytoneme annular neighborhood
%   Integer value, specifies the inner and outer radius of long-range
%   itneraction by cytonemes.
%
% BIRTH-DEATH PARAMETERS
%
%mu, nu, psi, p_d - 
% alpha, beta, eta, d_crowd - 
% phi_1, phi_2, kappa, d_rand - 
%
% INITIAL PATTERN PARAMETERS
% hmwidth - 
%

% Base on code by GOU Jia, and modified from ODESET on 26/01/23, which was by:
%   Mark W. Reichelt and Lawrence F. Shampine, 5/6/94
%   Copyright 1984-2016 The MathWorks, Inc.

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
  fprintf('Agent parameters... \n')
  fprintf('           ntype:  [ 2 ]\n');
  fprintf('            domx: [ 2000 ]\n');
  fprintf('            domy: [ 1000 ]\n');
  fprintf('           domxt: [ 130 ]\n');
  fprintf('           domyt: [ 130 ]\n');
  fprintf('              rm: [ 20 ]\n');
  fprintf('              rx: [ 20 ]\n');
  fprintf('       gamma_loc: [ 75 ]\n');
  fprintf('gamma_long_inner: [ 75 ]\n');
  fprintf('gamma_long_outer: [ 75 ]\n');
  fprintf('\n');
  fprintf('Birth-death parameters... \n');
  fprintf('              mu: [ 1] \n');
  fprintf('              nu: [ 1] \n');
  fprintf('             psi: [ 1.2] \n');
  fprintf('             p_d: [ 1/30] \n');
  fprintf('           alpha: [ 1] \n');
  fprintf('            beta: [ 3.5] \n');
  fprintf('             eta: [ 6] \n');
  fprintf('         d_crowd: [ 82] \n');
  fprintf('           phi_1: [ 1.3] \n');
  fprintf('           phi_2: [ 1.2] \n');
  fprintf('           kappa: [ 10] \n');
  fprintf('          d_rand: [ 100] \n');
  fprintf('               a: [ 3000] \n');
  fprintf('               b: [ 300] \n');
  fprintf('\n');
  fprintf('Initial pattern parameters... \n');
  fprintf('         hmwidth: [ 100] \n');
  fprintf('\n');
  fprintf('Simulation control parameters... \n');
  fprintf('       plot_init: [ false] \n');
  fprintf('\n');
  return;
end

Names = [
    'ntype           '
    'domx            '
    'domy            '
    'domxt           '
    'domyt           '
    'rm              '
    'rx              '
    'gamma_loc       '
    'gamma_long_inner'
    'gamma_long_outer'
    'mu              '
    'nu              '
    'psi             '
    'p_d             '
    'alpha           '
    'beta            '
    'eta             '
    'd_crowd         '
    'phi_1           '
    'phi_2           '
    'kappa           '
    'd_rand          '
    'a               '
    'b               '
    'hmwidth         '
    'plot_init       '
    ];
m = size(Names,1);
names = lower(Names);


%% The following is an almost exact copy from odeset's corresponding section.

% Combine all leading options structures o1, o2, ... in abmset(o1,o2,...).
options = [];
for j = 1:m
  options.(deblank(Names(j,:))) = [];
end
i = 1;

while i <= nargin
  arg = varargin{i};
  if ischar(arg) || (isstring(arg) && isscalar(arg)) % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(message('MATLAB:odeset:NoPropNameOrStruct', i));
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
        val = arg.(deblank(Names(j,:)));
      else
        val = [];
      end
      if ~isempty(val)
        options.(deblank(Names(j,:))) = val;
      end
    end
  end
  i = i + 1;
end
% Convert string arguments and options.
for ii = 1:nargin
    if isstring(varargin{ii}) && isscalar(varargin{ii})
        varargin{ii} = char(varargin{ii});
    end
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error(message('MATLAB:odeset:ArgNameValueMismatch'));
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~ischar(arg)
      error(message('MATLAB:odeset:NoPropName', i));
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(message('MATLAB:odeset:InvalidPropName', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
            matches = deblank(Names(j(1),:));
        for k = j(2:length(j))'
                matches = [matches ', ' deblank(Names(k,:))]; %#ok<AGROW>
        end
            error(message('MATLAB:odeset:AmbiguousPropName',arg,matches));
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error(message('MATLAB:odeset:NoValueForProp', arg));
end
