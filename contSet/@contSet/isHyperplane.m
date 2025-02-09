function res = isHyperplane(S)
% isHyperplane - (DEPRECATED -> representsa)
%
% Syntax:
%    res = isHyperplane(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

funcname = mfilename;
warning(sprintf(['The function ''' funcname ''' is deprecated (since CORA 2024) and has been replaced by ''representsa''.\n' ...
    '         When updating the code, please rename every function call ''' funcname '(S)'' -> ''representsa(S,''hyperplane'')''.\n' ...
    '         Note that the function ''' funcname ''' will be removed in a future release.']));
res = representsa(S,'hyperplane');

% ------------------------------ END OF CODE ------------------------------
