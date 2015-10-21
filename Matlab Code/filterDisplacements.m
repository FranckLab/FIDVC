function u = filterDisplacements(u0,filterSize,z)
% I = filterDisplacements(I0,filterSize,z) applies a low-pass convolution
% filter to the displacement field to mitigate divergence based on 
%
% F. F. J. Schrijer and F. Scarano. Effect of predictor corrector filtering
% on the stability and spatial resolution of iterative PIV interrogation. 
% Exp. Fluids, 45(5):927{941, May 2008. doi: 10.1007/s00348-008-0511-7
% 
% INPUTS
% -------------------------------------------------------------------------
%   u0: displacement field vector defined at every meshgrid point with 
%      spacing dm. Format: cell array, each containing a 3D matrix 
%         (components in x,y,z)
%         u0{1} = displacement in x-direction
%         u0{2} = displacement in y-direction
%         u0{3} = displacement in z-direction
%         u0{4} = magnitude
%   filterSize: size of the filter
%   z: filter strength
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: cell containing the filtered displacement field
%
% NOTES
% -------------------------------------------------------------------------
% none
% 
% For more information please see section 2.2.
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2
% 

% Parse inputs and set defaults
if nargin < 3,  z = 0.0075; end
if ~iscell(u0), u0 = {u0}; end
u = cell(size(u0));

if z == 0,
    u = cellfun(@double, u0, 'UniformOutput',0); % no filter
else
    rf = generateFilter(filterSize,z);
 
    % apply filter using convolution
    for i = 1:length(u0), u{i} = double(convn(u0{i}, rf,'same')); end
end

end

%% ========================================================================
function rf = generateFilter(filterSize,z)
% generates the filter convolution filter
l = filterSize;

[m{1}, m{2}, m{3}] = ndgrid(-l(1)/2:l(1)/2,-l(2)/2:l(2)/2,-l(3)/2:l(3)/2);
m = cellfun(@abs, m, 'UniformOutput', 0);


f1 = (l(1)/2)^z - (m{1}).^z;
f2 = (l(2)/2)^z - (m{2}).^z;
f3 = (l(3)/2)^z - (m{3}).^z;

i{1} = (m{1} >= m{2} & m{1} >= m{3});
i{2} = (m{1} < m{2} & m{2} >= m{3});
i{3} = (m{3} > m{2} & m{3} > m{1});

rf0 = f1.*i{1}+f2.*i{2}+f3.*i{3};

rf = rf0/sum(rf0(:));

end