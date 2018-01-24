function [u, du, cc, m] = addDisplacements(u0,du0,cc0,m0,dm)
% u = addDisplacements(u,thr,epsilon) removes outliers using the universal
% outlier test based on
%
% J. Westerweel and F. Scarano. Universal outlier detection for PIV data.
% Exp. Fluids, 39(6):1096{1100, August 2005. doi: 10.1007/s00348-005-0016-6
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
%   thr: theshold for passing residiual (default = 2)
%   epsilon: fluctuation level due to cross-correlation (default = 0.1)
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: cell containing the displacement field with outliers removed
%   normFluctValues: normalized fluctuation values based on the universal
%   outier test.
%
% NOTES
% -------------------------------------------------------------------------
% needs medFilt3 and John D'Errico's inpaint_nans3
% (http://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint-nans)function.
%
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2

for i = 1:3, du0{i} = inpaint_nans3(du0{i}); end % remove NaNs if present

idx = cell(1,3);
for i = 1:3, idx{i} = m0{i}(1):dm:m0{i}(end); end % construct new meshgrid

[m0_{1}, m0_{2}, m0_{3}] = ndgrid(m0{1},m0{2},m0{3});
[m{1}, m{2}, m{3}] = ndgrid(idx{1},idx{2},idx{3});
% sample to desired mesh spacing

du = cell(1,3);
for i = 1:3, 
    F = griddedInterpolant(m0_{1}, m0_{2}, m0_{3}, du0{i}, 'linear');
    du{i} = F(m{1},m{2},m{3});  
end

F = griddedInterpolant(m0_{1}, m0_{2}, m0_{3}, cc0, 'linear');
cc = F(m{1},m{2},m{3});

if  sum(cellfun(@numel, u0)) == 3, u = du; % on first iteration u = du
else u = cellfun(@plus,u0,du,'UniformOutput', 0); % else u^(k) = sum(u^(k-1)) + du (see eq. 7)
end

end