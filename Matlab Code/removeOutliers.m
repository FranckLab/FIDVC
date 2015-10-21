function [u,normFluctValues] = removeOutliers(u,thr,epsilon)
% u = removeOutliers(u,thr,epsilon) removes outliers using the universal
% outlier test based on
%
% J. Westerweel and F. Scarano. Universal outlier detection for PIV data.
% Exp. Fluids, 39(6):1096{1100, August 2005. doi: 10.1007/s00348-005-0016-6
%
% INPUTS
% -------------------------------------------------------------------------
%   u: cell containing the input displacement field. (u{1:3} = {u_x, u_y,
%   	u_z})
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
% For more information please see section 2.2.
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2

% set default values
if nargin < 3, epsilon = 0.1; end
if nargin < 2, thr = 2; end
if ~iscell(u), u = {u}; end


medianU = cell(size(u));
normFluct = cell(size(u));
normFluctMag = zeros(size(u{1}));

for i = 1:length(u)

    [medianU{i}, normFluct{i}] = funRemoveOutliers(u{i},epsilon);
    normFluctMag = normFluctMag + normFluct{i}.^2;
end

normFluctMag = sqrt(normFluctMag);
outlierIdx = find(normFluctMag > thr);
normFluctValues = normFluctMag(outlierIdx);

for i = 1:length(u),
    u{i}(outlierIdx) = nan;
    u{i} = inpaint_nans3(double(u{i}),0);
end

end

%% ========================================================================
function [medianU, normFluct] = funRemoveOutliers(u,epsilon)

nSize = 3*[1 1 1];
skipIdx = ceil(prod(nSize)/2);
padOption = 'symmetric';

u =  inpaint_nans3(double(u),0);

medianU = medFilt3(u,nSize,padOption,skipIdx);
fluct = u - medianU;
medianRes = medFilt3(abs(fluct),nSize,padOption,skipIdx);
normFluct = abs(fluct./(medianRes + epsilon));

end

%% ========================================================================
function Vr = medFilt3(V0,nSize, padoption, skipIdx)
% fast median filter for 3D data with extra options.

if nargin < 4, skipIdx = 0; end
if nargin < 3, padoption = 'symmetric'; end
if nargin < 2, nSize = [3 3 3]; end

nLength = prod(nSize);
if mod(nLength,2) == 1, padSize = floor(nSize/2);
elseif mod(nLength,2) == 0, padSize = [nSize(1)/2-1,nSize(2)/2];
end

if strcmpi(padoption,'none')
    V = V0;
else
    V = (padarray(V0,padSize(1)*[1,1,1],padoption,'pre'));
    V = (padarray(V,padSize(2)*[1,1,1],padoption,'post'));
end

S = size(V);
nLength = prod(nSize)-sum(skipIdx>1);
Vn = single(zeros(S(1)-(nSize(1)-1),S(2)-(nSize(2)-1),S(3)-(nSize(3)-1),nLength));  % all the neighbor

%%
% build the neighboor

i = cell(1,nSize(1)); j = cell(1,nSize(2)); k = cell(1,nSize(3));
for m = 1:nSize(1), i{m} = m:(S(1)-(nSize(1)-m)); end
for m = 1:nSize(2), j{m} = m:(S(2)-(nSize(2)-m)); end
for m = 1:nSize(3), k{m} = m:(S(3)-(nSize(3)-m)); end

p = 1;
for m = 1:nSize(1)
    for n = 1:nSize(2)
        for o = 1:nSize(3)
            if p ~= skipIdx || skipIdx == 0
                Vn(:,:,:,p) = V(i{m},j{n},k{o});
            end
            p = p + 1;
        end
    end
end

if skipIdx ~= 0, Vn(:,:,:,skipIdx) = []; end
% perform the processing
Vn = sort(Vn,4);

if mod(nLength,2) == 1 % if odd get the middle element
    Vr = Vn(:,:,:,ceil(nLength/2));
else % if even get the mean of the two middle elements
    Vr = mean(cat(4,Vn(:,:,:,nLength/2),Vn(:,:,:,nLength/2+1)),4);
end

end
