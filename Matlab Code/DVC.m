function [u, cc] = DVC(varargin)
% [du, cc] = DVC(I,sSize,sSpacing,ccThreshold) estimates
% displacements between two volumetric images through digital volume
% correlation.
%
% INPUTS
% -------------------------------------------------------------------------
%   I: cell containing the undeformed, I{1}, and deformed, I{2} 3-D images
%   sSize: interrogation window (subset) size
%   sSpacing: interrogation window (subset) spacing.  Determines window
%             overlap factor
%   ccThreshold: threshold value that defines a bad cross-correlation
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: cell containing the displacement field (u{1:3} = {u_x, u_y, u_z})
%   cc: peak values of the cross-correlation for each interrogation
%
% NOTES
% -------------------------------------------------------------------------
% all functions are self contained
%
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2

% Parse inputs and create meshgrid
[I,m,mSize,sSize,MTF,M,ccThreshold] = parseInputs(varargin{:});

% Initialize variables
mSize_ = prod(mSize);
u123 = zeros(mSize_,3);
cc = zeros(mSize_,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wb = findall(0,'Tag','TMWWaitbar'); wb = wb(1);

waitbar(1/7,wb,'Estimating Displacements (Time Remaining: )');

for k = 1:mSize_
    
    
    tStart = tic; % begin timer
    %--------------------------------------------------------------------------
    % grab the moving subset from the images
    A = I{1}(m{1}(k,:),m{2}(k,:),m{3}(k,:));
    B = I{2}(m{1}(k,:),m{2}(k,:),m{3}(k,:));
    
    % multiply by the modular transfer function to alter frequency content
    A = MTF.*A; B = MTF.*B;
    
    % run cross-correlation
    A = xCorr3(A,B,sSize);
    
    % find maximum index of the cross-correlaiton
    [cc(k), maxIdx] = max(A(:));
    
    % compute voxel resolution displacements
    [u1, u2, u3] = ind2sub(sSize,maxIdx);
    
    % gather the 3x3x3 voxel neighborhood around the peak
    try xCorrPeak = reshape(A(u1 + (-1:1), u2 + (-1:1), u3 + (-1:1)),27,1);        
        % last squares fitting of the peak to calculate sub-voxel displacements
        du123 = lsqPolyFit3(xCorrPeak, M{1}, M{2});       
        u123(k,:) = [u1 u2 u3] + du123' - sSize/2 - 1;
        %--------------------------------------------------------------------------
    catch
        u123(k,:) = nan;
    end    
%     xCorrPeak = reshape(A(u1 + (-1:1), u2 + (-1:1), u3 + (-1:1)),27,1);
%     
%     % least squares fitting of the peak to calculate sub-voxel displacements
%     du123 = lsqPolyFit3(xCorrPeak, M{1}, M{2});
%     u123(k,:) = [u1 u2 u3] + du123' - (sSize/2) - 1;
    %--------------------------------------------------------------------------
    
    % waitbar calculations (update only every 100 iterations)
    if rem(k,100) == 0
        tRemaining = (toc(tStart)*(mSize_ - k)); % Time remaining for waitbar
        waitbar(1/7*(k/mSize_ + 1),wb,['Estimating Displacements (Time Remaining: ', datestr(datenum(0,0,0,0,0,tRemaining),'MM:SS'),')'])
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reshape displacements and set bad correlations to zero
waitbar(2/7,wb,'Removing Bad Correlations')

cc = reshape(double(cc),mSize);
% cc = permute(cc,[2 1 3]);
[cc, ccMask] = removeBadCorrelations(I,cc,ccThreshold);

% u = cell(1,3);
% for i = 1:3
%     u{i} = reshape(double(u123(:,i)),mSize).*ccMask;
%     u{i} = permute(u{i},[2 1 3]);
% end
u{1} = reshape(double(u123(:,2)),mSize).*ccMask;
u{2} = reshape(double(u123(:,1)),mSize).*ccMask;
u{3} = reshape(double(u123(:,3)),mSize).*ccMask;

end

%% ========================================================================
function varargout = parseInputs(varargin)
% Parse inputs and create meshgrid

I{1} = varargin{1}{1};
I{2} = varargin{1}{2};
sSize = varargin{2};
sSpacing = varargin{3};
padSize = varargin{4};
ccThreshold = varargin{5};

% pad images with zeros so that we don't grab any subset outside of the image
% domain. This would produce an error
I{1} = padarray(I{1},padSize,0,'both');
I{2} = padarray(I{2},padSize,0,'both');
sizeV = size(I{1});

% Initialize Mesh Variables
idx = cell(1,3);
for i = 1:3, idx{i} = (1+padSize(i)) : sSpacing(i) : (sizeV(i)-sSize(i)-padSize(i)+1); end
[m{1},m{2},m{3}] = ndgrid(idx{:});


% sSize = [sSize(2) sSize(1) sSize(3)];
mSize = size(m{1});
mSize_ = prod(mSize);

m_ = cell(1,3);
for k = 1:3,
    m_{k} = zeros([mSize_,sSize(k)],'uint16');
    repmat_ = repmat((1:sSize(k))-1,mSize_,1);
    m_{k} = bsxfun(@plus, repmat_,m{k}(:));
end

% Initialize quadratic least squares fitting coefficients
[mx, my, mz] = meshgrid((-1:1),(-1:1),(-1:1));
m = [mx(:), my(:), mz(:)];

M{1} = zeros(size(m,1),10);
for i = 1:size(m,1)
    x = m(i,1); y = m(i,2); z = m(i,3);
    M{1}(i,:) = [1,x,y,z,x^2,x*y,x*z,y^2,y*z,z^2];
end

M{2} = M{1}'*M{1};

% Generate Moduluar transfer function (see eq. 3)
[~,~,MTF] = generateMTF(sSize);

%% Parse outputs

varargout{    1} = I;
varargout{end+1} = m_;
varargout{end+1} = mSize;
varargout{end+1} = sSize;
varargout{end+1} = MTF;
varargout{end+1} = M;
varargout{end+1} = ccThreshold;

end

%% ========================================================================
function A = xCorr3(A,B,sSize)
% performs fft based cross correlation of A and B (see equation 2)

A = fftn(A,sSize);
B = fftn(B,sSize);
B = conj(B);
A = A.*B;
A = ifftn(A);
A = real(A);
A = fftshift(A);
end

%% ========================================================================
function    duvw = lsqPolyFit3(b, M, trMM)
% LeastSqPoly performs a 3D polynomial fit in the least squares sense
% Solves M*x = b,
% trMM = transpose(M)*M
% trMb = tranpose(M)*b
%
% If you need to generate the coefficients then uncomment the following
% [mx, my, mz] = meshgrid(-1:1,-1:1,-1:1);
% m = [mx(:), my(:), mz(:)];
%
% for i = 1:size(m,1)
%    x = m(i,1); y = m(i,2); z = m(i,3);
%    M1(i,:) = [1,x,y,z,x^2,x*y,x*z,y^2,y*z,z^2];
% end
%
% trMM1 = M'*M;

% b = log(b);
trMb = sum(bsxfun(@times, M, b));

x = trMM\trMb'; %solve for unknown coefficients

A = [x(6), 2*x(5), x(7);
    2*x(8),  x(6) x(9)
    x(9),    x(7), 2*x(10)];

duvw = (A\(-x([2 3 4])));
end

%% ========================================================================
function varargout = generateMTF(sSize)
% MTF functions taken from
% J. Nogueira, A Lecuona, P. A. Rodriguez, J. A. Alfaro, and A. Acosta.
% Limits on the resolution of correlation PIV iterative methods. Practical
% implementation and design of weighting functions. Exp. Fluids,
% 39(2):314{321, July 2005. doi: 10.1007/s00348-005-1017-1

%% equation 4

if prod(single(sSize == 32)) || prod(single(sSize == 16))
    sSize = sSize(1);
    
    x = cell(1,3);
    [x{1}, x{2}, x{3}] = meshgrid(1:sSize,1:sSize,1:sSize);
    
    nu{1} = 1;
    for i = 1:3
        x{i} = x{i} - sSize/2 - 0.5;
        x{i} = abs(x{i}/sSize);
        nu{1} = nu{1}.*(3*(4*x{i}.^2-4*x{i}+1));
    end
    
    %% equation 5
    [x{1}, x{2}, x{3}] = meshgrid(1:sSize,1:sSize,1:sSize);
    
    for i = 1:3, x{i} = x{i} - sSize/2 - 0.5; end
    
    r = abs(sqrt(x{1}.^2 + x{2}.^2 + x{3}.^2)/sSize);
    nu{2}  = zeros(size(r));
    nu{2}(r < 0.5) = 24/pi*(4*r(r < 0.5).^2-4*r(r < 0.5)+1);
    
    %% equation 6
    [x{1}, x{2}, x{3}] = meshgrid(1:sSize,1:sSize,1:sSize);
    
    nu{3} = 1;
    for i = 1:3,
        x{i} = x{i} - sSize/2 - 0.5;
        x{i} = (x{i}/sSize);
        
        nu{3} = nu{3}.*(12*abs(x{i}).^2 - 12*abs(x{i}) + 3 + ...
            0.15*cos(4*pi*x{i}) + 0.20*cos(6*pi*x{i}) + ...
            0.10*cos(8*pi*x{i}) + 0.05*cos(10*pi*x{i}));
        
    end
    nu{3}(nu{3} < 0) = 0;
    
else
    
    nu{1} = ones(sSize(1),sSize(2),sSize(3));
    nu{2} = nu{1};
    nu{3} = nu{1};
    
end

nu = cellfun(@(x) x/sum(x(:)), nu, 'UniformOutput',0);
nu = cellfun(@sqrt, nu, 'UniformOutput',0);

varargout = nu;

end

%% ========================================================================
function [cc, ccMask] = removeBadCorrelations(I,cc,ccThreshold)
% removes bad correlations.  You can insert your own method here.
minOS = 1;
for i = 1:2
    zeroIdx = I{i} == 0;
    threshold = mean2(I{i}(~zeroIdx));
    I_ = I{i}.*(I{i} < threshold);
    minOS = minOS*sum(I_(:))/sum(~zeroIdx(:));
end
cc = cc - minOS;
cc = cc/(max(I{1}(:))*max(I{2}(:)));
ccMask = double(cc >= ccThreshold);

CC = bwconncomp(~ccMask);
[~,idx] = max(cellfun(@numel,CC.PixelIdxList));
if ~isempty(idx)
    ccMask(CC.PixelIdxList{idx}) = inf;
end
ccMask(cc == 0) = nan;
ccMask(~isfinite(ccMask)) = 0;
cc = cc.*ccMask;

end
