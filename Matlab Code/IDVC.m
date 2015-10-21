function [u, cc, dm] = IDVC(varargin)
% [u, cc] = IDVC(I,sSize,u0,className);
% I = filterDisplacements(I0,filterSize,z) applies a low-pass convolution
% filter to the displacement field to mitigate divergence based on 
%
% F. F. J. Schrijer and F. Scarano. Effect of predictor corrector filtering
% on the stability and spatial resolution of iterative PIV interrogation. 
% Exp. Fluids, 45(5):927{941, May 2008. doi: 10.1007/s00348-008-0511-7
% 
% INPUTS
% -------------------------------------------------------------------------
%   I0: cell containing the undeformed, I0{1}, and deformed, I0{2} 3-D
%       images
%   sSize: interrogation window (subset) size
%   u0: pre-estimated displacement field (typically zeros)
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: displacement field vector defined at every meshgrid point with 
%      spacing dm. Format: cell array, each containing a 3D matrix 
%         (components in x,y,z)
%         u{1} = displacement in x-direction
%         u{2} = displacement in y-direction
%         u{3} = displacement in z-direction
%         u{4} = magnitude
%   cc: peak values of the cross-correlation for each interrogation
%   dm: meshgrid spacing (8 by default)
%
% NOTES
% -------------------------------------------------------------------------
% To run you need a compatible C compiler. Please see
% (http://www.mathworks.com/support/compilers/R2014a/index.html)
% 
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2


wb = waitbar(0,'Parsing Inputs','Name','Running IDVC');

% PRESET CONSTANTS
maxIterations = 20; % maximum number of iterations
dm = 8; % desired output mesh spacing
convergenceCrit = [0.25, 0.5, 0.0625]; % convergence criteria
ccThreshold = 1e-4; % bad cross-correlation threshold

[I0, sSize, sSpacing, padSize, DVCPadSize, u] = parseInputs(varargin{:});

% START ITERATING
i = 2; converged01 = 0; SSE = []; I = I0;

t0 = tic;
while ~converged01 && i - 1 < maxIterations
    ti = tic;
    
    set(wb,'name',['Running IDVC (Iteration ',num2str(i-1),')']);
    waitbar(0/7,wb,'Checking Convergence');
    % Check for convergence
     [converged01, SSE(i-1) , sSize(i,:), sSpacing(i,:)] = checkConvergenceSSD(I,SSE,sSize,sSpacing,convergenceCrit);

    if ~converged01
        [I, m] = parseImages(I,sSize(i,:),sSpacing(i,:));
        
        % run cross-correlation to get an estimate of the displacements
        [du, cc] = DVC(I,sSize(i,:),sSpacing(i,:),DVCPadSize,ccThreshold);
        
        % add the displacements from previous iteration to current
        waitbar(3/7,wb,'Adding displacements from previous iteration');
        [u, ~, cc] = addDisplacements(u,du,cc,m,dm);
        
        % filter the  displacements using a predictor filter
        waitbar(4/7,wb,'Filtering Displacements');
        u = filterDisplacements(u,sSize(i,:)/dm);
        
        % remove outliers in displacement field
        waitbar(5/7,wb,'Removing Outliers');
        u = removeOutliers(u);

        % mesh and pad images based on new subset size and spacing
        [I, m] = parseImages(I0,sSize(i,:),sSpacing(i,:));
        
        % map volumes based on displacment field
        waitbar(6/7,wb,'Warping Images');
        I = volumeMapping(I,m,u);
        
        disp(['Elapsed time (iteration ',num2str(i-1),'): ',num2str(toc(ti))]);
        i = i + 1;
    end
    
end

[u,cc] = parseOutputs(u,cc,dm,padSize);

disp(['Convergence at iteration ',num2str(i)]);
disp(['Title time: ',num2str(toc(t0))]);
end



%% ========================================================================
function varargout = parseImages(varargin)
% pads images and creates meshgrid

I{1} = single(varargin{1}{1});
I{2} = single(varargin{1}{2});
sSize = varargin{2};
sSpacing = varargin{3};


prePad = sSize/2;
postPad = sSize/2;

sizeI = size(I{1});
I{1} = padarray(I{1},prePad,0,'pre');
I{1} = padarray(I{1},postPad,0,'post');

I{2} = padarray(I{2},prePad,0,'pre');
I{2} = padarray(I{2},postPad,0,'post');


idx = cell(1,3);
for i = 1:3, idx{i} = (1:sSpacing(i):(sizeI(i) + 1)) + sSize(i)/2; end

% [m{1},m{2},m{3}] = meshgrid(idx{:});

varargout{    1} = I;
varargout{end+1} = idx;

end

%% ========================================================================
function varargout = parseInputs(varargin)
% parses inputs and pads images so that there is an divisable meshgrid number.

I0{1} = single(varargin{1}{1});
I0{2} = single(varargin{1}{2});

% I0{1} = permute(I0{1},[2 1 3]);
% I0{2} = permute(I0{2},[2 1 3]);

sSize = varargin{2}; 
sSize = [sSize(2), sSize(1), sSize(3)];

sSpacing = sSize/2; 
u0 = varargin{3};

DVCPadSize = sSpacing/2;

sizeI0 = size(I0{1});
sizeI = ceil(sizeI0./sSpacing).*sSpacing;
prePad = ceil((sizeI - sizeI0)/2);
postPad = floor((sizeI - sizeI0)/2);

I{1} = padarray(I0{1},prePad,0,'pre');
I{1} = padarray(I{1},postPad,0,'post');

I{2} = padarray(I0{2},prePad,0,'pre');
I{2} = padarray(I{2},postPad,0,'post');

varargout{    1} = I;
varargout{end+1} = sSize;
varargout{end+1} = sSpacing;
varargout{end+1} = [prePad; postPad];
varargout{end+1} = DVCPadSize;
varargout{end+1} = u0;
end


function [u,cc] = parseOutputs(u,cc,filterSpacing,padSize)
% parses outputs and unpads the displacment field and cc.

unpadSize(1,:) = ceil(padSize(1,:)/filterSpacing);
unpadSize(2,:) = floor((padSize(2,:)+1)/filterSpacing); 
% +1 from the extra meshgrid point during the meshing of the DVC algorithm. EBK (10-23-2013)

for i = 1:3
    u{i} = u{i}(1+unpadSize(1,1):end-unpadSize(2,1),...
        1+unpadSize(1,2):end-unpadSize(2,2),...
        1+unpadSize(1,3):end-unpadSize(2,3));
end
u{4} = sqrt(u{1}.^2 + u{2}.^2 + u{3}.^2);

cc = cc(1+unpadSize(1,1):end-unpadSize(2,1),...
    1+unpadSize(1,2):end-unpadSize(2,2),...
    1+unpadSize(1,3):end-unpadSize(2,3));

end
