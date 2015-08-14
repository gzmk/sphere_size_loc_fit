%fitting brdf to images

% Idea: I might be able to store the params in an array and write them into
% the conditions file at each iteration - take a look at this
function costIm = renderIm_3params_locsize(var)

% check rho_s+rho_d <1, if this is bigger than 1, renderer will scale both.
% We don't want the renderer to scale as that messes up fminsearch
% parameters. 

% persistent past_params; 
% persistent cost_arr;
% if isempty(past_params) 
%     past_params = []; 
% end
% if isempty(cost_arr)
%     cost_arr = [];
% end

% print the new values of parameters for every fminsearch iteration
% sprintf('The new variables are: ro_s: %f ro_d: %f alphau: %f', var(1), var(2), var(3))
% sprintf('The new variables are: ro_s: %f ro_d: %f', var(1), var(2))
% sprintf('The new variables are: alphau: %f', var(1))


% % empty = isempty(past_params);
% if ~isempty(past_params)
%     sprintf('BBBBBBB')
% %     idx = find(ismember(past_params,var','rows'))
%     idx = find(ismember(past_params,var,'rows'))
%     if ~isempty(idx)
%         costIm = cost_arr(idx(1));
%         return; 
%     end
% end

%% to write new conditions file with replaced params
% write to file in tabular form
% var = XBest; % this is to re-render the best fits
% var = [0.0760; 0.2168; 0.0472]; % this is for test

% THIS IS FOR MONOCHROMATIC RENDERING
fixedalpha = getGlobalalpha;


% sprintf('Fixed alpha is: %f', fixedalpha)
% ro_s = var(1)/(var(1)+var(2));
% ro_d = var(2)/(var(1)+var(2));
% % alphau = var(3);
% alphau = fixedalpha;
% light = (var(1)+var(2));

% rho_s = 0.04;
% rho_d = 0.1;
% fixedalpha = 0.2;
% ro_s = ['300:',num2str(rho_s),' 800:',num2str(rho_s)];
% ro_d = ['300:', num2str(rho_d), ' 800:', num2str(rho_d)];

ro_s = ['300:',num2str(getGlobalros),' 800:',num2str(getGlobalros)];
ro_d = ['300:', num2str(getGlobalrod), ' 800:', num2str(getGlobalrod)];
alphau = fixedalpha; % alphau and alphav should always be the same value for isotropic brdf
locx = var(1);
locy = var(2);
scalex = var(3);
mycell = {ro_s, ro_d, alphau, locx, locy, scalex};


T = cell2table(mycell, 'VariableNames', {'ro_s' 'ro_d' 'alphau' 'locx' 'locy' 'scalex'});
writetable(T,'/scratch/gk925/sphere_size_loc_fit/sphere_3params_Conditions.txt','Delimiter','\t')
% writetable(T,'/Local/Users/gizem/Documents/Research/GlossBump/sphere_size_loc_fit/sphere_3params_Conditions.txt','Delimiter','\t')

%% Rendering bit

% Set preferences
setpref('RenderToolbox3', 'workingFolder', '/scratch/gk925/sphere_size_loc_fit');
% setpref('RenderToolbox3', 'workingFolder', '/Local/Users/gizem/Documents/Research/GlossBump/sphere_size_loc_fit/');

% use this scene and condition file. 
parentSceneFile = 'test_sphere.dae';
% WriteDefaultMappingsFile(parentSceneFile); % After this step you need to edit the mappings file

conditionsFile = 'sphere_3params_Conditions.txt';
mappingsFile = 'sphere_3params_DefaultMappings.txt';

% Make sure all illuminants are added to the path. 
addpath(genpath(pwd))

% Choose batch renderer options.

% hints.imageWidth = 4012;
% hints.imageHeight = 6034;
% hints.imageWidth = 600;% these are for quick rendering
% hints.imageHeight = 800;
hints.imageWidth = 668;% this is isotropic scaling (orig. size divided by 4)
hints.imageHeight = 1005;
hints.renderer = 'Mitsuba';

datetime=datestr(now);
datetime=strrep(datetime,':','_'); %Replace colon with underscore
datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore
%hints.recipeName = ['Test-SphereFit' datetime];
hints.recipeName = ['Test-SphereFit' date];

ChangeToWorkingFolder(hints);

% nativeSceneFiles = MakeSceneFiles(parentSceneFile, '', mappingsFile, hints);
nativeSceneFiles = MakeSceneFiles(parentSceneFile, conditionsFile, mappingsFile, hints);
radianceDataFiles = BatchRender(nativeSceneFiles, hints);

%comment all this out
toneMapFactor = 10;
isScale = true;
montageName = sprintf('%s (%s)', 'Test-SphereFit', hints.renderer);
montageFile = [montageName '.png'];
[SRGBMontage, XYZMontage] = ...
    MakeMontage(radianceDataFiles, montageFile, toneMapFactor, isScale, hints);

% load the monochromatic image and display it
% imPath = ['/Local/Users/gizem/Documents/Research/GlossBump/sphere_size_loc_fit/', hints.recipeName, '/renderings/Mitsuba/test_sphere-001.mat']
imPath = ['/scratch/gk925/sphere_size_loc_fit/', hints.recipeName, '/renderings/Mitsuba/test_sphere-001.mat']
load(imPath, 'multispectralImage');
im2 = multispectralImage;
% figure;imshow(im2(:,:,1))

%% calculate the ssd (error) between two images
% dcraw command: -4 -d -v -w -b 3.0 DSC_0111_70gloss.pgm
% -b 3.0 makes it 3 times brighter
% gloss40 = imread('registered_photo.pgm','pgm');
% gloss = imread('registered40.pgm','pgm'); % turn this into a variable

% prepare a mask image for %40
mask = zeros(1005,668);
mask(382:670,256:444)=1;

% here instead of the registered image, we should start using the original pgm
pgm_name = ['40gloss.pgm'];
glossIm = imread(pgm_name,'pgm');
% for 0percent we need to rotate and then flip the image
rescaledIm = double(glossIm)/65535;
photo = imresize(rescaledIm, [1005,668]);

% load('registered40.mat') % make this a variable
% photo = renderRegisteredAdjusted;
masked_photo = mask.*photo;

% black = imread('DSC_0112.pgm')';
% imblack = imresize(black, [1005,668]);
% imblack2 = double(imblack)/65535;
% image1 = photo-imblack2;

renderedIm = im2(:,:,1); %for multispectral rendering
% renderedIm = im2;


diff = masked_photo-renderedIm;
costIm = sum(sum(diff.^2));

% cost_arr = [cost_arr;costIm];
% % past_params = [past_params;var'];
% past_params = [past_params;var]; % this for grid search as it takes in row arrays

return;



