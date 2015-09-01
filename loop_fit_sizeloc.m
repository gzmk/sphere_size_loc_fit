%% Loop through two fits 
% Step 1: Start from alphau = 0.2, fit ro_s and ro_d
% Step 2: Get the XBest from above fit and fix them, fit only alphau
function loop_fit_sizeloc()

% iter = 5;
setGlobalalpha(0.2);
setGlobalrod(0.1);
setGlobalros(0.04); %40 gloss

% init 2 param fitting
LB = [-0.03; -0.03; 2.0];
UB = [0.03; 0.03; 2.4];
NumDiv = [10 10 10];
MinDeltaX = [1e-5 1e-5 1e-5];

bestfit_3pr = [];

bestLocSize = [];

fitname1 = '40percent_sizeloc.mat';

[LocSizeBest, BestF2] = gridsearch(@renderIm_3params_locsize, LB, UB, 10, 0.8, 1e-7, 1000, 1, 1);
% [LocSizeBest,BestF2,Iters2] = Grid_Search(3, LB', UB', NumDiv, MinDeltaX, 1e-7, 1000, 'renderIm_3params_locsize');
sprintf('This is LocSizeBest:');
LocSizeBest;
% setGloballightx(LocSizeBest(1))
% setGloballighty(LocSizeBest(2))
bestLocSize = [bestLocSize;LocSizeBest];
bestfit_3pr = [bestfit_3pr;BestF2];
% delta1 = bestLight(i,1)*0.01;
% delta2 = bestLight(i,2)*0.01;
% 
% if i>1
%     converge1 = abs(bestLight(i,1)-bestLight(i-1,1));
%     converge2 = abs(bestLight(i,2)-bestLight(i-1,2));
% end

imname = strcat('/scratch/gk925/sphere_size_loc_fit/fit_results/multispectral/', fitname1);
save(imname, 'bestLocSize','bestfit_3pr');

sprintf('Done!')

% if i>1 && (converge1<=delta1) && (converge2<=delta2)
%     sprintf('Parameters converged! Exiting the loop')
%     break
% end

return;
