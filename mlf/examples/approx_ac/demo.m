clearvars; close all; clc;
set(groot,'DefaultFigurePosition', [100 150 1300 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',14)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for iloe1 = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(iloe1)},'factory','default'),'latex'); end
%
addpath('/Users/charles/Documents/GIT/mLF')
%
SAVEIT = false;
% Load frequency evaluations of the CFD model
load('dataONERA_FlexibleAircraft')
WW          = union(W,logspace(log10(min(W)),log10(max(W)),200));
subSetOut   = [23:2:31 33:44];
wSpace      = 1:length(W)-200;
%H           = H(subSetOut,:,:);
H           = H(subSetOut,:,wSpace); W = W(wSpace);
% clear H
% Hf = @(i1,i2) freqresp(tf(i2,[.1*i2/100 .2/10 1]),i1);
% for i2 = 1:10
%     for i1 = 1:numel(W)
%         H(i2,1,i1) = Hf(W(i1),i2);
%         tab(i1,i2) = H(i2,1,i1);
%     end
% end
for i2 = 1:size(H,1)
    for i1 = 1:size(H,3)
        tab(i1,i2) = H(i2,1,i1);
    end
end

% Plot the frequency response 
ny = size(H,1);
figure, hold on
for ii = 1:ny
    mdspack.bodemag(H(ii,:,:),W,'.')
end
title('CFD data')
if SAVEIT; drawnow, pause(.5), figSaveJPG('bodemag1',.5); end
% figure, hold on
% for ii = 1:ny
%     mdspack.bodephase(H(ii,:,:),W,'.')
% end
% title('CFD data')
% if SAVEIT; drawnow, pause(.5), figSaveJPG('bodephase1',.5); end
%%
%%% From the data identify a parametric model 
% location = 1:size(H,1);%linspace(1,2,size(H,1));
% for i2 = 1:size(H,1)
%     for i1 = 1:size(H,3)
%         tab(i1,i2) = H(i2,1,i1);
%         %tab_(i1,i2) = Hf(W(i1),i2);
%     end
% end
ip = {1i*W(:); 1:ny};
tab = tab([2:2:end 1:2:end], ...
          [2:2:end 1:2:end]);
for i = 1:length(ip)
    p_c{i} = ip{i}(2:2:end);
    p_r{i} = ip{i}(1:2:end);
    %p_c{i} = ip{i}(1:floor(length(ip{i})/2));
    %p_r{i} = ip{i}(1+floor(length(ip{i})/2):end);
end

opt.ord_tol     = 1e-9;
opt.ord_show    = true;
opt.method      = 'full';
[g,iloe]        = mlf.alg1(tab,p_c,p_r,opt);

%%
%%% Along first and second variables 
x1      = WW+rand(1)/100;%imag(linspace(min(ip{1}),max(ip{1}),40)+rand(1)/10);
x2      = (1:ny)+rand(1)/100;%linspace(min(ip{2}),max(ip{2}),41)+rand(1)/100;
[X,Y]   = meshgrid(x1,x2);
for i1 = 1:numel(x1)
    for i2 = 1:numel(x2)
        param           = [1i*x1(i1) x2(i2)];
        %Gref(i2,1,i1)   = Hf(x1(i1),x2(i2));
        Gr2(i2,1,i1)    = g(param);
    end
end
%
figure, hold on
for ii = 1:ny
    mdspack.bodemag(H(ii,:,:),W,'.')
    %mdspack.bodemag(Gref(ii,:,:),WW,'-')
    mdspack.bodemag(Gr2(ii,:,:),WW,'k--')
end
% title(['Approximation with $r=' num2str(length(Hr2.a)) '$'])
% if SAVEIT; drawnow, pause(.5), figSaveJPG('bodemag3',.5); end
% figure, hold on
% for ii = 1:ny
%     mdspack.bodephase(H(ii,:,:),W,'.')
%     mdspack.bodephase(Gr2(ii,:,:),WW,'k--')
% end
% title(['Approximation with $r=' num2str(length(Hr2.a)) '$'])
% if SAVEIT; drawnow, pause(.5), figSaveJPG('bodephase3',.5); end
% 
% 
% %%
% % Plot models and errors for the last case, r = 100
% 
% %%% 
% r = length(Hr2.a);
% funPlotGLRAfreqWing(W,H,Hr2,['anim_r' num2str(r)])
