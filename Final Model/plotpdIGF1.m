close all;
clear all;

%[success, variant_out, model_out, exitInfo]= sbiosteadystate(ans);
m = createPAPPA2modelt2();

cs = m.getconfigset;
cs.SolverOptions.AbsoluteTolerance = 1e-8;
cs.SolverOptions.RelativeTolerance = 1e-8;
set(cs,'StopTime',24*100);
set(cs,'TimeUnits','hour');

% this is only for debug purposes
% make all parameter variable so you would get their value as part of the
% simulation result
p=sbioselect(m,'type','parameter');
for i=1:length(p) p(i).Constant=0; end

R = sbiosimulate(m,cs);

figure;
tl=tiledlayout("flow",'TileSpacing','compact','Padding','compact');
% Add layout title
title(tl,'IGF1 plots');

plotvars = {'IGF1','IGFBPc','PAPPA2','IGFBP','PAPPA','totIGFBP','mAb','PGE2','EC'};

for i=1:length(plotvars)
    nexttile;
    [t,x]=selectbyname(R,plotvars{i}); 
    plot(t,x);
    title(plotvars{i}); 
    xlabel(sprintf('Time(%s)',cs.TimeUnits));
    ylabel('Conc [nM]');
end
