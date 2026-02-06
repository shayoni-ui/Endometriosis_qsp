close all;
clear all;

%[success, variant_out, model_out, exitInfo]= sbiosteadystate(ans);
m = createEndometriosisModel2();

% d1 = adddose(m,'testdose');
% d1.AmountUnits = 'nanomole';
% d1.TargetName = 'IGF1';
% d1.TimeUnits = 'hour';
% d1.StartTime = 40;
% d1.Amount = 1000;

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
m.addparameter('h',1, 'Units','hour');
 m.addparameter('DeltaNO',50, 'Units','nanomolarity');
 on1 = addevent(m,'time/h>500','NO=NO+DeltaNO');
 off1 = addevent(m,'time/h>1000','NO=NO-DeltaNO');
%---------------------------------
m.addparameter('DeltaIGF1',50, 'Units','nanomolarity');
on = addevent(m,'time/h>500','IGF1=IGF1+DeltaIGF1');
off = addevent(m,'time/h>1000','IGF1=IGF1-DeltaIGF1');
% %---------------------------------
%---------------------------------
m.addparameter('DeltaEP3',50, 'Units','nanomolarity');
on = addevent(m,'time/h>500','EP3=EP3+DeltaEP3');
off = addevent(m,'time/h>1000','EP3=EP3-DeltaEP3');
% %---------------------------------
dose = 100;
for jj = 1:length(dose)

IGF1 = sbioselect(m, 'Type', 'species', 'Name', 'IGF1');IGF1.BoundaryCondition = true;
% IGF1.InitialAmount = dose(jj);
% argi = sbioselect(m, 'Type', 'species', 'Name', 'Arginine_i');argi.BoundaryCondition = true;
% argi.InitialAmount = dose(jj);
NO = sbioselect(m, 'Type', 'species', 'Name', 'NO');NO.BoundaryCondition = true;
% NO.InitialAmount = 50;
EP3 = sbioselect(m, 'Type', 'species', 'Name', 'EP3');EP3.BoundaryCondition = true;
% EP3.InitialAmount = 10;
 
sd(jj) = sbiosimulate(m,cs);
%sdnew(jj) = resample(sd(jj),[1:16],'linear');
end


figure;
tl=tiledlayout("flow",'TileSpacing','compact','Padding','compact');
% Add layout title
title(tl,'Endometriosis plots');

%plotvars = {'IGF1','IGFBPc','PAPPA2','IGFBP','PAPPA','Arginine_i','NO','PGE2','EP3','EC','Ca','JFlux','PainFlux','PainF'};
plotvars = {'IGF1','NO','EP3','EC','Ca','JFlux','Nchannels','Kaff','PainFlux'};
plotvars2 = {'IGF1','NO','EP3','EC','Ca','JFlux','Nchannels','Kaff','PainFlux'};
for i=1:length(plotvars)
    nexttile;
    [t,x]=selectbyname(R,plotvars{i}); 
    plot(t,x);
    hold on;
    [t1,x1]=selectbyname(sd,plotvars2{i}); 
    plot(t1,x1);
    title(plotvars{i}); 
    xlabel(sprintf('Time(%s)',cs.TimeUnits));
    ylabel('Conc [nM]');
end
