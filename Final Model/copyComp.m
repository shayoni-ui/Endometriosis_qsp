function m = copyEDmodel(m1,m2)


m  = copyobj(m1);
mm = copyobj(m2);
rename(mm.Compartments, 'ED1');

% Query object for all species
sp  = m.species;
sp2 = mm.Species;
%save reaction ans species in variable to parent m for deleting at later



% %add species to the model of interest by contacatenating comp.species
for i=1:length(sp),spnm{i} = sp(i).Name; end 
     for j=1:length(sp2)
         spnm2{j} = sp2(j).Name; 
     end 
 %spnm = unique(spnm);
 spnmOrg = strcat(m.Compartments.Name,'.',spnm);
 spnmOrg2 = strcat(m.Compartments.Name,'.',spnm2);
 %spnmOrg = unique(spnmOrg);
% spnm1 = cell(1,6);
% 
 for i=1:length(m.Compartments)
 
     spnm1(i,:)= strcat(m.Compartments(i).Name,'.',spnm);
     spnm2(i,:)= strcat(m.Compartments(i).Name,'.',spnm2);
 end

% % copy all species in m2 to m1 now m
 for k=1:length(m.Compartments)
     for j=1:length(sp2)-3
         copyobj(sp2(j), m.Compartments(k))
 
     end

 end
% % copy all  reaction in m1
 rx = get(mm,'Reactions');
 rx2 = cell(1,length(rx));
% 

 for i=1:length(rx)

% %for j=1:length(rxnew)
%     for i=1:6
%     %for i=1:(length(m1.Compartments)*3)
         rx2{i}= wordrep(rx(i).Reaction,spnm2(1,:),spnmOrg2(1,:));       
         rx4(i) = addreaction(m,rx2{i});
         rx4(i).ReactionRate = wordrep(rx(i).ReactionRate,spnm2(1,:),spnmOrg2(1,:));
% 
%     end
%    
end
% 
% %% delete species and reaction from m1
% %rmreactant(ans, 'all');
% delete(Mrx2del);
% delete(Msp2del)
% 
% 
