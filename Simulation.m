clear;

load(['\\unimaas.nl\users\Students\i6144836\data\My Documents\MATLAB\Results TESTABLE' '\' 'results_tab.mat']);


%Simulation

% load(['\\unimaas.nl\users\Students\i6144836\data\My Documents\MATLAB\Results TESTABLE\Final' '\' 'p_basel.mat']);
% load(['\\unimaas.nl\users\Students\i6144836\data\My Documents\MATLAB\Results TESTABLE\Final' '\' 'p_honest.mat']);
% load(['\\unimaas.nl\users\Students\i6144836\data\My Documents\MATLAB\Results TESTABLE\Final' '\' 'p_dishonest.mat']);
load(['\\unimaas.nl\users\Students\i6144836\data\My Documents\MATLAB\Results TESTABLE\Final' '\' 'parameters.mat']);

exclu = []
  for ii = 1:82
        if ~any(results_tab.lie(results_tab.subj_num == ii))
            results_tab(results_tab.subj_num == ii, :) = [];
            n = n - 1;
            exclu = [exclu, ii];
        elseif all(results_tab.lie(results_tab.subj_num == ii))
            results_tab(results_tab.subj_num == ii, :) = [];
            n = n - 1;
            exclu = [exclu, ii];
        end
    end
%%
dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 4,'n_t',50);
    
    u = [results_tab.soloTrueV(1:50,:)'; results_tab.soloFalseV(1:50,:)'; ...
        results_tab.soloTrueV(51:100,:)'; results_tab.soloFalseV(51:100,:)';...
        results_tab.soloTrueV(101:150,:)'; results_tab.soloFalseV(101:150,:)'];
    
        g_fname = @gModels2;
        options.sources = struct('type',1 ,'out', 1);
        options.inG = 'conformity_sim';
        options.DisplayWin = 0; options.verbose = 0;
        options.dim = dim;

final = []; 
subj = [];
basel = [];
honest = [];
dishonest = [];

param(:,1) = exp(param(:,1));

for ii = 1:70
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + 149;
    data = results_tab(nrows_init:nrows_final,:);
     u = [data.soloTrueV(1:50,:)'; data.soloFalseV(1:50,:)'; ...
        data.soloTrueV(51:100,:)'; data.soloFalseV(51:100,:)';...
        data.soloTrueV(101:150,:)'; data.soloFalseV(101:150,:)'];
    
  for block = 1:3
      switch block
          case 1
              p = param(:,1:4);
          case 2
              p = [param(:,1:3), param(:,5)];
          case 3
              p = [param(:,1:3), param(:,6)];   
      end
      
    [y] = VBA_simulate(1,[],g_fname,[], p(ii,:), u((1+2*(block-1)):(2+2*(block-1)),:), [],[],options);
    
    subj = [subj; y' ,u((1+2*(block-1)):(2+2*(block-1)),:)'];
    
    switch block
        case 1
            basel = [basel; y' ,u(1:2,:)'];
        case 2
            honest = [honest; y' ,u(3:4,:)'];
        case 3
             dishonest = [dishonest; y' ,u(5:6,:)'];
    end 
  end
    final = [final; subj];
    subj = [];
end
writematrix(basel,'sim_basel.csv');
writematrix(honest, 'sim_honest.csv');
writematrix(dishonest,'sim_dishonest.csv');
%%
figure
tab = [final(:,2), final(:,3), final(:,1)];
tab = array2table(tab, 'VariableNames', {'soloTrueV', 'soloFalseV', 'lie'});

c = heatmap(tab, 'soloTrueV', 'soloFalseV', 'ColorVariable', 3);
c.Title = 'Percentage of lying per parameters';
c.XLabel = 'True draw';
c.YLabel = 'Other draw';
caxis([0,0.5]);

%%
figure
sgtitle('Heatmap of overall lying, per parameters and block');

tabb.basel = [basel(:,2), basel(:,3), basel(:,1)];
tabb.basel = array2table(tabb.basel, 'VariableNames', {'soloTrueV', 'soloFalseV', 'lie'});
tabb.goodb =  [honest(:,2), honest(:,3), honest(:,1)];
tabb.goodb = array2table(tabb.goodb, 'VariableNames', {'soloTrueV', 'soloFalseV', 'lie'});
tabb.badb = [dishonest(:,2), dishonest(:,3), dishonest(:,1)];
tabb.badb = array2table(tabb.badb, 'VariableNames', {'soloTrueV', 'soloFalseV', 'lie'});

subplot(1,3,1)
c = heatmap(tabb.basel, 'soloTrueV', 'soloFalseV', 'ColorVariable', 3);
c.Title = 'Baseline';
c.XLabel = 'True draw';
c.YLabel = 'Other draw';
caxis([0,0.52]);

subplot(1,3,2)
c = heatmap(tabb.goodb, 'soloTrueV', 'soloFalseV', 'ColorVariable', 3);
c.Title = 'Honest group';
c.XLabel = 'True draw';
c.YLabel = 'Other draw';
caxis([0,0.52]);

subplot(1,3,3)
c = heatmap(tabb.badb, 'soloTrueV', 'soloFalseV', 'ColorVariable', 3);
c.Title = 'Dishonest group';
c.XLabel = 'True draw';
c.YLabel = 'Other draw';
caxis([0,0.52]);

%% DIFFERENCE

figure
tab = [final(:,2), final(:,3), final(:,1)-results_tab.lie];
tab = array2table(tab, 'VariableNames', {'soloTrueV', 'soloFalseV', 'lie'});

c = heatmap(tab, 'soloTrueV', 'soloFalseV', 'ColorVariable', 3);
c.Title = 'Percentage of lying per parameters';
c.XLabel = 'True draw';
c.YLabel = 'Other draw';
caxis([-0.2,0.2]);

