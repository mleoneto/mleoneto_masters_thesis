%% General code

clear;
load(['\\unimaas.nl\users\Students\i6144836\data\My Documents\MATLAB\Results TESTABLE' '\' 'results_tab.mat']);
%load('C:\Users\roolio\Dropbox\Research\Moral contagion\Results TESTABLE\results_tab.mat');

path = '\\unimaas.nl\users\Students\i6144836\data\My Documents\MATLAB\Results TESTABLE';
%path='C:\Users\roolio\Dropbox\Research\Moral contagion\Results TESTABLE';
files = dir(fullfile(path,'Results_Subj_*.csv'));

n_tot = length(files); % n subjects
ntrials = 150;
nmodels = 26; %n models being tested


every_participant = 1; %1 removes participants that always lied or always told the truth. 
exclu = [];n=n_tot;
if every_participant == 1
    for ii = 1:n_tot
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
end

fprintf('%d subjects excluded out of %d',length(exclu),n);

%L = NaN(nmodels, n);
%load(['\\unimaas.nl\users\Students\i6144836\data\My Documents\MATLAB\Results TESTABLE\Julien 2' '\' 'L_70.mat']);
model_count = 1;

%% Bias Model
% 
% %Initializing mu and sigma matrices
% 
% for ii = 1:n %loop that runs through n participants
%     fprintf('Processing subject No %d of %d... for Bias Model (%d/%d) \n',ii,n,model_count,nmodels);
%     
%     % Defining number of trials we're looking at, starting number of row and
%     % final number of row.
%     nrows_init = 1 + 150*(ii-1);
%     nrows_final = nrows_init + ntrials - 1;
%     
%     dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 2,'n_t',ntrials);
%     
%     y = results_tab.lie(nrows_init:nrows_final)';
%     u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
%     
%     g_fname = @gModels2;
%     
%     %Defining priors and options
%     options.sources = struct('type',1 ,'out', 1);
%     options.inG = 'bias';
%     options.DisplayWin = 0; options.verbose = 0;
%     
%     %Inversion/estimation
%   %  [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
%     
% %     mu_final = [mu_final;  posterior.muPhi'];
% %     sigma_final = [sigma_final; posterior.SigmaPhi];
%   %  L(model_count,ii) = model(model_count,ii).out.F;
%   
%   [bias(ii).posterior, bias(ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
%   L_bias(ii) = bias(ii).out.F;
%     clearvars -except results_tab model_count L ii n ntrials nmodels model exclu bias L_bias
% end
%     model_count = model_count + 1;
% 
% %% Bias Multisession Model
% 
% for ii = 1:n %loop that runs through n participants
%     fprintf('Processing subject No %d of %d... for Bias Multi Model (%d/%d) \n',ii,n,model_count,nmodels);
%     
%     % Defining number of trials we're looking at, starting number of row and
%     % final number of row.
%     nrows_init = 1 + 150*(ii-1);
%     nrows_final = nrows_init + ntrials - 1;
%     
%     dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 2,'n_t',ntrials);
%     
%     y = results_tab.lie(nrows_init:nrows_final)';
%     u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
%     
%     g_fname = @gModels2;
%     
%     %Defining priors and options
%     options.sources = struct('type',1 ,'out', 1);
%     options.inG = 'bias';
%     options.DisplayWin = 0; options.verbose = 0;
%     options.multisession.split = [50 50 50];
%     options.multisession.fixed.phi = [2]; %Temperature fix
%     
%     %Inversion/estimation
%     [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
% 
%     L(model_count,ii) = model(model_count,ii).out.F;
%     clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
% end
%     model_count = model_count + 1;
% 

%% Basic Model

%Initializing mu and sigma matrices

for ii = 1:n %loop that runs through n participants
    fprintf('Processing subject No %d of %d... for Basic Model (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 3,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    g_fname = @gModels2;
    
    %Defining priors and options
    options.sources = struct('type',1 ,'out', 1);
    options.inG = 'basic';
    options.DisplayWin = 0; options.verbose = 0;
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
    
%     mu_final = [mu_final;  posterior.muPhi'];
%     sigma_final = [sigma_final; posterior.SigmaPhi];
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
end
    model_count = model_count + 1;

%% Basic Multisession Model

for ii = 1:n %loop that runs through n participants
    fprintf('Processing subject No %d of %d... for Basic Multi Model (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 3,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    g_fname = @gModels2;
    
    %Defining priors and options
    options.sources = struct('type',1 ,'out', 1);
    options.inG = 'basic';
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [3]; %Temperature fix
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);

    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
end
    model_count = model_count + 1;
    
%% Basic Multisession alpha Model

for ii = 1:n %loop that runs through n participants
    fprintf('Processing subject No %d of %d... for Basic Multi Alpha Model (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 3,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    g_fname = @gModels2;
    
    %Defining priors and options
    options.sources = struct('type',1 ,'out', 1);
    options.inG = 'basic';
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [2 3]; %Temperature fix
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);

    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
end
    model_count = model_count + 1;

    %% Basic Multisession delta Model

for ii = 1:n %loop that runs through n participants
    fprintf('Processing subject No %d of %d... for Basic Multi delta Model (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 3,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    g_fname = @gModels2;
    
    %Defining priors and options
    options.sources = struct('type',1 ,'out', 1);
    options.inG = 'basic';
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 3]; %Temperature fix
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);

    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
end
    model_count = model_count + 1;


%% Fixed Cost  Model

for ii = 1:n %loop that runs through n participants
    fprintf('Processing subject No %d of %d... for Fixed Cost Model (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 3,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    g_fname = @gModels2;
    
    %Defining priors and options
    options.sources = struct('type',1 ,'out', 1);
    options.inG = 'fixed_cost';
    options.DisplayWin = 0; options.verbose = 0;
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
    
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
end
    model_count = model_count + 1;

%% Fixed Cost Model Multisession

for ii = 1:n %loop that runs through n participants
    fprintf('Processing subject No %d of %d... for Fixed Cost Model Multisession (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 3,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    % block = results_tab.block(nrows_init:nrows_final)';
    % bad_group = results_tab.bad_group(nrows_init:nrows_final)';
    
    u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    g_fname = @gModels2;
    
    %Defining priors and options
    options.sources = struct('type',1 ,'out', 1);
    options.inG = 'fixed_cost';
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [3]; %Temperature fix
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
    
end
    model_count = model_count + 1;

%% Fixed Cost alpha Model Multisession

for ii = 1:n %loop that runs through n participants
    fprintf('Processing subject No %d of %d... for Fixed Cost alpha Model Multisession (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 3,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    % block = results_tab.block(nrows_init:nrows_final)';
    % bad_group = results_tab.bad_group(nrows_init:nrows_final)';
    
    u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    g_fname = @gModels2;
    
    %Defining priors and options
    options.sources = struct('type',1 ,'out', 1);
    options.inG = 'fixed_cost';
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [2 3]; %Temperature fix
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
    
end
    model_count = model_count + 1;
 
 %% Fixed Cost delta Model Multisession

for ii = 1:n %loop that runs through n participants
    fprintf('Processing subject No %d of %d... for Fixed Cost delta Model Multisession (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 3,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    % block = results_tab.block(nrows_init:nrows_final)';
    % bad_group = results_tab.bad_group(nrows_init:nrows_final)';
    
    u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    g_fname = @gModels2;
    
    %Defining priors and options
    options.sources = struct('type',1 ,'out', 1);
    options.inG = 'fixed_cost';
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 3]; %Temperature fix
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
    
end
    model_count = model_count + 1;
    
 %% basic conformity
    
for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Basic Conformity (%d/%d) \n',ii,n,model_count,nmodels);
% Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 4,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    basel = [ones(1,50) zeros(1,50)  zeros(1,50)];
    u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'; basel];
    
    g_fname = @gModels2;
    
    %Defining priors and options
    options.sources = struct('type',1 ,'out', 1);
    options.inG = 'conformity';
    options.DisplayWin = 0; options.verbose = 0;
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
    
%     mu_final = [mu_final;  posterior.muPhi'];
%     sigma_final = [sigma_final; posterior.SigmaPhi];
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
    
end
    model_count = model_count + 1;

%% basic conformity multisession

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Conformity Multi Model (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 0,'p',1, 'n_theta', 0, 'n_phi', 4,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    basel = [ones(1,50) zeros(1,50)  zeros(1,50)];
    u = [results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'; basel];
        
    g_fname = @gModels2;
    
    %Defining priors and options
    options.sources = struct('type',1 ,'out', 1);
    options.inG = 'conformity';
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 2 3]; %Temperature fix
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
%[extra(ii).posterior, extra(ii).out] = VBA_NLStateSpaceModel(y, u,[], g_fname, dim, options);
%extraL(ii) = extra(ii).out.F;
    L(model_count,ii) = model(model_count,ii).out.F;
    

    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu extra extraL
    
   
end
    model_count = model_count + 1;
%% Action RL alpha Model

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Action Alpha Model (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 2,'p',1, 'n_theta', 1, 'n_phi', 3,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    diff_liars = [nan(50,1) ; 10 - 2.*results_tab.pred_true(nrows_init+50:nrows_final)]; % n  of liars - truth-tellers
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    u = [diff_liars' ; bad_group' ; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    
    g_fname = @gModels2;
    f_fname = @fModels2;
    
    
    %Defining priors and options
    options.inG = 'action_alpha';
    options.inF = 'action_alpha';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
    options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;   
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
end
    model_count = model_count + 1;
    
%% Action RL alpha Model Multi

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Action Alpha Model Multi (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 2,'p',1, 'n_theta', 1, 'n_phi', 3,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    diff_liars = [nan(50,1) ; 10 - 2.*results_tab.pred_true(nrows_init+50:nrows_final)]; % n  of liars - truth-tellers
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    u = [diff_liars' ; bad_group' ; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    
    g_fname = @gModels2;
    f_fname = @fModels2;
    
    
    %Defining priors and options
    options.inG = 'action_alpha';
    options.inF = 'action_alpha';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
        options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 2 3]; %Temperature is fixed
        
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
end
    model_count = model_count + 1;    
%% Action RL delta Model

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Action Delta Model (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 2,'p',1, 'n_theta', 1, 'n_phi', 3,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    diff_liars = [nan(50,1) ; 10 - 2.*results_tab.pred_true(nrows_init+50:nrows_final)]; % n  of liars - truth-tellers
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    u = [diff_liars' ; bad_group' ; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    
    g_fname = @gModels2;
    f_fname = @fModels2;
    
    
    %Defining priors and options
    options.inG = 'action_delta';
    options.inF = 'action_delta';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
        options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [];
    
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
end
    model_count = model_count + 1;
    
%% Action RL delta Model multi

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Action Delta Model Multi (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 2,'p',1, 'n_theta', 1, 'n_phi', 3,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    diff_liars = [nan(50,1) ; 10 - 2.*results_tab.pred_true(nrows_init+50:nrows_final)]; % n  of liars - truth-tellers
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    u = [diff_liars' ; bad_group' ; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    
    g_fname = @gModels2;
    f_fname = @fModels2;
    
    
    %Defining priors and options
    options.inG = 'action_delta';
    options.inF = 'action_delta';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
        options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 2 3]; %Temperature is fixed
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu
end
    model_count = model_count + 1;    
%% Action alpha and delta Model

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Action alpha and delta Model (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 4,'p',1, 'n_theta', 2, 'n_phi', 3,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    diff_liars = [nan(50,1) ; 10 - 2.*results_tab.pred_true(nrows_init+50:nrows_final)]; % n  of liars - truth-tellers
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    u = [diff_liars' ; bad_group' ; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];

    g_fname = @gModels2;
    f_fname = @fModels2;
    
    %Defining priors and options
    options.inG = 'action_alpha_delta';
    options.inF = 'action_alpha_delta';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
        options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [] ;
    
    %Inversion/estimation
    [model(15,ii).posterior, model(15,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(15,ii) = model(15,ii).out.F;
    
   %  [extra(ii).posterior, extra(ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
  %  L_extra(ii) = extra(ii).out.F;
    
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu extra L_extra
end
    model_count = model_count + 1;
    
%% Action alpha and delta Model Multi

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Action alpha and delta Model Multi (%d/%d) \n',ii,n,model_count,nmodels);
    
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 4,'p',1, 'n_theta', 2, 'n_phi', 3,'n_t',ntrials);
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    diff_liars = [nan(50,1) ; 10 - 2.*results_tab.pred_true(nrows_init+50:nrows_final)]; % n  of liars - truth-tellers
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    u = [diff_liars' ; bad_group' ; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];

    g_fname = @gModels2;
    f_fname = @fModels2;
    
    %Defining priors and options
    options.inG = 'action_alpha_delta';
    options.inF = 'action_alpha_delta';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
        options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 2 3]; %Temperature is fixed
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
 %[model(14,ii).posterior, model(14,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
  %  L(14,ii) = model(14,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu 
end
    model_count = model_count + 1;  

%% Outcome alpha model

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Outcome alpha Model (%d/%d) \n',ii,n,model_count,nmodels);
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 2,'p',1, 'n_theta', 1, 'n_phi', 3,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    truthtellers = results_tab.pred_true(nrows_init:nrows_final);
    predTrueV = results_tab.predTrueV(nrows_init:nrows_final);
    predFalseV = results_tab.predFalseV(nrows_init:nrows_final);
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    
    u = [truthtellers' ; predTrueV' ; predFalseV'; bad_group'; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    

    g_fname = @gModels2;
    f_fname = @fModels2;
    
    
    %Defining priors and options
    options.inG = 'outcome_alpha';
    options.inF = 'outcome_alpha';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
    options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model   
end
    model_count = model_count + 1;
    
%% Outcome alpha model Multi

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Outcome alpha Model multi (%d/%d) \n',ii,n,model_count,nmodels);
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 2,'p',1, 'n_theta', 1, 'n_phi', 3,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    truthtellers = results_tab.pred_true(nrows_init:nrows_final);
    predTrueV = results_tab.predTrueV(nrows_init:nrows_final);
    predFalseV = results_tab.predFalseV(nrows_init:nrows_final);
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    
    u = [truthtellers' ; predTrueV' ; predFalseV'; bad_group'; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    

    g_fname = @gModels2;
    f_fname = @fModels2;
    
    
    %Defining priors and options
    options.inG = 'outcome_alpha';
    options.inF = 'outcome_alpha';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
    options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 2 3]; %Temperature is fixed
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model   
end
    model_count = model_count + 1;
%% Outcome delta model

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Outcome delta Model (%d/%d) \n',ii,n,model_count,nmodels);
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 2,'p',1, 'n_theta', 1, 'n_phi', 3,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    truthtellers = results_tab.pred_true(nrows_init:nrows_final);
    predTrueV = results_tab.predTrueV(nrows_init:nrows_final);
    predFalseV = results_tab.predFalseV(nrows_init:nrows_final);
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    u = [truthtellers' ; predTrueV' ; predFalseV'; bad_group'; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    g_fname = @gModels2;
    f_fname = @fModels2;
    
    %Defining priors and options
    options.inG = 'outcome_delta';
    options.inF = 'outcome_delta';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
    options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model   
end
    model_count = model_count + 1;    
    
%% Outcome delta model multi

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Outcome delta Model multi (%d/%d) \n',ii,n,model_count,nmodels);
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 2,'p',1, 'n_theta', 1, 'n_phi', 3,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    truthtellers = results_tab.pred_true(nrows_init:nrows_final);
    predTrueV = results_tab.predTrueV(nrows_init:nrows_final);
    predFalseV = results_tab.predFalseV(nrows_init:nrows_final);
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    u = [truthtellers' ; predTrueV' ; predFalseV'; bad_group'; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    
    g_fname = @gModels2;
    f_fname = @fModels2;
    
    
    %Defining priors and options
    options.inG = 'outcome_delta';
    options.inF = 'outcome_delta';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
    options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 2 3]; %Temperature is fixed
    
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model   
end
    model_count = model_count + 1;        
%% Outcome alpha and delta Model

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Outcome alpha and delta Model (%d/%d) \n',ii,n,model_count,nmodels);
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 4,'p',1, 'n_theta', 2, 'n_phi', 3,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    truthtellers = results_tab.pred_true(nrows_init:nrows_final);
    predTrueV = results_tab.predTrueV(nrows_init:nrows_final);
    predFalseV = results_tab.predFalseV(nrows_init:nrows_final);
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    
    u = [truthtellers' ; predTrueV' ; predFalseV'; bad_group'; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    

    g_fname = @gModels2;
    f_fname = @fModels2;
    
    
    %Defining priors and options
    options.inG = 'outcome_alpha_delta';
    options.inF = 'outcome_alpha_delta';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
    options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model   
end
    model_count = model_count + 1;

%% Outcome alpha and delta Model multi

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Outcome alpha and delta Model Multi (%d/%d) \n',ii,n,model_count,nmodels);
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 4,'p',1, 'n_theta', 2, 'n_phi', 3,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    truthtellers = results_tab.pred_true(nrows_init:nrows_final);
    predTrueV = results_tab.predTrueV(nrows_init:nrows_final);
    predFalseV = results_tab.predFalseV(nrows_init:nrows_final);
    bad_group = [nan(50,1); zeros(50,1) ; ones(50,1)];
    
    u = [truthtellers' ; predTrueV' ; predFalseV'; bad_group'; results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'];
    

    g_fname = @gModels2;
    f_fname = @fModels2;
    
    
    %Defining priors and options
    options.inG = 'outcome_alpha_delta';
    options.inF = 'outcome_alpha_delta';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
        options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 2 3]; %Temperature is fixed
    
    %Inversion/estimation
    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
    L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model   
end
    model_count = model_count + 1;
    
%%  conformity action

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Conformity Action (%d/%d) \n',ii,n,model_count,nmodels);

    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 1,'p',1, 'n_theta', 1, 'n_phi', 4,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    diff_liars = [nan(50,1) ; 10 - 2.*results_tab.pred_true(nrows_init+50:nrows_final)];
   bad_group = [nan(50,1); ones(1,1); zeros(49,1) ; ones(1,1); zeros(49,1)];    
       basel = [ones(1,50) zeros(1,50)  zeros(1,50)];
    u = [diff_liars'; bad_group';results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'; basel];
       
    
    g_fname = @gModels2;
    f_fname = @fModels2;
    
    %Defining priors and options
    options.inG = 'conformity_action';
    options.inF = 'conformity_action';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
        options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    %Inversion/estimation
%     [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
%     L(model_count,ii) = model(model_count,ii).out.F;
    
    
  [ m_temp(1,ii).posterior, m_temp(1,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
 L_temp(1,ii) = m_temp(1,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu m_temp L_temp
end
    model_count = model_count + 1;

%  conformity action multisession

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Conformity Action multi (%d/%d) \n',ii,n,model_count,nmodels);
   
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 1,'p',1, 'n_theta', 1, 'n_phi', 4,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
   bad_group = [nan(50,1); ones(1,1); zeros(49,1) ; ones(1,1); zeros(49,1)];    
    diff_liars = [nan(50,1) ; 10 - 2.*results_tab.pred_true(nrows_init+50:nrows_final)];
        basel = [ones(1,50) zeros(1,50)  zeros(1,50)];
    u = [diff_liars';bad_group';results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'; basel];

    g_fname = @gModels2;
    f_fname = @fModels2;


    
    %Defining priors and options
    options.inG = 'conformity_action';
    options.inF = 'conformity_action';
    
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
        options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 2 3]; %Gamma and hidden states are estimated per block
    
    %Inversion/estimation
%  [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
%    L(model_count,ii) = model(model_count,ii).out.F;

   [ m_temp(2,ii).posterior, m_temp(2,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
 L_temp(2,ii) = m_temp(2,ii).out.F;

    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu m_temp L_temp
    
end
   model_count = model_count + 1;
 
 % conformity outcome

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Conformity outcome (%d/%d) \n',ii,n,model_count,nmodels);
   
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 1,'p',1, 'n_theta', 1, 'n_phi', 4,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    diff_liars = [nan(50,1) ; 10 - 2.*results_tab.pred_true(nrows_init+50:nrows_final)];
    predFalseV = results_tab.predFalseV(nrows_init:nrows_final);  
    predTrueV = results_tab.predTrueV(nrows_init:nrows_final);  
   bad_group = [nan(50,1);  ones(1,1); zeros(49,1) ; ones(1,1); zeros(49,1)];
       basel = [ones(1,50) zeros(1,50)  zeros(1,50)];
    u = [diff_liars' ;bad_group';predTrueV';predFalseV';results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'; basel];
    
    g_fname = @gModels2;
    f_fname = @fModels2;
    
    %Defining priors and options
    options.inG = 'conformity_outcome';
    options.inF = 'conformity_outcome';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
        options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    
         [ m_temp(3,ii).posterior, m_temp(3,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
 L_temp(3,ii) = m_temp(3,ii).out.F;
 
%     [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
%     L(model_count,ii) = model(model_count,ii).out.F;
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu L_temp m_temp
    

    
end
    model_count = model_count + 1;    
    
% basic conformity outcome multi

for ii = 1:n %loop that runs through n participants
    
    fprintf('Processing subject No %d of %d... for Conformity outcome multi (%d/%d) \n',ii,n,model_count,nmodels);
   
    % Defining number of trials we're looking at, starting number of row and
    % final number of row.
    nrows_init = 1 + 150*(ii-1);
    nrows_final = nrows_init + ntrials - 1;
    
    dim = struct('n', 1,'p',1, 'n_theta', 1, 'n_phi', 4,'n_t',ntrials); % Dimension
    
    y = results_tab.lie(nrows_init:nrows_final)';
    
    
    diff_liars = [nan(50,1) ; 10 - 2.*results_tab.pred_true(nrows_init+50:nrows_final)];
    predFalseV = results_tab.predFalseV(nrows_init:nrows_final);  
    predTrueV = results_tab.predTrueV(nrows_init:nrows_final);  
   bad_group = [nan(50,1); ones(1,1); zeros(49,1) ; ones(1,1); zeros(49,1)];
       basel = [ones(1,50) zeros(1,50)  zeros(1,50)];

    u = [diff_liars' ;bad_group';predTrueV';predFalseV';results_tab.soloTrueV(nrows_init:nrows_final)'; results_tab.soloFalseV(nrows_init:nrows_final)'; basel];
   
    
    g_fname = @gModels2;
    f_fname = @fModels2;
    
    %Defining priors and options
    options.inG = 'conformity_outcome';
    options.inF = 'conformity_outcome';
    options.priors.muPhi = zeros(dim.n_phi,1);
    options.priors.SigmaPhi = eye(dim.n_phi);
    options.sources = struct('type',1 ,'out', 1);
        options.skipf = [ones(1,50) zeros(1,ntrials-50)];
    options.DisplayWin = 0; options.verbose = 0;
    options.multisession.split = [50 50 50];
    options.multisession.fixed.phi = [1 2 3]; % Gamma and hidden states are estimated per block
        
%    [model(model_count,ii).posterior, model(model_count,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
%     L(model_count,ii) = model(model_count,ii).out.F;
     [ m_temp(4,ii).posterior, m_temp(4,ii).out] = VBA_NLStateSpaceModel(y, u,f_fname, g_fname, dim, options);
 L_temp(4,ii) = m_temp(4,ii).out.F;
 
    clearvars -except results_tab model_count L ii n ntrials nmodels model exclu L_temp m_temp
    
end

%% Bayesian model
% ntrials = 50;
% 
% for ii = 1:n %loop that runs through n participants
%     if ii == 5 || ii == 46 || ii == 48 || ii == 25 || ii == 23 || ii == 67
%         continue
%     end
%     for part = 1:2
%         for prior_diff = 1:2
% %             for repet = 1:10
% 
%                 fprintf('Processing subject No %d of %d, Part %d repetition %d \n',ii,n,part,prior_diff);
%                 nrows_init = 51 + 150*(ii-1) + 50*(part-1);
%                 nrows_final = nrows_init + ntrials - 1;
% 
%                 dim = struct('n', 9,'p',11, 'n_theta', 0, 'n_phi', 0,'n_t',ntrials);
% 
%                 % Multinomial
%                 tmp = results_tab.pred_dec(nrows_init:nrows_final);
%                 y = [];
%                 for gg = 1:50
%                     tmp(gg) = 10- tmp(gg);
%                     if part == 2 && ii == 5 && gg == 48
%                     y = [y, [zeros(1,3) 1 zeros(1,10 - 3)]'];  
%                     else
%                     y = [y, [zeros(1,tmp(gg)) 1 zeros(1,10 - tmp(gg))]'];
%                     end
%                 end
% 
%                 % Decision into binary
% %                 tmp = results_tab.pred_dec(nrows_init:nrows_final);   
% %                 tmp2 = results_tab.pred_true(nrows_init:nrows_final);
% %                 for gg = 1:ntrials
% %                     y(gg) = VBA_random('Binomial',1,tmp(gg)/10);
% %                     u3(gg) = VBA_random('Binomial',1,tmp2(gg)/10);
% %                 end
%                 tmp2 = results_tab.pred_true(nrows_init:nrows_final);
%                 u = [(results_tab.predFalseV(nrows_init:nrows_final) - results_tab.predTrueV(nrows_init:nrows_final))';...
%                     results_tab.predFalseV(nrows_init:nrows_final)';((10-tmp2)./10)'];
% 
%                 g_fname = @gBayesian;
%                 f_fname = @fBayesian;
% 
%                 %Defining priors and options
%                 options.sources = struct('type',2);
%                 options.DisplayWin = 0;
%                 options.verbose = 0;
%                 if prior_diff == 1
%                     priors.muX0 = [1 0 1 1 0 1 -1 0 1]';
%                 elseif prior_diff == 2 && part == 2
%                     priors.muX0 = [0.72 0 1 1.10 0 1 -1.26 0 1]';                     
%                 else
%                     priors.muX0 = [0.11 0 1 1.12 0 1 -1.24 0 1]';                    
%                 end
%                 priors.SigmaX0 = eye(dim.n); % prior covariance (hidden states)
%                 options.priors = priors;
% 
%                 %Inversion/estimation
%                 [model_b(ii,part,prior_diff).posterior, model_b(ii,part,prior_diff).out] = ...
%                     VBA_NLStateSpaceModel(y, u, f_fname, g_fname, dim, options);
%                 clearvars -except results_tab L ii n ntrials model part Part1 Part2 aa prior_diff exclu
% %             end
%         end
%     end
% 
% end
%       
%% Comparing models

%[p, o] = VBA_groupBMC(L_70);
[p, o] = VBA_groupBMC(L);

% %Should we be using families?
% options.families = {[1,2], [3,4], [5,6],[7,8] [9, 10], [11, 12], [13, 14], [15,16],[17,18]} ;
% [p, o] = VBA_groupBMC(L_70,options);

options.families = {[1,2,3,4], [5,6,7,8],[9, 10], [11, 12], [13, 14], [15,16], [17,18], [19,20], [21,22], [23,24]};
[p, o] = VBA_groupBMC(L,options);

L(17:22,:) = [];
options.families = {[1,2,3,4], [5,6,7,8],[9, 10], [11, 12], [13, 14], [15,16], [17,18], [19,20]};
[p, o] = VBA_groupBMC(L,options);


%% Plotting models results

for md = 1:26
        mdl(md).Phi = [];mdl(md).muX=[];
        mdl(md).Theta = [];
        
    for ii = 1:70
        mdl(md).Phi = [mdl(md).Phi, model(md,ii).posterior.muPhi];
        if md > 4
            mdl(md).Theta = [mdl(md).Theta, model(md,ii).posterior.muTheta];
            mdl(md).muX = [mdl(md).muX; model(md,ii).posterior.muX];
        end
    end
end

group = results_tab.bad_group_first(1:150:end);
bgf = find(group);
ggf = find(~group);

nt = 50;

idx.bias = [];idx.biasm = [];
idx.basic = 1;idx.basicm = 2;idx.basicma = 3;idx.basicmd = 4;
idx.fixed = 5;idx.fixedm = 6;idx.fixedma = 7;idx.fixedmd = 8;
idx.conf = 9; idx.confm = 10;
idx.acta = 11; idx.actam = 12; idx.actd = 13; idx.actdm = 14; idx.actad = 15; idx.actadm = 16;
idx.outa = 17; idx.outam = 18; idx.outd = 19; idx.outdm = 20; idx.outad = 21; idx.outadm = 22;
idx.confa = 23; idx.confam = 24;
idx.confo = 25; idx.confom = 26;

%% Parameters 
% 
% %Bias
% figure;
% bar(1:4,  [mean(mdl(idx.bias).Phi(1,:)), mean(mdl(idx.biasm).Phi(1,:)), mean(mdl(idx.biasm).Phi(3,:)), mean(mdl(idx.biasm).Phi(4,:))]);
% xticks(1:4);xticklabels({'Bias', 'Bias baseline', 'Bias good group', 'Bias bad group'});
% hold on
% e1 = errorbar(1:4,[mean(mdl(idx.bias).Phi(1,:)), mean(mdl(idx.biasm).Phi(1,:)), mean(mdl(idx.biasm).Phi(3,:)), mean(mdl(idx.biasm).Phi(4,:))],...
%     [std(mdl(idx.bias).Phi(1,:)), std(mdl(idx.biasm).Phi(1,:)), std(mdl(idx.biasm).Phi(3,:)), std(mdl(idx.biasm).Phi(4,:))]/sqrt(n));
% e1.Color = 'black';
% e1.LineStyle = 'none';   
% hold off
% title('Bias model');
% 
% histogram(mdl(1).Phi(1,:), 10);

%% Basic model
figure;

subplot(2,4,1)
bar(1:2,[mean(exp(mdl(idx.basic).Phi(1,:))),mean(mdl(idx.basic).Phi(2,:))],'blue');
xticks([1 2]);xticklabels({'Alpha' 'Delta'});
hold on
e1 = errorbar(1:2,[mean(exp(mdl(idx.basic).Phi(1,:))),mean(mdl(idx.basic).Phi(2,:))], ...
    [std(exp(mdl(idx.basic).Phi(1,:))),std(mdl(idx.basic).Phi(2,:))]/sqrt(n));
e1.Color = 'black';
e1.LineStyle = 'none';   
hold off
title('Basic model');

% Basic model multi
subplot(2,4,2)
e1 = errorbar(1:3,mean(exp(mdl(idx.basicm).Phi([1 4 6],:)),2),std(exp(mdl(idx.basicm).Phi([1 4 6],:)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(mdl(idx.basicm).Phi([2 5 7],:),2),std(mdl(idx.basicm).Phi([2 5 7],:),0,2)/sqrt(n));
e2.Color = 'green';
yline(0,'k');
xlim([0.5 3.5]); xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha','Delta', 'Location', 'best');title('Basic model multi');

%Basic model multi alpha
subplot(2,4,3)
e1 = errorbar(1:3,mean(exp(mdl(idx.basicma).Phi([1 4 5],:)),2),std(exp(mdl(idx.basicma).Phi([1 4 5],:)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(mdl(idx.basicma).Phi([2 2 2],:),2),std(mdl(idx.basicma).Phi([2 2 2],:),0,2)/sqrt(n));
e2.Color = 'green';
yline(0,'k');
xlim([0.5 3.5]); xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha','Delta', 'Location', 'best');title('Basic model multi alpha');

%Basic model multi alpha
subplot(2,4,4)
e1 = errorbar(1:3,mean(exp(mdl(idx.basicmd).Phi([1 1 1],:)),2),std(exp(mdl(idx.basicmd).Phi([1 1 1],:)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(mdl(idx.basicmd).Phi([2 4 5],:),2),std(mdl(idx.basicmd).Phi([2 4 5],:),0,2)/sqrt(n));
e2.Color = 'green';
yline(0,'k');
xlim([0.5 3.5]); xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha','Delta', 'Location', 'best');title('Basic model multi alpha');

% which group first
subplot(2,4,5)
bar(1:4,[mean(exp(mdl(idx.basic).Phi(1,ggf))),...
    mean(mdl(idx.basic).Phi(2,ggf)),...
    mean(exp(mdl(idx.basic).Phi(1,bgf))),...
    mean(mdl(idx.basic).Phi(2,bgf))],'blue');
xticks([1 4]);xticklabels({'Alpha GG First' 'Delta GG First' 'Alpha BG First' 'Delta BG First'});
hold on
e1 = errorbar(1:4,[mean(exp(mdl(idx.basic).Phi(1,ggf))),...
    mean(mdl(idx.basic).Phi(2,ggf)),...
    mean(exp(mdl(idx.basic).Phi(1,bgf))),...
    mean(mdl(idx.basic).Phi(2,bgf))],...
    [std(exp(mdl(idx.basic).Phi(1,ggf)),0,2)/sqrt(length(ggf)),...
    std(mdl(idx.basic).Phi(2,ggf),0,2)/sqrt(length(ggf)),...
    std(exp(mdl(idx.basic).Phi(1,bgf)),0,2)/sqrt(length(bgf)),...
    std(mdl(idx.basic).Phi(2,bgf),0,2)/sqrt(length(bgf))]);
e1.Color = 'black';
e1.LineStyle = 'none';   
hold off
title('Basic model per group type first');

subplot(2,4,6)
e1 = errorbar(1:3,mean(exp(mdl(idx.basicm).Phi([1 4 6],ggf)),2),std(exp(mdl(idx.basicm).Phi([1 4 6],ggf)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(exp(mdl(idx.basicm).Phi([1 4 6],bgf)),2),std(exp(mdl(idx.basicm).Phi([1 4 6],bgf)),0,2)/sqrt(n));
e2.Color = 'cyan';
e3 = errorbar(1:3,mean(mdl(idx.basicm).Phi([2 5 7],ggf),2),std(mdl(idx.basicm).Phi([2 5 7],ggf),0,2)/sqrt(n));
e3.Color = 'green';
e4 = errorbar(1:3,mean(mdl(idx.basicm).Phi([2 5 7],bgf),2),std(mdl(idx.basicm).Phi([2 5 7],bgf),0,2)/sqrt(n));
e4.Color = 'magenta';
yline(0,'k');
%xlim([0.5 3.5]); 
xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha GG First','Alpha BG First','Delta GG First','Delta BG First', 'Location', 'best');title('Basic model multi per group type first');


subplot(2,4,7)
e1 = errorbar(1:3,mean(exp(mdl(idx.basicma).Phi([1 4 5],ggf)),2),std(exp(mdl(idx.basicma).Phi([1 4 5],ggf)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(exp(mdl(idx.basicma).Phi([1 4 5],bgf)),2),std(exp(mdl(idx.basicma).Phi([1 4 5],bgf)),0,2)/sqrt(n));
e2.Color = 'cyan';
e3 = errorbar(1:3,mean(mdl(idx.basicma).Phi([2 2 2],ggf),2),std(mdl(idx.basicma).Phi([2 2 2],ggf),0,2)/sqrt(n));
e3.Color = 'green';
e4 = errorbar(1:3,mean(mdl(idx.basicma).Phi([2 2 2],bgf),2),std(mdl(idx.basicma).Phi([2 2 2],bgf),0,2)/sqrt(n));
e4.Color = 'magenta';
yline(0,'k');
%xlim([0.5 3.5]); 
xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha GG First','Alpha BG First','Delta GG First','Delta BG First', 'Location', 'best');title('Basic model multi alpha per group type first');


subplot(2,4,8)
e1 = errorbar(1:3,mean(exp(mdl(idx.basicmd).Phi([1 1 1],ggf)),2),std(exp(mdl(idx.basicmd).Phi([1 1 1],ggf)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(exp(mdl(idx.basicmd).Phi([1 1 1],bgf)),2),std(exp(mdl(idx.basicmd).Phi([1 1 1],bgf)),0,2)/sqrt(n));
e2.Color = 'cyan';
e3 = errorbar(1:3,mean(mdl(idx.basicmd).Phi([2 4 5],ggf),2),std(mdl(idx.basicmd).Phi([2 4 5],ggf),0,2)/sqrt(n));
e3.Color = 'green';
e4 = errorbar(1:3,mean(mdl(idx.basicmd).Phi([2 4 5],bgf),2),std(mdl(idx.basicmd).Phi([2 4 5],bgf),0,2)/sqrt(n));
e4.Color = 'magenta';
yline(0,'k');
%xlim([0.5 3.5]); 
xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha GG First','Alpha BG First','Delta GG First','Delta BG First', 'Location', 'best');title('Basic model multi delta per group type first');


%% Fixed cost one

figure;
subplot(2,4,1)
bar(1:2,[mean(exp(mdl(idx.fixed).Phi(1,:))),mean(mdl(idx.fixed).Phi(2,:))],'blue');
xticks([1 2]);xticklabels({'Alpha' 'Delta'});
hold on
e1 = errorbar(1:2,[mean(exp(mdl(idx.fixed).Phi(1,:))),mean(mdl(idx.fixed).Phi(2,:))], [std(exp(mdl(idx.fixed).Phi(1,:))),std(mdl(idx.fixed).Phi(2,:))]/sqrt(n));
e1.Color = 'black';
e1.LineStyle = 'none';   
hold off
title('Fixed cost model');

subplot(2,4,2)
e1 = errorbar(1:3,mean(exp(mdl(idx.fixedm).Phi([1 4 6],:)),2),std(exp(mdl(idx.fixedm).Phi([1 4 6],:)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(mdl(idx.fixedm).Phi([2 5 7],:),2),std(mdl(idx.fixedm).Phi([2 5 7],:),0,2)/sqrt(n));
e2.Color = 'green';
yline(0,'k');
%xlim([0.5 3.5]); 
xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha','Delta', 'Location', 'best');title('Fixed cost model multi');


%Fixed model multi alpha
subplot(2,4,3)
e1 = errorbar(1:3,mean(exp(mdl(idx.fixedma).Phi([1 4 5],:)),2),std(exp(mdl(idx.fixedma).Phi([1 4 5],:)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(mdl(idx.fixedma).Phi([2 2 2],:),2),std(mdl(idx.fixedma).Phi([2 2 2],:),0,2)/sqrt(n));
e2.Color = 'green';
yline(0,'k');
xlim([0.5 3.5]); xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha','Delta', 'Location', 'best');title('Fixed model multi alpha');

%Fixed model multi alpha
subplot(2,4,4)
e1 = errorbar(1:3,mean(exp(mdl(idx.fixedmd).Phi([1 1 1],:)),2),std(exp(mdl(idx.fixedmd).Phi([1 1 1],:)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(mdl(idx.fixedmd).Phi([2 4 5],:),2),std(mdl(idx.fixedmd).Phi([2 4 5],:),0,2)/sqrt(n));
e2.Color = 'green';
yline(0,'k');
xlim([0.5 3.5]); xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha','Delta', 'Location', 'best');title('Fixed model multi alpha');

% Which group was seen first 
subplot(2,4,5)
bar(1:4,[mean(exp(mdl(idx.fixed).Phi(1,ggf))),...
    mean(mdl(idx.fixed).Phi(2,ggf)),...
    mean(exp(mdl(idx.fixed).Phi(1,bgf))),...
    mean(mdl(idx.fixed).Phi(2,bgf))],'blue');
xticks([1 4]);xticklabels({'Alpha GG First' 'Delta GG First' 'Alpha BG First' 'Delta BG First'});
hold on
e1 = errorbar(1:4,[mean(exp(mdl(idx.fixed).Phi(1,ggf))),...
    mean(mdl(idx.fixed).Phi(2,ggf)),...
    mean(exp(mdl(idx.fixed).Phi(1,bgf))),...
    mean(mdl(idx.fixed).Phi(2,bgf))],...
    [std(exp(mdl(idx.fixed).Phi(1,ggf)),0,2)/sqrt(length(ggf)),...
    std(mdl(idx.fixed).Phi(2,ggf),0,2)/sqrt(length(ggf)),...
    std(exp(mdl(idx.fixed).Phi(1,bgf)),0,2)/sqrt(length(bgf)),...
    std(mdl(idx.fixed).Phi(2,bgf),0,2)/sqrt(length(bgf))]);
e1.Color = 'black';
e1.LineStyle = 'none';   
hold off
title('Fixed cost model per group type first');

subplot(2,4,6)
e1 = errorbar(1:3,mean(exp(mdl(idx.fixedm).Phi([1 4 6],ggf)),2),std(exp(mdl(idx.fixedm).Phi([1 4 6],ggf)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(exp(mdl(idx.fixedm).Phi([1 4 6],bgf)),2),std(exp(mdl(idx.fixedm).Phi([1 4 6],bgf)),0,2)/sqrt(n));
e2.Color = 'cyan';
e3 = errorbar(1:3,mean(mdl(idx.fixedm).Phi([2 5 7],ggf),2),std(mdl(idx.fixedm).Phi([2 5 7],ggf),0,2)/sqrt(n));
e3.Color = 'green';
e4 = errorbar(1:3,mean(mdl(idx.fixedm).Phi([2 5 7],bgf),2),std(mdl(idx.fixedm).Phi([2 5 7],bgf),0,2)/sqrt(n));
e4.Color = 'magenta';
yline(0,'k');
xlim([0.5 3.5]); xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha GG First','Alpha BG First','Delta GG First','Delta BG First', 'Location', 'best');title('Fixed cost model multi per group type first');

subplot(2,4,7)
e1 = errorbar(1:3,mean(exp(mdl(idx.fixedma).Phi([1 4 5],ggf)),2),std(exp(mdl(idx.fixedma).Phi([1 4 5],ggf)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(exp(mdl(idx.fixedma).Phi([1 4 5],bgf)),2),std(exp(mdl(idx.fixedma).Phi([1 4 5],bgf)),0,2)/sqrt(n));
e2.Color = 'cyan';
e3 = errorbar(1:3,mean(mdl(idx.fixedma).Phi([2 2 2],ggf),2),std(mdl(idx.fixedma).Phi([2 2 2],ggf),0,2)/sqrt(n));
e3.Color = 'green';
e4 = errorbar(1:3,mean(mdl(idx.fixedma).Phi([2 2 2],bgf),2),std(mdl(idx.fixedma).Phi([2 2 2],bgf),0,2)/sqrt(n));
e4.Color = 'magenta';
yline(0,'k');
%xlim([0.5 3.5]); 
xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha GG First','Alpha BG First','Delta GG First','Delta BG First', 'Location', 'best');title('Basic model multi alpha per group type first');


subplot(2,4,8)
e1 = errorbar(1:3,mean(exp(mdl(idx.fixedmd).Phi([1 1 1],ggf)),2),std(exp(mdl(idx.fixedmd).Phi([1 1 1],ggf)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(exp(mdl(idx.fixedmd).Phi([1 1 1],bgf)),2),std(exp(mdl(idx.fixedmd).Phi([1 1 1],bgf)),0,2)/sqrt(n));
e2.Color = 'cyan';
e3 = errorbar(1:3,mean(mdl(idx.fixedmd).Phi([2 4 5],ggf),2),std(mdl(idx.fixedmd).Phi([2 4 5],ggf),0,2)/sqrt(n));
e3.Color = 'green';
e4 = errorbar(1:3,mean(mdl(idx.fixedmd).Phi([2 4 5],bgf),2),std(mdl(idx.fixedmd).Phi([2 4 5],bgf),0,2)/sqrt(n));
e4.Color = 'magenta';
yline(0,'k');
%xlim([0.5 3.5]); 
xticks([1 2 3]);xticklabels({'Baseline' 'GG' 'BG'});
legend('Alpha GG First','Alpha BG First','Delta GG First','Delta BG First', 'Location', 'best');title('Basic model multi delta per group type first');


%% Conformity model 
figure;
subplot(2,2,1)
bar(1:3,[mean(exp(mdl(idx.conf).Phi(1,:))),mean(mdl(idx.conf).Phi(2,:)), mean(mdl(idx.conf).Phi(4,:))],'blue');
xticks([1 2 3]);xticklabels({'Alpha' 'Delta', 'Gamma'});
hold on
e1 = errorbar(1:3,[mean(exp(mdl(idx.conf).Phi(1,:))),mean(mdl(idx.conf).Phi(2,:)), mean(mdl(idx.conf).Phi(4,:))], ...
    [std(exp(mdl(idx.conf).Phi(1,:))),std(mdl(idx.conf).Phi(2,:)), std(mdl(idx.conf).Phi(4,:))]/sqrt(n));
e1.Color = 'black';
e1.LineStyle = 'none';   
hold off
title('Conformity model');

subplot(2,2,2)
e1 = errorbar(1:3,mean(exp(mdl(idx.confm).Phi([1 1 1],:)),2),std(exp(mdl(idx.confm).Phi([1 1 1],:)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(mdl(idx.confm).Phi([2 2 2],:),2),std(mdl(idx.confm).Phi([2 2 2],:),0,2)/sqrt(n));
e2.Color = 'green';
e2 = errorbar(1:3,mean(mdl(idx.confm).Phi([4 5 6],:),2),std(mdl(idx.confm).Phi([4 5 6],:),0,2)/sqrt(n));
e2.Color = 'red';
yline(0,'k');
%xlim([0.5 3.5]); 
xticks([1 2 3]);xticklabels({'Baseline' 'HG' 'DG'});
legend('Alpha','Delta','Gamma');title('Conformity model multi');

subplot(2,2,3)
bar(1:6,[mean(exp(mdl(idx.conf).Phi(1,ggf))),...
    mean(mdl(idx.conf).Phi(2,ggf)),...
    mean(mdl(idx.conf).Phi(4,ggf)),...
    mean(exp(mdl(idx.conf).Phi(1,bgf))),...
    mean(mdl(idx.conf).Phi(2,bgf)),...
    mean(mdl(idx.conf).Phi(4,bgf))],'blue');
xticks(1:6);xticklabels({'Alpha GG First' 'Delta GG First' 'Gamma GG First' 'Alpha BG First' 'Delta BG First' 'Gamma BG First'});
hold on
e1 = errorbar(1:6,[mean(exp(mdl(idx.conf).Phi(1,ggf))),...
    mean(mdl(idx.conf).Phi(2,ggf)),...
    mean(mdl(idx.conf).Phi(4,ggf)),...
    mean(exp(mdl(idx.conf).Phi(1,bgf))),...
    mean(mdl(idx.conf).Phi(2,bgf)),...
    mean(mdl(idx.conf).Phi(4,bgf))],...
    [std(exp(mdl(idx.conf).Phi(1,ggf)),0,2)/sqrt(length(ggf)),...
    std(mdl(idx.conf).Phi(2,ggf),0,2)/sqrt(length(ggf)),...
    std(mdl(idx.conf).Phi(4,ggf),0,2)/sqrt(length(ggf)),...
    std(exp(mdl(idx.conf).Phi(1,bgf)),0,2)/sqrt(length(bgf)),...
    std(mdl(idx.conf).Phi(2,bgf),0,2)/sqrt(length(bgf)),...
    std(mdl(idx.conf).Phi(4,bgf),0,2)/sqrt(length(bgf))]);
e1.Color = 'black';
e1.LineStyle = 'none';   
hold off
title('Conformity model per group type first');

subplot(2,2,4)
e1 = errorbar(1:3,mean(exp(mdl(idx.confm).Phi([1 1 1],ggf)),2),std(exp(mdl(idx.confm).Phi([1 1 1],ggf)),0,2)/sqrt(n));
e1.Color = 'blue';
hold on
e2 = errorbar(1:3,mean(exp(mdl(idx.confm).Phi([1 1 1],bgf)),2),std(exp(mdl(idx.confm).Phi([1 1 1],bgf)),0,2)/sqrt(n));
e2.Color = 'cyan';
e3 = errorbar(1:3,mean(mdl(idx.confm).Phi([2 2 2],ggf),2),std(mdl(idx.confm).Phi([2 2 2],ggf),0,2)/sqrt(n));
e3.Color = 'green';
e4 = errorbar(1:3,mean(mdl(idx.confm).Phi([2 2 2],bgf),2),std(mdl(idx.confm).Phi([2 2 2],bgf),0,2)/sqrt(n));
e4.Color = 'magenta';
e5 = errorbar(1:3,mean(mdl(idx.confm).Phi([4 5 6],ggf),2),std(mdl(idx.confm).Phi([4 5 6],ggf),0,2)/sqrt(n));
e5.Color = 'yellow';
e6 = errorbar(1:3,mean(mdl(idx.confm).Phi([4 5 6],bgf),2),std(mdl(idx.confm).Phi([4 5 6],bgf),0,2)/sqrt(n));
e6.Color = 'black';
yline(0,'k');
hold off
%xlim([0.5 3.5]); 
xticks([1 2 3]);xticklabels({'Baseline' 'HG' 'DG'});
legend('Alpha HG First','Alpha DG First', 'Delta HG First','Delta DG First','Gamma HG First', 'Gamma DG First','Location', 'bestoutside');title('Conformity model multi per group type first');


%% Action models

% Hidden states (sanity check)
figure;
subplot(2,3,1)
e1 = errorbar(1:50,mean(mdl(idx.acta).muX(1:2:end,51:100)),std(mdl(idx.acta).muX(1:2:end,51:100))/sqrt(n));
e1.Color = 'green';
hold on
e2 = errorbar(1:50,mean(mdl(idx.acta).muX(2:2:end,101:end)),std(mdl(idx.acta).muX(2:2:end,101:end))/sqrt(n));
e2.Color = 'red';
legend('good group','bad group');
title('RL Alpha');
xlim([0 50]);
hold off

subplot(2,3,2)
e3 = errorbar(1:50,mean(mdl(idx.actam).muX(3:6:end,51:100)),std(mdl(idx.actam).muX(3:6:end,51:100))/sqrt(n),'--');
e3.Color = 'green';
hold on
e3 = errorbar(1:50,mean(mdl(idx.actam).muX(6:6:end,101:end)),std(mdl(idx.actam).muX(6:6:end,101:end))/sqrt(n),'--');
e3.Color = 'red';
legend('good group','bad group');
title('RL Alpha Multi');
xlim([0 50]);
hold off

subplot(2,3,3)
e1 = errorbar(1:50,mean(mdl(idx.actd).muX(1:2:end,51:100)),std(mdl(idx.actd).muX(1:2:end,51:100))/sqrt(n));
e1.Color = 'green';
hold on
e2 = errorbar(1:50,mean(mdl(idx.actd).muX(2:2:end,101:end)),std(mdl(idx.actd).muX(2:2:end,101:end))/sqrt(n));
e2.Color = 'red';
legend('good group','bad group');
title('RL Delta');
xlim([0 50]);
hold off

subplot(2,3,4)
e3 = errorbar(1:50,mean(mdl(idx.actdm).muX(3:6:end,51:100)),std(mdl(idx.actdm).muX(3:6:end,51:100))/sqrt(n),'--');
e3.Color = 'green';
hold on
e3 = errorbar(1:50,mean(mdl(idx.actdm).muX(6:6:end,101:end)),std(mdl(idx.actdm).muX(6:6:end,101:end))/sqrt(n),'--');
e3.Color = 'red';
legend('good group','bad group');
title('RL Delta Multi');
xlim([0 50]);
hold off

subplot(2,3,5)
e1 = errorbar(1:50,mean(mdl(idx.actad).muX(1:4:end,51:100)),std(mdl(idx.actad).muX(1:4:end,51:100))/sqrt(n));
e1.Color = 'green';
hold on
e2 = errorbar(1:50,mean(mdl(idx.actad).muX(2:4:end,101:end)),std(mdl(idx.actad).muX(2:4:end,101:end))/sqrt(n));
e2.Color = 'red';
hold off
legend('Alpha GG','Alpha BG');
xlim([0 50]);
title('RL Alpha-Delta');

subplot(2,3,6)
e3 = errorbar(1:50,mean(mdl(idx.actad).muX(3:4:end,51:100)),std(mdl(idx.actad).muX(3:4:end,51:100))/sqrt(n),'--');
e3.Color = 'green';
hold on
e3 = errorbar(1:50,mean(mdl(idx.actad).muX(4:4:end,101:end)),std(mdl(idx.actad).muX(4:4:end,101:end))/sqrt(n),'--');
e3.Color = 'red';
legend('Delta GG','Delta BG');
title('RL Alpha-Delta');
xlim([0 50]);
hold off


figure;
e1 = errorbar(1:50,mean(mdl(idx.actadm).muX(5:12:end,51:100)),std(mdl(idx.actadm).muX(5:12:end,51:100))/sqrt(n));
e1.Color = 'green';
hold on
e2 = errorbar(1:50,mean(mdl(idx.actadm).muX(10:12:end,101:end)),std(mdl(idx.actadm).muX(10:12:end,101:end))/sqrt(n));
e2.Color = 'red';
e3 = errorbar(1:50,mean(mdl(idx.actadm).muX(7:12:end,51:100)),std(mdl(idx.actadm).muX(7:12:end,51:100))/sqrt(n),'--');
e3.Color = 'blue';
e3 = errorbar(1:50,mean(mdl(idx.actadm).muX(12:12:end,101:end)),std(mdl(idx.actadm).muX(12:12:end,101:end))/sqrt(n),'--');
e3.Color = 'magenta';
legend('Alpha GG','Alpha BG','Delta GG','Delta BG');
title('RL Alpha-Delta Multi');
hold off


% Learning rate & parameters
figure
subplot(2,2,1)
plot(mdl(idx.actam).Theta(2,:), mdl(idx.actam).Theta(3,:), 'o');
hold on
plot(mdl(idx.acta).Theta, mdl(idx.acta).Theta,'r*');
hold off
xlabel('good group'); ylabel('bad group');
title('LR RL-Alpha Multi. In red, LR of basic.');

subplot(2,2,2)
plot(mdl(idx.actdm).Theta(2,:), mdl(idx.actdm).Theta(3,:), 'o');
hold on
plot(mdl(idx.actd).Theta, mdl(idx.actd).Theta,'r*');
hold off
xlabel('good group'); ylabel('bad group');
title('LR RL-Delta Multi. In red, LR of basic.');

subplot(2,2,3)
plot(mdl(idx.actad).Theta(1,:), mdl(idx.actad).Theta(2,:), 'o');
xlabel('alpha'); ylabel('delta');
title('LR RL-Alpha-Delta');

subplot(2,2,4)
plot(mdl(idx.actadm).Theta(3,:), mdl(idx.actadm).Theta(5,:), 'o');
hold on
plot(mdl(idx.actadm).Theta(4,:), mdl(idx.actadm).Theta(6,:), 'r*');
hold off
legend('Alpha', 'Delta');
xlabel('good group'); ylabel('bad group');
title('LR RL-Alpha-Delta Multi');


%% Outcome models

% Hidden states (sanity check)
figure;
subplot(2,3,1)
e1 = errorbar(1:50,mean(mdl(idx.outa).muX(1:2:end,51:100)),std(mdl(idx.outa).muX(1:2:end,51:100))/sqrt(n));
e1.Color = 'green';
hold on
e2 = errorbar(1:50,mean(mdl(idx.outa).muX(2:2:end,101:end)),std(mdl(idx.outa).muX(2:2:end,101:end))/sqrt(n));
e2.Color = 'red';
legend('good group','bad group');
title('RL Outcome Alpha');
xlim([0 50]);
hold off

subplot(2,3,2)
e3 = errorbar(1:50,mean(mdl(idx.outam).muX(3:6:end,51:100)),std(mdl(idx.outam).muX(3:6:end,51:100))/sqrt(n),'--');
e3.Color = 'green';
hold on
e3 = errorbar(1:50,mean(mdl(idx.outam).muX(5:6:end,101:end)),std(mdl(idx.outam).muX(5:6:end,101:end))/sqrt(n),'--');
e3.Color = 'red';
legend('good group','bad group');
title('RL Outcome Alpha Multi');
xlim([0 50]);
hold off

subplot(2,3,3)
e1 = errorbar(1:50,mean(mdl(idx.outd).muX(1:2:end,51:100)),std(mdl(idx.outd).muX(1:2:end,51:100))/sqrt(n));
e1.Color = 'green';
hold on
e2 = errorbar(1:50,mean(mdl(idx.outd).muX(2:2:end,101:end)),std(mdl(idx.outd).muX(2:2:end,101:end))/sqrt(n));
e2.Color = 'red';
legend('good group','bad group');
title('RL Outcome Delta');
xlim([0 50]);
hold off

subplot(2,3,4)
e3 = errorbar(1:50,mean(mdl(idx.outdm).muX(3:6:end,51:100)),std(mdl(idx.outdm).muX(3:6:end,51:100))/sqrt(n),'--');
e3.Color = 'green';
hold on
e3 = errorbar(1:50,mean(mdl(idx.outdm).muX(5:6:end,101:end)),std(mdl(idx.outdm).muX(5:6:end,101:end))/sqrt(n),'--');
e3.Color = 'red';
legend('good group','bad group');
title('RL Outcome Delta Multi');
xlim([0 50]);
hold off

subplot(2,3,5)
e1 = errorbar(1:50,mean(mdl(idx.outad).muX(1:4:end,51:100)),std(mdl(idx.outad).muX(1:4:end,51:100))/sqrt(n));
e1.Color = 'green';
hold on
e2 = errorbar(1:50,mean(mdl(idx.outad).muX(2:4:end,101:end)),std(mdl(idx.outad).muX(2:4:end,101:end))/sqrt(n));
e2.Color = 'red';
hold off
legend('Alpha GG','Alpha BG');
xlim([0 50]);
title('RL Outcome Alpha-Delta');

subplot(2,3,6)
e3 = errorbar(1:50,mean(mdl(idx.outad).muX(3:4:end,51:100)),std(mdl(idx.outad).muX(3:4:end,51:100))/sqrt(n),'--');
e3.Color = 'green';
hold on
e3 = errorbar(1:50,mean(mdl(idx.outad).muX(4:4:end,101:end)),std(mdl(idx.outad).muX(4:4:end,101:end))/sqrt(n),'--');
e3.Color = 'red';
legend('Delta GG','Delta BG');
title('RL Outcome Alpha-Delta');
xlim([0 50]);
hold off


figure;
e1 = errorbar(1:50,mean(mdl(idx.outadm).muX(5:12:end,51:100)),std(mdl(idx.outadm).muX(5:12:end,51:100))/sqrt(n));
e1.Color = 'green';
hold on
e2 = errorbar(1:50,mean(mdl(idx.outadm).muX(10:12:end,101:end)),std(mdl(idx.outadm).muX(10:12:end,101:end))/sqrt(n));
e2.Color = 'red';
e3 = errorbar(1:50,mean(mdl(idx.outadm).muX(7:12:end,51:100)),std(mdl(idx.outadm).muX(7:12:end,51:100))/sqrt(n),'--');
e3.Color = 'blue';
e3 = errorbar(1:50,mean(mdl(idx.outadm).muX(12:12:end,101:end)),std(mdl(idx.outadm).muX(12:12:end,101:end))/sqrt(n),'--');
e3.Color = 'magenta';
legend('Alpha GG','Alpha BG','Delta GG','Delta BG');
title('RL Outcome Alpha-Delta Multi');
hold off


% Learning rate & parameters
figure
subplot(2,2,1)
plot(mdl(idx.outam).Theta(2,:), mdl(idx.outam).Theta(3,:), 'o');
hold on
plot(mdl(idx.outa).Theta, mdl(idx.outa).Theta,'r*');
hold off
xlabel('good group'); ylabel('bad group');
title('LR Outcome RL-Alpha Multi. In red, LR of basic.');

subplot(2,2,2)
plot(mdl(idx.outdm).Theta(2,:), mdl(idx.outdm).Theta(3,:), 'o');
hold on
plot(mdl(idx.outd).Theta, mdl(idx.outd).Theta,'r*');
hold off
xlabel('good group'); ylabel('bad group');
title('LR Outcome RL-Delta Multi. In red, LR of basic.');

subplot(2,2,3)
plot(mdl(idx.outad).Theta(1,:), mdl(idx.outad).Theta(2,:), 'o');
xlabel('alpha'); ylabel('delta');
title('LR RL-Alpha-Delta');

subplot(2,2,4)
plot(mdl(idx.outadm).Theta(3,:), mdl(idx.outadm).Theta(5,:), 'o');
hold on
plot(mdl(idx.outadm).Theta(4,:), mdl(idx.outadm).Theta(6,:), 'r*');
hold off
legend('Alpha', 'Delta');
xlabel('good group'); ylabel('bad group');
title('LR RL-Alpha-Delta Multi');




%% Conformity Action
figure;

subplot(1,3,1)
bar(1:3,[mean(exp(mdl(idx.confa).Phi(1,:))),mean(mdl(15).Phi(2,:)), mean(mdl(idx.confa).Phi(4,:))],'blue');
xticks([1 2 3]);xticklabels({'Alpha' 'Delta', 'Gamma'});
hold on
e1 = errorbar(1:3,[mean(exp(mdl(idx.confa).Phi(1,:))),mean(mdl(idx.confa).Phi(2,:)), mean(mdl(idx.confa).Phi(4,:))], ...
    [std(exp(mdl(idx.confa).Phi(1,:))),std(mdl(idx.confa).Phi(2,:)), std(mdl(idx.confa).Phi(4,:))]/sqrt(n));
e1.Color = 'black';
e1.LineStyle = 'none';   
hold off
title('Parameters Conformity Action');

subplot(1,3,2)
e1 = errorbar(1:50, mean(mdl(idx.confa).muX(:,51:100)), std(mdl(idx.confa).muX(:,51:100))/sqrt(n));
e1.Color = 'blue';
hold on
e4 = errorbar(1:50, mean(mdl(idx.confa).muX(:,101:150)), std(mdl(idx.confa).muX(:,101:150))/sqrt(n));
e4.Color = 'yellow';
e2 = errorbar(1:50, mean(mdl(idx.confam).muX(2:3:end,51:100)), std(mdl(idx.confam).muX(2:3:end,51:100))/sqrt(n));
e2.Color = 'red';
e3 = errorbar(1:50, mean(mdl(idx.confam).muX(3:3:end,101:150)), std(mdl(idx.confam).muX(3:3:end,101:150))/sqrt(n));
e3.Color = 'green';
hold off
xlim([0 50]);
legend('basic GG', 'basic BG', 'multi GG', 'multi BG');
title('Hidden states Conformity Action');

subplot(1,3,3)
plot(mdl(idx.confam).Theta(2,:), mdl(idx.confam).Theta(3,:), 'o');
hold on
plot(mdl(idx.confa).Theta, mdl(idx.confa).Theta, '*');
hold off
title('Learning Rate Conformity Action. In red, basic model')



%% Conformity Outcome
figure;

subplot(1,3,1)
bar(1:3,[mean(exp(mdl(idx.confo).Phi(1,:))),mean(mdl(idx.confo).Phi(2,:)), mean(mdl(idx.confo).Phi(4,:))],'blue');
xticks([1 2 3]);xticklabels({'Alpha' 'Delta', 'Gamma'});
hold on
e1 = errorbar(1:3,[mean(exp(mdl(idx.confo).Phi(1,:))),mean(mdl(idx.confo).Phi(2,:)), mean(mdl(idx.confo).Phi(4,:))], ...
    [std(exp(mdl(idx.confo).Phi(1,:))),std(mdl(idx.confo).Phi(2,:)), std(mdl(idx.confo).Phi(4,:))]/sqrt(n));
e1.Color = 'black';
e1.LineStyle = 'none';   
hold off
title('Parameters Conformity Outcome');

subplot(1,3,2)
e1 = errorbar(1:50, mean(mdl(idx.confo).muX(:,51:100)), std(mdl(idx.confo).muX(:,51:100))/sqrt(n));
e1.Color = 'blue';
hold on
e4 = errorbar(1:50, mean(mdl(idx.confo).muX(:,101:150)), std(mdl(idx.confo).muX(:,101:150))/sqrt(n));
e4.Color = 'yellow';
e2 = errorbar(1:50, mean(mdl(idx.confom).muX(2:3:end,51:100)), std(mdl(idx.confom).muX(2:3:end,51:100))/sqrt(n));
e2.Color = 'red';
e3 = errorbar(1:50, mean(mdl(idx.confom).muX(3:3:end,101:150)), std(mdl(idx.confom).muX(3:3:end,101:150))/sqrt(n));
e3.Color = 'green';
hold off
xlim([0 50]);
legend('basic GG', 'basic BG', 'multi GG', 'multi BG');
title('Hidden states Conformity Outcome');

subplot(1,3,3)
plot(mdl(idx.confom).Theta(2,:), mdl(idx.confom).Theta(3,:), 'o');
hold on
plot(mdl(idx.confo).Theta, mdl(idx.confo).Theta, '*');
hold off
title('Learning Rate Conformity Outcome. In red, basic model')

