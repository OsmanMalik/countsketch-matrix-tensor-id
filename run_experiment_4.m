% RUN_EXPERIMENT_4 Run tensor ID experiment for sparse tensors

%% Settings

N = 5;
%Is = [10*1e+3 25*1e+3 50*1e+3 100*1e+3 250*1e+3 500*1e+3 1e+6];
%Is = [1e+3 2e+3];
Is = [1e+3 2.5e+3 5e+3 1e+4 2.5e+4 5e+4 1e+5 2.5e+5];
K = 1e+3;
L = K + 10;
mn = 8;
no_trials = 10;
K_mult = 10;
lambda = 10.^(-((0:K_mult*K-1)/(K_mult*K))*mn);
fac_mat_dens = .01;
s_norm_tol = 1e-12;
maxit = 1000;
results_matlab_file = 'matlab_output_exp_4';
verbosity = 1;
cnt = 1;

%% Main loop

% Create mat and set up matfile for saving results computed in Matlab
save_mat = matfile(results_matlab_file, 'Writable', true);
save_mat.I = zeros(1, length(Is)*no_trials);
save_mat.trial = zeros(1, length(Is)*no_trials);
save_mat.time = zeros(3, length(Is)*no_trials);
save_mat.error = zeros(3, length(Is)*no_trials);

for i = 1:length(Is)

    I = Is(i)*ones(1, N);
    
    for tr = 1:no_trials
        
        if verbosity >= 1
            fprintf('\nStarting experiments for I = %.1e, trial = %d\n', Is(i), tr);
        end
        
        % Generate sparse tensor X
        if verbosity >= 1
            fprintf('Generating sparse tensor... ');
        end
        X = generate_tensor(N, I, K_mult*K, 'lambda_type', 'custom', 'lambda', lambda, 'density', fac_mat_dens);
        if verbosity >= 1
            fprintf('Done!\n')
        end
        
        % Compute gram tensor ID
        if verbosity >= 1
            fprintf('Running gram tensor ID... ');
        end
        gram_tic = tic;
        Xk_gram = gram_tensor_ID(X, K);
        gram_toc = toc(gram_tic);
        if verbosity >= 1
            fprintf('Done!\n');
        end
                
        % Compute Gaussian tensor ID
        if verbosity >= 1
            fprintf('Running Gaussian tensor ID... ');
        end
        gaussian_tic = tic;
        Xk_gaussian = gaussian_tensor_ID(X, K, L, 'qr');
        gaussian_toc = toc(gaussian_tic);
        if verbosity >= 1
            fprintf('Done!\n');
        end
        
        % Compute CountSketch tensor ID
        if verbosity >= 1
            fprintf('Running CountSketch tensor ID... ');
        end
        CS_tic = tic;
        Xk_CS = CS_tensor_ID(X, K, L, 'qr');
        CS_toc = toc(CS_tic);
        if verbosity >= 1
            fprintf('Done!\n');
        end        
        
        % Compute errors
        if verbosity >= 1
            fprintf('Computing errors...\n');
        end
        
        gram_error = s_norm(X-Xk_gram, s_norm_tol, 'verbosity', verbosity, 'init', 'mean', 'maxit', maxit);
        if verbosity >= 1
            fprintf('Gram tensor ID error: %.10e. Time: %.2f s.\n', gram_error, gram_toc);
        end
        
        gaussian_error = s_norm(X-Xk_gaussian, s_norm_tol, 'verbosity', verbosity, 'init', 'mean', 'maxit', maxit);
        if verbosity >= 1
            fprintf('Gaussian tensor ID error: %.10e. Time: %.2f s.\n', gaussian_error, gaussian_toc);
        end
        
        CS_error = s_norm(X-Xk_CS, s_norm_tol, 'verbosity', verbosity, 'init', 'mean', 'maxit', maxit);
        if verbosity >= 1
            fprintf('CountSketch tensor ID error: %.10e. Time: %.2f s.\n', CS_error, CS_toc);
        end
        
        % Save errors and times
        if verbosity >= 1
            fprintf('Done computing errors. \nSaving to file... ')
        end        
        save_mat.I(1, cnt)      = Is(i);
        save_mat.trial(1, cnt)  = tr;   
        save_mat.error(1, cnt)  = gram_error; 
        save_mat.error(2, cnt)  = gaussian_error; 
        save_mat.error(3, cnt)  = CS_error;
        save_mat.time(1, cnt)   = gram_toc;
        save_mat.time(2, cnt)   = gaussian_toc;
        save_mat.time(3, cnt)   = CS_toc;
        if verbosity >= 1
            fprintf('Done!\n')
        end
        cnt = cnt + 1; 
        
    end
    
end