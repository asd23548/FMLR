function [output, opt_model] = FMEM_v9(data,options)
%% Prepare the data
X = data.X;
Y = data.Y;
mask = data.mask;
imagesize = size(mask);
index = mask(:);
aug_mat = sparse(find(index), 1:sum(index), 1, length(index), sum(index));

%% Prepare the options for algorithm
q_max = options.q_max;
band_width_s = options.band_width_s;
lambda_s = options.lambda_s;
%% Start of the algorithm
K = cell(length(band_width_s),1);
for i = 1:length(band_width_s)
    band_width = band_width_s(i);
    K{i} = Gaussian_Kernal_Sparse(mask, 4, band_width);
end
output = cell(1);
for q = 1:q_max
    % Initial values
    n = length(Y);
    C_old = zeros(n,q);
    C = C_old;
    Temp = rand(n,q);
    Tau = Temp;
    for i = 1:n
        Tau(i,:) = Temp(i,:)./sum(Temp(i,:));
    end
    % Tau = ones(n,q)/q;
    Pi = sum(Tau)'/n;
    Phi = Tau;
    Sigma2 = var(Y);
    loss = zeros(q,1);
    loss_total = zeros(q,1);
    %% EM algorithm
    if q == 1
        iter_max =1;
        Tau_opt = ones(n,1);
    else
        iter_max = options.iter_max;
    end
    fprintf('Max Iteration: %5d\n',iter_max)
    fprintf('Number of Componets: %5d\n ',q)
    fprintf('========================EM algorithm started========================\n')
    tic;
    t_old = 0;
    K_old = 0;
    lambda_old = 0;
    opt_level = 1;
    fprintf('Iteration |    BIC    |   Inner Time  |  Total Time  |  Proportion\n')
    for iter = 1:iter_max
        %% M-Step
        BIC_o = inf;
        if opt_level<10
            for i_ker = 1:length(band_width_s)
                K_i = K{i_ker};
                W_tp = X*K_i*X';
                for i_lam = 1:length(lambda_s)
                    lambda_i = lambda_s(i_lam);
                    C_tp = C;
                    for j = 1:q
                        T = diag(Tau(:,j));
                        C_tp(:,j) = pinv(W_tp'*T*W_tp+n*lambda_i*W_tp)...
                            *(W_tp'*T*Y);
                    end
                    for j = 1:q
                        loss(j) = Tau(:,j)'*(Y - W_tp*C_tp(:,j)).^2;
                        loss_total(j) = loss(j) +n*lambda_i*C_tp(:,j)'*W_tp*C_tp(:,j);
                    end
                    Sigma2_tp = sum(loss_total)/n;
                    % Check selection creteria
                    H_tp = 0;
                    for i = 1:q
                        T = diag(Tau(:,i));
                        H_tp = H_tp + T*W_tp*pinv(W_tp'*T*W_tp+n*lambda_i*W_tp)...
                            *(W_tp'*T);
                    end
                    df = trace(H_tp);   y_hat = H_tp*Y;
                    MSE = mean((Y-y_hat).^2);
                    %     AIC = 2 * df + n * log(MSE);
                    BIC = log(n) * df + n * log(MSE);
                    %     GCV = MSE / (1 - df / n)^2;
                    if BIC<BIC_o
                        BIC_o = BIC;
                        C = C_tp;
                        Sigma2 = Sigma2_tp;
                        W = W_tp;
                        H = H_tp;
                        lambda_opt = lambda_i;
                        K_opt = K_i;
                    end
                end
            end
        else
            K_i = K_opt;
            W_tp = X*K_i*X';
            lambda_i = lambda_opt;
            C_tp = C;
            for j = 1:q
                T = diag(Tau(:,j));
                C_tp(:,j) = pinv(W_tp'*T*W_tp+n*lambda_i*W_tp)...
                    *(W_tp'*T*Y);
            end
            for j = 1:q
                loss(j) = Tau(:,j)'*(Y - W_tp*C_tp(:,j)).^2;
                loss_total(j) = loss(j) +n*lambda_i*C_tp(:,j)'*W_tp*C_tp(:,j);
            end
            Sigma2_tp = sum(loss_total)/n;
            % Check selection creteria
            H_tp = 0;
            for i = 1:q
                T = diag(Tau(:,i));
                H_tp = H_tp + T*W_tp*pinv(W_tp'*T*W_tp+n*lambda_i*W_tp)...
                    *(W_tp'*T);
            end
            df = trace(H_tp);   y_hat = H_tp*Y;
            MSE = mean((Y-y_hat).^2);
            %     AIC = 2 * df + n * log(MSE);
            BIC = log(n) * df + n * log(MSE);
            BIC_o = BIC;
            C = C_tp;
            Sigma2 = Sigma2_tp;
            W = W_tp;
            H = H_tp;
            lambda_opt = lambda_i;
            K_opt = K_i;
        end
        if isequal(K_old , K_opt)&& (lambda_old == lambda_opt)
            opt_level = opt_level+1;
        else
            opt_level = 1;
            K_old = K_opt;
            lambda_old = lambda_opt;
        end
        %% E-Step
        
        for j = 1:q
            c_j = C(:,j);
            Phi(:,j) = Pi(j)/sqrt(Sigma2)*exp(-(Y-W*c_j).^2/(2*Sigma2));
        end
        
        for i = 1 : n
            Sum_j = sum(Phi(i,:));
            for j = 1 : q
                Tau (i,j) = Phi(i,j)/Sum_j;
            end
        end
        Pi = sum(Tau)'/n;
        t_new = toc;
        
        %% Stop Check
        if iter>10
            if ((sum((C(:)-C_old(:)).^2)/sum(C_old(:).^2))<10^-6) || ...
                    min(Pi)<0.001
                iter_max = iter-1;
                break
            else
                C_old = C;
                Tau_opt = Tau;
            end
        end
        %         fprintf('Iteration |    BIC    |   Inner Time  |  Total Time  |  Proportion\n')
        fprintf(' %5d     |',iter)
        fprintf('%10.2f |',BIC)
        fprintf('%10.1f     |',t_new-t_old);
        fprintf('%9.1f     |  ',t_new);
        if q>1
            for prop_i = 1:q-1
                fprintf('%2.0f/',Pi(prop_i)*100)
            end
        end
        fprintf('%2.0f\n',Pi(q)*100)
        t_old = t_new;
        
    end
    fprintf('Training completed\n');
    
    %Prune the model
    col_index = Pi'>0.01;
    C = C(:,col_index);
    Pi = normalize_to_one(Pi',col_index);
    Tau_opt = normalize_to_one(Tau_opt,col_index);
    q_e = sum(col_index);
    
    model.df = trace(H);
    model.y_hat = H*Y;
    
    MSE = mean((Y-model.y_hat).^2);
    
    model.C = C;
    model.beta_est=K_opt*X'*C;
    model.beta_img = cell(q_e,1);
    for i_q = 1:q_e
        model.beta_img{i_q} = reshape(aug_mat*model.beta_est(:,i_q),imagesize);
    end
    model.Pi = Pi;
    model.Sigma2 = Sigma2;
    model.Tau = Tau_opt;
    model.MSE = MSE;
    model.iter_max = iter_max;
    model.AIC = 2 * model.df + n * log(MSE);
    model.BIC = log(n) * model.df + n * log(MSE);
    model.GCV = MSE / (1 - model.df / n)^2;
    alpha = 0.2;
    model.order_BIC(q) = sum(loss_total) + q_e^alpha*n^(1-alpha)*log(n);
    
    output{q} = model;

    if q>1
        if (model.order_BIC(q)>model.order_BIC(q-1)) ...
            || (min(model.Pi)<0.01)
            break
        end
    end
end

opt_model = output{q-1}