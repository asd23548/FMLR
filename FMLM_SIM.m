function [output, opt_model]=FMLM_SIM(num_comp)
% num_comp specifies the number of mixing components in the simulation
% study.

% Generate Data
n_each = 80;
n = n_each*num_comp;
D1 = 32;
D2 = 32; 

imagesize = [D1,D2];
switch num_comp
    case 2
        c1 = [D1/3,D1/3*2];
        c2 = [D2/3,D2/3*2];
    case 3
        c1 = [D1/4,D1/4,D1/4*3];
        c2 = [D2/4,D2/4*3,D2/2];
    case 4
        c1 = [D1/4,D1/4,D1/4*3,D1/4*3];
        c2 = [D2/4,D2/4*3,D2/4,D2/4*3];
    case 5
        c1 = [D1/4,D1/4,D1/4*3,D1/4*3,D1/2];
        c2 = [D2/4,D2/4*3,D2/4,D2/4*3,D2/2];
end
Beta = cell(num_comp,1);
for k =1:num_comp
    for  i =  1:32
        for j =1:32
            Beta{k}(i,j) = exp((-(i-c1(k))^2-(j-c2(k))^2)/10);
        end
    end
end

sigma = 1;
d1 = 0;
d2 = 0.5;
Corr.name = 'exp';
Corr.c0 = 1;
Corr.sigma = sigma;
x = linspace(d1, d2, D1)';
[X1, X2] = meshgrid(x, x);
Mesh = [X1(:), X2(:)]; % 2-D mesh
data.x = [d1 * ones(D2, 1), x; d2 * ones(D2, 1), x; x(2:end-1), ...
    d1 * ones(D1-2, 1); x(2:end-1), d2 * ones(D1-2, 1)];

[F, ~] = randomfield(Corr, Mesh, 'data', data, 'nsamples', n);

F = F / max(abs(F(:)));

X = F';
Y = zeros(n,1);

for i = 1:num_comp
    start_index = (i-1)*n_each;
    Y((1:n_each)+start_index) = X((1:n_each)+start_index,:)*Beta{i}(:)+randn(n_each,1);
end

% Prepare the data input

data.imagesize = imagesize;
data.mask = ones(imagesize);
data.mask(1,:) = 0;
data.mask(end,:) = 0;
data.mask(:, 1) = 0;
data.mask(:, end) = 0;
data.X = X(:, data.mask(:) > 0);
data.Y = Y;

% Prepare the options input

options.q_max = 10;
options.band_width_s = 1:3;
options.lambda_s = 2.^(-5:5);
options.iter_max = 500;

[output, opt_model] = FMEM_v9(data,options);

