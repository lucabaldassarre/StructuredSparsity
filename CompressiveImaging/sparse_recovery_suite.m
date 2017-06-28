%% CHOOSE ALGORITHMS

% Problem size - image side
% m = 2^15;
% m = 256;
m = 256;
% m = 256;
N = m^2;
image_choice = 4;
% subsamplings [6.25%, 12.5%, 25%, 0.1%]
subsampling = 1;
% DCT Sampling = 1, Hadamard Sampling = 2, Expander = 3;
sampling = 1;
% Representation - Wavelet = 1, Shearlet = 2;
representation = 1;

% Algorithms selection
bp_decomp = 1;
bp_pos_upper_bound_decomp = 0; % With positivity constraint

hgl_bp_decomp = 1;
hgl_bp_pos_decomp = 0; % With positivity constraint

family_lgl_decomp = 0;

bp_spgl = 0;

tv_decomp = 0;

display_results = 1;

if ~exist('Results', 'dir')
    mkdir('Results');
end


%%
addpath(genpath('Toolboxes')) % Add all toolboxes: RWT, FFST, Fast Hadamard transform, Spams
addpath('Helpers');
addpath('Operators');

%% Select images
switch m
   case {256, 512, 1024}
      switch image_choice
         case 1
            image = 'berlin_0256.png';
            image_str = 'berlin_0256';
         case 2
            image = 'earth_0256.png';
         case 3
            image = 'woman.jpg';
         case 4
            image = 'epfl_campus.jpg';
            image_str = 'epfl';
         case 5
            image = 'siddhartha.jpg';
         case 6
            image = 'mountains.jpg';
         otherwise
            error('Wrong selection of image');
      end
   case 2048
      switch image_choice
         case 1
            image = 'berlin_2048.png';
         case 2
            image = 'earth_2048.png';
         case 3
            image = 'siddhartha.jpg';
         case 4
            image = 'mountains.jpg';
         otherwise
            error('Wrong selection of image');
      end
   case 4096
      switch image_choice
         case 1
            image = 'berlin_4096.png';
            image_str = 'berlin_4096';
         case 2
            image = 'earth_4096.png';
            image_str = 'earth_4096';
         case 3
            image = 'mountains.jpg';
         otherwise
            error('Wrong selection of image');
      end
   case 8192
      switch image_choice
         case 1
            image = 'berlin_8192.png';
         case 2
            image = 'earth_8192.png';
         otherwise
            error('Wrong selection of image');
      end
   case 2^15
      image = 'world.png';
      image_str = 'world';
end

f = imread(['Images/',image]);
if numel(size(f)) == 3
   f = rgb2gray(f);
end
f = imresize(f,[m, m]);
f = double(f);

% figure, imagesc(f), colormap gray, axis image

%% REPRESENTATION

switch representation
   case 1 % WAVELET
      % define the function handles that compute 
      % the products by W (inverse DWT) and W' (DWT)
      wav = daubcqf(4);
      % Maximum level
      level = log2(m);
      W = @(x) midwt(x,wav,level); 
      WT= @(x) mdwt(x,wav,level);
      x_true = WT(f);
      x_true = x_true(:);
      
      representation_operator = @(x) reshape(W(reshape(x,[m,m])),[N,1]);
      representation_operator_trans = @(x) reshape(WT(reshape(x,[m,m])),[N,1]);
      
      str_repr = 'wav';
   case 2 % SHEARLETS
      tic;
      [ST, Psi] = shearletTransformSpect(f);
      toc
      x_true = ST(:);
      shearlet_levels = size(Psi,3);
      
      representation_operator = @(x) reshape(inverseShearletTransformSpect(reshape(x,[m,m,shearlet_levels]),Psi),[N,1]);
      representation_operator_trans = @(x) reshape(shearletTransformSpect(reshape(x,[m,m]), Psi),[N*shearlet_levels,1]);
      
      str_repr = 'shear';
end

N_effective = numel(x_true);

%% MEASUREMENTS on WAVELET BK

switch subsampling
   case 1
      str_sub = '06p25';
   case 2
      str_sub = '12p5';
   case 3
      str_sub = '25p';
end
switch m
   case 256
      str_size = '0256';
   case 512
      str_size = '0512';
   otherwise
      str_size = num2str(m);
end

load(sprintf('Indexes/inds_%s', str_sub));
str_ind = sprintf('ind_%s_%s', str_size,str_sub);
ind = eval(str_ind);
M = length(ind);

% Define operators
% Use Anders's operators
switch sampling
   case 1 % DCT SAMPLING
      measurement_forward = @(x) dct2fwd(reshape(x,size(f)),ind);
      measurement_backward = @(x) reshape(dct2adj(x,ind,m),[numel(f),1]);
      str_sam = 'dct';
   case 2 % HADAMARD SAMPLING
      measurement_forward = @(x) uhadtr2fwd(reshape(x,size(f)),ind);
      measurement_backward = @(x) reshape(uhadtr2adj(x,ind,m),[numel(f),1]);
      str_sam = 'had';
   case 3 % EXPANDER
      % Expander degree
      d = 8;
      disp('Creating Expander');
      time_A = tic;
      results.A = bin_mat(d,M,N); 
      results.time_A = toc(time_A)
      measurement_forward = @(x) 1/d*results.A*x;
      measurement_backward = @(x) 1/d*results.A'*x;
      str_sam = 'exp';
end

Phi.A  = @(x) measurement_forward(representation_operator(x));
Phi.At = @(x) representation_operator_trans(measurement_backward(x));

Phi.forward = Phi.A;
Phi.backward = Phi.At;

% Generate measurements
y = Phi.A(x_true);

fsave = sprintf('Results/%s_%s_%s_%s.mat',image_str,str_sub,str_sam, str_repr);

%% DECOMPOSITION COMMON PARAMETERS

decomp_RelTolX         = 1e-3;
results.decomp_RelTolX = decomp_RelTolX;
MaxIters               = 1e5;
Verbosity              = 3;
PrintStep              = 10;
saveHistMode           = 0;

%% BP - DECOMP

if bp_decomp
   disp('Running BP via the decomposition method');
   param_bp.MaxIters      = MaxIters;
   param_bp.Verbosity     = Verbosity; 
   param_bp.PrintStep     = PrintStep;
   param_bp.RelTolX       = decomp_RelTolX;
   param_bp.saveHistMode  = saveHistMode;
   time_bp = tic;
   [results.x_bp, results.out_bp] = lassoSolver('BP', Phi.forward, y, param_bp, 'AT', Phi.backward, 'nx', N_effective);
   results.time_bp = toc(time_bp)
   results.err_bp_l2 = norm(results.x_bp-x_true)/norm(x_true)
   % figure,
   % loglog(out_lasso.hist.rf_gap,'linewidth',2);
   % title('LASSO - Quoc''s','fontsize',16','interpreter','latex');
   image_bp = reshape(representation_operator(results.x_bp),[m,m]);
   results.psnr_bp = 20*log10(max(f(:))*sqrt(N)/norm(image_bp - f,'fro'));
%    figure, imshow(image_bp_decomp, [0 255]), axis image, colormap gray,  title('BP recovery')
   save(fsave,'results','-v7.3');
end

%% BP POS + UPPER BOUND - QUOC
if bp_pos_upper_bound_decomp
   disp('Running BP with positivity and upper bound via the decomposition method');
   param_bp_pos_decomp.MaxIters      = 1e5;
   param_bp_pos_decomp.Verbosity     = 3; 
   param_bp_pos_decomp.RelTolX       = decomp_RelTolX;
   param_bp_pos_decomp.saveHistMode  = 0;
   param_bp_pos.prox = @(x, t) sign(x).*max(abs(x) - 1/t, 0);
   param_bp_pos.pos_flag = true;
   param_bp_pos.upper_bound = 255;
   time_bp_pos_decomp = tic;
   bp_pos_forward = @(x) bp_pos_op_forward(x, measurement_forward, representation_operator, N);
   bp_pos_backward = @(x) bp_pos_op_backward(x, measurement_backward, representation_operator_trans, M);
   [x_bp_pos_decomp, out_bp_pos_decomp] = lassoSolver('BP', bp_pos_forward, [y; zeros(N,1)], param_bp_pos, 'AT', bp_pos_backward, 'nx', 2*N);
   x_bp_pos_decomp = x_bp_pos_decomp(1:N);
   time_bp_pos_decomp = toc(time_bp_pos_decomp)
   err_bp_pos_decomp_l2 = norm(x_bp_pos_decomp-x_true)/norm(x_true)
   image_bp_pos_decomp = W(reshape(x_bp_pos_decomp,size(f)));
end

%% BP - SPGL
if bp_spgl
   disp('Running BP via SPGL');
   spgl_opts = spgSetParms('optTol',1e-4,'iterations',1e4,'verbosity',1);
   spgl_OP = @(x,mode) create_OP(Phi, M, N, x, mode);
   sigma = 0;
   time_bp_spgl = tic;
   [x_bp_spgl, R_bp_spgl, G_bp_spgl, info_bp_spgl] = spg_bp( spgl_OP, y, spgl_opts);
   time_bp_spgl = toc(time_bp_spgl)
   nnz(x_bp_spgl)
   err_bp_spgl_l2 = norm(x_bp_spgl-x_true)/norm(x_true)
   image_bp_spgl = W(reshape(x_bp_spgl,size(f)));
end

%% Hierarchical Group Lasso - Define tree structure for running Julien's code for computing the prox
if hgl_bp_decomp || hgl_bp_pos_decomp || family_lgl_decomp
   % Load tree structure
   tree_structure_file_name = sprintf('tree_structure_%d.mat',m);
   if exist(tree_structure_file_name,'file')
      load(tree_structure_file_name);
   else
      D = 4;
      tree = create_tree_structure_HGL(m, D);
   end
   
   % Set spams parameters
   param_hgl_prox.regul = 'tree-l2';
   param_hgl_prox.verbose = false;
   param_hgl_prox.pos = false;
   param_hgl_prox.intercept = false;
   param_hgl_prox.num_threads = -1;
   param_hgl_prox.lambda = 1;

   setenv('MKL_NUM_THREADS','1')
   setenv('MKL_SERIAL','YES')
   setenv('MKL_DYNAMIC','NO')
   
end

%% HGL BP - DECOMP
if hgl_bp_decomp
   disp('Running HGL via the decomposition method');
   param_hgl.MaxIters      = 1e4;
   param_hgl.Verbosity     = 3; 
   param_hgl.RelTolX       = decomp_RelTolX;
%    param_hgl_bp.true_x        = x_true_wav;
   param_hgl.saveHistMode  = 1;
   param_hgl.prox = @(x, t) hgl_prox(x, t, tree, param_hgl_prox);
   time_hgl = tic;
   [results.x_hgl, results.out_hgl] = lassoSolver('BP', Phi.forward, y, param_hgl, 'AT', Phi.backward, 'nx', N);
   results.time_hgl = toc(time_hgl)
   results.err_hgl_l2 = norm(results.x_hgl - x_true)/norm(x_true)
   image_hgl = W(reshape(results.x_hgl,size(f)));
   results.psnr_hgl = 20*log10(max(f(:))*sqrt(N)/norm(image_hgl - f,'fro'));
%    figure; imshow(image_hgl_bp_decomp,[0 255]); axis image, colormap gray, title('Structured recovery')
   save(fsave,'results','-v7.3');
end

%% HGL BP POSITIVITY - DECOMP

if hgl_bp_pos_decomp
   disp('Running HGL with positivity and upper bound via the decomposition method');
   param_hgl_bp_pos.MaxIters      = 1e4;
   param_hgl_bp_pos.Verbosity     = 3; 
   param_hgl_bp_pos.RelTolX       = decomp_RelTolX;
   param_hgl_bp_pos.saveHistMode  = 1;
   param_hgl_bp_pos.prox = @(x, t) hgl_prox(x, t, tree, param_hgl_prox);
   param_hgl_bp_pos.pos_flag = true;
   param_hgl_bp_pos.upper_bound = 255;
   time_hgl_bp_pos_decomp = tic;
   hgl_bp_pos_forward = @(x) bp_pos_op_forward(x, measurement_forward, representation_operator, N);
   hgl_bp_pos_backward = @(x) bp_pos_op_backward(x, measurement_backward, representation_operator_trans, M);
   [x_hgl_bp_pos_decomp, out_hgl_bp_pos_decomp] = lassoSolver('BP', hgl_bp_pos_forward, [y; zeros(N,1)], param_hgl_bp_pos, 'AT', hgl_bp_pos_backward, 'nx', 2*N);
   x_hgl_bp_pos_decomp = x_hgl_bp_pos_decomp(1:N);
   time_hgl_bp_pos_decomp = toc(time_hgl_bp_pos_decomp)
   err_hgl_bp_pos_decomp_l2 = norm(x_hgl_bp_pos_decomp-x_true)/norm(x_true)
   image_hgl_bp_pos_decomp = W(reshape(x_hgl_bp_pos_decomp,size(f)));
end

%% FAMILY latent group lasso

if family_lgl_decomp
   disp('Running FAMILY LGL via the decomposition method');
   param_family_lgl.MaxIters      = 1e4;
   param_family_lgl.Verbosity     = 3; 
   param_family_lgl.PrintStep     = 50;
   param_family_lgl.RelTolX       = decomp_RelTolX;
   param_family_lgl.saveHistMode  = 1;
   C = create_family_matrix_for_computing_norms(m);
   B = create_family_structure(m);
   
   n_subtree = (m^2-1)/3;
   D = 4;
   depth = log(n_subtree*(D-1)+1)/log(D) - 1;
   group_weights = [];
   for kd = 1:depth
      group_weights = [group_weights; kd.^2*ones(4^(kd-1),1)];
   end
   group_weights = repmat(group_weights,3,1);
   group_weights = [1; group_weights];
   
   param_family_lgl.prox = @(x, t) prox_family_group(x, 1/t, C, group_weights);
   NN = 5/4*N - 1;
   
   forward = @(x) Phi.forward(B*x);
   backward = @(x) B'*Phi.backward(x);
   
   time_family_lgl = tic;
   [x_family_lgl, results.out_family_lgl] = lassoSolver('BP', forward, y, param_family_lgl, 'AT', backward, 'nx', NN);
   results.x_family_lgl = B*x_family_lgl;
   clear x_family_lgl;
   results.time_family_lgl = toc(time_family_lgl)
   results.err_family_lgl_l2 = norm(results.x_family_lgl - x_true)/norm(x_true)
   image_family_lgl = W(reshape(results.x_family_lgl,size(f)));
   
   results.psnr_family_lgl = 20*log10(max(f(:))*sqrt(N)/norm(image_family_lgl - f,'fro'));
%    figure, imshow(image_family_lgl,[0 255]), axis image, colormap gray, title('FAM recovery') 
   save(fsave,'results','-v7.3');
end

%% TV norm minimization - decomposition
if tv_decomp
   disp('Running the TV norm minimization via the decomposition method');
   xv0 = zeros(N,1);
   param_tv = [];
   param_tv.MaxIters      = 1e3;
   param_tv.Verbosity     = 3; 
   param_tv.RelTolX       = decomp_RelTolX;
   param_tv.saveHistMode  = 1;
   NN = (m-1)*m*2;
   time_tv_decomp = tic;
   [x_tv, results.out_tv] = tvNormSolver(N, M, NN, measurement_forward, measurement_backward, ...
                                       @(x) gradOperator_luca(x,m), @(x) divOperator_luca(x,m), y, param_tv, xv0, 0, 255);
   results.time_tv_decomp = toc(time_tv_decomp)
   results.x_tv = x_tv.x;
   results.err_tv_l2 = norm(results.x_tv-f(:))/norm(f(:))
   results.psnr_tv = 20*log10(max(f(:))*sqrt(N)/norm(results.x_tv-f(:),'fro'));
%    figure, imagesc(reshape(x_tv_decomp,[m,m])), axis image, colormap gray
   save(fsave,'results','-v7.3');
end

if false
if tv_decomp
   disp('Running the TV norm minimization via the decomposition method');
   xv0 = 125*ones(N,1);
   MW  = length(gradOperator_luca(xv0, m));
   param_tv = [];
   param_tv.MaxIters      = 1e3;
   param_tv.Verbosity     = 3; 
   param_tv.RelTolX       = decomp_RelTolX;
   param_tv.saveHistMode  = 1;
   time_tv_decomp = tic;
   %y1 = measurement_forward(f);
%    [x_tv_decomp, tv_output] = tvNormSolver(N, M, MW, measurement_forward, measurement_backward, ...
%                                        @(x) gradOperator_luca(x,m), @(x) divOperator_luca(x,m), y, param_tv, xv0, 0, 255);
%    x_tv_decomp = x_tv_decomp.x;
   
   [results.x_tv, results.out_tv] = tvNormProxSolver(measurement_forward, measurement_backward, y, 1, 0, 255, param_tv, xv0); 
   results.time_tv = toc(time_tv_decomp)
   
   results.err_tv_l2 = norm(results.x_tv-f(:))/norm(f(:))
   results.psnr_tv = 20*log10(max(f(:))*sqrt(N)/norm(results.x_tv-f(:),'fro'));
%    figure, imshow(reshape(results.x_tv,[m,m]),[0 255]), axis image, colormap gray, title('TV recovery') 
end

end

%%
if display_results
   cont = 0;
   clear time err
   labels = {};

   if hgl_bp_decomp
      cont = cont + 1;
      time(cont) = results.time_hgl;
      err(cont) = results.err_hgl_l2;
      labels = [labels; 'HGL'];
   end
   if hgl_bp_pos_decomp
      cont = cont + 1;
      time(cont) = time_hgl_bp_pos_decomp;
      err(cont) = err_hgl_bp_pos_decomp_l2;
      labels = [labels; 'HGL POS'];
   end
   if bp_decomp
      cont = cont + 1;
      time(cont) = results.time_bp;
      err(cont) = results.err_bp_l2;
      labels = [labels; 'BP'];
   end
   if bp_pos_upper_bound_decomp
      cont = cont + 1;
      time(cont) = time_bp_pos_decomp;
      err(cont) = err_bp_pos_decomp_l2;
      labels = [labels; 'BP POS UB'];
   end
   if bp_spgl
      cont = cont + 1;
      time(cont) = time_bp_spgl;
      err(cont) = err_bp_spgl_l2;
      labels = [labels; 'BP spgl'];
   end
   if tv_decomp
      cont = cont + 1;
      time(cont) = results.time_tv;
      err(cont) = results.err_tv_l2;
      labels = [labels; 'TV'];
   end
   if family_lgl_decomp
      cont = cont + 1;
      time(cont) = results.time_family_lgl;
      err(cont) = results.err_family_lgl_l2;
      labels = [labels; 'FAM LGL'];
   end

   ref_time = min(time);
   figure, 
   subplot(121),
   bar(time/ref_time);
   title(sprintf('Computation times - %s',image),'fontsize',20,'fontname','garamond', 'Interpreter', 'none');
   ylabel('Relative computation time to best','fontsize',16,'interpreter','latex');
   set(gca,'XTickLabel',labels,'fontsize',12);
   subplot(122),
   ref_err = min(err);

%    bar((err/ref_err-1)*100);
   bar(err);
   title(sprintf('l2 errors - %s',image),'fontsize',20,'fontname','garamond', 'Interpreter', 'none');
   ylabel('Difference to best value (in perc.)','fontsize',16,'interpreter','latex');
   set(gca,'XTickLabel',labels,'fontsize',12);

end

clear tree
clear x_true_wav
clear B C group_weights x_family_lgl_decomp
save(fsave,'results','-v7.3');