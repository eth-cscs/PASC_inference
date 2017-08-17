clear all

addpath('common')
addpath('compute')
addpath('compute/wavelet')
addpath('compute_optimal')

plot_data = true;
plot_solution = true;

compute_fourier = true;
compute_fourierl2 = true;
compute_fouriersobolev = true;
compute_fouriertvr = true;
compute_gauss = true;
compute_median = true;
compute_wiener = true;
compute_denS2D = true;
compute_denC2D = true;
compute_denR2D = true;


Sigmaid = 7;

Sigma_values = [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256];
sigma_value = Sigma_values(Sigmaid+1); 
n_mc=0;
name='usi';

C_true = loadbin(['data/image_' name '/solution.bin'])';
C_data = loadbin(['data/image_' name '/' name '_id' num2str(n_mc) '_idSigma' num2str(Sigmaid) '.bin'])';

width = 1024;
height = 512;

image_solution = reshape(C_true,width,height)';
image_data = reshape(C_data,width,height)';

if plot_data
    figure
    hold on
    %showimage(img, mytitle, subx, suby, subnmb, scaleimg)
    showimage(image_solution, 'solution', 2, 2, 1, false)
    showimage(image_data, 'data', 2, 2, 2, false)
    showimage(image_solution, 'solution scaled', 2, 2, 3, true)
    showimage(image_data, 'data scaled', 2, 2, 4, true)
    hold off
end

% compute snr
[snrvalue, snrvalue_sigma] = snr(image_data, image_solution, sigma_value);
disp(['SNR = std(solution)/std(data) = ' num2str(snrvalue) ])
disp(['SNR = std(solution)/sigma     = ' num2str(snrvalue_sigma) ])

used_methods_counter = 0;

% compute fourier
if compute_fourier 
    disp(' - computing fourier')
    S=[5 20 30 40 60 80];
    [F_fourier, Ferr_fourier, Sbest_fourier] = optimal_fourier(image_data,image_solution, S);
    disp(['   - best for S = ' num2str(Sbest_fourier)])
    disp(['   - error      = ' num2str(Ferr_fourier)])

    used_methods_counter = used_methods_counter + 1;
end
    
    
% compute fourierl2
if compute_fourierl2
    disp(' - computing fourier L2')
    Lambda=[0.1 1 10 50 100];
    S=[5 20 30 40 60 80];
    [F_fourierl2, Ferr_fourierl2, Sbest_fourierl2, Lambdabest_fourierl2] = optimal_fourierl2(image_data,image_solution, S, Lambda);
    disp(['   - best for S      = ' num2str(Sbest_fourierl2)])
    disp(['   - best for Lambda = ' num2str(Lambdabest_fourierl2)])
    disp(['   - error           = ' num2str(Ferr_fourierl2)])

    used_methods_counter = used_methods_counter + 1;
end

% compute fouriersobolev
if compute_fouriersobolev
    disp(' - computing fourier sobolev')

    Lambda=[0.1 1 10 50 100];
    S=[5 20 30 40 60 80];

    [F_fouriersobolev, Ferr_fouriersobolev, Sbest_fouriersobolev, Lambdabest_fouriersobolev] = optimal_fouriersobolev(image_data,image_solution, S, Lambda);
    disp(['   - best for S      = ' num2str(Sbest_fouriersobolev)])
    disp(['   - best for Lambda = ' num2str(Lambdabest_fouriersobolev)])
    disp(['   - error           = ' num2str(Ferr_fouriersobolev)])

    used_methods_counter = used_methods_counter + 1;
end

% compute fouriertvr
if compute_fouriertvr
    disp(' - computing fourier-vtr')    

    Lambda=[0.1 100 1000];
    S=[5 20 50];
    Epsilon = [0.001 0.01 0.1];

    [F_fouriertvr, Ferr_fouriertvr, Sbest_fouriertvr, Lambdabest_fouriertvr, Epsilonbest_fouriertvr] = optimal_fouriertvr(image_data,image_solution, S, Lambda, Epsilon);
    disp(['   - best for S       = ' num2str(Sbest_fouriertvr)])
    disp(['   - best for Lambda  = ' num2str(Lambdabest_fouriertvr)])
    disp(['   - best for Epsilon = ' num2str(Epsilonbest_fouriertvr)])
    disp(['   - error            = ' num2str(Ferr_fouriertvr)])

    used_methods_counter = used_methods_counter + 1;
end

% compute gauss
if compute_gauss 
    disp(' - computing gauss')

    N=[5 11 13 15];
    Sigma=[2 3 4 6 8];

    [F_gauss, Ferr_gauss, Nbest_gauss, Sigmabest_gauss] = optimal_gauss(image_data,image_solution, N, Sigma);
    disp(['   - best for N     = ' num2str(Nbest_gauss)])
    disp(['   - best for Sigma = ' num2str(Sigmabest_gauss)])
    disp(['   - error          = ' num2str(Ferr_gauss)])

    used_methods_counter = used_methods_counter + 1;
end

% compute median
if compute_median
    disp(' - computing median')

    N=[5 10 20 30];

    [F_median, Ferr_median, Nbest_median] = optimal_median(image_data,image_solution, N);
    disp(['   - best for N     = ' num2str(Nbest_median)])
    disp(['   - error          = ' num2str(Ferr_median)])

    used_methods_counter = used_methods_counter + 1;
end

% compute wiener
if compute_wiener
    disp(' - computing wiener')

    N=[5 10 20 30];

    [F_wiener, Ferr_wiener, Nbest_wiener] = optimal_wiener(image_data,image_solution, N);
    disp(['   - best for N     = ' num2str(Nbest_wiener)])
    disp(['   - error          = ' num2str(Ferr_wiener)])

    used_methods_counter = used_methods_counter + 1;
end

% compute denS2D
if compute_denS2D
    disp(' - computing denS2D')

    T=[0.01 0.05 0.1 0.5 1 2 5 10 15 20 25 30];

    [F_denS2D, Ferr_denS2D, Tbest_denS2D] = optimal_denS2D(image_data,image_solution, T);
    disp(['   - best for T     = ' num2str(Tbest_denS2D)])
    disp(['   - error          = ' num2str(Ferr_denS2D)])

    used_methods_counter = used_methods_counter + 1;
end

% compute denC2D
if compute_denC2D
    disp(' - computing denC2D')

    T=[0.01 0.05 0.1 0.5 1 2 5 10 15 20 25 30];

    [F_denC2D, Ferr_denC2D, Tbest_denC2D] = optimal_denC2D(image_data,image_solution, T);
    disp(['   - best for T     = ' num2str(Tbest_denC2D)])
    disp(['   - error          = ' num2str(Ferr_denC2D)])

    used_methods_counter = used_methods_counter + 1;
end

% compute denR2D
if compute_denR2D
    disp(' - computing compute_denR2D')

    T=[0.01 0.05 0.1 0.5 1 2 5 10 15 20 25 30];

    [F_denR2D, Ferr_denR2D, Tbest_denR2D] = optimal_denR2D(image_data,image_solution, T);
    disp(['   - best for T     = ' num2str(Tbest_denR2D)])
    disp(['   - error          = ' num2str(Ferr_denR2D)])

    used_methods_counter = used_methods_counter + 1;
end





% plot best results
if plot_solution
    figure
    hold on
    
    used_methods_counter_plot = 1;
    
    cols = 2;
    rows = ceil(used_methods_counter/cols);
    
    %showimage(img, mytitle, subx, suby, subnmb, scaleimg)
    if compute_fourier
        showimage(F_fourier,   ['fourier, ' num2str(Ferr_fourier)], rows, cols, used_methods_counter_plot, true)
        used_methods_counter_plot = used_methods_counter_plot + 1;
    end

    if compute_fourierl2
        showimage(F_fourierl2,   ['fourier l2, ' num2str(Ferr_fourierl2)], rows, cols, used_methods_counter_plot, true)
        used_methods_counter_plot = used_methods_counter_plot + 1;
    end
    
    if compute_fouriersobolev
        showimage(F_fouriersobolev, ['fourier Sobolev, ' num2str(Ferr_fouriersobolev)], rows, cols, used_methods_counter_plot, true)
        used_methods_counter_plot = used_methods_counter_plot + 1;
    end

    if compute_fouriertvr
        showimage(F_fouriertvr, ['fourier TVR, ' num2str(Ferr_fouriertvr)], rows, cols, used_methods_counter_plot, true)
        used_methods_counter_plot = used_methods_counter_plot + 1;
    end

    if compute_gauss
        showimage(F_gauss(5:end-5,5:end-5), ['gauss, ' num2str(Ferr_gauss)], rows, cols, used_methods_counter_plot, true)
        used_methods_counter_plot = used_methods_counter_plot + 1;
    end
    
    if compute_median
        showimage(F_median(5:end-5,5:end-5), ['median, ' num2str(Ferr_median)], rows, cols, used_methods_counter_plot, true)
        used_methods_counter_plot = used_methods_counter_plot + 1;
    end
    
    if compute_wiener
        showimage(F_wiener(5:end-5,5:end-5), ['wiener, ' num2str(Ferr_wiener)], rows, cols, used_methods_counter_plot, true)
        used_methods_counter_plot = used_methods_counter_plot + 1;
    end

    if compute_denS2D
        showimage(F_denS2D, ['denS2D, ' num2str(Ferr_denS2D)], rows, cols, used_methods_counter_plot, true)
        used_methods_counter_plot = used_methods_counter_plot + 1;
    end
    
    if compute_denC2D
        showimage(F_denC2D, ['denC2D, ' num2str(Ferr_denC2D)], rows, cols, used_methods_counter_plot, true)
        used_methods_counter_plot = used_methods_counter_plot + 1;
    end

    if compute_denR2D
        showimage(F_denR2D, ['denR2D, ' num2str(Ferr_denR2D)], rows, cols, used_methods_counter_plot, true)
        used_methods_counter_plot = used_methods_counter_plot + 1;
    end
    
    hold off
end



