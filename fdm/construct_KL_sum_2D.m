function [kl_sum] = construct_KL_sum_2D(x, y, rand_tensor_list, eta_x, eta_y)
    sigma_x = 1;
    sigma_y = 1;
    omega_x = calculate_omega_based_on_eta(eta_x);
    omega_y = calculate_omega_based_on_eta(eta_y);
    
    % omegaX omegaY are vectors
    lambda_x = 2.0 * eta_x * sigma_x ./ (1.0 + (eta_x * omega_x).^2);
    lambda_y = 2.0 * eta_y * sigma_y ./ (1.0 + (eta_y * omega_y).^2);
    
    kl_sum = 0;
    
    for i=1:6
        kl_sum = kl_sum + rand_tensor_list(i) * sqrt(lambda_x(i)) * sqrt(lambda_y(i)) * (eta_x * omega_x(i) * cos(omega_x(i) * x) + sin(omega_x(i) * x)) * (eta_y * omega_y(i) * cos(omega_y(i) * y) + sin(omega_y(i) * y));  
    end
end