function [M,G,MT] = Gradient_renewal(type,m,grad,accumu_grad,Mt,lr_initial,iteration)
    switch type
        case 1
            if iteration < 99
                M = m-grad*lr_initial*(0.99^iteration);
            else
                M = m-grad*lr_initial*(0.99^100);
            end
            G = 0;
            MT = 0;
        case 2 
            gg = grad.*grad;
            G = gg+accumu_grad;
            Theta = -lr_initial*grad./sqrt(G+0.00001);
            M = m+Theta;
            MT = 0;
        case 3 
            G = 0.9*accumu_grad+0.1*grad.*grad;
            Theta = -lr_initial*grad./sqrt(G+0.00001);
            M = m+Theta;
            MT = 0;
        case 4 
            G = 0.9*accumu_grad-lr_initial*grad;
            M = m+G;
            MT = 0;
        case 5 
            MT = (0.9*Mt+0.1*grad)/(1-0.9^iteration);
            G = (0.99*accumu_grad+0.01*grad.*grad)/(1-0.99^iteration);
            Theta = -lr_initial*MT./sqrt(G+0.00001);
            M = m+Theta;
        case 6 
            G = accumu_grad*0.999+0.001*grad.*grad;
            MT = 0.9*Mt+0.1*grad;
            alpha = lr_initial*sqrt(1-0.999^iteration)/(1-0.9^iteration);
            M = m-alpha*MT./(sqrt(G)+1e-8);
        case 7 
            G = accumu_grad*0.99+0.01*grad.*grad;
            MT = 0.9*Mt+0.1*grad;
            alpha = lr_initial*sqrt(1-0.99^iteration)/(1-0.9^iteration);
            M = m-alpha*MT./(sqrt(G)+1e-8);
        case 8 
            G = 0;
            MT = 0;
            M = m-lr_initial*grad;
    end
end