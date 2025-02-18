function [X] = Activate_grad(x,type)
    switch type
        case 1 
            X = ones(size(x));
            X(x<=0) = 0;
        case 2 
            X = ones(size(x));
            X(x<=0) = 0.01;
        case 3 
            tem = 1./(1+exp(-x));
            X = tem.*(1-tem);
        case 4 
            X = sech(x).^2;
    end
end