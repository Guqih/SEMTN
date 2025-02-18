function [X,Type] = Activate(x,type)
    switch type
        case 1 
            X = max(0,x);
        case 2 
            X = max(0.01*x,x);
        case 3 
            X = 1./(1+exp(-x));
        case 4 
            X = tanh(x);
    end
    Type=type;
end