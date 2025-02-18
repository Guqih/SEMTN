clc; clear;

Order_highest = 2;
lr_init = 0.03;
norm_sf = 0.01;   
iteration_max = 1000;
Acti_type = 4;
expe_num = 100;
k = 5; 
type_num = 7;

load Data_AWGN_2dB_70000
format long

total_num = size(feature,2); 
train_num = total_num*(k-1)/k;
len_data = size(feature,1);  
feature_size = len_data; 

expan_size = 1;
for i = 1:Order_highest
    expan_size = expan_size+prod(feature_size:(feature_size+i-1))/prod(1:i); 
end
hidden_size = feature_size; 
output_size = type_num; 

onehot = zeros(total_num,type_num);
identity_matrix = eye(type_num);
for i = 1:total_num
    onehot(i,:) = identity_matrix(Labels(i),:);
end
onehot = onehot';

confusion_matrices = cell(1,expe_num); 
TrainLabels_valid_Cell = cell(1,expe_num);
TestLabels_Cell = cell(1,expe_num);
train_predictedLabels_Cell = cell(1,expe_num);
predictedLabels_Cell = cell(1,expe_num);
time = zeros(1,expe_num);

for expe = 1:expe_num
        tic;
        cv = cvpartition(Labels','kfold',k); 
        i = 1;
        init_Train_Data = feature(:,cv.training(i));
        init_Test_Data = feature(:,cv.test(i));
        Tmin = min(init_Train_Data,[],2); 
        Tmax = max(init_Train_Data,[],2); 
        init_Train_Data = (init_Train_Data-Tmin)./(Tmax-Tmin)*(1-norm_sf)+norm_sf; 
        init_Test_Data = (init_Test_Data-Tmin)./(Tmax-Tmin)*(1-norm_sf)+norm_sf; 

        TrainData = (Taylor_expan(init_Train_Data,Order_highest));
        TrainData_valid = TrainData;
        TestData = (Taylor_expan(init_Test_Data,Order_highest));
        TrainLabels = (onehot(:,cv.training(i))); 
        TrainLabels_valid = (Labels(:,cv.training(i))); 
        TestLabels = (Labels(:,cv.test(i)));
         
        Wt = (Initialization(expan_size,expan_size,7)); 
        La = (Initialization(hidden_size,expan_size,7));
        K = (Initialization(output_size,hidden_size,7));
        MWt = zeros(size(Wt));
        MLa = zeros(size(La));
        MK = zeros(size(K));  
        GWt = zeros(size(Wt));
        GLa = zeros(size(La));
        GK = zeros(size(K)); 
        
        train_accuracy = 0;
        test_accuracy = 0; 
        
        for iteration = 1:iteration_max 
            We = Wt*TrainData; 
            S = 1./(1+exp(-We)); 
            Z = TrainData.*S; 
            Y = La*Z;  
            hidden = Activate(Y,Acti_type); 
            J = K*hidden; 
            output = softmax(J); 
            d_softmax = output-TrainLabels; 
            d_K = d_softmax*hidden'/train_num; 
            d_hidden = K'*d_softmax; 
            d_Y = d_hidden.*Activate_grad(Y,Acti_type); 
            d_La = d_Y*Z'/train_num; 
            d_Z = La'*d_Y; 
            d_S = d_Z.*TrainData; 
            d_We = d_S.*S.*(1-S); 
            d_Wt = d_We*TrainData'/train_num; 
            [Wt,GWt,MWt] = Gradient_renewal(7,Wt,d_Wt,GWt,MWt,lr_init,iteration);
            [La,GLa,MLa] = Gradient_renewal(7,La,d_La,GLa,MLa,lr_init,iteration); 
            [K,GK,MK] = Gradient_renewal(7,K,d_K,GK,MK,lr_init,iteration);

            if iteration == iteration_max
                We = Wt*TrainData_valid; 
                S = 1./(1+exp(-We)); 
                Z = TrainData_valid.*S; 
                Y = La*Z; 
                hidden = Activate(Y,Acti_type); 
                J = K*hidden; 
                output = softmax(J); 
                [~,train_predictedLabels] = max(output); 
                train_accuracy = mean(train_predictedLabels==TrainLabels_valid);
                We = Wt*TestData ; 
                S = 1./(1+exp(-We)); 
                Z = TestData.*S; 
                Y = La*Z; 
                hidden = Activate(Y,Acti_type); 
                J = K*hidden; 
                output = softmax(J);
                [~,predictedLabels] = max(output); 
                test_accuracy = mean(predictedLabels==TestLabels);
                TestLabels_onehot = (zeros(type_num,length(TestLabels))); 
                for j = 1:length(TestLabels)
                    TestLabels_onehot(TestLabels(j), j)=1; 
                end
                cross_entropy_loss_test = -sum(sum(TestLabels_onehot.*log(output + 0.00001)))/length(TestLabels); 
                fprintf('Train accuracy:%.6f,Test accuracy:%.6f\n',train_accuracy*100,test_accuracy*100);
            end
        end
    elapsed_time = toc;
    confusion_matrices{expe} = confusionmat(TestLabels,predictedLabels);
    time(expe) = elapsed_time;
    train_predictedLabels_Cell{expe} = train_predictedLabels;
    TrainLabels_valid_Cell{expe} = TrainLabels_valid;
    predictedLabels_Cell{expe} = predictedLabels;
    TestLabels_Cell{expe} = TestLabels;
end

accuracy_ten = zeros(1,expe_num);
precision_ten = zeros(1,expe_num);
specificity_ten = zeros(1,expe_num);
F1_ten = zeros(1,expe_num);
recall_ten = zeros(1,expe_num);

for idx = 1:expe_num
    matrix = confusion_matrices{idx};  
    numClasses = size(matrix,1);  
    TP = zeros(numClasses,1);
    TN = zeros(numClasses,1);
    FP = zeros(numClasses,1);
    FN = zeros(numClasses,1);
    accuracy = zeros(numClasses,1);
    specificity = zeros(numClasses,1);
    recall = zeros(numClasses,1);
    precision = zeros(numClasses,1);
    F1 = zeros(numClasses,1);
    
    for i = 1:numClasses
        TP(i) = matrix(i,i);
        FN(i) = sum(matrix(i,:))-TP(i);
        FP(i) = sum(matrix(:, i))-TP(i);
        TN(i) = sum(matrix(:))-TP(i)-FP(i)-FN(i);
        accuracy(i) = TP(i)/(TP(i)+FN(i));
        specificity(i) = TN(i)/(TN(i)+FP(i));
        recall(i) = accuracy(i); 
        precision(i) = TP(i)/(TP(i)+FP(i));
        F1(i) = (2 *precision(i)*recall(i))/(precision(i)+recall(i));
    end
   
    avg_accuracy = mean(accuracy);
    avg_specificity = mean(specificity);
    avg_recall = mean(recall);
    avg_precision = mean(precision);
    avg_F1 = mean(F1);
    accuracy_ten(idx) = avg_accuracy;
    precision_ten(idx) = avg_precision;
    specificity_ten(idx) = avg_specificity;
    F1_ten(idx) = avg_F1;
    recall_ten(idx) = avg_recall;
end

Accuracy_Mean = mean(accuracy_ten);Accuracy_STD = std(accuracy_ten);
Precision_Mean = mean(precision_ten);Precision_STD = std(precision_ten);
Specificity_Mean = mean(specificity_ten);Specificity_STD = std(specificity_ten);
F1_Mean = mean(F1_ten);F1_STD = std(F1_ten);
Recall_Mean = mean(recall_ten);Recall_STD = std(recall_ten);
Training_Time_Mean = mean(time); Training_Time_STD = std(time);
fprintf('Accuracy_Mean: %.6f ± %.6f\n', Accuracy_Mean, Accuracy_STD);
fprintf('Precision_Mean: %.6f ± %.6f\n', Precision_Mean, Precision_STD);
fprintf('Specificity_Mean: %.6f ± %.6f\n', Specificity_Mean, Specificity_STD);
fprintf('F1_Mean: %.6f ± %.6f\n', F1_Mean, F1_STD);
fprintf('Recall_Mean: %.6f ± %.6f\n', Recall_Mean, Recall_STD);
fprintf('Training_Time_Mean: %.6f ± %.6f\n', Training_Time_Mean, Training_Time_STD);