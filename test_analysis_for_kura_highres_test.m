---CODE SAMPLE. THIS SCRIPT PLOTS THE AVERAGE VARIANCE IN ORDER PARAMETER CURVES OF A LARGE NUMBER OF KURAMOTO MODEL SIMULATIONS.---

    %%find the maximum index of OR_STD as critical point of syncrhonization
    for m = 1:size(Or_std_data_reshaped,1)
        [maxValue(m), Index(m)] =max(Or_std_data_reshaped(m,:));
        disp(['maximum orstd value at ' num2str(Index(m)) ' value is ' ...
            num2str(max(Or_std_data_reshaped(m,:)))]);

    end
    %%align curves for figure
    shift = Index-max(Index);
    Extended_mat = [Or_std_data_reshaped, zeros(size(Or_std_data_reshaped,1), abs(min(shift)))];
    for o = 1:size(Or_std_data_reshaped,1)
        Extended_mat(o,:) = circshift(Extended_mat(o,:),abs(shift(o)));
        [maxValueShift(o), Index2(o)] = max(Extended_mat(o,:));
        disp(['maximum aligned orstd value at ' num2str(Index2(o))])
    end

 x = mean(Or_std_data_reshaped, 1);
 y = mean(Extended_mat,1);

 figure();
 plot(x)

 figure();
 plot(y(100:2000))


figure();
 for i = 1:size(Or_data_reshaped,1)
     hold on 
     plot(Or_data_reshaped(i,:))
 end
   