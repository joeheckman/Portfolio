%%%CODE SAMPLE. THIS CREATES A NETWORK OF COUPLED OSCILLATORS USING THE 'TARGET FREQUENCY GAP' PARAMETER


sizeMAT=100;
savepath = 'D:\Joe\Zauberbaum\AdjustFrequencyGapNetworks\ripser_test_4.5\';
cores = 16;
num_iterations=cores * 5;  %one set of oscillators
TargetNumberofEdges = 20 * sizeMAT;
%  target_FG = [0.38:.01:.5]
target_FG = [.34:.01:.45];
for FG = 1:length(target_FG)

    W_save = zeros(num_iterations,sizeMAT);
    MAT = zeros(sizeMAT,sizeMAT,num_iterations);
    for iter_no=1:num_iterations
        L = 0;
        total = 0;
        w0=0;     
        W = 1*rand(1,sizeMAT);
%         W = w0*ones(1,sizeMAT)+range*rand(1,sizeMAT)
        W_save(iter_no,:)=W; %%1x500 with natural freq [1 2]
        while L < TargetNumberofEdges
%             disp(L)
            node1 = (randi([1 sizeMAT]));
            node2 = (randi([1 sizeMAT]));
            if node1 == node2||MAT(node1,node2,iter_no)>0
                continue
            end
            NodeDifference = abs(W_save(iter_no,node2)-W_save(iter_no,node1));
            total_temp  = total + NodeDifference;
            if total_temp / L > target_FG(FG)
                MAT(node1,node2,iter_no) = 1;
                MAT(node2,node1,iter_no) = 1;
                L = L + 1;
                total = total_temp;
            end
            
        end
        total_temp / L
        x = sum(MAT(:,:,iter_no),1);
        disp(['mean degree for iteration ' num2str(iter_no) ' is ' num2str(mean(x))]);
   
    end

   
        save(sprintf('%sNaturalFrequencies_and_MAT_FreqGap=%.3f.mat', savepath, target_FG(FG)),'W_save', 'MAT')



    %%figures for the last iteration.
    %%so you know every iteration works lol
    figure();
    title(['FIGURE FOR FG = ' num2str(target_FG(FG))])
    for iter_no = 1:cores
        subaxis(2,cores,iter_no,'SpacingHoriz',0,'Margin',0)
        scatter(W_save(iter_no,:),sum(MAT(:,:,iter_no)),Marker = '.', MarkerFaceColor="g")
        hold on
        xlim([w0+.0001 w0+1])
        title([' CDF =' num2str(target_FG(FG))])
        % xlabel('Natural Frequency')
        % ylabel('Degree')
        % set(gca,'fontname', 'times new roman')
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);

        %subplot(2,num_iterations,num_iterations+iter_no)
        subaxis(2,cores,iter_no+cores)
        %make frequency to average frequency plot
        d = zeros(size(MAT,1),1);
        for a = 1:size(MAT,1)
            b = MAT(a,:,iter_no);
            c = b .* W_save(iter_no,:);
            d(a) = sum(c) / sum(b);
        end
        scatter(W_save(iter_no,:),d, marker ='.')
        xlim([w0+.0001 w0+1])
        % title(['FA ' num2str(target_FG(FG))])
        % xlabel('Natural Frequency')
        % ylabel('Neighborhood Frequency')
        % set(gca, 'fontname', 'times new roman')
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);

    end


%%MAKE FIGURES TO CHECK IF IT LOOKS LIKE THE PAPER'S FIGURES 


     saveas(gcf,sprintf('%sfigure_target_FG=%.3f.png',savepath, target_FG(FG)))
end

     title(['FIGURE FOR FG = ' num2str(target_FG(FG))])
  
         subplot(length(target_FG),2,(2*FG)-1)
         scatter(W_save(iter_no,:),sum(MAT(:,:,iter_no)),Marker = '.', MarkerFaceColor="g")
         hold on
         xlim([1.001 2])
         title([' CDF =' num2str(target_FG(FG))])
         xlabel('Natural Frequency')
         ylabel('Degree')
         set(gca,'fontname', 'times new roman')
         set(gca, 'XTick', []);
         set(gca, 'YTick', []);
         set(gca, 'XTickLabel', []);
         set(gca, 'YTickLabel', []);
 
         %subplot(2,num_iterations,num_iterations+iter_no)
         subplot(length(target_FG),2,2*FG)
         %make frequency to average frequency plot
         d = zeros(size(MAT,1),1);
         for a = 1:size(MAT,1)
             b = MAT(a,:,iter_no);
             c = b .* W_save(iter_no,:);
             d(a) = sum(c) / sum(b);
         end
         scatter(W_save(iter_no,:),d, marker ='.')
         xlim([1.001 2])
         title(['FA ' num2str(target_FG(FG))])
         xlabel('Natural Frequency')
         ylabel('Neighborhood Frequency')
         set(gca, 'fontname', 'times new roman')
         set(gca, 'XTick', []);
         set(gca, 'YTick', []);
         set(gca, 'XTickLabel', []);
         set(gca, 'YTickLabel', []);

     saveas(gcf,sprintf('%sfigure_target_FG=%.3f.png',savepath, target_FG(FG)))

