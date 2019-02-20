function SadSeeBlue(X)

if nargin==0
    load SadSeeBlue
    X = Study1;
end

% Bar plot
figure
barploterr([mean(X(X(:,1)==0,2:3)); mean(X(X(:,1)==1,2:3))]', [sem(X(X(:,1)==0,2:3)); sem(X(X(:,1)==1,2:3))]',2); 
colormap([1 .5 0; 0 .5 1]); 
if nargin==0 
    ylim([.8 .94]);
end
set(gca, 'fontsize', 12, 'xtick', 1:2, 'xticklabel', {'Blue-Yellow' 'Red-Green'}); 
ylabel('Accuracy'); 
legend({'Control' 'Sadness'});

% Scatter plot
figure; hold on
scatter(randn(sum(X(:,1)==0),1)*.02+0.9, X(X(:,1)==0,2), 80, [1 .5 0], 'filled'); 
scatter(randn(sum(X(:,1)==1),1)*.02+1.1, X(X(:,1)==1,2), 80, [0 .5 1], 'filled'); 
scatter(randn(sum(X(:,1)==0),1)*.02+1.9, X(X(:,1)==0,3), 80, [1 .5 0], 'filled'); 
scatter(randn(sum(X(:,1)==1),1)*.02+2.1, X(X(:,1)==1,3), 80, [0 .5 1], 'filled'); 
xlim([0.5 2.5]);
set(gca, 'fontsize', 12, 'xtick', 1:2, 'xticklabel', {'Blue-Yellow' 'Red-Green'}); 
ylabel('Accuracy'); 
legend({'Control' 'Sadness'});

% Box & whisker plot
figure; hold on
boxplot([X(X(:,1)==0,2); X(X(:,1)==1,2); X(X(:,1)==0,3); X(X(:,1)==1,3)], [1*ones(sum(X(:,1)==0),1); 2*ones(sum(X(:,1)==1),1); 3*ones(sum(X(:,1)==0),1); 4*ones(sum(X(:,1)==1),1)],'notch','on'); 
set(gca, 'fontsize', 12, 'xtick', 1:4, 'xticklabel', {'BY-Ctrl' 'BY-Sad' 'RG-Ctrl' 'RG-Sad'}); 
ylabel('Accuracy'); 

% Cat-eye plot
figure; hold on
h1 = cateye(X(X(:,1)==0,2:3), [0.7 2.7], [1 .5 0]);
h2 = cateye(X(X(:,1)==1,2:3), [1.3 3.3], [0 .5 1]);
% h1 = cateye(X(X(:,1)==0,2:3), [0.7 2.7], [1 .5 0], 2, true);
% h2 = cateye(X(X(:,1)==1,2:3), [1.3 3.3], [0 .5 1], 2, true);
xlim([0 4]);
set(gca, 'fontsize', 12, 'xtick', [1 3], 'xticklabel', {'Blue-Yellow' 'Red-Green'}); 
ylabel('Accuracy'); 
legend([h1(1) h2(1)],{'Control' 'Sadness'});
