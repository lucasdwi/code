% go to where output data files exist
cd('G:\GreenLab\data\angela\output\')
% pre-allocate arrays
[thisAUC,thisAUCR] = deal(zeros(1,100));
[x,y,xR,yR] = deal(cell(1,100));
% cycle through all 100 files, load them, apply model to test set to get
% AUC
for ii = 1:100
   % Load iith file
    load(['miaData',num2str(ii),'.mat'])
   % Apply mdl to test set
   pred = cvglmnetPredict(acc{1}.mdl{1},hist.cfg.naive.testX);
   % Compare prediction to actual Y; generates AUC; uses thisAUC to not
   % overwrite auc variable
   [x{ii},y{ii},~,thisAUC(ii)] = perfcurve(hist.cfg.naive.testY,...
       pred,1);
   % Do the same as above with randomized data
   predR = cvglmnetPredict(accR{1}.mdl{1},histR.cfg.naive.testX);
   [xR{ii},yR{ii},~,thisAUCR(ii)] = perfcurve(histR.cfg.naive.testY,...
       predR,1);
end