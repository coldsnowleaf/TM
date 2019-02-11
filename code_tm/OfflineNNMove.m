clear
load OfflineTrain.mat
global ErrorNN;
global ErrorArrayNN;
global PathNum;
global NNnum;
global ReScale;

% RSSScaleMax = 170;
% RSSScaleMin = 40;

AbRSSNorm = [];
AbRSSFitNorm = [];
RePoNorm = [];
AbPoNorm = [];
AbRSSNormTran = [];
AbRSSFitNormTran = [];
RePoNormTran = [];
AbPoNormTran = [];
for i = 1 : PathNum
     %AbRSSNorm = [AbRSSNorm; (TrainRSS(i).AbRSS - RSSScaleMin)/ (RSSScaleMax - RSSScaleMin)];
     AbRSSNorm = [AbRSSNorm; TrainRSS(i).AbRSS];
     %AbRSSNorm = [AbRSSNorm; (TrainRSS(i).AbRSSPCA - repmat(RSSScaleMin, size(TrainRSS(i).AbRSSPCA, 1), size(TrainRSS(i).AbRSSPCA, 2))) ./ (repmat(RSSScaleMax - RSSScaleMin, size(TrainRSS(i).AbRSSPCA, 1), size(TrainRSS(i).AbRSSPCA, 2)))];
     %AbRSSFitNorm = [AbRSSFitNorm; (TrainRSS(i).AbRSSFit - repmat(RSSScaleMin, size(TrainRSS(i).AbRSS, 1), size(TrainRSS(i).AbRSS, 2))) ./ (repmat(RSSScaleMax - RSSScaleMin, size(TrainRSS(i).AbRSS, 1), size(TrainRSS(i).AbRSS, 2)))];
     %AbRSSNorm = [AbRSSNorm; TrainRSS(i).AbRSSPCA];
     RePoNorm = [RePoNorm; TrainRSS(i).RePo];
     %RePoNorm = [RePoNorm; TrainRSS(i).RePo];
     %AbPoNorm = [AbPoNorm; TrainRSS(i).AbPo ./ repmat(StepSize, size(TrainRSS(i).AbPo, 1), size(TrainRSS(i).AbPo, 2))];
end
%AbRSSNorm = AbRSSNorm - repmat(AllRSSMean(1:Nap, :)', size(AbRSSNorm, 1), 2);
RePoMax = max(max(RePoNorm));
RePoMin = min(min(RePoNorm));
%RePoNorm = (RePoNorm- RePoMin)/(RePoMax- RePoMin);

%AbRSSNorm = (AbRSSNorm - repmat(RSSMean, size(AbRSSNorm, 1), 1)) ./ repmat(RSSVar, size(AbRSSNorm, 1), 1);

%RePoNorm = (RePoNorm - repmat(NagMean, size(RePoNorm, 1), 1)) ./ repmat(NagVar, size(RePoNorm, 1), 1);

AbRSSNormTran(:, 1:NapPca) = AbRSSNorm(:, NapPca+1:2*NapPca);
AbRSSNormTran(:, NapPca+1:2*NapPca) = AbRSSNorm(:, 1:NapPca);
% AbRSSFitNormTran(:, 1:NapPca) = AbRSSFitNorm(:, NapPca+1:2*NapPca);
% AbRSSFitNormTran(:, NapPca+1:2*NapPca) = AbRSSFitNorm(:, 1:NapPca);
% AbPoNormTran(:, 1:2) = AbPoNorm(:, 3:4);
% AbPoNormTran(:, 3:4) = AbPoNorm(:, 1:2);
RePoNormTran = - RePoNorm;
% RePoNorm = (RePoNorm- RePoMin)/(RePoMax- RePoMin);
% RePoNormTran = (RePoNormTran- RePoMin)/(RePoMax- RePoMin);

AbRSSDB = [AbRSSNorm; AbRSSNormTran];
AbRSSMinusDB = AbRSSDB(:, NapPca+1:2*NapPca) - AbRSSDB(:, 1:NapPca);
RePoDB = [RePoNorm; RePoNormTran];

% if Restart == 1
%     RSSMean = mean(AbRSSDB(:, 1:NapPca));
%     RSSVar = sqrt(var(AbRSSDB(:, 1:NapPca)));
%     NagMean = mean(RePoDB);
%     NagVar = sqrt(var(RePoDB));
% end

ReScale = 1;
%RePoDBE = repmat((RePoDB - repmat(NagMean,size(RePoDB,1), 1))./repmat(NagVar,size(RePoDB,1), 1), 1, size(AbRSSDB, 2) / 4);
%RePoDBE = repmat(RePoDB, 1, size(AbRSSDB, 2) / 4);
% 
% 
% Abtm1ReDB = [(AbRSSDB(:, 1:NapPca) - repmat(RSSMean, size(AbRSSDB, 1), 1))./repmat(RSSVar, size(AbRSSDB, 1), 1), RePoDBE * ReScale];
% AbtDB = (AbRSSDB(:, NapPca+1:2*NapPca) - repmat(RSSMean, size(AbRSSDB, 1), 1))./repmat(RSSVar, size(AbRSSDB, 1), 1);


Abtm1ReDB = [AbRSSDB(:, 1:NapPca), RePoDB * ReScale];
AbtDB = AbRSSDB(:, NapPca+1:2*NapPca);

% InputData = [];
% InputDataFit = [];
% InputTarget = [];
% InputTargetFit = [];
% InputTargetPos = [];
% InputData = [InputData; AbRSSNorm(:, 1:NapPca), RePoNorm];
% InputData = [InputData; AbRSSNormTran(:, 1:NapPca), RePoNormTran];
% % InputDataFit = [InputDataFit; AbRSSFitNorm(:, 1:NapPca), RePoNorm];
% % InputDataFit = [InputDataFit; AbRSSFitNormTran(:, 1:NapPca), RePoNormTran];
% InputTarget = [InputTarget; AbRSSNorm(:, NapPca+1:2*NapPca)];
% InputTarget = [InputTarget; AbRSSNormTran(:, NapPca+1:2*NapPca)];
% % InputTargetFit = [InputTargetFit; AbRSSFitNorm(:, NapPca+1:2*NapPca)];
% % InputTargetFit = [InputTargetFit; AbRSSFitNormTran(:, NapPca+1:2*NapPca)];
% InputTargetPos = [InputTargetPos; AbPoNorm(:, 3:4)];
% InputTargetPos = [InputTargetPos; AbPoNormTran(:, 3:4)];
% 
% % BatchSize = 10;
% % DataSize = size(InputData, 1);
% % BatchNum = fix(DataSize/BatchSize);
% 
% net = [];
% net = newff(InputData', InputTarget',[10],{'logsig','logsig'},'trainlm');%'purelin'
% net.trainParam.max_fail = 6;  
% net.trainParam.min_grad = 10^(-15);
% net.trainparam.epochs = 100;
% 
% % netFit = [];
% % netFit = newff(InputDataFit', InputTargetFit',[20],{'logsig','purelin'},'trainlm');%'purelin'
% % netFit.trainParam.max_fail = 6;  
%     
% [net,~] = train(net, InputData', InputTarget');
% OutputTarget = net(InputData')';
% ErrorNN = norm(OutputTarget - InputTarget)
% 
% % [netFit,~] = train(netFit, InputDataFit', InputTargetFit');
% % OutputTargetFit = netFit(InputDataFit')';
% % ErrorNNFit = norm(OutputTargetFit - InputTargetFit)
% 
% % [zIT, mapping] = compute_mapping(InputTarget, 'GPLVM', 1);
% % [zOT, mapping] = compute_mapping(OutputTarget, 'GPLVM', 1);
% % PlotGriddata(InputTargetPos(:,1),InputTargetPos(:,2),zIT);
% % PlotGriddata(InputTargetPos(:,1),InputTargetPos(:,2),zOT);
% % 
% % for i = 1: size(InputTarget, 2)
% %         figure
% %         PlotGriddata(InputTargetPos(:,1),InputTargetPos(:,2),InputTarget(:,i),i);
% %         figure
% %         PlotGriddata(InputTargetPos(:,1),InputTargetPos(:,2),OutputTarget(:,i),i);
% % %         figure
% % %         PlotGriddata(InputTargetPos(:,1),InputTargetPos(:,2),OutputTargetFit(:,i),i);
% % %     axis([0 5 0 5 0 1])
% % end
% 
% % %for i = 1: size(InputTarget, 1)
% %         figure
% %         PlotGriddata(InputData(:,1),InputData(:,2),InputTarget(:,1),i);
% %         figure
% %         PlotGriddata(InputData(:,1),InputData(:,2),OutputTarget(:,1),i);
% % %end

save OfflineNNMove.mat