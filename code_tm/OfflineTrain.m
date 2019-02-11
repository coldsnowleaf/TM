clear
%load RSSMap.mat
global Snag;
global Svar;
global ErrorMeanArrayNg;
global ErrorMeanArrayEKF;
global PathNum;
global Nap;
global Restart;
% global TrainRSS;
global SnagTrain;
global PathTrain;
global ApLoc;
global NapPca;
global tt;
global ErrorNgTt;
global ErrorEKFTt;
global ErrorFPTt;
global mapping;
global emp;
global trainrandp;

Width = 200;
Height = 200;
St = 150;
S0 = 1;
Svar = 5; 
StepSize = 20;
Nap = 30; 
NapPca = Nap; %fix(Nap/10)+1;
PathNum = 50;
Snag = 1;
Restart = 1;
% 
% 
emp = 0;

% ApLoc = [];
% ApLoc(:,1) = rand(Nap,1) * Width;
% ApLoc(:,2) = rand(Nap,1) * Height;

    
if Restart == 1
    
    ApLoc = [];
    ApLoc(:,1) = rand(Nap,1) * Width;
    ApLoc(:,2) = rand(Nap,1) * Height;
    
    PathTrain = [];
%     ApLoc = [];
%     ApLoc(:,1) = rand(Nap,1) * Width;
%     ApLoc(:,2) = rand(Nap,1) * Height;
    for i = 1 : PathNum
%         TrainRSS(i).Path(1,:) = rand(1,2) * Width;
%         TrainRSS(i).Path(2,:) = rand(1,2) * Height;
        PathTrain(i,:) = [rand(1,2) * Width, rand(1,2) * Height];
        %PathTrain(i,:) = [[0,0] * Width, rand(1,2) * Height];
    end
end
%PathTrain(PathNum,:) = [rand(1,2) * Width, rand(1,2) * Height];

for i = 1 : PathNum
    %PathTrain(i,:) = [rand(1,2) * Width, rand(1,2) * Height];
    TrainRSS(i).Path(1,:) = PathTrain(i,1:2);
    TrainRSS(i).Path(2,:) = PathTrain(i,3:4);       
end

% for i = 1 : size(SegRSS, 1)
%     TrainRSS(i).Path = SegRSS(i).Segs;
% end

% PathNum = size(TrainRSS, 2);
% 
% inc = 1;
% Snag = 0.75;
% for k = 2 : inc
%     for i = 1 : PathNum
%         for j = 1 : size(TrainRSS(i).Path, 2)
%             TrainRSS((k-1)*PathNum+i).Path(1,j) = TrainRSS(i).Path(1,j) + normrnd(0,Snag*StepSize);
%             TrainRSS((k-1)*PathNum+i).Path(2,j) = TrainRSS(i).Path(2,j) + normrnd(0,Snag*StepSize);
%         end
%     end
% end
% 
% PathNumInc = PathNum * inc;


for i = 1 : PathNum
    temps = 0;
    for j = 2 : size(TrainRSS(i).Path, 2)
        dis = norm(TrainRSS(i).Path(:,j) - TrainRSS(i).Path(:,j-1));
        samplePoints = round(dis/StepSize) + 1;
        for s = 1 : samplePoints
            temps = temps +1;
            offsetX = (TrainRSS(i).Path(1, j) - TrainRSS(i).Path(1, j-1)) / samplePoints * (s-1);
            offsetY = (TrainRSS(i).Path(2, j) - TrainRSS(i).Path(2, j-1)) / samplePoints * (s-1);
            
            for k = 1 : Nap
                d = pdist2(ApLoc(k,:), [TrainRSS(i).Path(1, j-1) + offsetX, TrainRSS(i).Path(2, j-1) + offsetY]);
                TrainRSS(i).RSSs(k,temps) = St-S0-10*2*log(d/1)+normrnd(0,Svar);
                %TrainRSS(i).RSSs(k,temps) = normrnd(0,Svar);
                TrainRSS(i).RSSsExp(k,temps) = St-S0-10*2*log(d/1);
            end
            TrainRSS(i).RSSs(k+1, temps) = TrainRSS(i).Path(1, j-1) + offsetX;
            TrainRSS(i).RSSs(k+2, temps) = TrainRSS(i).Path(2, j-1) + offsetY;
        end
    end
end
if Restart == 1
    trainrandp = [];
    for i = 1 : PathNum
        trainrandp(i).p = randperm(size(TrainRSS(i).RSSs, 2)*Nap);
    end
end

for i = 1 : PathNum
    for j = 1 : fix(Nap*size(TrainRSS(i).RSSs, 2)*emp)
%         w = fix(rand*size(TrainRSS(i).RSSs, 2)) +1;
%         h = fix(rand*Nap)+1;
        h = fix((trainrandp(i).p(j)-1)/size(TrainRSS(i).RSSs, 2))+1;
        w = trainrandp(i).p(j)-(h-1)*size(TrainRSS(i).RSSs, 2);
        if w > size(TrainRSS(i).RSSs, 2)
            w = size(TrainRSS(i).RSSs, 2);
        end
        if h > Nap
            h = Nap;
        end
        TrainRSS(i).RSSs(h,w) = -20;
    end
%     TrainRSS(i).RSSs(1:Nap, 1:size(TrainRSS(i).RSSs, 2)) = 0;
%     TrainRSS(i).RSSs(1,fix(size(TrainRSS(i).RSSs, 2)/2)+1) = 1000;
end

AllRSS = [];
AllRSSExp = [];
for i = 1 : PathNum
    AllRSS = [AllRSS, TrainRSS(i).RSSs];
    AllRSSExp = [AllRSSExp, TrainRSS(i).RSSsExp];
end
%AllRSSMean = mean(AllRSS(:, 1:Nap), 2);

if Restart == 1
%     coeff = pca(AllRSS(1:Nap,:)');
%     mapping.M = coeff(:,1:NapPca);%eye(Nap);%coeff(:,1:NapPca);   
     mapping.M = eye(Nap);
     NapPca = Nap;
end

AllRSSPCA = [(AllRSS(1:Nap,:)'*mapping.M)'; AllRSS(Nap+1:Nap+2,:)];
% for i = 1 : PathNum
%     csize = size(TrainRSS(i).RSSs, 2);
%     for j = 1 : Nap
%         A = polyfit([1:csize],TrainRSS(i).RSSs(j,:),1);
%         TrainRSS(i).RSSsFit(j,:) = polyval(A,[1:csize]);
%     end
% end

 if Restart == 1 %IMU
    for i = 1 : PathNum
        for j = 1 : 100%size(TrainRSS(i).RSSs, 2) - 1 
            SnagTrain(i,j).Norm = [normrnd(0, Snag * StepSize), normrnd(0, Snag * StepSize)];
        end
    end
 end


RSSScaleMax = - 10000000;
RSSScaleMin = 10000000;

for i = 1 : PathNum
    for j = 1 : size(TrainRSS(i).RSSs, 2) - 1 
        TrainRSS(i).AbRSS(j,:) = [TrainRSS(i).RSSs(1:Nap, j)', TrainRSS(i).RSSs(1:Nap, j+1)'];
        TrainRSS(i).AbRSSExp(j,:) = [TrainRSS(i).RSSsExp(1:Nap, j)', TrainRSS(i).RSSsExp(1:Nap, j+1)'];
        %TrainRSS(i).AbRSSFit(j,:) = [TrainRSS(i).RSSsFit(1:Nap, j)', TrainRSS(i).RSSsFit(1:Nap, j+1)'];
        TrainRSS(i).AbRSSPCA(j,:) = [TrainRSS(i).RSSs(1:Nap, j)'*mapping.M, TrainRSS(i).RSSs(1:Nap, j+1)'*mapping.M];
        RSSScaleMax = max(RSSScaleMax, max(max(TrainRSS(i).AbRSSPCA(j, 1:NapPca)), max(TrainRSS(i).AbRSSPCA(j, NapPca+1:2*NapPca))));
        RSSScaleMin = min(RSSScaleMin, min(min(TrainRSS(i).AbRSSPCA(j, 1:NapPca)), min(TrainRSS(i).AbRSSPCA(j, NapPca+1:2*NapPca))));
        %RSSScaleMin = min(RSSScaleMin, min(min(TrainRSS(i).RSSs(1:Nap, j)), min(TrainRSS(i).RSSs(1:Nap, j+1))));
        %TrainRSS(i).RePo(j,:) = [TrainRSS(i).RSSs(Nap+1:Nap+2, j+1)' - TrainRSS(i).RSSs(Nap+1:Nap+2, j)'] + [normrnd(0, Snag * StepSize), normrnd(0, Snag * StepSize)];
        TrainRSS(i).RePo(j,:) = [TrainRSS(i).RSSs(Nap+1:Nap+2, j+1)' - TrainRSS(i).RSSs(Nap+1:Nap+2, j)'] + SnagTrain(i,j).Norm;
        TrainRSS(i).RePoTrue(j,:) = [TrainRSS(i).RSSs(Nap+1:Nap+2, j+1)' - TrainRSS(i).RSSs(Nap+1:Nap+2, j)'];
        TrainRSS(i).AbPo(j,:) = [TrainRSS(i).RSSs(Nap+1:Nap+2, j)', TrainRSS(i).RSSs(Nap+1:Nap+2, j+1)'];    
    end
%     for j = 1 : size(TrainRSS(i).RSSs, 2)  
%         TrainRSS(i).AbRSS(size(TrainRSS(i).RSSs, 2)+j,:) = [TrainRSS(i).RSSs(1:Nap, j)', TrainRSS(i).RSSs(1:Nap, j)'];
%         TrainRSS(i).RePo(size(TrainRSS(i).RSSs, 2)+j,:) = [TrainRSS(i).RSSs(Nap+1:Nap+2, j)' - TrainRSS(i).RSSs(Nap+1:Nap+2, j)'];
%         TrainRSS(i).AbPo(size(TrainRSS(i).RSSs, 2)+j,:) = [TrainRSS(i).RSSs(Nap+1:Nap+2, j)', TrainRSS(i).RSSs(Nap+1:Nap+2, j)'];    
%     end
end

save OfflineTrain.mat
