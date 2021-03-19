function ARNSGAIII(Global)
% <algorithm> <A>
% 基于参考点选择策略的改进型NSGA-III算法

%------------------------------- Reference --------------------------------
% 耿焕同,戴中斌,王天雷,许可, 等.基于参考点选择策略的改进型NSGA-III算法[J].模式识别与人工智能,2020.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the reference points and random population

    [Z,numOfReferencePoint] = DzbUniformPoint(Global.N * 1.19,Global.M);	%初始化参考点数目>=N*1.20,亦可以写作>N*1.19
    
    Population   = Global.Initialization();     %初始化种群
    
    M = Global.M;                               %目标数
    N = Global.N;                               %种群规模
    D = Global.D;                               %决策变量维度
    ub =  Global.upper;                         %决策变量上界
    lb = Global.lower;                          %决策变量下界
    Gmax = Global.maxgen;                       %种群最大进化代数
    
    
    threshold = D*((0.5+1/N)*log10(0.5+1/N)-0.5*log10(0.5));    %种群进化阶段阈值
	fprintf('阈值：%f\n',threshold);

    
    [dec, ~] = getPopulationDecAndObj(Population);              %提取种群的决策向量
    X = sort(dec);                                              %对决策向量每个维度进行排序
    mid = X(int32(N*2/4),:);                                    %决策向量每个维度的中位数
   
    i=1;
    entropyArr = zeros(1,Gmax);                                 %用来存储每代种群的熵
    entropyArr2 = zeros(1,Gmax);                                %用来存储相邻两代种群的熵差
    entropy = getEntropy(Population, ub, lb, mid);              %计算熵
    entropyArr(i) = entropy;

    Zmin = min(Population(all(Population.cons<=0,2)).objs,[],1);
    

    totalOfRPointAssociation = [0];                             %每个参考点关联的个体数统计值
    Dsum = 0;                                                   %统计熵差小于阈值出现的次数
    flag = 0;                                                   %判定是否已经进行参考点删除操作的标记
    fprintf('Dzb种群数：%d\t原始参考点数：%d\t保留参考点数：%d\n',Global.N,numOfReferencePoint,Global.N);
    
    
    %% Optimization
    while Global.NotTermination(Population)

        MatingPool = TournamentSelection(2,Global.N,sum(max(0,Population.cons),2));
        Offspring  = GA(Population(MatingPool));
        Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
        
        
        [Population,rho] = EnvironmentalSelection([Population,Offspring],Global.N,Z,Zmin);
        totalOfRPointAssociation = totalOfRPointAssociation + rho;  %统计参考点关联的个体数，其中rho为每一代参考点关联的个体数目
        
        
        i=i+1;
        entropy = getEntropy(Population, ub, lb, mid);              %计算该种群决策变量的熵
        entropyArr(i) = entropy;
        entropyArr2(i) = abs(entropyArr(i-1) - entropyArr(i));      %计算相邻两代种群的熵差
        
        %如果熵差小于阈值
        if entropyArr2(i) < threshold
            Dsum = Dsum + 1;
        end

        
        [dec, ~] = getPopulationDecAndObj(Population);              %提取种群的决策向量
        X = sort(dec);                                              %对决策向量每个维度进行排序
        mid = X(int32(N*2/4),:);                                    %决策向量每个维度的中位数
        
             
%         if Global.evaluated >= Global.evaluation / 2 && flag == 0
        %参考点筛选操作
        if Dsum >= Gmax * 0.1 && flag == 0
            fprintf('第%d代进入“探究阶”段\n',i);
            flag = flag + 1;
            numOfDel = numOfReferencePoint - N;                     %计算删除的参考点数目
            k=0;
            while k < numOfDel
                [~,n]=min(totalOfRPointAssociation);                %获取关联个体数目最少的参考点
                totalOfRPointAssociation(n) = [];
                Z(n,:)=[];                                          %删除参考点
                k = k + 1;
            end
%             save('Zn.mat','Z');
%             fprintf('保留参考点数：%d\n',size(Z,1));
        end
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
end