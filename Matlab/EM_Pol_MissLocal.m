function [Priors, Mu, Sigma, Pix,Miss_Data] = EM_Pol_MissLocal(Data,Miss_Data,P_Dim, Priors0, Mu0, Sigma0, Miss_Dim, Obs_Dim)
%
% Expectation-Maximization estimation of GMM parameters.
% This source code is the implementation of the algorithms described in 
% Section 2.6.1, p.47 of the book "Robot Programming by Demonstration: A 
% Probabilistic Approach".
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
%
% This function learns the parameters of a Gaussian Mixture Model 
% (GMM) using a recursive Expectation-Maximization (EM) algorithm, starting 
% from an initial estimation of the parameters.
%
%
% Inputs -----------------------------------------------------------------
%   o Data:    D x N array representing N datapoints of D dimensions.
%   o Priors0: 1 x K array representing the initial prior probabilities 
%              of the K GMM components.
%   o Mu0:     D x K array representing the initial centers of the K GMM 
%              components.
%   o Sigma0:  D x D x K array representing the initial covariance matrices 
%              of the K GMM components.
% Outputs ----------------------------------------------------------------
%   o Priors:  1 x K array representing the prior probabilities of the K GMM 
%              components.
%   o Mu:      D x K array representing the centers of the K GMM components.
%   o Sigma:   D x D x K array representing the covariance matrices of the 
%              K GMM components.
%
% This source code is given for free! However, I would be grateful if you refer 
% to the book (or corresponding article) in any academic publication that uses 
% this code or part of it. Here are the corresponding BibTex references: 
%
% @book{Calinon09book,
%   author="S. Calinon",
%   title="Robot Programming by Demonstration: A Probabilistic Approach",
%   publisher="EPFL/CRC Press",
%   year="2009",
%   note="EPFL Press ISBN 978-2-940222-31-5, CRC Press ISBN 978-1-4398-0867-2"
% }
%
% @article{Calinon07,
%   title="On Learning, Representing and Generalizing a Task in a Humanoid Robot",
%   author="S. Calinon and F. Guenter and A. Billard",
%   journal="IEEE Transactions on Systems, Man and Cybernetics, Part B",
%   year="2007",
%   volume="37",
%   number="2",
%   pages="286--298",
% }

%% Criterion to stop the EM iterative update
loglik_threshold = 1e-10; %likelihood threshold for termination of further computation

%% Initialization of the parameters
[nbVar, nbData] = size(Data);
[~, nbData2] = size(Miss_Data);
nbStates = size(Sigma0,3);
loglik_old = -realmax;
nbStep = 0;

Mu = Mu0;
Sigma = Sigma0;
Priors = Priors0;
PredictMiss = zeros(nbStates,nbData2); %initializiation of matrices, for code optimization
PredictMissh = zeros(nbStates,nbData2);
PredictMissCov = zeros(length(Miss_Dim),length(Miss_Dim),nbData2);
TotalMat=zeros(nbVar,nbVar);

%% EM fast matrix computation (see the commented code for a version 
%% involving one-by-one computation, which is easier to understand)
for Beta=0.5:0.1:1 % not used anywhere 
    Beta
if(Beta==1)
    Iterations = 1000;
else
    Iterations = 1000;
end
    
for MaxSteps=1:Iterations
  %% E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i=1:nbStates
    %Compute probability p(x|i)
    Data_Tmp = Data(P_Dim,:); % for reference data
    Data_Tmp = FixPolarRange(Data_Tmp,Mu(P_Dim,i)-pi,Mu(P_Dim,i)+pi);
    Data2(:,:,i)=Data;
    Data2(P_Dim,:,i)=Data_Tmp;
    
    Data_Tmp2 = Miss_Data(P_Dim,:); %for missing data
    Data_Tmp2 = FixPolarRange(Data_Tmp2,Mu(P_Dim,i)-pi,Mu(P_Dim,i)+pi);
    Data3(:,:,i)=Miss_Data;
    Data3(P_Dim,:,i)=Data_Tmp2;
    
%     figure;
%     plot(Data2(1,:),Data2(2,:),'*','color',[0.5 0 0.5]);
    Pxi(:,i) = gaussPDF(Data2(:,:,i), Mu(:,i), Sigma(:,:,i)); %calculating gaussian value for this data point, for missing and complete data
    Pxi2(:,i)= gaussPDF(Data3(:,:,i), Mu(:,i), Sigma(:,:,i));
  end 
  %Compute posterior probability p(i|x)
  Pix_tmp = (repmat(Priors,[nbData 1]).*Pxi).^Beta; 
  Pix = Pix_tmp ./ repmat(sum(Pix_tmp,2),[1 nbStates]); %pim
  
  Pix_tmp2 = (repmat(Priors,[nbData2 1]).*Pxi2).^Beta; 
  Pix2 = Pix_tmp2 ./ repmat(sum(Pix_tmp2,2),[1 nbStates]); %qjm
  
  %Compute cumulated posterior probability
  E = sum(Pix); %hm in paper, for complete and missing data
  E2 = sum(Pix2); %rm in paper
  Priors_Old = Priors;
  Mu_Old = Mu;
  Sigma_Old = Sigma;
  %% M-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updating values
  for i=1:nbStates
    %Update the priors
    Priors(i) = (E(i)+E2(i)) / (nbData+nbData2); %contributuion of this gaussian on the whole data, including missing anf complete data
%     tic
%     for j=1:nbData2
%         PredictMiss(i,j) = [Mu(Miss_Dim,i) + Sigma(Miss_Dim,Obs_Dim,i)*inv(Sigma(Obs_Dim,Obs_Dim,i))*(Miss_Data(:,j)-Mu(Obs_Dim,i))];
%     end
%     toc
    PredictMiss(i,:) = [repmat(Mu(Miss_Dim,i),1,nbData2) + Sigma(Miss_Dim,Obs_Dim,i)/Sigma(Obs_Dim,Obs_Dim,i)*(Miss_Data(Obs_Dim,:)-repmat(Mu(Obs_Dim,i),1,nbData2))]; % GMR, :) in paper, clock signal predicted using current state of Gaussian 
%     [val ind]=max(Pix2(:,i));
%     PredictMiss(i,:) = PredictMiss(i,ind)*ones(size(PredictMiss(i,:)));
    
%     c_mean = sum(cos(PredictMiss(i,:)).*Pix2(:,i)')/sum(Pix2(:,i));
%     s_mean = sum(sin(PredictMiss(i,:)).*Pix2(:,i)')/sum(Pix2(:,i));
    
    c_mean = sum(cos([Data(1,:) PredictMiss(i,:)]).*[Pix(:,i)' Pix2(:,i)'])/sum([Pix(:,i)' Pix2(:,i)']); % cosine terms
    s_mean = sum(sin([Data(1,:) PredictMiss(i,:)]).*[Pix(:,i)' Pix2(:,i)'])/sum([Pix(:,i)' Pix2(:,i)']); % sine terms
    
       for i2=1:length(c_mean) 
            if(c_mean(i2)<0)
                mean_angle(i2) = atan(s_mean(i2)/c_mean(i2))+pi;
            elseif(s_mean(i2)>0)
                mean_angle(i2) = atan(s_mean(i2)/c_mean(i2));
            else
                mean_angle(i2) = atan(s_mean(i2)/c_mean(i2))+2*pi;
            end
       end
     %Alternate option
%      mean_angle = sum(PredictMiss(i,:).*Pix2(:,i)')/sum(Pix2(:,i));
%      mean_angle = FixPolarRange(mean_angle,Mu(P_Dim,i)-pi,Mu(P_Dim,i)+pi);
%         if(mean_angle<0)
%            abcd = 120 
%         end
%     end
    
    %Update the centers
%     Mu(:,i) = ( (Data2(:,:,i)*Pix(:,i)) + [repmat(mean_angle,1,nbData2);Miss_Data]*Pix2(:,i) ) / (E(i)+E2(i));
    Mu(2:end,i) = ( (Data2(2:end,:,i)*Pix(:,i)) + [Miss_Data(Obs_Dim,:)]*Pix2(:,i) ) / (E(i)+E2(i));
    Mu(P_Dim,i) = mean_angle;
    Mu(P_Dim,i) = FixPolarRange(Mu(P_Dim,i),0,2*pi);
    
    
    Data_Tmp = Data(P_Dim,:); % now that mean is updated and we transform the data again to place it in the -pi and +pi range
    Data_Tmp = FixPolarRange(Data_Tmp,Mu(P_Dim,i)-pi,Mu(P_Dim,i)+pi);
    PredictMiss(i,:) = FixPolarRange(PredictMiss(i,:),Mu(P_Dim,i)-pi,Mu(P_Dim,i)+pi);
    
%     tic
%     for j=1:nbData2
%         PredictMissh(i,j) = Pix2(j,i)*PredictMiss(i,j); 
%         PredictMissCov(:,:,j) = Pix2(j,i)*(Sigma(Miss_Dim,Miss_Dim,i)-Sigma(Miss_Dim,Obs_Dim,i)*inv(Sigma(Obs_Dim,Obs_Dim,i))*Sigma(Miss_Dim,Obs_Dim,i)'+PredictMiss(i,j)*PredictMiss(i,j)');
%         TotalMat(Miss_Dim,Miss_Dim,j) = PredictMissCov(:,:,j) + Pix2(j,i)*Mu(Miss_Dim,i)*Mu(Miss_Dim,i)' - 2*PredictMissh(i,j)*Mu(Miss_Dim,i)';
%         TotalMat(Obs_Dim,Obs_Dim,j) = Pix2(j,i)*(Miss_Data(:,j)-Mu(Obs_Dim,i))*(Miss_Data(:,j)-Mu(Obs_Dim,i))';
%         TotalMat(Obs_Dim,Miss_Dim,j) = (PredictMissh(i,j) + Pix2(j,i)*Mu(Miss_Dim,i))*(Miss_Data(:,j)-Mu(Obs_Dim,i))';
%         TotalMat(Miss_Dim,Obs_Dim,j) = TotalMat(Obs_Dim,Miss_Dim,j)';
%     end
%     toc
%     TotalMat2=zeros(4,4);

    PredictMissh(i,:) = Pix2(:,i)'.*PredictMiss(i,:);
    PredictMissCov = reshape(Pix2(:,i)'.*(Sigma(Miss_Dim,Miss_Dim,i)-Sigma(Miss_Dim,Obs_Dim,i)/Sigma(Obs_Dim,Obs_Dim,i)*Sigma(Miss_Dim,Obs_Dim,i)'+PredictMiss(i,:).*PredictMiss(i,:)),1,1,[]);
    TotalMat(Miss_Dim,Miss_Dim) = sum(PredictMissCov + reshape(Pix2(:,i)'*Mu(Miss_Dim,i)*Mu(Miss_Dim,i)' - 2*PredictMissh(i,:)*Mu(Miss_Dim,i)',1,1,nbData2),3);% A22 calculated
    Data_tmp2 = Miss_Data(Obs_Dim,:) - repmat(Mu(Obs_Dim,i),1,nbData2);
    TotalMat(Obs_Dim,Obs_Dim) = (repmat(Pix2(:,i)',length(Obs_Dim), 1) .* Data_tmp2*Data_tmp2'); % A11 calculated
    TotalMat(Obs_Dim,Miss_Dim) = (PredictMissh(i,:) - Pix2(:,i)'*Mu(Miss_Dim,i))*(Miss_Data(Obs_Dim,:)-repmat(Mu(Obs_Dim,i),1,nbData2))';% A21 calculated
    TotalMat(Miss_Dim,Obs_Dim) = TotalMat(Obs_Dim,Miss_Dim)'; % A12 calculated
    
    Data2(:,:,i)=Data; %not transformed
    Data2(P_Dim,:,i)=Data_Tmp; %copying transformed polar region into Data2
    
    %Update the covariance matrices
    Data_tmp1 = Data2(:,:,i) - repmat(Mu(:,i),1,nbData); %intermediate step
%     Sigma(:,:,i) = ((repmat(Pix(:,i)',nbVar, 1) .* Data_tmp1*Data_tmp1') + sum(TotalMat,3)) / (E(i)+E2(i));
    Sigma(:,:,i) = ((repmat(Pix(:,i)',nbVar, 1) .* Data_tmp1*Data_tmp1') + TotalMat) / (E(i)+E2(i));
    %% Add a tiny variance to avoid numerical instability % to avoid singularity of Gaussian, to avoid variance becoming zero, where likelihood becomes finite (not good) 
    Sigma(:,:,i) = Sigma(:,:,i) + 1E-3.*diag(ones(nbVar,1));
    
%     [V D] = eig(Sigma(:,:,i));
%     D=diag(D);
%     D(D<1E-3) = 1E-3;
%     D=diag(D);
%     Sigma(:,:,i) = V*D*V';
%     
  end
  Mu(P_Dim,:) = FixPolarRange(Mu(P_Dim,:),0,2*pi); %Mu transformed to be between 0 and 2pi
%   [Miss_Data(Miss_Dim,:), Sigma_y] = GMR_Polar_max2(Priors, Mu, Sigma, Miss_Data, [2 3 4], 1, 1);
%   [Miss_Data(Miss_Dim,:), Sigma_y] = GMR_Polar_max2(Priors, Mu, Sigma, Miss_Data, Obs_Dim, 1, 1);

% Till now gaussian updated, clock signal to be updated now, missing data
% updated, first updation after one iteration
  [Temp_Miss_Data, Sigma_y] = GMR_Polar_max2(Priors, Mu, Sigma, Miss_Data, Obs_Dim, 1, 1, Beta);
   Miss_Data(Miss_Dim,:) = Temp_Miss_Data;
  Miss_Data(P_Dim,:) = FixPolarRange(Miss_Data(P_Dim,:),0,2*pi);
  
  %% Stopping criterion %%%%%%%%%%%%%%%%%%%%
%   for i=1:nbStates
%     %Compute the new probability p(x|i)
%     Pxi(:,i) = gaussPDF(Data, Mu(:,i), Sigma(:,:,i));
%   end
%   %Compute the log likelihood
%   F = Pxi*Priors';
%   F(find(F<realmin)) = realmin;
%   loglik = mean(log(F));
loglik = sum((Priors_Old-Priors).^2) + sum(sum((Mu_Old-Mu).^2)) + sum(sum(sum((Sigma_Old-Sigma).^2))); %for monitoring the parameters to see if theare close to the threshold value for termination for iteration

  %Stop the process depending on the increase of the log likelihood 
  if abs((loglik/loglik_old)-1) < loglik_threshold
    MaxSteps  
    break;
  end
  loglik_old = loglik;
%   nbStep = nbStep+1;
end
end
% %% EM slow one-by-one computation (better suited to understand the
% %% algorithm) 
% while 1
%   %% E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   for i=1:nbStates
%     %Compute probability p(x|i)
%     Pxi(:,i) = gaussPDF(Data, Mu(:,i), Sigma(:,:,i));
%   end
%   %Compute posterior probability p(i|x)
%   for j=1:nbData
%     Pix(j,:) = (Priors.*Pxi(j,:))./(sum(Priors.*Pxi(j,:))+realmin);
%   end
%   %Compute cumulated posterior probability
%   E = sum(Pix);
%   %% M-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   for i=1:nbStates
%     %Update the priors
%     Priors(i) = E(i) / nbData;
%     %Update the centers
%     Mu(:,i) = Data*Pix(:,i) / E(i);
%     %Update the covariance matrices 
%     covtmp = zeros(nbVar,nbVar);
%     for j=1:nbData
%       covtmp = covtmp + (Data(:,j)-Mu(:,i))*(Data(:,j)-Mu(:,i))'.*Pix(j,i);
%     end
%     Sigma(:,:,i) = covtmp / E(i);
%   end
%   %% Stopping criterion %%%%%%%%%%%%%%%%%%%%
%   for i=1:nbStates
%     %Compute the new probability p(x|i)
%     Pxi(:,i) = gaussPDF(Data, Mu(:,i), Sigma(:,:,i));
%   end
%   %Compute the log likelihood
%   F = Pxi*Priors';
%   F(find(F<realmin)) = realmin;
%   loglik = mean(log(F));
%   %Stop the process depending on the increase of the log likelihood 
%   if abs((loglik/loglik_old)-1) < loglik_threshold
%     break;
%   end
%   loglik_old = loglik;
%   nbStep = nbStep+1;
% end

%% Add a tiny variance to avoid numerical instability
% for i=1:nbStates
%   Sigma(:,:,i) = Sigma(:,:,i) + 1E-4.*diag(ones(nbVar,1));
% end
end

