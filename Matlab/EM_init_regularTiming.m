function [Priors, Mu, Sigma] = EM_init_regularTiming(Data, nbStates)
%
% This function initializes the parameters of a Gaussian Mixture Model 
% (GMM) by using k-means clustering algorithm.
%
% Inputs -----------------------------------------------------------------
%   o Data:     D x N array representing N datapoints of D dimensions.
%   o nbStates: Number K of GMM components.
% Outputs ----------------------------------------------------------------
%   o Priors:   1 x K array representing the prior probabilities of the
%               K GMM components.
%   o Mu:       D x K array representing the centers of the K GMM components.
%   o Sigma:    D x D x K array representing the covariance matrices of the 
%               K GMM components.
% Comments ---------------------------------------------------------------
%   o This function uses the 'kmeans' function from the MATLAB Statistics 
%     toolbox. If you are using a version of the 'netlab' toolbox that also
%     uses a function named 'kmeans', please rename the netlab function to
%     'kmeans_netlab.m' to avoid conflicts. 
%
% Copyright (c) 2006 Sylvain Calinon, LASA Lab, EPFL, CH-1015 Lausanne,
%               Switzerland, http://lasa.epfl.ch

% TimingSep = linspace(min(Data(1,:)), max(Data(1,:)), nbStates+1);
TimingSep = linspace(1, length(Data(1,:)), nbStates+1); %states are the no. gaussians
% TimingSep = linspace(1, length(Data(1,:)), nbStates+2);
% TimingSep = TimingSep(2:end-1)
for i=1:nbStates
  idtmp = find( [1: length(Data(1,:))]>=TimingSep(i) & [1: length(Data(1,:))]<TimingSep(i+1));
  Priors(i) = length(idtmp);
  Mu(:,i) = mean(Data(:,idtmp)');
  Sigma(:,:,i) = cov(Data(:,idtmp)');% + 1E-10.*diag(ones(size(Mu(:,i),1),1));
end
Priors = Priors ./ sum(Priors);
nbData = length(Data);
  for i=1:nbStates
    %Compute probability p(x|i)
    Pxi(:,i) = gaussPDF(Data, Mu(:,i), Sigma(:,:,i));
  end
  Pix_tmp = repmat(Priors,[nbData 1]).*Pxi;
  Pix = Pix_tmp ./ repmat(sum(Pix_tmp,2),[1 nbStates]);
  %Compute cumulated posterior probability
  E = sum(Pix);
  Priors = E./sum(E);
end


