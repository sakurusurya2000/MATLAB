clc
clear all
%% Principal Components of a Data Set  

%% 
% Load the sample data set. 
load hald 

%%
% The ingredients data has 13 observations for 4 variables.  

%% 
% Find the principal components for the ingredients data. 
 coeff = pca(ingredients) 

%%
% The rows of |coeff| contain the coefficients for the four ingredient variables,
% and its columns correspond to four principal components.   

