function [outmat] = interpmmlarson(inmat,numinterp);
%this function interpolates in the mm-plane for a single set of data with
%one concentration (each csv file)

outmat = [linspace(inmat(1,2),inmat(end,2),numinterp)' zeros(numinterp,1)];
outmat(:,2) = interp1(inmat(:,2),inmat(:,3),outmat(:,1),'spline');
outmat = [ones(numinterp,1)*inmat(1,1) outmat];

