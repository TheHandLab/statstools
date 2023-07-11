function [M,SD,N,SE,Z,T,P]=describe(data)
    %% process the inputs
    if nargin==0
        error('Missing required input - 1xn data array');
    end
    if ndims(data)>2
        error('input data must be 1xn array');
    end
    if size(data,2)>size(data,1)
        data=data';
        warning('Rotating data');
    end
    if size(data,2)>1
        warning('Using first column of data only');
        data=data(:,1);
    end
    M=nanmean(data);
    SD=nanstd(data);
    N=size(data,1);
    SE=SD./sqrt(N);
    Z=M./SD;
    T=M./SE;
    P=tcdf(T,N-1,'Upper').*2;
end