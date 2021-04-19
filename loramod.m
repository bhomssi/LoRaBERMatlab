function [y] = loramod(x,SF,BW,fs,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if (nargin < 4)
    error(message('comm:pskmod:numarg1'));
end

if (nargin > 5)
    error(message('comm:pskmod:numarg2'));
end

% Check that x is a positive integer
if (~isreal(x) || any(any(ceil(x) ~= x)) || ~isnumeric(x))
    error(message('comm:pskmod:xreal1'));
end

M = 2^SF ;

% Check that M is a positive integer
if (~isreal(M) || ~isscalar(M) || M<=0 || (ceil(M)~=M) || ~isnumeric(M))
    error(message('comm:pskmod:Mreal'));
end

% Check that x is within range
if ((min(min(x)) < 0) || (max(max(x)) > (M-1)))
    error(message('comm:pskmod:xreal2'));
end

if nargin == 4
    Inv = 1 ;
elseif nargin == 5
    Inv = varargin{1} ;
end

Ts = 2^SF/BW ;
beta = BW/(2*Ts) ;
n_symbol = fs.*M/BW ;
t_symbol = (0:n_symbol-1).*1/fs ;

y = [] ;
for ctr = 1 : length(x)
    gamma = (x(ctr) - M/2)*BW/M ;
    lambda = 1 - x(ctr)/M ;
    t1 = t_symbol(1:end*lambda) ;
    t2 = t_symbol(end*lambda+1:end) ;
    y = [y; exp(-j.*2.*pi.*(t1'*gamma + beta*t1'.^2)*Inv); exp(-j.*2.*pi.*(t2'*(-BW + gamma) + beta*t2'.^2)*Inv)] ;
end
y = reshape(y,1,numel(y))' ;
end


