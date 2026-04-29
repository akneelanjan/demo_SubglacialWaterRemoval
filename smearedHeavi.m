function Y = smearedHeavi(X,eps)
%Builds smeared Heaviside function (0-to-1) with smearing band eps

Y = zeros(size(X));
Y(X >= 0) = 0.5 + (X(X >= 0)-(eps/2))/(eps) + ...
               sin(2*pi*(X(X >= 0)-(eps/2))/eps)/(2*pi);
Y(X > eps) = 1;