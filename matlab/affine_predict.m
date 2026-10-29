function target=affine_predict(measurement,V,testWhite,targetWhite)
% Explicit side information: white responses for the actual normalized operators.
% Rows of V follow FL3 XYZ, FL8 XYZ, FL12 XYZ; columns follow A,D65,FL2,FL7,FL11 XYZ.
% This is a floating evaluation. Exact coefficients and guarantees are in R177-predictor-exact.json.
assert(isequal(size(V),[9 15]));
assert(numel(measurement)==15 && numel(testWhite)==15 && numel(targetWhite)==9);
assert(all(isfinite([measurement(:);V(:);testWhite(:);targetWhite(:)])));
target=V*(measurement(:)-testWhite(:)/2)+targetWhite(:)/2;
end
