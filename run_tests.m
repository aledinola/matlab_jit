%% Run utility tests followed by the full-model validation in main.m.
clear;
clc;

% Scalar interpolation is exact for affine data, including endpoints.
x_test = (0:4).';
y_test = 2*x_test+3;
query_test = [0,0.25,2.75,4];
for iq = 1:length(query_test)
    expected = 2*query_test(iq)+3;
    actual = interp1_scal(x_test,y_test,query_test(iq));
    assert(abs(actual-expected) < 1.0e-13,'interp1_scal test failed.');
end

% Golden search recovers the unique maximum of a concave quadratic.
[x_star,f_star] = golden(@(x) -(x-1.25)^2+4,-2,5,1.0e-10);
assert(abs(x_star-1.25) < 5.0e-8,'golden maximizer test failed.');
assert(abs(f_star-4) < 1.0e-12,'golden objective-value test failed.');

fprintf('Utility tests passed. Running the full-model tests in main.m.\n');
main;
