prob = testProblem(15);
step = 10;
eps1 = 10^(-3);
eps2 = 10^(-12);
eps3 = 10^(-6);
max_calls = 10000;

DQL(prob.Objective, prob.x_0, step, eps1, eps2, eps3, max_calls, 2, 3, 5, 1)


