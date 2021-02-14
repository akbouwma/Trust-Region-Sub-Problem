A1 = load('A1.mat');
A2 = load('A2.mat');


[D,V] = eig(A1.Expression1, A2.Expression1) % should this be -1 * A2.Expression1 ?  I'm not sure lol