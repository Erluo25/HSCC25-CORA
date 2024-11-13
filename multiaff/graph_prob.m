clc
clear

load struct_vardubins.mat
B = pz.expMat;

A = B*B';

figure(1)
clf
subplot(2,2,1)
spy(A)
title('no reordering')


subplot(2,2,2)
p = amd(A);
p = p(end:-1:1);
spy(A(p,p))
title('amd')


subplot(2,2,3)
p = dissect(A);
spy(A(p,p))
title('dissect')



subplot(2,2,4)
p = symrcm(A);
spy(A(p,p))
title('symrcm')