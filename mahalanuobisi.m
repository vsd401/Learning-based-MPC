%% pdist2函数的理解
x1 =[0.0873         0
    0.1417   -0.2514
    0.1785   -0.4825
    0.1991   -0.6645
    0.2071   -0.7893
    0.2017   -0.9572]';
x2=[0.0873         0
    0.1417   -0.2514
    0.1785   -0.4825
    0.1991   -0.6645
    0.2071   -0.7893]';
M= [0.0100         0
         0    0.0100];
kkk=pdist2(x1',x2','mahalanobis',M)

K = mahanuobisi(x1',x2',M)



function MM =mahanuobisi(X1,X2, M)

for i = 1:size(X1,1)
    for j=1:size(X2,1)
        
        MM(i,j)=sqrt((X1(i,:)-X2(j,:))*M^(-1)*(X1(i,:)-X2(j,:))');
    end
    det =0.01;
    MM
    kernel=det* exp( -0.5 * MM)
end




end