function [rrefM,solution] = solveMatrix(A,b)

% 参数A,b构成一个增广矩阵；返回值rrefA求A的最简行阶梯形，rank求A的秩

r = rank(A);

Y = [A,b];
e = 1e-5;
[m,n] = size(Y);

%只适用于有唯一解的线性代数方程
if(r~=n-1)
    disp('方程非唯一解！！！');
else
    %从前往后消
    for k=1:(m-1)
        if(abs(Y(k,k))<e)
            [M,I] = max(abs(Y(k+1:m,k)));
            temp = Y(k,:);
            Y(k,:) = Y(k+I,:);
            Y(k+I) = temp;
        end
        for i = k+1:m
            factor = Y(i,k)/Y(k,k);
            %若ai,k元素已经为0则不用相减
            if(abs(factor)>e)
                Y(i,:) = Y(i,:) - factor* Y(k,:);
            end
        end
    
    end
    
    %从后往前消：
    for k=r:-1:2
        Y(k,:) = Y(k,:)/Y(k,k);
        for j = k-1:-1:1
            factor = Y(j,k)/Y(k,k);
            if(abs(factor)>e)
                Y(j,:) = Y(j,:) - factor* Y(k,:);
            end
        end
    end
    Y(1,:) = Y(1,:)/Y(1,1);
    rrefM = Y;
    solution = Y(:,end);
    disp('输出的第一个参数为最简行列式，第二个参数是线性方程组的解');
end

end
