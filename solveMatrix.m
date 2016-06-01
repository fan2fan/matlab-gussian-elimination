function [rrefM,solution] = solveMatrix(A,b)

% ����A,b����һ��������󣻷���ֵrrefA��A������н����Σ�rank��A����

r = rank(A);

Y = [A,b];
e = 1e-5;
[m,n] = size(Y);

%ֻ��������Ψһ������Դ�������
if(r~=n-1)
    disp('���̷�Ψһ�⣡����');
else
    %��ǰ������
    for k=1:(m-1)
        if(abs(Y(k,k))<e)
            [M,I] = max(abs(Y(k+1:m,k)));
            temp = Y(k,:);
            Y(k,:) = Y(k+I,:);
            Y(k+I) = temp;
        end
        for i = k+1:m
            factor = Y(i,k)/Y(k,k);
            %��ai,kԪ���Ѿ�Ϊ0�������
            if(abs(factor)>e)
                Y(i,:) = Y(i,:) - factor* Y(k,:);
            end
        end
    
    end
    
    %�Ӻ���ǰ����
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
    disp('����ĵ�һ������Ϊ�������ʽ���ڶ������������Է�����Ľ�');
end

end
