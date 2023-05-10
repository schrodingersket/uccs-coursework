close all
clear all

load('Data2D_Beethoven0')

x = dsites(:,1);
y = dsites(:,2);

figure

tes = delaunayn(dsites);
triplot(tes,x,y,'g')

for l=1:5
    for j=1:1500
        n = size(dsites,1);
        d = zeros(1,n);
        for i=1:n
            temp = dsites;
            temp(i,:) = [];
            [k,d(i)] = dsearchn(temp,dsites(i,:));
            if (k >= i)
                k=k+1;
            end
        end
        r = min(d);
        idx = find(d==r);
    
        dsites(idx(1),:) = [];
        x(idx(1)) = [];
        y(idx(1)) = [];
        rhs(idx(1)) = [];
    end
    figure
    tes = delaunayn(dsites);
    triplot(tes,x,y,'r')
    name = sprintf('Data2D_Beethoven%d', l);
    save(name, 'dsites', 'rhs')
end
