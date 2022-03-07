%ICM Refugee Problem Model v4.0
clear
n = 14;

node(1).incoming = [];

node(2).incoming = [];

node(3).incoming = [];

node(4).incoming = [1];

node(5).incoming = [4];

node(6).incoming = [5];

node(7).incoming = [5 6];

node(8).incoming = [7 9 10];

node(9).incoming = [5];

node(10).incoming = [8 9 11];

node(11).incoming = [10];

node(12).incoming = [8];

node(13).incoming = [8];

node(14).incoming = [8];

A = xlsread('MigrationMatrix.xlsx');
A(isnan(A)) = 0;
Output = zeros(n,12);

capviolation = zeros(11, 10);
burdenvar = zeros(1, 10);
burden = zeros(11, 10);
reptot = 1000;
Dist = zeros(n, 12, reptot);
MeanDist = zeros(n, 12);

for rep=1:reptot
    Q = A(2:15, 2:15);
    R = .95 * Q;                                                            %Lower bound of flow
    R(:, :, 2) = Q;                                                         %Upper bound
    
    D = diag(rand(n,1));                                                    %Simuluates uncertainty at each node
    Q = (R(:, :, 1) + R(:, :, 2))/2 + (R(:, :, 2) - R(:, :, 1))/2 * D ;     %Flow matrix
    
    for i = 1:n
        if sum(Q(i, :)) >= 1
            rowsum = sum(Q(i, :));
            Q(i, :) = Q(i, :)/(rowsum + .01);
        end
    end
    N = (eye(n) - Q)^(-1);
    
refmax = [5000000, 300000, 120000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
mm = sum(refmax);

m = [1000000, 300000, 120000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
Flow = Q;
sitrep = m;
mcount = 1;
flowcount = 0;
flowaltered = zeros(n, 2);
for time = 1:12
    cap = [0, 0, 0, mm*0.041453131, mm*0.009543048, mm*0.007404276, mm*0.022300053, mm*0.457184165, mm*0.216475406, mm*0.167412967,mm*0.066818511, mm*0.003112962, mm*0.00218508,mm*0.0061104];
    nations = randi([5, 14], 1, 3);
    change = randn(1, 3);
    change = .08 * change + .75;
    cap(nations) =  cap(nations) .* change;
    flowaltered(:, 1) = flowaltered(:, 2);
    for i = 4:n
        
        if sitrep(i) > cap(i)
            mcount = mcount + 1;
            m(mcount, :) = m(mcount - 1, :)*Flow^flowcount;
            flowcount = 0;
            
            flowaltered(i, 2) = 1;
            
            Flow(i, node(i).incoming) = (sitrep(i) - cap(i))/sitrep(i) / length(node(i).incoming);
            
            for j = 4:n
                if any(node(i).incoming == j) && all(node(j).incoming ~= i)
                    Flow(j, i) = 0;
                end
                if j == i
                    Flow(i, j) = .5;
                end
            end
            
        elseif flowaltered(i, 2) == 1
            
            flowaltered(i, 2) = 0;
            
            for j = 4:n
                if any(node(i).incoming == j) && all(node(j).incoming ~= i)
                    Flow(i, j) = 0;
                end
            end
            
            for j = 1:n
                Flow(j, :) = (1-Q(j, i))*Flow(j, :);
            end
            Flow(:, i) = Q(:, i);
        end
    end
      
    for j = 1:n
        if sum(Flow(j, :)) > 0
            Flow(j, :) = Flow(j, :)/sum(Flow(j, :));
        end
    end
    
    if any(flowaltered(:, 1) ~= flowaltered(:, 2))
        mcount = mcount + 1;
        m(mcount, :) = m(mcount - 1, :)*Flow^flowcount;
        flowcount = 0;
    end
    flowcount = flowcount + 1;
    sitrep = m(mcount, :)*Flow^flowcount;
    Dist(:, time, rep) = sitrep;
end


a = Dist(:, 12, rep)./cap';
b = Dist(:, 12, rep) > cap';

burden(:, rep) = a(4:n);
capviolation(:, rep) = b(4:n);
burdenvar(rep) = var(burden(:, rep));

end
for i = 1:reptot
    MeanDist = MeanDist + Dist(:, :, i);
end
MeanDist = MeanDist/reptot;
display(mean(burdenvar));