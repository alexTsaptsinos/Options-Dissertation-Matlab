
exponentialTimes = [];
pairs = [];
indexes = [1:10];

for i = 10:-1:2;
    parameter = nchoosek(i,2);
    expRV = exprnd(parameter);
    exponentialTimes = [exponentialTimes expRV];
    pair = randperm(i,2);
    index1 = pair(1);
    index2 = pair(2);
    pair1 = indexes(index1);
    pair2 = indexes(index2);
    pairs = [pairs; pair1,pair2];
    
    if pair1 > pair2
        indexes(index1) = [];
    else
        indexes(index2) = [];
    end
    
end

exponentialTimes
pairs