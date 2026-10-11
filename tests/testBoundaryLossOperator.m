function tests = testBoundaryLossOperator
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.oldPath=path;
addpath(root);
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testBothBoundariesAndSectors(testCase)
for L=1:6
    for N=0:L
        for site=unique([1,L])
            actual=BoundaryLossOperator(L,N,0.7,site);
            expected=legacyBoundary(L,N,0.7,site);
            verifyEqual(testCase,actual,expected);
        end
    end
end
end

function result=legacyBoundary(L,N,gamma,site)
% Frozen reference construction from the original La/Lb routines.
basis=[];
for k1=1:N
    N1=N-k1+1;
    key=nchoosek(1:L,N1);
    basis1=zeros(size(key,1),L);
    for i=1:size(key,1)
        for j=1:N1
            basis1(i,key(i,j))=1;
        end
    end
    basis=[basis;basis1];
end
basis=[basis;zeros(1,L)];
dim=size(basis,1);
result=zeros(dim,dim);
bz=zeros(1,L);
bz(1,site)=1;
for i=1:dim-1
    for j=i+1:dim
        if basis(i,:)-basis(j,:)==bz
            result(j,i)=gamma;
        end
    end
end
end

