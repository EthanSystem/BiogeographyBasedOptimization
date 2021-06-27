clear;
clc;

% A是10行3列的随机数矩阵。其中，10行表示1000个点，每行的3列表示每个点的x,y,z坐标值
A=round(100*rand(10,3));

% B是根据A求出的值。B有10行，每一行表示一个点的函数值。
B=sum(A,2);

% 初始化结构体E和F。
% E有两个字段。字段1：位置数组position 包含x,y,z值；字段2：值cost。E由10个这样的字段构成的元组组成。
% F和E一样。
for i=1:length(A)
    E(i).position=zeros(1,length(A));
    E(i).cost=0;
end
F=E;

% 现在要把A赋值给E和F的position，B赋值给E和F的cost。
% 以下有两种方法：
% 方法一：结果表明赋值成功，但是用了for循环，影响程序运行速度。
tic;
for i=1:length(A)
    E(i).position=A(i,:);
    E(i).cost=B(i);
end
time01=toc

% 方法二：不用for循环，以下是论坛“伏戈”的提供的回答。
tic;
tpt=cellstr(num2str(A));
p=cellfun(@str2num,tpt,'uni',false);
c=num2cell(B);
v=[p,c];
f={'position','cost'};
F=cell2struct(v,f,2);
time02=toc

% 结果表明，方法二在数据量大的情况下，程序运行速度远快于方法一。







