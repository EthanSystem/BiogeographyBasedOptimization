clear;
clc;

% A��10��3�е�������������У�10�б�ʾ1000���㣬ÿ�е�3�б�ʾÿ�����x,y,z����ֵ
A=round(100*rand(10,3));

% B�Ǹ���A�����ֵ��B��10�У�ÿһ�б�ʾһ����ĺ���ֵ��
B=sum(A,2);

% ��ʼ���ṹ��E��F��
% E�������ֶΡ��ֶ�1��λ������position ����x,y,zֵ���ֶ�2��ֵcost��E��10���������ֶι��ɵ�Ԫ����ɡ�
% F��Eһ����
for i=1:length(A)
    E(i).position=zeros(1,length(A));
    E(i).cost=0;
end
F=E;

% ����Ҫ��A��ֵ��E��F��position��B��ֵ��E��F��cost��
% ���������ַ�����
% ����һ�����������ֵ�ɹ�����������forѭ����Ӱ����������ٶȡ�
tic;
for i=1:length(A)
    E(i).position=A(i,:);
    E(i).cost=B(i);
end
time01=toc

% ������������forѭ������������̳�����ꡱ���ṩ�Ļش�
tic;
tpt=cellstr(num2str(A));
p=cellfun(@str2num,tpt,'uni',false);
c=num2cell(B);
v=[p,c];
f={'position','cost'};
F=cell2struct(v,f,2);
time02=toc

% ������������������������������£����������ٶ�Զ���ڷ���һ��







