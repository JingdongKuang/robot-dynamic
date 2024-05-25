%本函数用于获取零相位滤波信号 输入为列向量
function filteredData=ZeroPhaseAverageFilter(width,data)
    filter=ones(width,1)/width;
    filteredData_=data(1:floor(width/2));
    for i=floor(width/2)+1:length(data)-floor(width/2)
        temp=data(i-floor(width/2):i+floor(width/2))'*filter;
        filteredData_=[filteredData_;temp];
    end
    filteredData_=[filteredData_;data(length(data)-floor(width/2)+1:end)];
    data_=flip(filteredData_);
    filteredData=data_(1:floor(width/2));
    for i=floor(width/2)+1:length(data_)-floor(width/2)
        temp=data_(i-floor(width/2):i+floor(width/2))'*filter;
        filteredData=[filteredData;temp];
    end
    filteredData=[filteredData;data_(length(data)-floor(width/2)+1:end)];
    filteredData=flip(filteredData);
end