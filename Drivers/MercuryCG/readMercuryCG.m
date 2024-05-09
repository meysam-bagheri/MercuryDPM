function data = MercuryCG(name)
if ~exist('name','var')
    name = 'Chain.stat';
end
[raw, header1, header2] = getRawData(name);
data = readHeader1(header1);
data.name = name;
disp(name)
disp(header1)
%disp(header2)
dim = getDimension(data,raw);
data = getCoordinates(data,raw,dim);
data = getVariables(data,raw,header2,dim);
data = extendVariables(data);
end

function [raw, header1, header2] = getRawData(name)
raw = importdata(name,' ',2);
header1 = raw.textdata{1};
if size(raw.textdata,2)==1
    header2 = textscan(raw.textdata{2},'%s');
    header2 = header2{1};
else
    header2 = raw.textdata(2,:);
end
raw = raw.data;
end

function dim = getDimension(data,raw)
%order: z,y,x,t (such that you can plot surface plots with coordinates ordered as t,x,y,z)
dim = [data.n(end:-1:1)' size(raw,1)/prod(data.n)];
dim(dim==1)=[];
dim = [dim 1];
end

function data = getCoordinates(data,raw,dim)
data.t = reshape(raw(:,1),dim);
coord='xyz';
coord(data.n==1)=[];
for i=1:length(coord)
    data.(coord(i)) = reshape(raw(:,i+1),dim);
end
end

function data = getVariables(data,raw,names,dim)
contains = @(str, pattern) ~cellfun('isempty', strfind(str, pattern));
isVariable = contains(names,':');
for i=1:length(names)
    if (isVariable(i))
        colon = strfind(names{i},':');
        hyphen = strfind(names{i},'-');
        a = str2double(names{i}(1:hyphen-1));
        b = str2double(names{i}(hyphen+1:colon-1));
        n = names{i}(colon+1:end);
        if isempty(hyphen)
            a = str2double(names{i}(1:colon-1));
            data.(n) = reshape(raw(:,a),dim);
        elseif b-a==2
            data.([n 'X']) = reshape(raw(:,a),dim);
            data.([n 'Y']) = reshape(raw(:,a+1),dim);
            data.([n 'Z']) = reshape(raw(:,a+2),dim);
        elseif b-a==5 && ~strcmp(n,'ParticleSize')
            data.([n 'XX']) = reshape(raw(:,a),dim);
            data.([n 'XY']) = reshape(raw(:,a+1),dim);
            data.([n 'XZ']) = reshape(raw(:,a+2),dim);
            data.([n 'YY']) = reshape(raw(:,a+3),dim);
            data.([n 'YZ']) = reshape(raw(:,a+4),dim);
            data.([n 'ZZ']) = reshape(raw(:,a+5),dim);
        elseif b-a==8
            data.([n 'XX']) = reshape(raw(:,a),dim);
            data.([n 'XY']) = reshape(raw(:,a+1),dim);
            data.([n 'XZ']) = reshape(raw(:,a+2),dim);
            data.([n 'YX']) = reshape(raw(:,a+3),dim);
            data.([n 'YY']) = reshape(raw(:,a+4),dim);
            data.([n 'YZ']) = reshape(raw(:,a+5),dim);
            data.([n 'ZX']) = reshape(raw(:,a+6),dim);
            data.([n 'ZY']) = reshape(raw(:,a+7),dim);
            data.([n 'ZZ']) = reshape(raw(:,a+8),dim);
        else
            for j=a:b
                data.([n num2str(j-a)]) = reshape(raw(:,j),dim);
            end
        end
    end
end
end

function data = extendVariables(data)
if isfield(data,'ParticleSize0')
    data.ParticleNumber = data.ParticleSize0;
    mean = data.ParticleSize1./data.ParticleNumber;
    % get momenta
    mom2 = data.ParticleSize2./data.ParticleNumber;
    mom3 = data.ParticleSize3./data.ParticleNumber;
    mom4 = data.ParticleSize4./data.ParticleNumber;
    mom5 = data.ParticleSize5./data.ParticleNumber;
    % get central momenta
    mom5 = mom5 -5*mean.*mom4 +10*mean.^2.*mom3 -10*mean.^3.*mom2 +4*mean.^5;
    mom4 = mom4 -4*mean.*mom3 + 6*mean.^2.*mom2 - 3*mean.^4;
    mom3 = mom3 -3*mean.*mom2 + 2*mean.^3;
    mom2 = mom2 -mean.^2;
    std = sqrt(mom2);
    % get standardised momenta
    data.ParticleSize.Mean = mean;
    data.ParticleSize.Variance = mom2;
    data.ParticleSize.Skewness = mom3./std.^3;
    data.ParticleSize.Kurtosis = mom4./std.^4;
    data.ParticleSize.Hyperskewness = mom5./std.^5;
    % reset inf to zero
%     data.ParticleSize.Mean(isnan(data.ParticleSize.Mean)) = 0;
%     data.ParticleSize.Std(isnan(data.ParticleSize.Std)) = 0;
%     data.ParticleSize.Skewness(isnan(data.ParticleSize.Skewness)) = 0;
%     data.ParticleSize.Kurtosis(isnan(data.ParticleSize.Kurtosis)) = 0;
%     data.ParticleSize.Hyperskewness(isnan(data.ParticleSize.Hyperskewness)) = 0;
    %remove original fields
%     data = rmfield(data,'ParticleSize0');
%     data = rmfield(data,'ParticleSize1');
%     data = rmfield(data,'ParticleSize2');
%     data = rmfield(data,'ParticleSize3');
%     data = rmfield(data,'ParticleSize4');
%     data = rmfield(data,'ParticleSize5');
end
if isfield(data,'Density')
    data.VelocityX = data.MomentumX./data.Density;
    data.VelocityY = data.MomentumY./data.Density;
    data.VelocityZ = data.MomentumZ./data.Density;
    data.VelocityX(data.Density==0)=0;
    data.VelocityY(data.Density==0)=0;
    data.VelocityZ(data.Density==0)=0;
    data.KineticStressXX = data.MomentumFluxXX - data.MomentumX.*data.VelocityX;
    data.KineticStressXY = data.MomentumFluxXY - data.MomentumX.*data.VelocityY;
    data.KineticStressXZ = data.MomentumFluxXZ - data.MomentumX.*data.VelocityZ;
    data.KineticStressYY = data.MomentumFluxYY - data.MomentumY.*data.VelocityY;
    data.KineticStressYZ = data.MomentumFluxYZ - data.MomentumY.*data.VelocityZ;
    data.KineticStressZZ = data.MomentumFluxZZ - data.MomentumZ.*data.VelocityZ;
    data.StressXX = data.ContactStressXX + data.KineticStressXX;
    data.StressXY = data.ContactStressXY + data.KineticStressXY; 
    data.StressXZ = data.ContactStressXZ + data.KineticStressXZ;
    data.StressYX = data.ContactStressYX + data.KineticStressXY;
    data.StressYY = data.ContactStressYY + data.KineticStressYY;
    data.StressYZ = data.ContactStressYZ + data.KineticStressYZ;
    data.StressZX = data.ContactStressZX + data.KineticStressXZ;
    data.StressZY = data.ContactStressZY + data.KineticStressYZ;
    data.StressZZ = data.ContactStressZZ + data.KineticStressZZ;
    data.ContactPressure = (data.ContactStressXX + data.ContactStressYY + data.ContactStressZZ)/3;
    data.KineticPressure = (data.KineticStressXX + data.KineticStressYY + data.KineticStressZZ)/3;
    data.FluctuatingKineticEnergy = data.KineticStressXX + data.KineticStressYY + data.KineticStressZZ;
    data.KineticEnergy = data.MomentumFluxXX + data.MomentumFluxYY + data.MomentumFluxZZ;
    data.Pressure = (data.StressXX + data.StressYY + data.StressZZ)/3;
    data.Temperature = data.KineticPressure./data.Density;
end
end

function data = readHeader1(header1)
opt = textscan(header1,'%s');
opt = opt{1};
i = find(strcmp(opt,'n'));
data.n = str2double(opt(i+1:i+3));
i = find(strcmp(opt,'width'));
data.width = str2double(opt(i+1));
i = find(strcmp(opt,'timeMin'));
data.timeMin = str2double(opt(i+1));
i = find(strcmp(opt,'min'));
data.min = str2double(opt(i+1:i+3));
i = find(strcmp(opt,'max'));
data.max = str2double(opt(i+1:i+3));
end
