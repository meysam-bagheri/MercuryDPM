function data = read_fstat(filename)

%if no filename is given, pick the first fstat file in the directory
if (~exist('filename','var'))
  filenames = ls('*.fstat -t');
  filename = strtok(filenames);
  disp(filename);
end

%write data into a big array
fid=fopen(filename);
rawdata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
fclose(fid);

%write it into structure
data.t=rawdata{1};
data.PI=rawdata{2};
data.PJ=rawdata{3};
data.centre=[rawdata{4} rawdata{5} rawdata{6}];
data.deltan=rawdata{7};
data.deltat=rawdata{8};
data.forcen=rawdata{9};
data.forcet=rawdata{10};
data.normal=[rawdata{11} rawdata{12} rawdata{13}];
data.tangential=[rawdata{14} rawdata{15} rawdata{16}];

return
data