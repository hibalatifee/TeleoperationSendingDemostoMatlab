%% Matlab Subscriber for C++ demos using ZMQ 
%%
clc
clear workspace
clear Demo
clear all
close all
addpath(genpath('C:\Users\TeleOperation\Downloads\matlab-zmq-master'))
% addpath('R:\Python_C++_Development\MATLAB\matlab-zmq')
% 
% addpath('C:\Users\Arslan Ali\Desktop')
% load('Traj.mat')
% 

%%
Demo2=zeros(3,100000);
topicfilter = '10001';
port = 5556;
counter1=0;
%counter2=1;

% Socket to talk to server
context = zmq.core.ctx_new();
socket = zmq.core.socket(context, 'ZMQ_SUB');

% Subscribe to the server
fprintf('Collecting updates from Phantom slave...\n');

address = sprintf('tcp://localhost:%d', port);
zmq.core.connect(socket, address);


zmq.core.setsockopt(socket, 'ZMQ_SUBSCRIBE', topicfilter);


%%
% for i=1:3
%     i=i
%     counter1=0;
%     
%     while(1)
%     message=char(zmq.core.recv(socket));
%     parts=strsplit(message);
%     [topic,flag,px,py,pz]=parts{:};
%         
%         if (flag=='2')
%             tic;
%                 while toc < 4;
%                 end
%             break;
%         end
%         
%     counter1=counter1+1;
%     
%     topic=str2double(topic);
%     flag=str2double(flag);
%     Px=str2double(px);
%     Py=str2double(py);
%     Pz=str2double(pz);
% 
%     %fprintf('Recieved Positions (Px, py, pz) %0.4f %0.4f %0.4f on Topic %d\n', Px, Py, Pz,topic);
%     fprintf('Recieved Positions (Px, py, pz) %0.4f %0.4f %0.4f with Flag %d on Topic %d\n', Px, Py, Pz,flag,topic);  
%         
%     Demo(counter2:(counter2)+2,counter1)=[Px Py Pz];
%     end
% flag='1'   
% counter2=counter2+3
%  
% end  

%%    
while(1)
    message=char(zmq.core.recv(socket));
    parts=strsplit(message);
    [topic,flag,px,py,pz]=parts{:};
    
        if (flag=='2')  
        break;
        end
    counter1=counter1+1;
    
    topic=str2double(topic);
    flag=str2double(flag);
    Px=str2double(px);
    Py=str2double(py);
    Pz=str2double(pz);

    %fprintf('Recieved Positions (Px, py, pz) %0.4f %0.4f %0.4f on Topic %d\n', Px, Py, Pz,topic);
    fprintf('Received Positions (px, py, pz) %0.4f %0.4f %0.4f with Flag %d on Topic %d\n', Px, Py, Pz,flag,topic);  
    Demo2(:,counter1)=[Px Py Pz];
end

zmq.core.disconnect(socket, address);

zmq.core.close(socket);

zmq.core.ctx_shutdown(context);
zmq.core.ctx_term(context);

%%
save demo2.mat Demo2;
clear Demo2
load demo2.mat

%%

