%% 3-DOF Peg-in-Hole task GMR/DMP

%function Teleoperation_3DOF_DMP
addpath(genpath('C:\Users\TeleOperation\Downloads\matlab-zmq-master'))
clc
clear all
close all

%% ZMQ Subscriber

topicfilter1='10001';
port1=5556;

Samples=zeros(3,100000);
Button=zeros (1,100000);

counter1=0;

% Socket to talk to server
context1 = zmq.core.ctx_new();
subscriber = zmq.core.socket(context1, 'ZMQ_SUB');  

% Subscribe to the server
fprintf('Collecting updates from Phantom...\n');

address1 = sprintf('tcp://localhost:%d', port1);
zmq.core.connect(subscriber, address1);

zmq.core.setsockopt(subscriber, 'ZMQ_SUBSCRIBE', topicfilter1);

while (1)
    
        message1=zmq.core.recv(subscriber);
        
        parts=strsplit(char(message1));
        
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
        
        fprintf('Recieved Positions (Px, py, pz) %0.4f %0.4f %0.4f with Flag %d on Topic %d\n', Px, Py, Pz,flag,topic);  
        
        Samples(:,counter1)=[Px Py Pz];
        Button(:,counter1)=flag;

end

% while (1)
%     
%         message1=zmq.core.recv(subscriber);
%         
%         parts=strsplit(char(message1));
%         
%         [topic,num]=parts{:};
% 
%         if (strcmp(num, '1'))  
%             break;
%         end
%         
%         counter1=counter1+1;
%         
%         topic=str2double(topic);
%         num=str2double(num);
%         
%         fprintf('Recieved Topic %d\n', topic);  
%         
%         Samples(:,counter1)=[1 2 3];
%         Button(:,counter1)=num;
% 
% end
   
zmq.core.disconnect(subscriber, address1);

zmq.core.close(subscriber);

zmq.core.ctx_shutdown(context1);
zmq.core.ctx_term(context1); 
 


%% Positions
Position3D_Trial1=Samples(:,1:counter1);
Position3D_Trial1=Position3D_Trial1';
save Position3D_Trial1.mat Position3D_Trial1;


Position3D_Trial2=Samples(:,1:counter1);
Position3D_Trial2=Position3D_Trial2';
save Position3D_Trial2.mat Position3D_Trial2;

%%
Position3D_Trial3=Samples(:,1:counter1);
Position3D_Trial3=Position3D_Trial3';
save Position3D_Trial3.mat Position3D_Trial3;
% 
%% Initilizations

time_duration = 10;
tao = 1;          
Omega = 1/18174;
wn=200;
zeta=0.707;
kv = 2*zeta*wn; 
kp = (wn^2)/kv; 

components=20;

%%
figure;

plot3(Position3D_Trial1(:,1),Position3D_Trial1(:,2), Position3D_Trial1(:,3), 'Color',[0,0,0.9])
hold on;
plot3(Position3D_Trial2(:,1),Position3D_Trial2(:,2), Position3D_Trial2(:,3), 'Color',[0,0.9,0])
hold on
plot3(Position3D_Trial3(:,1),Position3D_Trial3(:,2), Position3D_Trial3(:,3), 'Color',[0.9,0,0])
hold off;
%%
length1=size(Position3D_Trial1,1);   
length2=size(Position3D_Trial2,1);  
length3=size(Position3D_Trial3,1);   

Lengths=[length1 length2 length3];  

SIZE=min(Lengths);  

yy11 = spline(linspace(0,1,length1),Position3D_Trial1(1:end,[1 2 3])',linspace(0,1,SIZE));
yy3=yy11';
yy3_v = (yy3(1:end-1,:)-yy3(2:end,:))/(time_duration/SIZE);
yy3_a = (yy3_v(1:end-1,:)-yy3_v(2:end,:))/(time_duration/SIZE);


yy22 = spline(linspace(0,1,length2),Position3D_Trial2(1:end,[1 2 3])',linspace(0,1,SIZE));
yy4=yy22';
yy4_v = (yy4(1:end-1,:)-yy4(2:end,:))/(time_duration/SIZE);
yy4_a = (yy4_v(1:end-1,:)-yy4_v(2:end,:))/(time_duration/SIZE);


yy33 = spline(linspace(0,1,length3),Position3D_Trial3(1:end,[1 2 3])',linspace(0,1,SIZE));
yy5=yy33';
yy5_v = (yy5(1:end-1,:)-yy5(2:end,:))/(time_duration/SIZE);
yy5_a = (yy5_v(1:end-1,:)-yy5_v(2:end,:))/(time_duration/SIZE);

figure;
grid on;

plot3(yy11(1,:),yy11(2,:),yy11(3,:), 'Color',[0,0,0.9])
hold on;
plot3(yy22(1,:),yy22(2,:), yy22(3,:), 'Color',[0,0.9,0])
plot3(yy33(1,:),yy33(2,:), yy33(3,:), 'Color',[0.9,0,0])
hold off;

%%

length= size(yy3(1:end-2,:),1);  

g = mean([yy3; yy4; yy5],1);

yy3r = (yy3_a - tao*kv*(kp*( repmat(g,length,1) - yy3(1:end-2,:) ) - yy3_v(1:end-1,:) ))./(tao*(repmat(g,length,1)-repmat(yy3(1,:),length,1)) );
yy4r = (yy4_a - tao*kv*(kp*( repmat(g,length,1) - yy4(1:end-2,:) ) - yy4_v(1:end-1,:) ))./(tao*(repmat(g,length,1)-repmat(yy4(1,:),length,1)));
yy5r = (yy5_a - tao*kv*(kp*( repmat(g,length,1) - yy5(1:end-2,:) ) - yy5_v(1:end-1,:) ))./(tao*(repmat(g,length,1)-repmat(yy5(1,:),length,1)));

PositionData = [[linspace(0,(2*pi/SIZE)*length,length)' yy3r];[linspace(0,(2*pi/SIZE)*length,length)' yy4r];[linspace(0,(2*pi/SIZE)*length,length)' yy5r]];

PositionData = PositionData';

Priors = ones(1,components)/components;
ind = ceil(rand(components,1) * size(PositionData,2));
Mu = PositionData(:,ind);
Sigma = repmat(cov(PositionData')/10,1,1,components);

[Priors, Mu, Sigma] = EM_Pol(PositionData, 1, Priors, Mu, Sigma);
%%
figure;
hold on;
for abc=1:size(Mu,2)
    plot_gaussian_ellipsoid(Mu([2 3 4],abc),Sigma([2 3 4],[2 3 4],abc),1);

end 

Clock = 0.0;

%predict clock value for random starting point from the trajectories
% RandNum = randi(3);

% if(RandNum==1)
%     Current_Point = yy3(1,:)';
%     Current_Velocity = (yy3(2,:)-yy3(3,:))'/time_duration;    
% elseif(RandNum==2)
%     Current_Point = yy4(1,:)';
%     Current_Velocity = (yy4(2,:)-yy4(1,:))'/time_duration;
% else
%     Current_Point = yy5(1,:)';
%     Current_Velocity = (yy5(2,:)-yy5(1,:))'/time_duration;
% end

Current_Point = yy5(1,:)';
Current_Velocity = (yy3(2,:)-yy3(3,:))'/time_duration;    
 
Start_Point = Current_Point;

Traj=[];
Traj_Length = 20000;

for j=1:Traj_Length

    [y, Sigma_y] = GMR_Polar(Priors, Mu, Sigma, Clock, 1, [2 3 4],1);    
    
    Acc =  tao*kv*(kp*( g' - Current_Point ) - Current_Velocity ) + tao*(g'-Start_Point).*y;
    Current_Velocity = Current_Velocity + Acc*(time_duration/SIZE);
    Current_Point = Current_Point + Current_Velocity*(time_duration/SIZE);
    
    Traj = [Traj Current_Point];
    Clock = Clock + (2*pi/SIZE);
end

Clock
Traj
figure;
plot3(Traj(1,:),Traj(2,:),Traj(3,:));
hold on;

plot3(Start_Point(1),Start_Point(2),Start_Point(3),'*r','MarkerSize',20);
plot3(Position3D_Trial1(:,1),Position3D_Trial1(:,2),Position3D_Trial1(:,3));
plot3(Position3D_Trial2(:,1),Position3D_Trial2(:,2),Position3D_Trial2(:,3));
plot3(Position3D_Trial3(:,1),Position3D_Trial3(:,2),Position3D_Trial3(:,3));


% Dont use the below code for publishing to C
% ZMQ (REPLY)
%%
% topicfilter2='10002';
% port2=5552;
% 
% %Socket talk to Server
% context2 = zmq.core.ctx_new();
% client = zmq.core.socket(context2, 'ZMQ_REP');   % socket for client
% %client = zmq.core.socket(context2, 'ZMQ_PUB');   % Hiba did this
% 
% % Server talk to client
% fprintf('Waiting for client requests...\n');
% 
% address2=sprintf('tcp://*:%d',port2);
% zmq.core.bind(client, address2);
%     
% 
% for i=1:15
% %for i=1:length(expData)
%     
%     Px=sprintf('%0.4f ',Traj(1,i));
%     Py=sprintf('%0.4f ',Traj(2,i));
%     Pz=sprintf('%0.4f ',Traj(3,i));
%    
%     data=[Px Py Pz];
%     
%     message2=char(zmq.core.recv(client));
%     
%     fprintf('Recieved request from client %s\n',message2);
%     fprintf('Sending Position (Px Py Pz) %0.4f % 0.4f % 0.4f to client\n', Traj(1,i), Traj(2,i), Traj(3,i));
%     
%     zmq.core.send(client, uint8(data));  
%     
% end
% %
% % topicfilter2='10002';
% % port2=5556;
% % 
% % %Socket talk to Server
% % context2 = zmq.core.ctx_new();
% % client = zmq.core.socket(context2, 'ZMQ_REP');   % socket for client
% % 
% % % Server talk to client
% % fprintf('Waiting for client requests...\n');
% % 
% % address2=sprintf('tcp://*:%d',port1);
% % zmq.core.bind(client, address2);
% %     
% % %for i=1:length(expData)
% % for i=1:15
% %     
% %     Px=sprintf('%0.4f ',Traj(1,i));
% %     Py=sprintf('%0.4f ',Traj(2,i));
% %     Pz=sprintf('%0.4f ',Traj(3,i));
% %    
% %     data=[Px Py Pz];
% %     
% %     message2=char(zmq.core.recv(client));
% %     
% %     fprintf('Recieved request from client %s\n',message2);
% %     fprintf('Sending Position (Px Py Pz) %0.4f % 0.4f % 0.4f to client\n', Traj(1,i), Traj(2,i), Traj(3,i));
% %     
% %     zmq.core.send(client, uint8(data));  
% %     
% % end
% %
% zmq.core.disconnect(cleint, address2);
% zmq.core.close(client);
% zmq.core.ctx_shutdown(context2);
% zmq.core.ctx_term(context2); 

%end





