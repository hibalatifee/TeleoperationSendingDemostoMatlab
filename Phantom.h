#ifndef PHANTOM_H
#define PHANTOM_H

#define _CRT_SECURE_NO_WARNINGS

#include "Global.h"

HDCallbackCode HDCALLBACK DeviceCalibrate(void *pUserData)
{
	hdBeginFrame(Master);
	hdGetDoublev(HD_CURRENT_POSITION,init_position[0]);
	hdEndFrame(hdGetCurrentDevice());

	hdBeginFrame(Slave);
	hdGetDoublev(HD_CURRENT_POSITION,init_position[1]);
	hdEndFrame(hdGetCurrentDevice());
		
	return HD_CALLBACK_DONE;
}


HDCallbackCode HDCALLBACK start_demonstrations(void *pUserData)

{
	if (count1 <2000)
	{
	
		count1++;
	
	}

	getDemonstrations();

	return HD_CALLBACK_CONTINUE;

}

HDCallbackCode HDCALLBACK Start_teleoperation(void *pUserData)

{

	agent();

	number+=1; //Hiba commented this

	if (number>=20000) //Hiba commented this
	{ //Hiba commented this
	
	close(); //Hiba commented this
	} //Hiba commented this

	saveData(); //Hiba commented this
	
	return HD_CALLBACK_CONTINUE;
}


void initiliaze()
{

#ifdef DEVICE_NAME_1
	Master = hdInitDevice(DEVICE_NAME_1);
	hdEnable(HD_FORCE_OUTPUT);
#endif

#ifdef DEVICE_NAME_2
	Slave=hdInitDevice(DEVICE_NAME_2);
	hdEnable(HD_FORCE_OUTPUT);
#endif

	hdStartScheduler();
	
	for (int i = 0; i<3; i++)
	{
		for (int j = 0; j<=2; j++)
		
		{
			force[j][i] = 0.0;
			position[j][i] = 0.0;
			init_position[j][i] = 0.0;
		}
		
		mst1.prev[i] = 0.0;
		mst1.position[i] = 0.0;
		mst1.force[i] = 0.0;
		mst1.disp[i] = 0.0;

		slv.prev[i] = 0.0;
		slv.position[i] = 0.0;
		slv.force[i] = 0.0;
		slv.disp[i] = 0.0;

		diff_e1[i] = 0.0;   prev_e1[i] = 0.0;	
		diff_e2[i] = 0.0;	prev_e2[i] = 0.0;
	}

	count1 = 0;
	count2 = 0;

	Phantom_Started = false;

}

void rotation()

{

	theta1_y*= DEGTORAD;    
	theta2_y*= DEGTORAD;    

	Rotation1(0,0)=cos(theta1_y);      Rotation1(0,1)= 0;	  Rotation1(0,2)=sin(theta1_y);
	Rotation1(1,0)=0;				   Rotation1(1,1)=1;	  Rotation1(1,2)=0;
	Rotation1(2,0)=-sin(theta1_y);     Rotation1(2,1)=0;      Rotation1(2,2)=cos(theta1_y);

	Rotation2(0,0)=cos(theta2_y);      Rotation2(0,1)=0;	  Rotation2(0,2)=sin(theta2_y);
	Rotation2(1,0)=0;				   Rotation2(1,1)=1;	  Rotation2(1,2)=0;
	Rotation2(2,0)=-sin(theta2_y);     Rotation2(2,1)=0;      Rotation2(2,2)=cos(theta2_y);

	R1_inverse=Rotation1.transpose();
	R2_inverse=Rotation2.transpose();

}

void camera1()
{
	
	//VideoCapture cap(0); // open the video camera no. 0

	//int width=1200;
	//int height=800;

	//char* windowName = "Webcam Feed";
	//namedWindow(windowName,0); //create a window to display our webcam feed
	//resizeWindow(windowName,width,height);

	//while (1)
	//
	//{

	//	Mat frame;

	//	bool bSuccess = cap.read(frame); // read a new frame from camera feed

	//	if (!bSuccess) //test if frame successfully read
	//	{
	//		cout << "ERROR READING FRAME FROM CAMERA FEED" << endl;
	//		break;
	//	}

	//}


}


void camera2()
{





}


void calibarate()
{
	if (!Phantom_Started) 
	{
		hdScheduleSynchronous(DeviceCalibrate, 0, HD_MIN_SCHEDULER_PRIORITY);
	}
}

void demonstrations()
{

		
	if (!Phantom_Started)
	
	{
		zmq_bind (publisher, "tcp://*:5556");   
		Phantom_Started = true;
        teleoperation1= hdScheduleAsynchronous(start_demonstrations, 0, HD_DEFAULT_SCHEDULER_PRIORITY);
	}

}

void start()
{
	if (!Phantom_Started)
	{
		zmq_connect (requester, "tcp://localhost:5552");  
		Phantom_Started = true;
        teleoperation2= hdScheduleAsynchronous(Start_teleoperation, 0, HD_DEFAULT_SCHEDULER_PRIORITY);
	}
}


void stopDemonstrations()

{

	//ZMQ STOP REQUESTER
	zmq_close (publisher);
    zmq_ctx_destroy (context2);

	hdStopScheduler();
	hdUnschedule(teleoperation1);
	hdDisableDevice(Master);
	hdDisableDevice(Slave);


}

void close ()
{
	//ZMQ STOP REQUESTER
	zmq_close (requester);
    zmq_ctx_destroy (context1);

	hdStopScheduler();
	hdUnschedule(teleoperation2);
	hdDisableDevice(Master);
	hdDisableDevice(Slave);

	filePosition.close();
	fileForce.close();
	
}

//*****************************AGENT******************************************

inline void agent()
{
	
	hdBeginFrame(Master);

    mst1.prev[0] = mst1.position[0];
	mst1.prev[1] = mst1.position[1];
	mst1.prev[2] = mst1.position[2];

	hdGetDoublev(HD_CURRENT_POSITION, mst1.position);

	mst1.position[0] = mst1.position[0] - init_position[0][0];
	mst1.position[1] = mst1.position[1] - init_position[0][1];
	mst1.position[2] = mst1.position[2] - init_position[0][2];

	mst1.disp[0] =  mst1.position[0] - mst1.prev[0];	
	mst1.disp[1] =  mst1.position[1] - mst1.prev[1];
	mst1.disp[2] =  mst1.position[2] - mst1.prev[2];	

/*******************************************************************************
SLAVE
*******************************************************************************/

	hdBeginFrame(Slave);

	slv.prev[0] = slv.position[0];
	slv.prev[1] = slv.position[1];
	slv.prev[2] = slv.position[2];
	 
	hdGetDoublev(HD_CURRENT_POSITION, slv.position);

	slv.position[0] = slv.position[0] - init_position[1][0];
	slv.position[1] = slv.position[1] - init_position[1][1];
	slv.position[2] = slv.position[2] - init_position[1][2];
	
	
	//ZMQ CLIENT (REQUESTER)
	sprintf_s(REQUEST,"%0.3f %0.3f %0.3f ", mst1.position[0],mst1.position[1],mst1.position[2]);
    zmq_send (requester, REQUEST, 30, 0);
    zmq_recv (requester, REPLY, 30, 0);

	char *token;
	char *rest=REPLY;

	token=strtok_s(rest, " ", &rest);
	px=atof(token);

	token=strtok_s(rest, " ", &rest);
	py=atof(token);
	    
	token=strtok_s(rest, " ", &rest);
	pz=atof(token);

    cout<<"Recieved Position (Px, Py, Pz)"<< px << py << pz<<  endl;


/*******************************************************************************
Forces
*******************************************************************************/
	
	slv.force[0] =  omniSaturation(Kp*(px - slv.position[0]));
	slv.force[1] =  omniSaturation(Kp*(py - slv.position[1]));
	slv.force[2] =  omniSaturation(Kp*(pz - slv.position[2]));
	
    mst1.force[0] =  -1 * omniSaturation(Kp*(mst1.position[0] - slv.position[0]));
	mst1.force[1] =  -1 * omniSaturation(Kp*(mst1.position[1] - slv.position[1]));
	mst1.force[2] =  -1 * omniSaturation(Kp*(mst1.position[2] - slv.position[2]));
	
	
	hdMakeCurrentDevice(Master);
	hdSetDoublev(HD_CURRENT_FORCE, mst1.force);

	hdMakeCurrentDevice(Slave);
	hdSetDoublev(HD_CURRENT_FORCE, slv.force);
	
	hdEndFrame(Master);
	hdEndFrame(Slave);

}

//*************************************demonstrations***************************

inline void  getDemonstrations()

{

/*******************************************************************************
Master 1
*******************************************************************************/
	
	hdBeginFrame(Master);

    mst1.prev[0] = mst1.position[0];
	mst1.prev[1] = mst1.position[1];
	mst1.prev[2] = mst1.position[2];

	hdGetDoublev(HD_CURRENT_POSITION, mst1.position);
	hdGetIntegerv(HD_CURRENT_BUTTONS,&buttonState);

	mst1.position[0] = mst1.position[0] - init_position[0][0];
	mst1.position[1] = mst1.position[1] - init_position[0][1];
	mst1.position[2] = mst1.position[2] - init_position[0][2];

	mst1.disp[0] =  mst1.position[0] - mst1.prev[0];	
	mst1.disp[1] =  mst1.position[1] - mst1.prev[1];
	mst1.disp[2] =  mst1.position[2] - mst1.prev[2];	
	/*
	master1_position_vector(0)=mst1.position[0]; 
	master1_position_vector(1)=mst1.position[1];
	master1_position_vector(2)=mst1.position[2];

	m1_transform_position=master1_position_vector;
	*/
	
	/*
	//ZMQ PUBLISHER
	sprintf(PUBLISH,"%d %d %0.4f %0.4f %0.4f ", topic,buttonState,slv.position[0],slv.position[1],slv.position[2]);
	//sprintf(PUBLISH,"%d %s ",topic,"1");
	zmq_send(publisher,PUBLISH,30,0);
	*/


	/*cout<<slv.position<<endl;*/
/*******************************************************************************
SLAVE
*******************************************************************************/

	hdBeginFrame(Slave);

	slv.prev[0] = slv.position[0];
	slv.prev[1] = slv.position[1];
	slv.prev[2] = slv.position[2];
	 
	hdGetDoublev(HD_CURRENT_POSITION, slv.position);

	slv.position[0] = slv.position[0] - init_position[1][0];
	slv.position[1] = slv.position[1] - init_position[1][1];
	slv.position[2] = slv.position[2] - init_position[1][2];
	
	slv.disp[0] = slv.position[0] - slv.prev[0]; 
	slv.disp[1] = slv.position[1] - slv.prev[1];
	slv.disp[2] = slv.position[2] - slv.prev[2];
/*
	slave_position_vector(0)=slv.position[0]; 
	slave_position_vector(1)=slv.position[1];
	slave_position_vector(2)=slv.position[2];

	slave_transform_position=slave_position_vector;
	*/
/*******************************************************************************
Forces
*******************************************************************************/
	
	slv.force[0] =  omniSaturation(Kp*(mst1.position[0] - slv.position[0]));
	slv.force[1] =  omniSaturation(Kp*(mst1.position[1] - slv.position[1]));
	slv.force[2] =  omniSaturation(Kp*(mst1.position[2] - slv.position[2]));

	//ZMQ PUBLISHER
	sprintf(PUBLISH,"%d %d %0.4f %0.4f %0.4f ", topic,buttonState,slv.position[0],slv.position[1],slv.position[2]);
	//sprintf(PUBLISH,"%d %s ",topic,"1");
	zmq_send(publisher,PUBLISH,30,0);

    mst1.force[0] =  -1 * omniSaturation(Kp*(mst1.position[0] - slv.position[0]));
	mst1.force[1] =  -1 * omniSaturation(Kp*(mst1.position[1] - slv.position[1]));
	mst1.force[2] =  -1 * omniSaturation(Kp*(mst1.position[2] - slv.position[2]));

	hdMakeCurrentDevice(Slave);
	hdSetDoublev(HD_CURRENT_FORCE, slv.force);

	hdMakeCurrentDevice(Master);
	hdSetDoublev(HD_CURRENT_FORCE, mst1.force);

	hdEndFrame(Slave);
	hdEndFrame(Master);
	
}

//*******************************************************************


inline void saveData()

    {

	filePosition<<master1_position_vector(0)<<" "<<master1_position_vector(1)<<" "<<master1_position_vector(2)<<endl;

	fileForce<<mst1.force<<" "<<slv.force<<endl;


	}

double omniSaturation(double force)
{
	if (force<-LIMIT_FORCE_Omni)

		force = -LIMIT_FORCE_Omni;

	if (force>LIMIT_FORCE_Omni)
		
		force = LIMIT_FORCE_Omni;

	return force;
}


double premiumSaturation(double force)
{
	if (force<-LIMIT_FORCE_Premium)

		force = -LIMIT_FORCE_Premium;

	if (force>LIMIT_FORCE_Premium)
		
		force = LIMIT_FORCE_Premium;

	return force;
}



#endif