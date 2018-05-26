#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
//#include <opencv\cv.h>
//#include <opencv\highgui.h>
#include <HD/hd.h>
#include <HDU/hduVector.h>
#include <HDU/hduError.h>
#include <zmq.h>
#include <Eigen/Dense>
#include <stdlib.h>


using namespace std;
using namespace Eigen;
//using namespace cv;

#define DEVICE_NAME_1  "Omni"
#define DEVICE_NAME_2  "slave"
#define LIMIT_FORCE_Omni   2.0
#define LIMIT_FORCE_Premium   2.5
#define Kp  0.1
#define RADTODEG	180/3.142
#define DEGTORAD	3.142/180
#define BUFSIZE 256

HHD Master,Slave;

hduVector3Dd position[3], force[3];
hduVector3Dd init_position[3];

ofstream filePosition("Position.txt");
ofstream fileForce("Force.txt");

int count1,count2=0;

bool Phantom_Started=false;

typedef struct _PACK {
	
	hduVector3Dd position;
	hduVector3Dd force;
	hduVector3Dd prev;
	hduVector3Dd disp;
	hduVector3Dd velocity;
	hduVector3Dd angle;
	hduVector3Dd maxforce;

}PACK;

PACK mst1,slv;

void initiliaze();
void calibarate();
void demonstrations();
void start();
void close();
void rotation();
void stopDemonstrations();
void camera1();
void camera2();
void DAD();
void compareEnergy();

double omniSaturation(double force);
double premiumSaturation(double force);

inline void saveData();
inline void getDemonstrations();
inline void agent();

HDSchedulerHandle teleoperation1, teleoperation2;

HDint button;

 // ZMQ Client
 void *context1 = zmq_ctx_new ();
 void *requester = zmq_socket (context1, ZMQ_REQ);

 //ZMQ Publisher
 void *context2 = zmq_ctx_new ();
 void *publisher= zmq_socket (context2,ZMQ_PUB);


 //char str[5]="STOP";

 char REQUEST [30]={};
 char REPLY[30]={};

 char PUBLISH [30]={};
 char BUFFER[10];

 double px=0;
 double py=0;
 double pz=0;
 double theta1_y=270;
 double theta2_y=0;

 double diff_e1[3], diff_e2[3], prev_e1[3], prev_e2[3];

 int number=0;
 int topic=10001;

Matrix3d Rotation1,Rotation2;
Matrix3d R1_inverse, R2_inverse;

Vector3d master1_position_vector,slave_position_vector,m1_transform_position;
Vector3d slave_transform_position;

HDint buttonState;


	//VideoCapture cap(0); // open the video camera no. 0

	//int width=1200;
	//int height=800;

	//char* windowName = "Webcam Feed";
	//namedWindow(windowName,0); //create a window to display our webcam feed
	//resizeWindow(windowName,width,height);

#endif