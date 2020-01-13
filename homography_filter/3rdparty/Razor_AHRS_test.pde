
/// include library
import java.awt.Frame; 
import gab.opencv.*;
import java.nio.*;
import java.util.*;
import org.opencv.core.Core;
import org.opencv.core.Mat;
import org.opencv.core.CvType;
import org.opencv.core.Scalar;
import org.opencv.core.Point;
import org.opencv.calib3d.Calib3d;
import org.opencv.core.MatOfPoint2f;
import processing.video.*;
import java.io.*;
import  java.lang.Math;
//import processing.opengl.*;
import processing.serial.*;
/////////////////////////////////////////////////
/// parameter settings for the I2A optical flow algorithm
float FPS = 30;
int REF1 = 15;
int REF2 = 7;
int SAMP1 = 4;
int SAMP2 = 3;
int PSIZE = 16;
int WIDTH = 320;
int HEIGHT = 240;
int XPATCHES = 10;
int YPATCHES = 8;
int XSPACING = 14;
int YSPACING = 14;
int XCENTRE;
int YCENTRE;
int PATCHSIZE_X;  
int PATCHSIZE_Y;  
int NUM_PATCHES;
int XMAX1, YMAX1;
int XMAX2, YMAX2;
////////////////
float[] ave_of;
float[] med_of;
float[] ave_ofwp;
float[] camAngular;
float[] rotGyr2cam;
//
int[][] imgmid; 
int[][] fltpre31; // filtered image using average filter
int[][] fltcur31;
int[][] fltpre15;
int[][] fltcur15;
int[][] gWarpimg15; // warped image using gyro measurements
int[][] gWarpimg31;
int[][] lutGwarp1x; // look up table for gyro warp
int[][] lutGwarp2x;
int[][] lutGwarp3x;
int[][] lutGwarp4x;
int[][] lutGwarp1y;
int[][] lutGwarp2y;
int[][] lutGwarp3y;
int[][] lutGwarp4y;
float[][] imgpts0x; // image point locations
float[][] imgpts0y;
float[][] imgpts1x;
float[][] imgpts1y;
float[][] imgpts2x;
float[][] imgpts2y;
float[][] imgpts3x;
float[][] imgpts3y;
float[][] imgpts4x;
float[][] imgpts4y;
///sensor data
float accx = 0.0f;
float accy = 0.0f;
float accz = 0.0f;
float magx = 0.0f;
float magy = 0.0f;
float magz = 0.0f;
float gyrx = 0.0f;
float gyry = 0.0f;
float gyrz = 0.0f;
float gxoff = 0.0f;
float gyoff = 0.0f;
float gzoff = 0.0f;
float gxc, gyc, gzc;
/////////////////////////
Mat Hest;
PFrame f;
Mat goodMatch;
Capture frmcap;
secondApplet sndWin;
File flowGyroFile;
File medofGyroFile;
FileWriter flowGyroWrite;
FileWriter medofGyroWrite;
float outLierMask = 1111.1111;
float gyrtoflowScale = 11.5;
float angarrow;
float arrowx1, arrowy1;
float arrowx2, arrowy2;
int inNum = 0;
int keyState = -1;
int mofClosed = 1;
int flowClosed = 1;
int mofGyrLogNum = 0;
int flowGyrLogNum = 0;
int offsetcount = 0;
int frameratetest = 1;
float gyrCaliNum = 30;
int saveFlowImg = 0;
float gzDisct;
/////////////////////
PFont font;
Serial serial;
boolean synched = false;
final static int SERIAL_PORT_NUM = 1;
final static int SERIAL_PORT_BAUD_RATE = 57600;
///////////////////////////////////////////////////////////////////////////

void setup() 
{  
  System.loadLibrary(Core.NATIVE_LIBRARY_NAME);// for java opencv
  //PFrame f = new PFrame();//////////////////////////////
  XCENTRE = WIDTH/2;
  YCENTRE = HEIGHT/2;
  NUM_PATCHES = XPATCHES*YPATCHES; 
  ave_of = new float[3];
  med_of = new float[3];
  ave_ofwp = new float[3];
  camAngular = new float[3];
  camAngular[0] = 0;
  camAngular[1] = 0;
  camAngular[2] = 0;
  rotGyr2cam = new float[9];
  imgmid = new int[WIDTH][HEIGHT];
  fltpre31 = new int[WIDTH][HEIGHT];
  fltpre15 = new int[WIDTH][HEIGHT];
  fltcur31 = new int[WIDTH][HEIGHT];
  fltcur15 = new int[WIDTH][HEIGHT]; 
  gWarpimg15 = new int[WIDTH][HEIGHT];
  gWarpimg31 = new int[WIDTH][HEIGHT];
  lutGwarp1x = new int[WIDTH][HEIGHT]; 
  lutGwarp2x = new int[WIDTH][HEIGHT]; 
  lutGwarp3x = new int[WIDTH][HEIGHT]; 
  lutGwarp4x = new int[WIDTH][HEIGHT]; 
  lutGwarp1y = new int[WIDTH][HEIGHT]; 
  lutGwarp2y = new int[WIDTH][HEIGHT]; 
  lutGwarp3y = new int[WIDTH][HEIGHT]; 
  lutGwarp4y = new int[WIDTH][HEIGHT];
  imgpts0x = new float[XPATCHES][YPATCHES];
  imgpts0y = new float[XPATCHES][YPATCHES];
  imgpts1x = new float[XPATCHES][YPATCHES];
  imgpts1y = new float[XPATCHES][YPATCHES];
  imgpts2x = new float[XPATCHES][YPATCHES];
  imgpts2y = new float[XPATCHES][YPATCHES];
  imgpts3x = new float[XPATCHES][YPATCHES];
  imgpts3y = new float[XPATCHES][YPATCHES];
  imgpts4x = new float[XPATCHES][YPATCHES];
  imgpts4y = new float[XPATCHES][YPATCHES];
  Hest = new Mat(3, 3, CvType.CV_32FC1);
  goodMatch = new Mat(NUM_PATCHES, 1, CvType.CV_32FC1);
  initSystem(imgpts0x, imgpts0y, rotGyr2cam);// generage feature locations
  constuctLUT(lutGwarp1x, lutGwarp1y, -8);
  constuctLUT(lutGwarp2x, lutGwarp2y, -4);
  constuctLUT(lutGwarp3x, lutGwarp3y,  4);
  constuctLUT(lutGwarp4x, lutGwarp4y,  8);
  copyPoints(imgpts0x, imgpts0y, imgpts1x, imgpts1y);
  copyPoints(imgpts0x, imgpts0y, imgpts2x, imgpts2y);
  copyPoints(imgpts0x, imgpts0y, imgpts3x, imgpts3y);
  copyPoints(imgpts0x, imgpts0y, imgpts4x, imgpts4y);
  /// open camera
  //frmcap = new Capture(this, WIDTH, HEIGHT, 30);
  //frmcap.start(); 
  String[] cameras = Capture.list(); 
  if (cameras.length == 0) {
    println("There are no cameras available for capture.");
    exit();
  } else {
    println("Available cameras:");
    for (int i = 0; i < cameras.length; i++) {
      println(cameras[i]);
    }
    frmcap = new Capture(this, cameras[2]);//15,2
    frmcap.start();
  }
  ////////////////////////////////////////////
  /// for inertial sensor
  size(WIDTH*2+0, HEIGHT);//////
  //frameRate(30); ///////////////////////
  font = loadFont("Univers-66.vlw");
  textFont(font);
  println("AVAILABLE SERIAL PORTS:");
  println(Serial.list());
  String portName = Serial.list()[SERIAL_PORT_NUM];
  println();
  println("HAVE A LOOK AT THE LIST ABOVE AND SET THE RIGHT SERIAL PORT NUMBER IN THE CODE!");
  println("  -> Using port " + SERIAL_PORT_NUM + ": " + portName);
  serial = new Serial(this, portName, SERIAL_PORT_BAUD_RATE);
}
/*
////////////////////////////////////////////////
// process image stream
void captureEvent(Capture frmcap) {
  int i, j, k;
  if (frmcap.available() == true) { 
    frmcap.read();
    frmcap.filter(GRAY); // convert color to gray image       
  } //frmcap.available() == true
}
*/
///////////////////////
// get raw IMU data
void draw() {
  //background(10); 
  if (!synched) {
    textFont(font, 20);
    textAlign(CENTER);
    fill(255);
    text("Connecting to Razor...", width/2, height/2, -200);
    if (frameCount == 2)
      setupRazor();  // Set ouput params and request synch token
    else if (frameCount > 2)
      synched = readToken(serial, "#SYNCH00\r\n");  // Look for synch token
    return;
  }
  // read raw sensor data
  while (serial.available () >= 36) {
    accx = readFloat(serial); 
    accy = readFloat(serial); 
    accz = readFloat(serial); 
    magx = readFloat(serial); 
    magy = readFloat(serial); 
    magz = readFloat(serial); 
    gyrx = readFloat(serial); 
    gyry = readFloat(serial); 
    gyrz = readFloat(serial); 
    accx = accx/256;
    accy = accy/256;
    accz = accz/256;
    magx = magx/256;
    magy = magy/256;
    magz = magz/256;
    gyrx = gyrx/256;
    gyry = gyry/256;
    gyrz = gyrz/256;
  } 
  /*
 // display sensor data 
  textFont(font, 15);
  textAlign(LEFT);
  text("ax: " + String.format("%4.2f", accx), WIDTH/4-30, HEIGHT/5);
  text("ay: " + String.format("%4.2f", accy), WIDTH*2/4-30, HEIGHT/5);
  text("az: " + String.format("%4.2f", accz), WIDTH*3/4-30, HEIGHT/5);
  text("gx: " + String.format("%4.2f", gyrx-gxoff), WIDTH/4-30, HEIGHT*2/5);
  text("gy: " + String.format("%4.2f", gyry-gyoff), WIDTH*2/4-30, HEIGHT*2/5);
  text("gz: " + String.format("%4.2f", gyrz-gzoff), WIDTH*3/4-30, HEIGHT*2/5);
  text("medianOfx: " + String.format("%4.2f", med_of[0]), WIDTH/4-30, HEIGHT*3/5);
  text("medianOfy: " + String.format("%4.2f", med_of[1]), WIDTH/4-30, HEIGHT*4/5);*/
  ////////////////////////////////////////////////////////////////////////////
  int i, j, k;
  int curtime;
  float elaptime;  
  if (frmcap.available() == true) { 
    frmcap.read();
    frmcap.filter(GRAY); // convert color to gray image 
    curtime = millis();   
    image(frmcap, 0, 0);  
    k = 0;   
    float sft = 0;
    strokeWeight(2); 
    textSize(20);
    fill(0, 250, 0);
    text("standalone OF", 80, 25);
    for (i=0; i<XPATCHES; i++) {
      for (j=0; j<YPATCHES; j++) {            
        stroke(0, 250, 0);
        line(imgpts3x[i][j]+sft, imgpts3y[i][j]+sft, imgpts4x[i][j]+sft, imgpts4y[i][j]+sft);
        angarrow = (float) Math.atan2(imgpts3y[i][j]-imgpts4y[i][j], imgpts3x[i][j]-imgpts4x[i][j]);
        arrowx1 = (float) (imgpts4x[i][j] + 3*Math.cos(angarrow + PI/4));
        arrowy1 = (float) (imgpts4y[i][j] + 3*Math.sin(angarrow + PI/4));
        line(arrowx1+sft, arrowy1+sft, imgpts4x[i][j]+sft, imgpts4y[i][j]+sft);
        arrowx2 = (float) (imgpts4x[i][j] + 3*Math.cos(angarrow - PI/4));
        arrowy2 = (float) (imgpts4y[i][j] + 3*Math.sin(angarrow - PI/4));
        line(arrowx2+sft, arrowy2+sft, imgpts4x[i][j]+sft, imgpts4y[i][j]+sft);
        k++;                
      }
    }      
    k = 0;  
    translate(WIDTH,0);   
    image(frmcap, 0, 0);   
    textSize(20);
    fill(0, 0, 200);
    text("gyro-aided OF", 80, 25);    
    for (i=0; i<XPATCHES; i++) {
      for (j=0; j<YPATCHES; j++) {    
        stroke(0, 0, 200);
        line(imgpts0x[i][j], imgpts0y[i][j], imgpts2x[i][j], imgpts2y[i][j]);
        angarrow = (float) Math.atan2(imgpts0y[i][j]-imgpts2y[i][j], imgpts0x[i][j]-imgpts2x[i][j]);
        arrowx1 = (float) (imgpts2x[i][j] + 3*Math.cos(angarrow + PI/4));
        arrowy1 = (float) (imgpts2y[i][j] + 3*Math.sin(angarrow + PI/4));
        line(arrowx1, arrowy1, imgpts2x[i][j], imgpts2y[i][j]);
        arrowx2 = (float) (imgpts2x[i][j] + 3*Math.cos(angarrow - PI/4));
        arrowy2 = (float) (imgpts2y[i][j] + 3*Math.sin(angarrow - PI/4));
        line(arrowx2, arrowy2, imgpts2x[i][j], imgpts2y[i][j]);
        k++;                
      }
    }
    //sndWin.redraw(); 
    ////////////////////////////////////////////////////////////
    preFilter3(frmcap, fltcur31, fltcur15);
    aveFilter(fltcur31, fltcur31, 29); // image filtering
    aveFilter(fltcur15, fltcur15, 13);
    copyPoints(imgpts0x, imgpts0y, imgpts1x, imgpts1y);   
    copyPoints(imgpts0x, imgpts0y, imgpts2x, imgpts2y);   
    copyPoints(imgpts0x, imgpts0y, imgpts3x, imgpts3y);
    copyPoints(imgpts0x, imgpts0y, imgpts4x, imgpts4y);
    if (abs(camAngular[2]) < 2){
      gyroWarpPoints(imgpts1x, imgpts1y, imgpts2x, imgpts2y, camAngular, 0);///
      I2A_OF(fltpre31, fltcur31, imgpts1x, imgpts1y, imgpts2x, imgpts2y, REF1, SAMP1);      
      I2A_OF(fltpre15, fltcur15, imgpts1x, imgpts1y, imgpts2x, imgpts2y, REF2, SAMP2); 
    }    
    else{
      if(camAngular[2]<-6) gzDisct = -8;
      if((camAngular[2]>=-6)&&(camAngular[2]<=-2)) gzDisct = -4;
      if((camAngular[2]>=2)&&(camAngular[2]<=6)) gzDisct = 4;
      if(camAngular[2]>6) gzDisct = 8;
      gyroWarpPoints(imgpts1x, imgpts1y, imgpts2x, imgpts2y, camAngular, gzDisct);
      if(gzDisct==-8) yawWarpimg(gWarpimg31, gWarpimg15, fltpre31, fltpre15, lutGwarp1x, lutGwarp1y);
      if(gzDisct==-4) yawWarpimg(gWarpimg31, gWarpimg15, fltpre31, fltpre15, lutGwarp2x, lutGwarp2y);
      if(gzDisct== 4) yawWarpimg(gWarpimg31, gWarpimg15, fltpre31, fltpre15, lutGwarp3x, lutGwarp3y);
      if(gzDisct== 8) yawWarpimg(gWarpimg31, gWarpimg15, fltpre31, fltpre15, lutGwarp4x, lutGwarp4y);
      I2A_OF(gWarpimg31, fltcur31, imgpts1x, imgpts1y, imgpts2x, imgpts2y, REF1, SAMP1);      
      I2A_OF(gWarpimg15, fltcur15, imgpts1x, imgpts1y, imgpts2x, imgpts2y, REF2, SAMP2);   
    }
    //elaptime = (float)(millis() - curtime);
    //println(elaptime);
    copyPoints(imgpts3x, imgpts3y, imgpts4x, imgpts4y);     
    I2A_OF(fltpre31, fltcur31, imgpts3x, imgpts3y, imgpts4x, imgpts4y, REF1, SAMP1);      
    I2A_OF(fltpre15, fltcur15, imgpts3x, imgpts3y, imgpts4x, imgpts4y, REF2, SAMP2); 
    //Hest = calHomography(imgpts1x, imgpts1y, imgpts2x, imgpts2y, goodMatch);/////
    //if (inNum>NUM_PATCHES/2) projectPoints(imgpts1x, imgpts1y, imgpts2x, imgpts2y, Hest);
    calAveFlow(imgpts0x, imgpts0y, imgpts2x, imgpts2y, ave_ofwp);
    calAveFlow(imgpts3x, imgpts3y, imgpts4x, imgpts4y, ave_of);
    //calMedFlow(imgpts1x, imgpts1y, imgpts2x, imgpts2y, med_of);
    imgmid = fltpre31; // sawp the pointer for previous and current image  
    fltpre31 = fltcur31;  
    fltcur31 = imgmid;
    imgmid = fltpre15;  
    fltpre15 = fltcur15;  
    fltcur15 = imgmid;
    if (keyState == 1) { // calibrate gyroscope offset
      offsetcount++;
      gxoff += gyrx;
      gyoff += gyry;
      gzoff += gyrz;
      if (offsetcount==gyrCaliNum) {
        keyState = 0;
        offsetcount = 0;
        gxoff = gxoff/gyrCaliNum;
        gyoff = gyoff/gyrCaliNum;
        gzoff = gzoff/gyrCaliNum;
        println("gyro offset calibrated");
      }
    }
    if (keyState == 2) { // record median flow and gyro data
      try {
        medofGyroWrite.write(String.valueOf(ave_of[0])); // elaptime
        medofGyroWrite.write("\r\n");   
        medofGyroWrite.write(String.valueOf(ave_of[1]));
        medofGyroWrite.write("\r\n"); 
        medofGyroWrite.write(String.valueOf(ave_of[2]));
        medofGyroWrite.write("\r\n"); 
        medofGyroWrite.write(String.valueOf(ave_ofwp[0])); // 
        medofGyroWrite.write("\r\n");   
        medofGyroWrite.write(String.valueOf(ave_ofwp[1]));
        medofGyroWrite.write("\r\n"); 
        medofGyroWrite.write(String.valueOf(ave_ofwp[2]));
        medofGyroWrite.write("\r\n"); 
        medofGyroWrite.write(String.valueOf(camAngular[0]));
        medofGyroWrite.write("\r\n"); 
        medofGyroWrite.write(String.valueOf(camAngular[1]));
        medofGyroWrite.write("\r\n"); 
        medofGyroWrite.write(String.valueOf(camAngular[2]));
        medofGyroWrite.write("\r\n");
      } 
      catch (Exception e) {
        println("Fail to Log median flow and gyro data");
      }
    }
    if ((keyState == 3) && (inNum>NUM_PATCHES*4/5)) { // record flow field and gyro data
      try {
        flowGyroWrite.write(String.valueOf(gyrx));
        flowGyroWrite.write("\r\n"); 
        flowGyroWrite.write(String.valueOf(gyry));
        flowGyroWrite.write("\r\n"); 
        flowGyroWrite.write(String.valueOf(gyrz));
        flowGyroWrite.write("\r\n");
        k = 0;
        for (i=0; i<XPATCHES; i++) {
          for (j=0; j<YPATCHES; j++) {
            if (goodMatch.get(k, 0)[0] == 1.0) {
              flowGyroWrite.write(String.valueOf(imgpts2x[i][j]-imgpts1x[i][j]));
              flowGyroWrite.write("\r\n");
              flowGyroWrite.write(String.valueOf(imgpts2y[i][j]-imgpts1y[i][j]));
              flowGyroWrite.write("\r\n");
            } else {
              flowGyroWrite.write(String.valueOf(outLierMask)); // outlier
              flowGyroWrite.write("\r\n");
              flowGyroWrite.write(String.valueOf(outLierMask));  
              flowGyroWrite.write("\r\n");
            }
            k++;
          }
        }
      }
      catch (Exception e) {
        println("Fail to Log median flow and gyro data");
      }
    }
    // calculate the angular rate of camera from gyro readings
    camAngular[0] = (rotGyr2cam[0]*(gyrx-gxoff) + rotGyr2cam[1]*(gyry-gyoff) + rotGyr2cam[2]*(gyrz-gzoff));
    camAngular[1] = (rotGyr2cam[3]*(gyrx-gxoff) + rotGyr2cam[4]*(gyry-gyoff) + rotGyr2cam[5]*(gyrz-gzoff));
    camAngular[2] = (rotGyr2cam[6]*(gyrx-gxoff) + rotGyr2cam[7]*(gyry-gyoff) + rotGyr2cam[8]*(gyrz-gzoff));
    if (keyState == 2) {  
      frameratetest++;
      if (frameratetest==8) {
        saveFlowImg++;
        frameratetest = 0;
        //save("flowimg"+String.format("%02d", saveFlowImg)+".jpg");
        //saveFrame();
      }
    }  
  } //frmcap.available() == true
}

////////////////////////////////////////////////
//////////////////////////////////////////////////////
// run when a key is pressed and released
void keyReleased() {
  if (key == 99) { // press c for calibrating gyro offset
    keyState = 1;
    gxoff = 0;
    gyoff = 0;
    gzoff = 0;
    println("calibrating gyro offset \r\n");
    println("please do not move the system!");
  }  
  if ((key == 101) && (mofClosed==1)) { // press e for extrinsic calibration
    keyState = 2;
    mofClosed = 0;
    mofGyrLogNum++;
    medofGyroFile = new File("medofgyrLog"+String.format("%02d", mofGyrLogNum)+".txt");
    try {
      medofGyroWrite = new FileWriter(medofGyroFile);
      println("FileWriter created, start logging median flow and gyro data");
    } 
    catch (Exception e) {
      println("Create FileWriter Failed");
    }
  }
  if ((key == 102) && (flowClosed==1)) { // press f to save whole optical flow field and gyro     
    keyState = 3;
    flowClosed = 0;
    flowGyrLogNum++;
    flowGyroFile = new File("flowgyrLog"+String.format("%02d", flowGyrLogNum)+".txt");
    try {
      flowGyroWrite = new FileWriter(flowGyroFile);
      println("FileWriter created, start logging flow field and gyro data");
    } 
    catch (Exception e) {
      println("Create FileWriter Failed");
    }
  }
  if (key == 116) { // press t to terminate data logging 
    keyState = 0;
    try {
      if (flowClosed==0) {        
        flowGyroWrite.write(String.valueOf(gxoff));
        flowGyroWrite.write("\r\n");
        flowGyroWrite.write(String.valueOf(gyoff));
        flowGyroWrite.write("\r\n");
        flowGyroWrite.write(String.valueOf(gzoff));
        flowGyroWrite.write("\r\n");
        flowGyroWrite.close();
        println("FileWrite Closed, Stop Logging");
      }
      if (mofClosed==0) {
        medofGyroWrite.write(String.valueOf(gxoff));
        medofGyroWrite.write("\r\n");
        medofGyroWrite.write(String.valueOf(gyoff));
        medofGyroWrite.write("\r\n");
        medofGyroWrite.write(String.valueOf(gzoff));
        medofGyroWrite.write("\r\n");
        medofGyroWrite.close();
        println("FileWrite Closed, Stop Logging");
      }
      mofClosed = 1;
      flowClosed = 1;
    } 
    catch (Exception e) {
      println("Close File Failed");
    }
  }
}

/////////////////////////////////////////////////////////
// calculate homography matrix using optica flow
Mat calHomography(float[][] p1x, float[][] p1y, float[][] p2x, float[][] p2y, Mat inlier)
{
  int row, col;
  ArrayList<Point> pts1 = new ArrayList<Point>();
  ArrayList<Point> pts2 = new ArrayList<Point>(); 

  int k = 0;
  for (col=0; col<XPATCHES; col++) {
    for (row=0; row<YPATCHES; row++) {      
      pts1.add(new Point(p1x[col][row], p1y[col][row]));
      pts2.add(new Point(p2x[col][row], p2y[col][row]));
      k++;
    }
  }  
  MatOfPoint2f mp1 = new MatOfPoint2f();
  mp1.fromList(pts1);
  MatOfPoint2f mp2 = new MatOfPoint2f(); 
  mp2.fromList(pts2);   
  Mat H = Calib3d.findHomography(mp1, mp2, Calib3d.RANSAC, 3, inlier);
  inNum = 0;
  for (k=0; k<NUM_PATCHES; k++) {
    if (inlier.get(k, 0)[0]==1.0) inNum++;
  }    
  return H;
}

////////////////////////////////////////
// project the image points in the first view to the second view using the homography matrix
void  projectPoints(float[][] p1x, float[][] p1y, float[][] p2x, float[][] p2y, Mat H)
{
  int row, col;
  double h1[]; 
  double h2[]; 
  double h3[];
  double h4[]; 
  double h5[]; 
  double h6[];
  double h7[]; 
  double h8[]; 
  double h9[];  
  double hx, hy, hz;
  for (col=0; col<XPATCHES; col++) {
    for (row=0; row<YPATCHES; row++) {      
      h1 = H.get(0, 0);  
      h2 = H.get(0, 1);  
      h3 = H.get(0, 2);
      h4 = H.get(1, 0);  
      h5 = H.get(1, 1);  
      h6 = H.get(1, 2);
      h7 = H.get(2, 0);  
      h8 = H.get(2, 1);  
      h9 = H.get(2, 2);
      hx = h1[0]*p1x[col][row] + h2[0]*p1y[col][row] + h3[0];
      hy = h4[0]*p1x[col][row] + h5[0]*p1y[col][row] + h6[0];
      hz = h7[0]*p1x[col][row] + h8[0]*p1y[col][row] + h9[0];
      p2x[col][row] = (float) (hx/hz);
      p2y[col][row] = (float) (hy/hz);
    }
  }
}

////////////////////
// construct look up table for gyro-aided image warping
void constuctLUT(int[][] lutx, int[][] luty, float gz)
{ 
  int i, j;
  int m, n;
  int row, col;
  float singz = (float) Math.sin((double) gz/FPS);
  float cosgz = (float) Math.cos((double) gz/FPS); 
  for (col=0; col<WIDTH; col++) {
    for (row=0; row<HEIGHT; row++) {
      j = col - XCENTRE;
      i = row - YCENTRE;
      n = (int) (cosgz*j - singz*i + XCENTRE);
      m = (int) (singz*j + cosgz*i + YCENTRE);
      if (m < 0) m = 0;
      if (m > (HEIGHT-1)) m = HEIGHT-1;
      if (n < 0) n = 0;
      if (n > (WIDTH-1)) n = WIDTH-1;
      lutx[col][row] = n;
      luty[col][row] = m;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void yawWarpimg(int[][] imgwp1, int[][] imgwp2, int[][] imgor1, int[][] imgor2, int[][] lutx, int[][] luty)
{
  int i, j;
  int row, col;
  for (col=0; col<WIDTH; col++) {
    for (row=0; row<HEIGHT; row++) {     
      j = lutx[col][row];
      i = luty[col][row];     
      imgwp1[col][row] = imgor1[j][i];
      imgwp2[col][row] = imgor2[j][i];
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void gyroWarpPoints(float[][] p1x, float[][] p1y, float[][] p2x, float[][] p2y, float[] camrot, float gz)
{  
  float i, j;
  float m, n;
  int row, col;
  float singz = (float) Math.sin((double) gz/FPS);
  float cosgz = (float) Math.cos((double) gz/FPS); 
  for (col=0; col<XPATCHES; col++) {
    for (row=0; row<YPATCHES; row++) {   
      if (gz==0.0){   
        p2x[col][row] = p1x[col][row] + camrot[1]*gyrtoflowScale;
        p2y[col][row] = p1y[col][row] - camrot[0]*gyrtoflowScale;
      }
      else{
        j = (float) (p1x[col][row] - XCENTRE);
        i = (float) (p1y[col][row] - YCENTRE);
        n =  cosgz*j + singz*i + XCENTRE;
        m = -singz*j + cosgz*i + YCENTRE;
        if (m < 0) m = 0;
        if (m > (HEIGHT-1)) m = HEIGHT-1;
        if (n < 0) n = 0;
        if (n > (WIDTH-1)) n = WIDTH-1;
        p1x[col][row] = n;
        p1y[col][row] = m;
        p2x[col][row] = n + camrot[1]*gyrtoflowScale;
        p2y[col][row] = m - camrot[0]*gyrtoflowScale;       
      }
    }
  }
}

////////////////////////////////////////////////////////////////
void initSystem(float[][] p1x, float[][] p1y, float[] rotGyr2cam)
{
  int k = 0;
  int rowp, colp;  
  //feature locations
  for (colp=0; colp<XPATCHES; colp++) {
    for (rowp=0; rowp<YPATCHES; rowp++) {
      goodMatch.get(k, 0)[0] = 1.0;
      k++;
      p1x[colp][rowp] = XCENTRE - XSPACING*(XPATCHES-1)/2 + colp*XSPACING;
      p1y[colp][rowp] = YCENTRE - YSPACING*(YPATCHES-1)/2 + rowp*YSPACING;
    }
  }   
  XMAX1 = XCENTRE - XSPACING*(XPATCHES-1)/2 - PSIZE*SAMP1/2 - REF1 - 3;    
  YMAX1 = YCENTRE - YSPACING*(YPATCHES-1)/2 - PSIZE*SAMP1/2 - REF1 - 3;  
  XMAX2 = XCENTRE - XSPACING*(XPATCHES-1)/2 - PSIZE*SAMP2/2 - REF2 - 3;    
  YMAX2 = YCENTRE - YSPACING*(YPATCHES-1)/2 - PSIZE*SAMP2/2 - REF2 - 3;
 
 // camera IMU relative orientation
  rotGyr2cam[0] = 0.7851;   rotGyr2cam[1] = 0.1515;   rotGyr2cam[2] = 0.6005;
  rotGyr2cam[3] = 0.6031;   rotGyr2cam[4] = 0.0337;   rotGyr2cam[5] =-0.7970;
  rotGyr2cam[6] =-0.1410;   rotGyr2cam[7] = 0.9879;   rotGyr2cam[8] =-0.0649;
}

//////////////////////////////////////////////////////////////
// calculate the average optical flow
void calAveFlow(float[][] p1x, float[][] p1y, float[][] p2x, float[][] p2y, float[] aveof)
{
  int row, col;
  aveof[0] = 0;
  aveof[1] = 0;  
  aveof[2] = 0;
  float uxAve = 0;
  for (col=0; col<XPATCHES; col++) {
    for (row=0; row<YPATCHES; row++) {
      aveof[0] += (p2x[col][row]-p1x[col][row]);
      aveof[1] += (p2y[col][row]-p1y[col][row]);
      if (col < XPATCHES/2) {        
        uxAve = uxAve - p1x[col][row];
        aveof[2] = aveof[2] - (p2y[col][row]-p1y[col][row]);
      } else {
        uxAve = uxAve + p1x[col][row];
        aveof[2] = aveof[2] + (p2y[col][row]-p1y[col][row]);
      }
    }
  } 
  aveof[0] = aveof[0] / NUM_PATCHES;
  aveof[1] = aveof[1] / NUM_PATCHES;
  aveof[2] = aveof[2] / uxAve * FPS;
}

//////////////////////////////////////////////////////////////
//calculate the median flow 
void calMedFlow(float[][] p1x, float[][] p1y, float[][] p2x, float[][] p2y, float[] medof)
{
  int i = 0;
  int j = 0;
  int k = 0;
  int row, col;
  float medRflow, medLflow;
  float[] flowx = new float[NUM_PATCHES];  
  float[] flowy = new float[NUM_PATCHES];
  float[] flowR = new float[NUM_PATCHES/2];  
  float[] flowL = new float[NUM_PATCHES/2];
  for (col=0; col<XPATCHES; col++) {
    for (row=0; row<YPATCHES; row++) {
      if (col < XPATCHES/2) {
        flowL[j] = (p2y[col][row] - p1y[col][row])/2;
        j++;
      } else {
        flowR[k] = (p2y[col][row] - p1y[col][row])/2;
        k++;
      }
      flowx[i] = p2x[col][row] - p1x[col][row];
      flowy[i] = p2y[col][row] - p1y[col][row];
      i++;
    }
  } 
  Arrays.sort(flowx);
  Arrays.sort(flowy);
  Arrays.sort(flowL);
  Arrays.sort(flowR);
  if (flowx.length % 2 == 0) {
    medof[0] = (flowx[flowx.length/2] + flowx[flowx.length/2 - 1])/2;
    medof[1] = (flowy[flowy.length/2] + flowy[flowy.length/2 - 1])/2;
  } else {
    medof[0] = flowx[flowx.length/2];
    medof[1] = flowy[flowy.length/2];
  }
  if (flowR.length % 2 == 0) {
    medLflow = (flowL[flowL.length/2] + flowL[flowL.length/2 - 1])/2;
    medRflow = (flowR[flowR.length/2] + flowR[flowR.length/2 - 1])/2;
  } else {
    medLflow = flowL[flowL.length/2];
    medRflow = flowR[flowR.length/2];
  }
  medof[2] = medRflow - medLflow;
}

////////////////////////////////////////////////////
void copyPoints(float[][] p1x, float[][] p1y, float[][] p2x, float[][] p2y)
{
  int row, col; 
  for (col=0; col<XPATCHES; col++) {
    for (row=0; row<YPATCHES; row++) {
      p2x[col][row] = p1x[col][row];
      p2y[col][row] = p1y[col][row];
    }
  }
}

////////////////////////////////////////////////////////
void preFilter3(Capture fcap, int[][] preflt31, int[][] preflt15)
{
  int row, col; 
  int index, Isum;
  int[][] imgtmp = new int[WIDTH][HEIGHT];
  frmcap.loadPixels();
  for (row=0; row<HEIGHT; row++) {
    Isum = 0 ;  
    for (col=0; col<3; col++) {
      index = row*WIDTH + col;
      Isum += fcap.pixels[index];
    }
    imgtmp[1][row] = (Isum>>2);////
    for (col=3; col<WIDTH; col++) {
      index = row*WIDTH + col;
      Isum += (fcap.pixels[index] - fcap.pixels[index-3]);
      imgtmp[col-1][row] = (Isum>>2);
    }
  }
  for (col=0; col<WIDTH; col++) {
    Isum = 0;
    for (row=0; row<3; row++) {
      Isum += imgtmp[col][row];
    }
    preflt31[col][1] = (Isum>>2);////
    preflt15[col][1] = (Isum>>2);
    for (row=3; row<HEIGHT; row++) {
      Isum += (imgtmp[col][row] - imgtmp[col][row-3]);
      preflt31[col][row-1] = (Isum>>2);
      preflt15[col][row-1] = (Isum>>2);
    }
  }
}

////////////////////////////////////////////////////////
// average (boxcar) filter
void aveFilter(int[][] preflt, int[][] fltimg, int win)
{
  int i, row, col; 
  int index, Isum;
  int pow2n = 2;
  int minindex = 1;
  int mindiff = 5000;
  int halfwin = (win-1)/2;
  int[][] imgtmp = new int[WIDTH][HEIGHT];

  for (i=1; i<10; i++) {
    if (abs(pow2n-win)<mindiff) {
      minindex = i;
      mindiff = abs(pow2n-win);
    }
    pow2n = pow2n*2;
  }  

  for (row=0; row<HEIGHT; row++) {
    Isum = 0 ;  
    for (col=0; col<win; col++) {  
      Isum += preflt[col][row];
    }
    imgtmp[win-halfwin-1][row] = (Isum>>minindex);
    for (col=win; col<WIDTH; col++) {
      Isum += (preflt[col][row] - preflt[col-win][row]);
      imgtmp[col-halfwin][row] = (Isum>>minindex);
    }
  }
  for (col=0; col<WIDTH; col++) {
    Isum = 0;
    for (row=0; row<win; row++) {
      Isum += imgtmp[col][row];
    }
    fltimg[col][win-halfwin-1] = (Isum>>minindex);////
    for (row=win; row<HEIGHT; row++) {
      Isum += (imgtmp[col][row] - imgtmp[col][row-win]);
      fltimg[col][row-halfwin] = (Isum>>minindex);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// calculate optica flow using the image interpolation algorithm
void I2A_OF(int[][] imgpre, int[][] imgcur, float[][] p1x, float[][] p1y, 
float[][] p2x, float[][] p2y, int rshift, int step) 
{
  int row, col;
  int colp, rowp;
  float scale = 100;
  float dx = 0, dy = 0;
  float f_f0, f2_f1, f4_f3;
  float A, B, C, D, E, denom;
  int prex, prey;
  int rowstart, rowend;
  int colstart, colend;
  int PATCHLEN = step*PSIZE;
  int rowmax1 = HEIGHT-PATCHLEN-rshift-1;
  int colmax1 = WIDTH -PATCHLEN-rshift-1;
  int rowmax2, colmax2 ;

  for (colp=0; colp<XPATCHES; colp++) {
    for (rowp=0; rowp<YPATCHES; rowp++) {  
      prex = (int) (p2x[colp][rowp]-p1x[colp][rowp]);   
      prey = (int) (p2y[colp][rowp]-p1y[colp][rowp]);  
      if (prex>100) prex = 100;
      if (prex<-100) prex = -100;
      if (prey>100) prey = 100;
      if (prey<-100) prey = -100;
      rowmax2 = HEIGHT-PATCHLEN-prey-1;
      colmax2 = WIDTH -PATCHLEN-prex-1;
      rowstart = (int) (p1y[colp][rowp]-PATCHLEN/2);
      if (rowstart<rshift) rowstart = rshift;
      if (rowstart<-prey) rowstart = -prey;
      if (rowstart>rowmax1) rowstart = rowmax1;      
      if (rowstart>rowmax2) rowstart = rowmax2;
      rowend = rowstart + PATCHLEN;
      colstart = (int) (p1x[colp][rowp]-PATCHLEN/2);
      if (colstart<rshift) colstart = rshift;
      if (colstart<-prex) colstart = -prex;
      if (colstart>colmax1) colstart = colmax1;       
      if (colstart>colmax2) colstart = colmax2;
      colend = colstart + PATCHLEN;
      A = 0;  
      B = 0; 
      C = 0; 
      D = 0;  
      E = 0;
      for (row=rowstart; row<=rowend; row+=step) {
        for (col=colstart; col<=colend; col+=step) {                      
          f2_f1 = imgpre[col-rshift][row] - imgpre[col+rshift][row];            
          f4_f3 = imgpre[col][row-rshift] - imgpre[col][row+rshift];
          f_f0  = imgcur[col+prex][row+prey] - imgpre[col][row];
          A += f2_f1 * f2_f1 / scale;
          B += f4_f3 * f2_f1 / scale;
          C += f_f0  * f2_f1 * 2 / scale;
          D += f4_f3 * f4_f3 / scale;
          E += f_f0  * f4_f3 * 2 / scale;
        }
      }
      denom = (A*D) - (B*B);
      if (denom != 0) {
        dx = ((C*D)-(B*E)) / denom;
        dy = ((A*E)-(B*C)) / denom;
      }      
      p2x[colp][rowp] = p1x[colp][rowp] + dx * rshift + prex; 
      p2y[colp][rowp] = p1y[colp][rowp] + dy * rshift + prey;
    }
  }
}

////////////////////////////////////////////////////////////
void setupRazor() {
  println("Trying to setup and synch Razor...");
  // On Mac OSX and Linux (Windows too?) the board will do a reset when we connect, which is really bad.
  // See "Automatic (Software) Reset" on http://www.arduino.cc/en/Main/ArduinoBoardProMini
  // So we have to wait until the bootloader is finished and the Razor firmware can receive commands.
  // To prevent this, disconnect/cut/unplug the DTR line going to the board. This also has the advantage,
  // that the angles you receive are stable right from the beginning. 
  delay(3000);  // 3 seconds should be enough
  // Set Razor output parameters
  serial.write("#osrb");  // Turn on binary output  ///////ob
  serial.write("#o1");  // Turn on continuous streaming output
  serial.write("#oe0"); // Disable error message output
  // Synch with Razor
  serial.clear();  // Clear input buffer up to here
  serial.write("#s00");  // Request synch token // #s00
}

float readFloat(Serial s) {
  // Convert from little endian (Razor) to big endian (Java) and interpret as float
  return Float.intBitsToFloat(s.read() + (s.read() << 8) + (s.read() << 16) + (s.read() << 24));
}

// Skip incoming serial stream data until token is found
boolean readToken(Serial serial, String token) {
  // Wait until enough bytes are available
  if (serial.available() < token.length())
    return false;
  // Check if incoming bytes match token
  for (int i = 0; i < token.length (); i++) {
    if (serial.read() != token.charAt(i))
      return false;
  }
  return true;
}

////////////////////////////////////////
// for a second window display
public class PFrame extends Frame {
  public PFrame() {
    setBounds(100, 100, 335, 278);
    sndWin = new secondApplet();
    add(sndWin);
    sndWin.init();
    show();
  }
}
public class secondApplet extends PApplet {
  public void setup() {
    size(320, 240);
    noLoop();
  }
  public void draw() {
  }
}

