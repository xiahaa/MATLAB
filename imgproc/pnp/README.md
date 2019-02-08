In this folder, I tried to learn the famous Perspective N Point algorithm which aims to estimate transformation information from a batch (>=3) 3D-2D correspondence points.
For some algorithms, I re-implemented them, while for others, I just write a warpper to call them.

## Already Imeplemented:
1. P3P Gao;
2. P3P Kneip;
3. P3P Pst;
4. PnP AK;
5. PnP QuanLong;
6. LHM;
7. EPnP (intermediate result differs with the EPnP, but final result approaches with numerical error);
8. POSIT;
9. softPOSIT (not very stable, in some case fail to converge);

## use only Warpper to call 3rdparty implementation:
### MLEstimator variant:
1. MLPnP: MLPnP may not be a good choice in my opinion: 1) extensive null vector computation for each 2D points; 2) ordinary case needs at least 6 points; 3) omit the constraint imposed by R; 4) use all 3D points for GN optimization.
2. CEPPnP: basically EEPnP + Prior covriance = setup likelihood function, and then use a Gauss-Newton-liked routine to iteratively optimize.
### 2nd class: seems a little bit tricky to select the basis and then do the elimination, so not so easy to reproduce.....
1. OPnP: I do implement a OPnP trial in order to see if I understand how to use the Automatic Minimal Solver correctly. It works well as in the OPnP_trial.m, however, it is pretty slow compared with the original OPnP solution. 
2. AsPnP: AsPnP seems to be predecessor of OPnP since they share a similar cost function and general framework. However, AsPnP divides into four different cases and extract answer from the 4 cases, then return them ranked by algebraic error. OPnP summaries those 4 cases as a uniform solver.
3. DLS
4. UPnP: 
### 3nd class: 
1. RPnP: .

## PnP tightly coupled with RANSAC:
1. REEPnP: EPnP variance with RANSAC coupled. Normal treatement of adding RANSAC is after the estimation. Here, RANSAC is placed just after the computation of the virtual 4 control points. So a speed-up can be achieved.


## Useful links:
PnP toolbox: 
1. https://github.com/haoyinzhou/PnP_Toolbox
2. https://github.com/urbste/MLPnP_matlab_toolbox

## Not yet read:
5. R1PnP;https://github.com/haoyinzhou/PnP_Toolbox
