## continuous estimation of homography matrix using SL3 filter

1. with known group velocity (trivial);

### use homography matrix as measurement

> Dynamic estimation of homography transformations on the special linear group for visual servo control

2. unknown group velocity: but assume it is constant (done);

### use point correspondence as measurement

> Hamel, Tarek, et al. "Homography estimation on the special linear group based on direct point correspondence." 2011 50th IEEE Conference on Decision and Control and European Control Conference. IEEE, 2011.

3. unknown group velocity: but assume it is constant (done);

4. unknown group velocity: assume the group undergoes rigid motion (done); choose the second form where the V/d is a constant (cyclical motion)

### use line correspondece as measurement
> Hua, Minh-Duc, et al. "Point and line feature-based observer design on SL (3) for Homography estimation and its application to image stabilization." (2017).

5. unknown group velocity: assume the group undergoes rigid motion (done);

**Note**: with 4 and 5, it is trivial to apply the combined point and line correspondece method.

**Note**: from the demo video made by the author, it can be inferred that the filter prediction will not be very accurate if the motion assumption doesn't hold (which is quite common). My personal opinion is that it can be used in a short time for prediction.

### extract optical flow from continuous homography
> Manerikar N, Hua M D, Hamel T. "Homography Observer Design on Special Linear Group SL (3) with Application to Optical Flow Estimation"[C], 2018 European Control Conference (ECC). IEEE, 2018: 1-5.

利用有限差分，log(H)获取群速度U，进而根据rigid motion下的U构成，去掉旋转向量之后，再提取optical flow。

**此文可以反过来用，利用相机获取连续H，差分得到U，motion时假设V极小或者d极大，则U和Omega对应，反过来可以用来做对齐和标定。**


### EKF version


### coupled with Optical flow
survey:
https://scholar.google.com.sg/citations?user=qpu09zcAAAAJ&hl=en
