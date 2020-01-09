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
