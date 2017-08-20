# Unscented Kalman Filter

---

The aim of this project is to implement an unscented Kalman filter using the Constant Turn Rate and Velocity (CTRV) motion model. We will be using the same bicycle simulation data set from the Extended Kalman filter [project](https://github.com/coldKnight/ExtendedKalmanFilter). That way we can compare our results with the EKF project.

Remember that all Kalman filters have the same three steps:

    Initialization
    Prediction
    Update

A standard Kalman filter can only handle linear equations. Both the extended Kalman filter and the unscented Kalman filter allow us to use non-linear equations; the difference between EKF and UKF is how they handle non-linear equations. But the basics are the same: initialize, predict, update.

Some additional resources: 
- https://balzer82.github.io/Kalman/
- www.mdpi.com/1424-8220/14/3/5239/pdf

## Results

### RMSE Values

* sample data 1 - **[0.0724284, 0.0781081, 0.57687, 0.568334]**
* sample data 2 - **[0.161822, 0.17826, 0.222124, 0.290671]**

### Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

### Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`

### Code Style

Please stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html) as much as possible.

### Generating Additional Data

If you'd like to generate your own radar and lidar data, see the
[utilities repo](https://github.com/udacity/CarND-Mercedes-SF-Utilities) for
Matlab scripts that can generate additional data.
