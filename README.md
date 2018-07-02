# Unscented Kalman Filter Project

---

## Overview

In this project I'he utilize an Unscented Kalman filter to estimate the state of a moving object of interest with noisy lidar and radar measurements.

## Compile

This project require [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) to be installed. This repository includes two files that can be used to set up and install it for either Linux or Mac systems.

Once the install for uWebSocketIO is complete, the main program can be built by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make

Dependencies:

* cmake >= 3.5
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
* gcc/g++ >= 5.4


## Run

This project needs Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

To run the UKF algorith estimation use

```

./build/UnscentedKF
```

command in project directory.


Here is the protocol that UKF estimator uses for uWebSocketIO in communicating with the simulator.


INPUT: values provided by the simulator to the c++ program:

* ["sensor_measurement"] - the measurement that the simulator observed (lidar and radar)


OUTPUT: values provided by the c++ program to the simulator:

* ["estimate_x"] - UKF filter estimated position x
* ["estimate_y"] - UKF filter estimated position y
* ["rmse_x"] - root mean square error between real and estimated position x
* ["rmse_y"] - root mean square error between real and estimated position y
* ["rmse_vx"] - root mean square error between real and estimated velocity x
* ["rmse_vy"] - root mean square error between real and estimated velocity y
