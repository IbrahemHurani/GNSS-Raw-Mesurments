# GNSS-Raw-Mesurments
The First assignment in Course Autonomous Robotics
## Source
https://www.johnsonmitchelld.com/2021/03/14/least-squares-gps.html


## How To Run:
* Download the files from github
* open the code (for example pyCharm)
* open the file integration.py
* in the main you need to do this before run the code
* first put the path file for what file you want(file fixed or walking or driving) in line 335 like this:
```python
ecef_List=process('C:\\Users\\user\Desktop\gnss-analysis-main\data\sample\walking.txt')
```
* Run the Code.

* the result for Q2:
GPS time, SatPRN (ID), Sat.X, Sat.Y, Sat.Z, Pseudo-Range, CN0, Doppler 
 you can found in file with name Measurements_Output.csv
* the result for Q3&Q4:
you can found in file android_position.csv and calculated_postion.csv
* the result for Q4:
you can found in file output_kml.kml




## Plot result:

# 1)File Fixed:

![Fixed_plot](https://github.com/IbrahemHurani/GNSS-Raw-Mesurments/assets/86603326/bf476885-72c9-42bb-830e-357179e702b5)

# 2)File Walking:

![walking_plot](https://github.com/IbrahemHurani/GNSS-Raw-Mesurments/assets/86603326/c4b02a4d-a143-456c-a85b-93f0d3287c7a)

# 3)File Driving:
![Driving_plot](https://github.com/IbrahemHurani/GNSS-Raw-Mesurments/assets/86603326/34115547-dcbd-4b72-9249-58fa7fb9c730)



