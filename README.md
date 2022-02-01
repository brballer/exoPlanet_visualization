# exoPlanet_visualization
A Python based tool for visualizing transit curves of exoplanets and eclipsing binary star systems. The major improvements of this version are 

* A star system animation
* Option to include a circumsecondary dust disk
* Speed improvements
* Option to bulk-load model inputs

### Intended audience
Citizen scientists who would like to understand what causes features of transit curves, for example those involved with PlanetHunters NGTS.

### Requirements
Users need to have Python3, matplotlib and numpy installed.

## Installation
1. Press "Code" (which is next to "Go to File and "Add to File"). Then choose "Download ZIP".
2. Open Terminal (Mac software pre-loaded). You can directly search for this in Finder. All the following steps are in the Terminal.
3. Type "/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)" and wait for completion.
4. Type "export PATH="/usr/local/opt/python/libexec/bin:$PATH" or if you have OS X 10.12 (Sierra) or older, type "export PATH=/usr/local/bin:/usr/local/sbin:$PATH". Wait for completion.
5. Type "brew install python3" and wait for download to complete.
6. Type "pip install numpy" and press enter. Wait for download to complete.
7. Type "pip install matplotlib" and press enter. Wait for download to complete.
8. Type "cd downloads" and press enter.
9. Type "cd exoPlanet_visualization-main" and press enter.

### Using the tool

The GUI is launched using the command: "python3 exoTransit_visualization.py". There are five planetary system model parameters that are user-defined within the GUI. The model allows specifying a limb-darkening factor that is applied to the primary star and a binary companion star. 

Detailed instructions can be found in ExoTransit_Visualization_Guide.pdf. The guide  includes screen shots of run results of planet transits and binary star eclipses, and technical details about how the transit curves are modeled. 
