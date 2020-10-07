# Scattering Process
The interaction of electromagnetic radiation with matter modifies to some extent the incident wave. The medium therefore produces a signature in the amplitude, phase or spectral composition which depends on composition and structure of the medium. In the atmosphere atoms and molecules are exposed to sun light and they absorb and re-emit it (according to Kirchhoff's Law) at any wavelength (or wavenumber) in the infrared or microwave regions in different directions with different intensity. This is a brief description of a scattering process. When radiation is only scattered by one localizated scattering center, this is called single scattering. Different regimes are identified by defining the Mie parameter:

![Mie Parameter](https://latex.codecogs.com/gif.latex?x%20%3D%20%5Cfrac%7B2%20%5Cpi%20a%7D%7B%20%5Clambda%20%7D)


which indicates the relationship between the dimensions of the scatterer (assumed spherical radius a) and the wavelength incident wave. If the scatterer is excessively small, it does not occur interactions with incident radiation. 

![Scattering regimes](https://images.slideplayer.com/26/8579527/slides/slide_2.jpg) 

Generally, in any volume of air, scatterers with different properties (dimensions, type of material) are found. The transfer of radiation through a volume is very dependent on the details of the type, dimension and number of these particles. Generally multiple scattering becomes more important as the number of particles increase. To simulate quantitatively the radiance in a dense media we need to know the properties of the scatterers considered as independent, that is the single scattering properties.
The project code is implemented to investigate the single scattering solutions (not considered multiple scattering) and the "colours" of the light reaching the observer from an horizontal path. In particular, the light reflected by a black object placed at a certain distance from the observer is analyzed: in the absence of molecules / aerosols between the object and the observer, the latter would observe the object completely black, while the situation changes when there are scatterers in the air. 
The atmosphere may contain molecules and aerosol.

# Structure of the code
1. The code has been implemented to give the opportunity to users to execute different cases which depends on different input values. The input values can be modified in the file [Configuration_maker.py](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/Configuration_maker.py) and the variables are shown below:
* **starting_wav** represents the starting value of wavelength vector (supposed to be measured in microns), essential for the correct functioning of the code. The value must be positive and greater then zero;
* **final_wav** represents the final value of wavelength vector (in microns). The value must be positive, greater then zero and also the starting value;
* **zenith_angle** *(0-90°)* is the angle between the zenith and the centre of the Sun's disc. It's important for the calculation of the irradiance;
* **scattering_angle** *(0-180°)* is the angle in the scattering plane (defined in turn as the plan formed by the direction of the incident radiation and that of scattered radiation) beetween the direction of propagation of incident e.m. field and the direction of the scattered wave;
*  **mixing_ratio** is the ratio between molecules and aerosols;
*  **wavelength_step** is the step with which the wavelength vector is constructed;
* **BrightnessTemperature_of_the_sun** represents the temperature of a black body, in this case the sun. Its value must be greater then zero.
* **Vertical_scale_Height** is vertical distance over which the density and pressure fall by a factor of 1/e and describe the atmosphere. Its value is expressed in kilometers and is fundamental for calculation of optical paths;
* **wavelength_norm** represents a specific value in the range of wavelength vector expressed in microns;
* **distance_from_black_wall** is the distance from observer used for scattering analysis.

2. The output of the code is summerized in a folder named [FiguresNew](https://github.com/franespo/Progetto-Software-and-Computing/tree/master/FiguresNew). Each figure is named according to its content (and has the title of the graphic inside it) and shows the results of implementation. [Irradiances.png](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/FiguresNew/Irradiances.png) and [Normalized Irradiances.png](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/FiguresNew/Normalized%20Irradiances.png) show how the irradiance change with the wavelength and respect to atmosphere level, represented by TOA (top of the atmosphere, E0) and by SL (surface level, E0S).  
[Optical Paths.png](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/FiguresNew/Optical%20Paths.png) shows how vertical optical path and horizontal optical path change with wavelength and scattering coefficient. The vertical path depends also on vertical scale height, whereas the horizontal path depends also on distance of an observer placed, in the current case, 25 meters away from a completely black object.
[Scattering coefficient trend.png](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/FiguresNew/Scattering%20coefficient%20trend.png) shows how molecular and aerosols scattering coefficients changes with the wavelength. Clearly, if in the atmosphere we have no aerosols, the scattering coefficient is zero.
[Scattering diagram.png](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/FiguresNew/Scattering%20diagram.png) shows molecular and aerosol scattering direction and how total scattering direction change at three different wavelengths (0.3, 0.64 and 1 microns).
Lastly, [Transmittance.png](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/FiguresNew/Transmittance.png) shows how the transmittance change in a vertical path and in an horizontal path with wavelength.

## The project is divided in different files:
* [FiguresNew](https://github.com/franespo/Progetto-Software-and-Computing/tree/master/FiguresNew) contains the figures output of the code;
* [Configuration_maker.py](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/Configuration_maker.py) contains the parameters used in the code and can be modified before running the code.
* [README.md](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/README.md) contains the instructions on how use the code;
* [ScatterProfile_Configuration.ini] contains the output of Configuration_maker.py and this file change respect to the modifies applied;
* [Scatter_Profile.py](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/Scatter_Profile.py) is the core of the code;
* file [functions](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/functions.py) contains the main functions used in the code;
* file [testing_functions.py](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/testing_functions.py) contains different tests for file functions to ensure that all of them work properly.

## How to run the code
In order to run properly the code, this steps can be followed:
1. The folder containing the code can be downloaded to this link: https://github.com/franespo/Progetto-Software-and-Computing.git 
2. The user can use the code running the file [Scatter_Profile.py](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/Scatter_Profile.py) with python via terminal;
3. If the user wants, can modify the input values in the [Configuration_maker.py](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/Configuration_maker.py) file and launch it before execution of step 2.

An example of figure output is reported below:


![Transmittance](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/FiguresNew/Transmittance.png)



