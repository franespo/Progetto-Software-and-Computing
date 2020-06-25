# Scattering Process
The interaction of electromagnetic radiation with matter modifies to some extent the incident wave. The medium therefore produces a signature in the amplitude, phase or spectral composition which depends on composition and structure of the medium. In the atmosphere atoms and molecules are exposed to sun light and they absorb and re-emit it (according to Kirchhoff's Law) at any wavelength (or wavenumber) in the infrared or microwave regions in different directions with different intensity. This is a brief description of a scattering process. When radiation is only scattered by one localizated scattering center, this is called single scattering. Different regimes are identified by defining the Mie parameter:

![Mie Parameter](https://latex.codecogs.com/gif.latex?x%20%3D%20%5Cfrac%7B2%20%5Cpi%20a%7D%7B%20%5Clambda%20%7D)


which indicates the relationship between the dimensions of the scatterer (assumed spherical radius a) and the wavelength incident wave. If the scatterer is excessively small, it does not occur interactions with incident radiation. 

![Scattering regimes](https://images.slideplayer.com/26/8579527/slides/slide_2.jpg) 

Generally, in any volume of air, scatterers with different properties (dimensions, type of material) are found. The transfer of radiation through a volume is very dependent on the details of the type, dimension and number of these particles. Generally multiple scattering becomes more important as the number of particles increase. To simulate quantitatively the radiance in a dense media we need to know the properties of the scatterers considered as independent, that is the single scattering properties.
The project code is implemented to investigate the single scattering solutions and the "colours" of the light reaching the observer from an horizontal path. The atmosphere may contain molecules and aerosol.

# Structure of the code
1. The code has been implemented to give the opportunity to users to execute different cases which depends on different input values as sun zenith angle, scattering angle, mixing ratio (between molecules and aerosols), Angstrom coefficient and wavelength. Some hints have been added to the first few lines of code but users are free to enter the values they want as long as they are consistent with the problem.  
2. The output of the code are three different figures saved in the folder [figures](https://github.com/franespo/Progetto-Software-and-Computing/tree/master/Figures). Each figure is composed by two plots and shows the results of implementation. In the first one (fig1.png) are showed how molecular and aerosols scattering coefficients changes with the wavelength and in the second one how the vertical and horizontal optical path change with wavelength; in the second figure (fig2.png) are showed in the upper plot how the transmittance change in a vertical path and in an horizontal path with wavelength, while in the plot below are showed molecular and aerosol scattering direction and how total scattering direction change at three different wavelengths (0.3, 0.64 and 1 microns); in the last figure (fig3.png) are showed sun radiance E0 at TOA, radiance reaching surface layer E0S and radiance reaching the observer L in the upper plot, while below are showed same curves but normalized.

## The project is divided in different files:
* file [single scatter code](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/SingleScatter.py) contains the code of single scatter implementation. It is essential to ensure that the code turns correctly insert the correct inputs. 
* file [functions](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/functions.py) contains the main functions used in the code.
* file [testing](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/testing_functions.py) contains different tests for file functions to ensure that all of them work properly, using hypothesis testing.
* file [figures](https://github.com/franespo/Progetto-Software-and-Computing/tree/master/Figures) contains the figures output of the code.

An example of figures output is reported below:
![fig1](https://github.com/franespo/Progetto-Software-and-Computing/blob/master/Figures/fig1.png)



