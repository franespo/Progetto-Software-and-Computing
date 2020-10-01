from configparser import ConfigParser

config = ConfigParser()

config['General_Variables'] = {
    'starting_wav' : '0.3',
    'final_wav' : '1.01',
    'zenith_angle' : '40',
    'scattering_angle' : '130',
    'mixing_ratio' : '0',
    'wavelength_step' : '0.01',
    'angstrom_coefficient_alpha' : '1',
    'Henyey_Greenstein_g_parameter' : '0.85',
    'BeightnessTemperature_of_the_sun' : '5800',
    'Vertical_scale_Height' : '7.7',
    'ith_wav' : '0.52',
    }

config['Output'] = {
    'output_figures' : './Figures/'}

with open('./Single_scatterConfiguration.ini','w') as file:
    config.write(file)

