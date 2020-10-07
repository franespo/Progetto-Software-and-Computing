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
    'BeightnessTemperature_of_the_sun' : '5800',
    'Vertical_scale_Height' : '7.7',
    'wavelength_norm' : '0.52',
    'distance_from_black_wall' :'25'
    }

config['Output'] = {
    'output_figures' : './FiguresNew/'}

with open('./ScatterProfile_Configuration.ini','w') as file:
    config.write(file)

