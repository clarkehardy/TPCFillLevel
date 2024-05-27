import yaml
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba


def rectangular_profile(properties):
    '''
    Returns a function to calculate the cumulative volume to a given
    height for an object or set of objects with a rectangular profile
    when viewed from the side (such as a vertical cylinder).
    '''
    if properties['profile'] != 'rectangular':
        return lambda y: 0
    
    height = properties['height']
    volume = properties['volume']
    cross_sectional_area = volume/height
    y_bottom = properties['y_position']
    number = properties['number']
    spacing = properties['spacing']
    sign = -1 if properties['solid'] else 1

    def func(y):
        vol = 0
        for i in range(number):
            y_low = y_bottom + i*spacing
            if y < y_low:
                pass
            elif y > y_low + height:
                vol += sign*volume
            else:
                vol += sign*cross_sectional_area*(y-y_low)
        return vol
    return func


def circular_profile(properties):
    '''
    Returns a function to calculate the cumulative volume to a given
    height for an object or set of objects with a circular profile
    when viewed from the side (such as a horizontal cylinder).
    '''
    if properties['profile'] != 'circular':
        return lambda y: 0
    
    diameter = properties['height']
    volume = properties['volume']
    length = 4*volume/(np.pi*diameter**2)
    y_bottom = properties['y_position']
    sign = -1 if properties['solid'] else 1

    def func(y):
        if y < y_bottom:
            vol = 0
        elif y > y_bottom + diameter:
            vol = sign*volume
        else:
            y_subtracted = y - y_bottom
            if y_subtracted < diameter/2:
                vol = sign*length*circular_segment_area(y_subtracted,diameter)
            else:
                vol = sign*length*(np.pi*diameter**2/4 - circular_segment_area(diameter-y_subtracted,diameter))
        return vol
    return func


def circular_segment_area(y,diameter):
    '''
    Returns the area of a circular segment given y, the perpendicular distance
    from a tangent in the direction of the center of the circle, and the diameter.
    '''
    return diameter**2*np.arccos((diameter - 2*y)/diameter)/4 - (diameter/2 - y)*np.sqrt(y*(diameter - y))


class FillLevel:

    def __init__(self,detector_yaml):
        '''
        Initialize the TPCObject with a yaml file containing the detector geometry.
        '''
        with open(detector_yaml,'r') as infile:
            self.geometry = yaml.safe_load(infile)

        # largest uncertainty (apart from approximations in geometry) is in the density
        # provide a lower, median, and upper value to get error bars if desired
        self.lxe_density = [2.89,2.95,2.97] # lower, median, upper in g/cm^3
        self.build_fill_funcs()
    

    def build_fill_funcs(self):
        '''
        Build functions that can be used to output the total mass filled
        given a height y and the reverse.
        '''
        volume_funcs = []
        for component in self.geometry:
            volume_funcs.append(rectangular_profile(self.geometry[component]))
            volume_funcs.append(circular_profile(self.geometry[component]))

        def mass_filled(y,return_errors=False):
            # if a single value is passed, return a single value
            return_num = False
            if np.shape(y) == ():
                return_num = True
                y = np.array([y])

            # integrate up to the given height to get the volume
            volume = np.zeros_like(y)
            for i in range(len(y)):
                volume[i] = np.sum([func(y[i]) for func in volume_funcs])

            if return_num:
                volume = volume[0]

            # return either the median value or the lower, median, and upper values
            if return_errors:
                return volume*self.lxe_density[0],volume*self.lxe_density[1],volume*self.lxe_density[2]
            else:
                return volume*self.lxe_density[1]
        
        self.mass_filled = mass_filled

        def fill_level(mass,return_errors=False):
            # if a single value is passed, return a single value
            return_num = False
            if np.shape(mass) == ():
                return_num = True
                mass = np.array([mass])

            # get the maximum height of the detector
            heights = []
            for component in self.geometry:
                comp = self.geometry[component]
                heights.append(comp['y_position']+comp['height']+(comp['number']-1)+comp['spacing'])

            # create an array of heights to interpolate the mass filled
            fill_pts = np.linspace(0,max(heights),100000)
            mass_pts = self.mass_filled(fill_pts,return_errors)

            if return_errors:
                lower_mass_pts,mass_pts,upper_mass_pts = mass_pts
                lower_levels = np.interp(mass,lower_mass_pts,fill_pts)
                upper_levels = np.interp(mass,upper_mass_pts,fill_pts)
                lower_levels[lower_levels < 0] = 0
                upper_levels[upper_levels < 0] = 0

            # interpolate the mass filled to get the fill level
            levels = np.interp(mass,mass_pts,fill_pts)
            levels[levels < 0] = 0

            if return_num:
                levels = levels[0]

            # return either the median value or the lower, median, and upper values
            if return_errors:
                return lower_levels,levels,upper_levels
            else:
                return levels
        
        self.fill_level = fill_level


    def plot_fill_level(self,masses,plot_errors=False):
        '''
        Plot the fill level in the detector as a function of the mass filled.
        '''

        fill_levels = self.fill_level(masses,return_errors=True)

        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        fig,ax = plt.subplots()
        ax.plot(masses,fill_levels[1],color=colors[0],label='Fill level')
        if plot_errors:
            ax.fill_between(masses,fill_levels[0],fill_levels[2],color=colors[0],lw=0,alpha=0.3)
        ax.set_xlabel('Mass filled [g]')
        ax.set_ylabel('Fill level [cm]')
        ax.set_title('Liquid xenon level in the detector')
        ax.set_xlim([np.amin(masses),np.amax(masses)])
        ax.set_ylim([np.amin(fill_levels),np.amax(fill_levels)])
        colors = colors[1:]

        for j,component in enumerate(self.geometry):
            if self.geometry[component]['y_position'] < 0:
                continue
            for i in range(self.geometry[component]['number']):
                y_lower = self.geometry[component]['y_position']+i*self.geometry[component]['spacing']
                y_upper = self.geometry[component]['y_position']+i*self.geometry[component]['spacing']\
                          + self.geometry[component]['height']
                mass_range = [self.mass_filled(y_lower),self.mass_filled(y_upper)]
                plt.fill_between(mass_range,y_lower,y_upper,edgecolor=colors[j%len(colors)],\
                                 facecolor=to_rgba(colors[j%len(colors)],alpha=0.1),\
                                 label=component if i == 0 else None)
                
        ax.legend(ncol=4,fontsize=8)

        return fig,ax


    def plot_mass_filled(self,fill_levels,plot_errors=False):
        '''
        Plot the mass filled as a function of the fill level in the detector.
        '''

        masses = self.mass_filled(fill_levels,return_errors=True)

        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        fig,ax = plt.subplots()
        ax.plot(fill_levels,masses[1],color=colors[0],label='Fill level')
        if plot_errors:
            ax.fill_between(fill_levels,masses[0],masses[2],color=colors[0],lw=0,alpha=0.3)
        ax.set_ylabel('Mass filled [g]')
        ax.set_xlabel('Fill level [cm]')
        ax.set_title('Xenon mass in the detector')
        ax.set_xlim([np.amin(fill_levels),np.amax(fill_levels)])
        ax.set_ylim([np.amin(masses),np.amax(masses)])
        colors = colors[1:]

        for j,component in enumerate(self.geometry):
            if self.geometry[component]['y_position'] < 0:
                continue
            for i in range(self.geometry[component]['number']):
                y_lower = self.geometry[component]['y_position']+i*self.geometry[component]['spacing']
                y_upper = self.geometry[component]['y_position']+i*self.geometry[component]['spacing']\
                          + self.geometry[component]['height']
                mass_range = [self.mass_filled(y_lower),self.mass_filled(y_upper)]
                plt.fill_betweenx(mass_range,y_lower,y_upper,edgecolor=colors[j%len(colors)],\
                                  facecolor=to_rgba(colors[j%len(colors)],alpha=0.1),\
                                  label=component if i == 0 else None)
                
        ax.legend(ncol=4,fontsize=8)

        return fig,ax
    

    def plot_geometry(self):
        '''
        Visualize the features of the detector geometry.
        '''
        
        fig,ax = plt.subplots()

        # get the top of the highest component
        heights = []
        for component in self.geometry:
            comp = self.geometry[component]
            heights.append(comp['y_position']+comp['height']+(comp['number']-1)+comp['spacing'])

        y_vals = np.linspace(0,np.amax(heights),1000)
        areas = np.diff(self.mass_filled(y_vals)/self.lxe_density[1])/(y_vals[1]-y_vals[0])

        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        ax.plot(areas,y_vals[:-1],color=colors[0],label='Area')
        ax.set_xlabel('Cross-sectional area [cm$^2$]')
        ax.set_ylabel('Height [cm]')
        ax.set_title('Detector geometry')
        ax.set_ylim([0,np.amax(y_vals)])
        lims = ax.get_xlim()
        colors = colors[1:]

        for j,component in enumerate(self.geometry):
            if self.geometry[component]['y_position'] < 0:
                continue
            for i in range(self.geometry[component]['number']):
                y_lower = self.geometry[component]['y_position']+i*self.geometry[component]['spacing']
                y_upper = self.geometry[component]['y_position']+i*self.geometry[component]['spacing']\
                          + self.geometry[component]['height']
                plt.fill_between(lims,y_lower,y_upper,edgecolor=colors[j%len(colors)],\
                                 facecolor=to_rgba(colors[j%len(colors)],alpha=0.1),\
                                 label=component if i == 0 else None)

        ax.legend(ncol=4,fontsize=8)

        return fig,ax