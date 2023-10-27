from montu import *
###############################################################
# Module constants
###############################################################
STELLAR_CATALOGUE = 'montu_stellar_catalogue_v37.csv'

"""References
    Styles: https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
    #sphx-glr-gallery-style-sheets-style-sheets-reference-py
"""
PLT_DEFAULT_STYLE = 'default' # others: ggplot, default, classic
SET_PLT_DEFAULT_STYLE = lambda:plt.style.use(PLT_DEFAULT_STYLE)

###############################################################
# Stars Class
###############################################################
class Stars(Sebau):
    """Stellar catalogue

    Initialization parameters:

        data: pandas.Dataframe, default = None:
            Pandas dataframe containing the stars.

        filename: string, default = None:
            File containing the database with stars.
            If None it uses the official MontuCatalogue, namely 
            montu_stellar_catalogue_x.csv (where x is the version of the catalogue)

    Attributes:
        
        data: pandas.DataFrame:
            Data containing the information on stars.

        number: int:
            Number of stars in data.

    Methods:
        get_stars:
            Get a subset of the stars in catalogue.

        get_stars_around:
            Get stars around a given position in the sky.

            
    """

    def __init__(self,data=None,filename=None):

        if data is not None:
            # Load data for stars from a dataframe already loaded
            self.data = copy.deepcopy(data)
            
        elif filename:
            # Load data from a file
            self.data = pd.read_csv(filename)

        else:
            # Load data from the database provided with package
            self.data = pd.read_csv(Util._data_path(STELLAR_CATALOGUE,check=True))

        self.number = len(self.data)

    def get_stars(self,**args):
        """Filter stars by criteria

        Examples:
            # Get a single stars
            aldebaran = allstars.get_stars(ProperName='Aldebaran')

            # All visible stars in the sky
            visible = allstars.get_stars(Mag=[-2,6.5])

            # All visible stars with declination less than 1 deg in absolute value
            equator = allstars.get_stars(Mag=[-3,6.5],Dec=[-1,1])
        """

        # If no args get all stars in data base
        if len(args)==0:
            return self
        
        # If args provided it will try to filter database according to conditions
        cond = np.array([True]*len(self.data))
        for key,item in args.items():
            if key == 'suffix':continue
            if isinstance(item,list):
                min = float(item[0])
                max = float(item[1])
                cond = (self.data[key]>=min)&(self.data[key]<=max)&(cond)
            elif isinstance(item,tuple):
                cond_or = np.array([False]*len(self.data))
                for it in item:
                    cond_or = (self.data[key]==it)|cond_or
                cond = (cond_or)&(cond)
            else:
                cond = (self.data[key]==item)&(cond)
    
        return Stars(self.data[cond])
    
    def get_stars_around(self,
                         center=[0,0],radius=10,
                         coords=['RAJ2000','DecJ2000'],**kwargs):
        """Get stars around a point in the sky

        Parameters:
            center: list, default = [0,0] 
                Center around which the region will be extracted.

            radius: float, default = 10 [same units as center]:
                Radius of the region

            coords: list of strings, default = ['RAJ2000','DecJ2000'];
                Name of coordinates on which the center is calculated.

        Return: 
            stars: Stars:
                Stars in area.

        Examples:
            # Get all stars around aldebaran in a radius of 5 degrees and with magnitudes between -1 and 4
            hyades = stars.get_stars_around(center=[aldebaran.data.RAJ2000,aldebaran.data.DecJ2000],radius=15,Vmag=[-1,4])
        """
        kwargs.update({
            coords[0]:[float(center[0]-radius/15),float(center[0]+radius/15)],
            coords[1]:[float(center[1]-radius),float(center[1]+radius)],
        })
        stars = self.get_stars(**kwargs)
        return stars
    
    def plot_stars(self,coords=['RAJ2000','DecJ2000'],
                   labels=True,pad=0,figargs=dict(),stargs=dict()):
        """Plot all stars in data.

        Parameters:
            coords: list of strings, default = ['RAJ2000','DecJ2000']
                Coordinates used in representation.

            labels: Boolean, default = True:
                Do you want to see stellar labels.

            figargs: dictionary:
                Additional options for the figure.
            
            starargs: dictionary:
                Additional options for the scatter.

        Return:
            fig,axs: 
                Figure and axes.
        """
        # Black background
        plt.style.use('dark_background')

        # Create figure
        dfigargs = dict(figsize=(8,8))
        dfigargs.update(figargs)
        fig,axs = plt.subplots(1,1,**dfigargs)

        # Axis
        axs.set_facecolor('black')

        # Scatter
        dstargs = dict(marker='*',color='y')
        dstargs.update(stargs)

        size_by_mag = Util._linear_map([-1.5,5],[200,1])
        axs.scatter(15*self.data[coords[0]],
                    self.data[coords[1]],
                    s=size_by_mag(self.data.Vmag),
                    **dstargs)
        
        # Labels
        fontsize = Util._linear_map([6,-2],[4,14])
        if labels:
            for index in self.data.index:
                star = self.data.loc[index]
                star.fillna('',inplace=True)
                name = star.Name
                """
                axs.text(star[coords[0]],star[coords[1]],f'{name}',
                         color='w',fontsize=fontsize(star.Vmag))
                         """
                axs.annotate(f'{name}',xy=[15*star[coords[0]],star[coords[1]]],xycoords='data',
                             xytext=[5,5],textcoords='offset points',
                             fontsize=fontsize(star.Vmag))

        # Decoration
        axs.set_xlabel(f'{coords[0]} [hh:mm]',fontsize=10)
        axs.set_ylabel(f'{coords[1]} [deg]',fontsize=10)
        
        # Range
        rang = max(((15*self.data[coords[0]]).max()-(15*self.data[coords[0]]).min()),
                   (self.data[coords[1]]).max()-(self.data[coords[1]]).min())
        axs.margins(pad*rang)
        
        axs.grid(alpha=0.2)
        axs.axis('equal')
        fig.tight_layout()

        # Change tick labels
        ra_ticks = axs.get_xticks()
        ra_tick_labels = []
        for ra in ra_ticks:
            comps = D2H(ra/15,string=False)
            ra_tick_labels += [f'{int(comps[0]):02d}:{comps[1]:02d}']
        axs.set_xticklabels(ra_tick_labels)

        dec_ticks = axs.get_yticks()
        dec_tick_labels = []
        for dec in dec_ticks:
            comps = D2H(dec,string=False)
            dec_tick_labels += [f'{int(comps[0]):02d}:{comps[1]:02d}']
        axs.set_yticklabels(dec_tick_labels,rotation=90)

        # Montu water mark
        Util.montu_mark(axs)

        SET_PLT_DEFAULT_STYLE()
        return fig,axs
    
    def __repr__(self):
        repr = f"Stars(number={len(self.data)})"
        return repr

    def __str__(self):
        desc = str(self.data)
        return desc