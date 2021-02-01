import matplotlib.pyplot as plt
import matplotlib.image as mgimg
from matplotlib import animation
import escher
from escher import Builder
import os


from skimpy.viz.utils import EMBEDD_CSS

DEFAULT_CHROME='/opt/google/chrome/chrome'

def plot_fluxes(flux_dict,
                escher_map,
                output_file='map.html',
                height=600,
                width=800,
                reaction_scale=None,
                min_flux = -10,
                max_flux = 10,):

    if min_flux is None:
        min_flux = min(flux_dict)

    if max_flux is None:
        max_flux = max(flux_dict)

    if reaction_scale is None:
        reaction_scale = [
            { 'type': 'value', 'value': min_flux, 'color': 'red', 'size': 32 },
            { 'type': 'value', 'value': 0, 'color': '#c8c8c8', 'size': 12 },
            { 'type':  'value', 'value': max_flux, 'color':'green', 'size': 32 }
        ]

    builder = Builder(
        height=height,
        width=width,
        map_json=escher_map,
        reaction_scale=reaction_scale
    )

    builder.reaction_data = flux_dict
    builder.save_html(output_file)
    builder.close()

def animate_fluxes(flux_time_data,
                   escher_map,
                   outputfile='animation.mp4',
                   height=600,
                   width=800,
                   time_interval_ms=100,
                   chrome=DEFAULT_CHROME,
                   reaction_scale=None,
                   min_flux=None,
                   max_flux=None,
                   time_size=12,
                   time_unit='h',
                   x_time=0.95,
                   y_time=0.9,
                   ):

    if min_flux is None:
        min_flux = flux_time_data.min().min()

    if max_flux is None:
        max_flux = flux_time_data.max().max()

    if reaction_scale is None:
        reaction_scale = [
            { 'type': 'value', 'value': min_flux, 'color': 'red', 'size': 32 },
            { 'type': 'value', 'value': 0, 'color': '#c8c8c8', 'size': 12 },
            { 'type':  'value', 'value': max_flux, 'color':'green', 'size': 32 }
        ]

    builder = Builder(
        height=height,
        width=width,
        map_json=escher_map,
        embedd_css=EMBEDD_CSS,
        menu='none',
        reaction_scale=reaction_scale
    )

    myimages = []

    fig = plt.figure()

    XVFB_DOCKER = '/usr/bin/xvfb-run -a -s \"-screen 0 {}x{}x24\"'.format(width,height)
    SCREENSHOT = "--headless --disable-gpu --no-sandbox  --virtual-time-budget=10000 --screenshot=\'{}\' {}"

    for t, fluxes in flux_time_data.iterrows():
        builder.reaction_data = fluxes
        builder.save_html('tmp.html',)

        # Hacky hack hack ...
        # Use chrome to make a screenshot
        cmd = "{} {} {}".format(XVFB_DOCKER, chrome, SCREENSHOT.format('tmp.png', 'tmp.html'))
        os.system(cmd)

        # Add time text
        ylim = plt.gca().get_ylim()
        xlim = plt.gca().get_xlim()
        y = (ylim[1] - ylim[0]) * y_time + ylim[0]
        x = (xlim[1] - xlim[0]) * x_time + xlim[0]
        text = plt.text(x, y, '{:.1f} {}'.format(t,time_unit),
                        horizontalalignment='right',
                        fontsize=time_size)


        img = mgimg.imread('tmp.png')
        imgplot = plt.imshow(img)

        # append AxesImage object to the list
        myimages.append([imgplot,text])

    plt.axis('off')
    fig.tight_layout()
    anim = animation.ArtistAnimation(fig, myimages, interval=time_interval_ms)
    anim.save(outputfile, dpi=300)

    builder.close()

