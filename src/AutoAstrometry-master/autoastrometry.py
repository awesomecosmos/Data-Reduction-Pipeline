"""
This script has many abilities. Using astroquery it can upload the image to astrometry.net and get a plate solution.
Along with this it can use the SIMBAD service to get coordinates for your target. Then using astropy it can convert
those coordinates to pixel coordinates within the image and get them back to you.
"""

import os
import webbrowser
from astropy.wcs.wcsapi.fitswcs import SlicedFITSWCS
from astroquery.astrometry_net import AstrometryNet
from astroquery.simbad import Simbad
from rich import print
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

# Helper functions ( basically everything outside the FITSU class )


def ask_for(prompt, error_msg=None, _type=None):
    """ While the desired prompt is not given, it repeats the prompt. """
    while True:
        inp = input(prompt).strip()
        if not inp:
            if error_msg:
                print(error_msg)
            continue

        if _type:
            try:
                inp = _type(inp)
            except ValueError:
                if error_msg:
                    print(error_msg)
                continue

        return inp


def redirect_to(url=str):
    """
    Takes in a URL as a string and asks the user if they would like to be redirected
    to that URL.

    Args:
        url (str, optional): URL that the user needs to be redirected to. Defaults to google.com.
    """
    invalid = True
    while invalid:
        # Asks the user if they want to redirected to the website
        print(
            f'\nWould you like to be [bold blue]redirected[/] to: {url} ?')
        redirect_to = ask_for(
            '(y/n): ', error_msg='\nPlease put in a letter. Y or N', _type=str).lower()

        #  If they do want to go to the website it opens it then breaks out of the loop.
        if redirect_to[0] == 'y':
            webbrowser.open_new_tab(f'{url}')
            invalid = False

        # They don't want to go to the website so it breaks out of the loop.
        elif redirect_to[0] == 'n':
            print('\nContinuing...')
            invalid = False


def simbad_query():
    """
    Looks up target and prints out the information such as RA and Dec.
    When using the function, it automatically asks the user for a target.
    """

    # Tries to query the target using astroquery if it fails it just opens the SIMBAD website
    try:
        # Tells user why they need this.
        print(
            '\nTo find the [bold blue]RA and Dec[/] of your target, please put it in here.')
        print("If your target can't be found, it will automatically redirect you to the website to put it in again.")

        #  Asks for target name and tries to look it up then if it can, prints it out.
        target = ask_for('\nTarget name: ')
        target = target.strip(' ')
        query = Simbad.query_object(f'{target}')
        print('\n*********************************************************************************************************************************************')
        print('[bold blue]Target info[/]:')
        query.pprint()
        print('*********************************************************************************************************************************************')

        #  Asks the user if they wanted to be redirected to the website.
        redirect_to('http://simbad.u-strasbg.fr/simbad/sim-fbasic')
    except:
        # * If theres an error it automatically just opens the website.
        webbrowser.open('http://simbad.u-strasbg.fr/simbad/sim-fbasic')


def find_image():
    look_for = True

    while look_for:
                # * Asks for the directory of the FITS file.
        print('\n-------------------------------------------------------')
        print('What is the path to the image file or directory?')
        print(
            'Please make sure the file [bold blue]ends in .FITS, JPEG, or .PNG[/]')
        print('-------------------------------------------------------')

        print(
            '\nYou can put in a path for a directory or file. \n( Directory EX: /Users/user/Pictures | File EX: /Users/user/Pictures/Kepler-1b.FITS )')
        print('\nTo [bold blue]copy a file path[/] on Mac right click on the file or directory \nthen hold alt, then select [bold blue]copy as pathname[/].')
        dir_or_file = ask_for('\n: ', _type=str)

        # * Supported file extensions (mainly those just supported by astrometry.net)
        file_extensions = ['.FITS', '.JPEG', '.PNG',
                           '.FIT', '.fits', '.fit', '.fts']

        #  File counter
        counter = 0

        # * Given path is directory
        if os.path.isdir(dir_or_file):
            #  For every file in the directory print out how many there are.
            for file_name in os.listdir(dir_or_file):
                counter += 1
                print(f'\nFile {counter}: ')
                #  Checks if the file is one of the supported file extensions.
                if file_name.endswith(tuple(file_extensions)):
                        #  Prints out file path
                    print(
                        f'{os.path.join(dir_or_file)}/{file_name}')

                    look_for = False
            #  Asks for the image they want to upload.
            print(
                '\n-------------------------------------------------------------------------------------------')
            print(
                'Which [bold blue]image[/] would you like to [bold blue]upload[/]? (One of the paths above)')
            not_file = True
            while not_file:
                upload_image = ask_for(': ', _type=str)

                # * Checks that the path given is a file
                if os.path.isfile(upload_image):
                    not_file = False
                else:
                    print(
                        '\nPlease [bold blue]re-enter[/] the path, the one you gave was invalid.')
                    print(
                        'If you accidentally put in / , just go into Finder and copy the file as a pathname, then paste it here.')
            print(
                '-------------------------------------------------------------------------------------------')
            # Unbound ( might not return anything )
            return upload_image

        # * Given path is file
        if os.path.isfile(dir_or_file):

            # * Finds out if the file is one of the supported file extensions.
            if dir_or_file.endswith(tuple(file_extensions)):

                # * Prints out path and breaks out of loop.
                print(f'\n{os.path.join(dir_or_file)}')
                look_for = False
                return dir_or_file


class FITSUploader():
    def __init__(self, image_path):
        self.image_path = image_path

    def upload_file(self):

        #  Creating instance of astrometry.net and API key
        ast = AstrometryNet()
        ast.api_key = 'bchkvzadjuswddhg'

        try_again = True
        submission_id = None

        while try_again:
            try:
                if not submission_id:
                    # Solves the image from the file path
                    wcs_header = ast.solve_from_image(f'{self.image_path}', force_image_upload=True,
                                                      submission_id=submission_id, solve_timeout=1000)
                else:
                    #  Time is in seconds.
                    wcs_header = ast.monitor_submission(
                        submission_id, solve_timeout=1000)
            except TimeoutError as e:
                #  Timeout error, never triggers. Basically useless code since it never triggers during timeout error
                submission_id = e.args[1]
                print('\nThere was a timeout error. ( Process took to long ).')
                print('Astometry.net could also be down at the moment.')
            else:
                #! got a result, so terminate
                try_again = False

        if wcs_header:
            #  Code to execute when solve succeeds
            print('\nSuccess! :thumbs_up:')
            print(
                '\nTo get the most possible information out of your image please visit the website below.')
            redirect_to('http://nova.astrometry.net/users/20995')

            print('\n*********************************************************************************************************************************************')

            # Looks up target with astroquery then can inp.redirect_to user to the website
            # to use the aladin lite view to find comp stars, look around, etc.
            simbad_query()

        else:
            #! Code to execute when solve fails
            print('\n[bold red]Failed[/bold red] to solve.')

        # do other shit here
    def find_px_coords(self, solved_image):
        not_radec = True
        while not_radec:
            try:
                # * Uses SkyCoord to verify correct RA and Declinaction values
                print('\n*********************************************************************************************************************************************')
                print(
                    '[bold blue]Right Ascension[/] and [bold blue]Declination[/] for your target.')
                print(
                    'Enter the values [bold]one at a time[/]. ( EX: 19 ( hit enter ), 07 ( hit enter ), 14 ( hit enter )')
                print('The same applies for [bold blue]Declination values[/].')
                print('*********************************************************************************************************************************************')

                #  RA
                print('\n*********************************************************************************************************************************************')

                ra1 = ask_for('\n: ', _type=int)
                ra2 = ask_for(': ', _type=int)
                ra3 = ask_for(': ', _type=int)

                #  Declination
                print("\n[bold blue]Declination[/], don't forget the + or -")
                dec1 = ask_for('\n: ', _type=int)
                dec2 = ask_for(': ', _type=int)
                dec3 = ask_for(': ', _type=int)
                print('*********************************************************************************************************************************************')

                #! Just keeping this here to ensure that the right values for RA and Dec are correct.
                c = SkyCoord(f'{ra1}h{ra2}m{ra3}s',
                             f'{dec1}d{dec2}m{dec3}s', frame='fk5', unit='deg')

            # * Value error, user didn't enter in the right values
            except ValueError:
                print('\n[red]Value Error ocurred[/]')
                print('Please [bold blue]re-enter your RA and Dec[/]')
            else:
                # Got a result so break out of loop
                not_radec = False

        def pixel_pos():
            not_image = True
            while not_image:
                try:
                    # * Asks user to put in the path to the plate solved image from https://nova.astrometry.net.
                    print('\n*********************************************************************************************************************************************')
                    print(
                        'Please put in the [bold blue]plate solved image[/] from https://nova.astrometry.net.')
                    print('It should be titled [bold blue]new-image.fits[/].')
                    print('*********************************************************************************************************************************************')

                    filename = solved_image

                    # Opens the file and looks at header data
                    hdu = fits.open(filename)
                    header = hdu[0].header

                    #  Applies WCS to header ( world coordinate system ). Also checks that the RA and Dec values are in fk5.
                    wcs = WCS(header)
                    coord = SkyCoord(
                        f'{ra1}h{ra2}m{ra3}s {dec1}d{dec2}m{dec3}s', frame='fk5')

                    # * Converts the RA and Dec values to pixel values within the image. It then also prints them out.
                    px = wcs.world_to_pixel(coord)
                    print('\n*********************************************************************************************************************************************')
                    print('Pixel coordinates:')
                    print(px)
                    print('*********************************************************************************************************************************************')
                    return px

                except:
                    #  File that was given was not the plate solved image from https://nova.astrometry.net.
                    print(
                        '\nPlease put in the [bold blue]plate solved image[/] from https://nova.astrometry.net.')
                else:
                    # Got correct result so break out of loop
                    not_image = False
        pixel_pos()

    def check_comp_stars(self):
        comp_stars = ask_for(
            '\nDo you have any comparison stars you would like to get the pixel coordinates from? (y/n): ', 'Not a yes or no', str).lower()
        if comp_stars[0] == 'y':
            num_comp_stars = ask_for(
                '\nHow many comparison stars do you have? ( Integer )', 'Not an integer.', int)
            counter = 0
            for _ in range(num_comp_stars):
                counter += 1
                print(f'\n[bold blue]Comparison star {counter}[/]:')
                fitsu.find_px_coords(find_image())
        elif comp_stars[0] == 'n':
            print('\nContinuing...')


if __name__ == "__main__":
    print('\n[bold blue]Following prompt is for the image to be uploaded to https://nova.astrometry.net[/]')
    fitsu = FITSUploader(find_image())
    fitsu.upload_file()

    find_coords = ask_for(
        'Would you like to find the pixel coordinates of your target? ( IMPORTANT: You must use the image from astrometry.net  y/n ): ', 'error, put in string', str).lower()
    print('*********************************************************************************************************************************************')

    if find_coords[0] == 'y':
        fitsu.find_px_coords(find_image())
        fitsu.check_comp_stars()
        print('\n**************\n[bold blue]End of program[/]\n**************')
    elif find_coords[0] == 'n':
        print('\n**************\n[bold blue]End of program[/]\n**************')
