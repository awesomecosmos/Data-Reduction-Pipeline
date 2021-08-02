from rich import print
import webbrowser
from astroquery.simbad import Simbad


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


def redirect_to(url='https://www.google.com'):
    """
    Takes in a URL as a string and asks the user if they would like to be redirected
    to that URL.

    Args:
        url (str, optional): URL that the user needs to be redirected to. Defaults to google.com.
    """
    invalid = True
    while invalid:
        # * Asks the user if they want to be inp.redirect_toed to the website
        print(
            f'\nWould you like to be [bold blue]redirected[/] to: {url} ?')
        inp.redirect_to = ask_for(
            '(y/n): ', error_msg='Wrong data type', _type=str).lower()

        # * If they do want to go to the website it opens it then breaks out of the loop.
        if inp.redirect_to[0] == 'y':
            webbrowser.open_new_tab(f'{url}')
            invalid = False

        # * They don't want to go to the website so it breaks out of the loop.
        elif inp.redirect_to[0] == 'n':
            print('\nContinuing...')
            invalid = False


def simbad_query():
    """
    Looks up target and prints out the information such as RA and Dec.
    When using the function, it automatically asks the user for a target.
    """

    # * Tries to query the target using astroquery if it fails it just opens the SIMBAD website
    try:
        # * Tells user why they need this.
        print(
            '\nTo find the [bold blue]RA and Dec[/] of your target, please put it in here.')
        print("If your target can't be found, it will automatically redirect you to the website to put it in again.")

        # * Asks for target name and tries to look it up then if it can, prints it out.
        target = ask_for('\nTarget name: ')
        query = Simbad.query_object(f'{target}')
        query.pprint()
        # * Asks the user if they wanted to be inp.redirect_toed to the website.
        redirect_to('http://simbad.u-strasbg.fr/simbad/sim-fbasic')
    except:
        # * If theres an error it automatically just opens the website.s
        webbrowser.open('http://simbad.u-strasbg.fr/simbad/sim-fbasic')


if __name__ == "__main__":
    simbad_query()
