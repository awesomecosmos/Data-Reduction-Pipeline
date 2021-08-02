from astroquery.simbad import Simbad
from time import sleep as s
import webbrowser
from rich import print


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


class Finder():
    """
    Using a target name, and the number of targets you would like to query
    it can use the SIMBAD query to look up those targets and provide you with the necessary information.
    """

    def __init__(self, target_name=str, multiple_queries=int):
        self.target_name = target_name
        self.multiple_queries = multiple_queries

    def find_target(self):
        """
        Finds target with SIMBAD query. It can also look up multiple objects.
        """

        # User is querying more than one object
        if self.multiple_queries > 1:
            for _ in range(self.multiple_queries):
                print(
                    '\nThese prompts [bold blue]will repeat for every target you have[/].')
                self.target_name = ask_for(
                    '\nTarget name: ', 'Not a string', str)
                print('\n[bold]Looking up target...[/]')
                try:
                    query = Simbad.query_object(self.target_name)
                    print('\n[bold blue]Target found[/]! Information below:')
                    s(0.25)
                    query.pprint()
                except:
                    print('\nTarget could not be found, opening website.')
                    s(1.5)
                    webbrowser.open(
                        'http://simbad.u-strasbg.fr/simbad/sim-fbasic')

        # User is only going to query one object
        elif self.multiple_queries == 1:
            print('\n[bold]Looking up target...[/]')
            try:
                query = Simbad.query_object(self.target_name)
                print('\n[bold blue]Target found[/]! Information below:')
                s(0.25)
                query.pprint()
            except:
                print('\nTarget could not be found, opening website.')
                s(1.5)
                webbrowser.open(
                    'http://simbad.u-strasbg.fr/simbad/sim-fbasic')


if __name__ == "__main__":
    # Getting target name
    print('\n[bold blue]Target name[/]')
    target = ask_for('\n: ', 'Not a string', str)

    # Number of targets to query
    print('\nDo you have [Bold blue]more than one target[/] to query? ( y/n )')
    multiple_queries_answer = ask_for('\n: ', 'Not an str', str).lower()

    # User has more than one target to query
    if multiple_queries_answer[0] == 'y':
        print('\nPlease input the number of targets')
        num_targets = ask_for(': ', 'Not an integer', int)
        f = Finder(target_name=target, multiple_queries=num_targets)
        f.find_target()

    # User has only one target to query
    elif multiple_queries_answer[0] == 'n':
        print('\nContinuing...')
        f = Finder(target, 1)
        f.find_target()
