"""
The aim of this script is to change configuration file from command line.
It takes as input the base config file from data.
"""

import argparse
import re
import os

def __cli():
    """
    Command line interface.
    """

    d = "Arguments to change the config file before running a Tree"
    parser = argparse.ArgumentParser(description=d)
    # Logs and saving information
    parser.add_argument("--DB_CACHE", type = lambda x: (str(x).lower() == 'true'), default=False)
    parser.add_argument("--DB_REPLACE", type = lambda x: (str(x).lower() == 'true'), default=False)
    parser.add_argument("--DB_time", default = 1, type=float)
    parser.add_argument("--biosensor", type = lambda x: (str(x).lower() == 'true'), default=False)
    parser.add_argument("--use_cache", type = lambda x: (str(x).lower() == 'true'), default=False)
    parser.add_argument("--add_Hs", type = lambda x: (str(x).lower() == 'true'), default=False)
    parser.add_argument("--use_transpositions", type = lambda x: (str(x).lower() == 'true'), default=False)
    parser.add_argument("--use_transpositions_depth", type = lambda x: (str(x).lower() == 'true'), default=False)
    parser.add_argument("--folder_to_save", default=os.path.dirname(os.path.abspath(__file__)))
    args = parser.parse_args()

    def change_dB_setting(DB_CACHE,
                          DB_REPLACE,
                          DB_time,
                          biosensor,
                          use_cache,
                          add_Hs,
                          use_transpositions,
                          use_transpositions_depth,
                          folder_to_save):
        with open("{}/data/base_config.py".format(os.path.dirname(os.path.abspath(__file__))), 'r') as file_original:
            whole_text = file_original.read()
        with open("{}/config.py".format(folder_to_save), "w") as replacement_text:
            # Changing DB_cache
            if DB_CACHE:
                if "DB_CACHE = True"  not in whole_text:
                    whole_text = whole_text.replace("DB_CACHE = False", "DB_CACHE = True")
            else:
                if "DB_CACHE = False"  not in whole_text:
                    whole_text = whole_text.replace("DB_CACHE = True", "DB_CACHE = False")
            # Changing DB replace
            if DB_REPLACE:
                if "DB_REPLACE = True" not in whole_text:
                    whole_text = whole_text.replace("DB_REPLACE = False", "DB_REPLACE = True")
            else:
                if "DB_REPLACE = False" not in whole_text:
                    whole_text = whole_text.replace("DB_REPLACE = True", "DB_REPLACE = False")
            # Changing DB_time:
            whole_text = re.sub("DB_time = \d+.\d+", 'DB_time = {}'.format(DB_time), whole_text)

            # Changing running mode from biosensor to retrosynthesis
            if biosensor:
                if "biosensor = True" not in whole_text:
                    whole_text = whole_text.replace("biosensor = False", "biosensor = True")
                    whole_text = whole_text.replace("retrosynthesis = True", "retrosynthesis = False")
            else:
                if "biosensor = False" not in whole_text:
                    whole_text = whole_text.replace("biosensor = True", "biosensor = False")
                    whole_text = whole_text.replace("retrosynthesis = False", "retrosynthesis = True")
            # Changing use_cache
            if use_cache:
                if "use_cache = True" not in whole_text:
                    whole_text = whole_text.replace("use_cache = False", "use_cache = True")
            else:
                if "use_cache = False"  not in whole_text:
                    whole_text = whole_text.replace("use_cache = True", "use_cache = False")

            # Hydrogen handling:
            if add_Hs:
                if "add_Hs = True" not in whole_text:
                    whole_text = whole_text.replace("add_Hs = False", "add_Hs = True")
            else:
                if "add_Hs = False"  not in whole_text:
                    whole_text = whole_text.replace("add_Hs = True", "add_Hs = False")

            # Changing use_transpositions
            if use_transpositions:
                if "use_transpositions = True" not in whole_text:
                    whole_text = whole_text.replace("use_transpositions = False", "use_transpositions = True")
            else:
                if "use_transpositions = False"  not in whole_text:
                    whole_text = whole_text.replace("use_transpositions = True", "use_transpositions = False")
            # Changing use_transpositions_depth
            if use_transpositions_depth:
                if "use_transpositions_depth = True" not in whole_text:
                    whole_text = whole_text.replace("use_transpositions_depth = False", "use_transpositions_depth = True")
            else:
                if "use_transpositions_depth = False"  not in whole_text:
                    whole_text = whole_text.replace("use_transpositions_depth = True", "use_transpositions_depth = False")
            replacement_text.write(whole_text)


    change_dB_setting(DB_CACHE = args.DB_CACHE,
                      DB_REPLACE = args.DB_REPLACE,
                      DB_time = args.DB_time,
                      biosensor = args.biosensor,
                      use_cache = args.use_cache,
                      add_Hs = args.add_Hs,
                      use_transpositions = args.use_transpositions,
                      use_transpositions_depth = args.use_transpositions_depth,
                      folder_to_save = args.folder_to_save)


if __name__ == "__main__":
    __cli()
