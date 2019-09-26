"""
The aim of this file is to define a representation class for tree printing.
It is useful to switch between the 2 for terminal or text file output.
"""

class Representation(object):
    """ Contains all things necessary for representing my nodes and trees"""
    def __init__(self, delimiter = "|", color = "red", printing_solved = "- solved"):
        self.delimiter = delimiter  # Delimiter between nodes
        if color == "red":
            self.color_begin = '\033[91m'
            self.color_end = '\033[0m'
        elif color == "":
             self.color_begin = ''
             self.color_end = ''
        else:
            raise NotImplementedError
        self.printing_solved = printing_solved

Test_representation = Representation(delimiter = "|", color = "red", printing_solved = "")
Test_to_file = Representation(delimiter = "|", color = "", printing_solved = "- solved")
