#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This script contains unit tests of the :mod:`arkane.main` module.
"""

import unittest
import os
import shutil

import rmgpy
from arkane import Arkane

################################################################################


class TestArkaneExamples(unittest.TestCase):
    """
    Run all of Arkane's examples, and report which one failed
    """
    @classmethod
    def setUp(self):
        """A function that is run ONCE before all unit tests in this class."""
        self.base_path = os.path.join(os.path.dirname(os.path.dirname(rmgpy.__file__)), 'examples', 'arkane')
        self.failed = []
        self.example_types = ['species', 'reactions', 'networks']

    def test_arkane_examples(self):
        for example_type in self.example_types:
            example_type_path = os.path.join(self.base_path, example_type)
            examples = []
            for (dirpath, dirnames, filenames) in os.walk(example_type_path):
                examples.extend(dirnames)
                break
            for example in examples:
                path = os.path.join(example_type_path, example)
                arkane = Arkane(inputFile=os.path.join(path, 'input.py'), outputDirectory=path)
                arkane.plot = True
                arkane.execute()
                with open(os.path.join(path, 'arkane.log'), 'r') as f:
                    log = f.readlines()
                for line in log[::-1]:
                    if 'execution terminated' in line:
                        break
                else:
                    self.failed.append([example_type, example + 'FreeRotor'])
        error_message = 'Arkane example(s) failed: '
        for type_name_tuple in self.failed:
            error_message += '{1} in {0}; '.format(type_name_tuple[0], type_name_tuple[1])
        self.assertTrue(len(self.failed) == 0, error_message)

    @classmethod
    def tearDown(self):
        """A function that is run ONCE after all unit tests in this class."""
        self.extensions_to_delete = ['pdf', 'csv', 'txt', 'inp']
        self.files_to_delete = ['arkane.log', 'output.py']
        self.files_to_keep = ['README.txt']  # files to keep that have extentions marked for deletion
        self.base_path = os.path.join(os.path.dirname(os.path.dirname(rmgpy.__file__)), 'examples', 'arkane')
        self.example_types = ['species', 'reactions', 'networks']
        for example_type in self.example_types:
            example_type_path = os.path.join(self.base_path, example_type)
            examples = []
            for (dirpath, dirnames, filenames) in os.walk(example_type_path):
                examples.extend(dirnames)
                break
            for example in examples:
                # clean working folder from all previous test output
                example_path = os.path.join(example_type_path, example)
                dirs = [d for d in os.listdir(example_path)
                        if not os.path.isfile(os.path.join(example_path, d))]
                for d in dirs:
                    shutil.rmtree(os.path.join(example_path, d, ''))
                files = [f for f in os.listdir(example_path) if os.path.isfile(os.path.join(example_path, f))]
                for f in files:
                    extension = f.split('.')[-1]
                    if f in self.files_to_delete or\
                            (extension in self.extensions_to_delete and f not in self.files_to_keep):
                        os.remove(os.path.join(example_path, f))


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
