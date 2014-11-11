"""
@package mi.dataset.parser.test.test_ctdmo_ghqr_imodem
@file mi/dataset/parser/test/test_ctdmo_ghqr_imodem.py
@author Maria Lutz
@brief test of the parser for ctdmo_ghqr_imodem dataset driver
"""

__author__ = 'Maria Lutz'
__license__ = 'Apache 2.0'

import os
from nose.plugins.attrib import attr

from mi.core.log import get_logger
log = get_logger()
from mi.core.exceptions import SampleException, RecoverableSampleException
from mi.dataset.test.test_parser import ParserUnitTestCase, BASE_RESOURCE_PATH
from mi.dataset.parser.ctdmo_ghqr_imodem import CtdmoGhqrImodemParser

RESOURCE_PATH = os.path.join(BASE_RESOURCE_PATH, 'ctdmo_ghqr', 'imodem', 'resource')


@attr('UNIT', group='mi')
class CtdmoGhqrImodemParserUnitTestCase(ParserUnitTestCase):


    def test_simple(self):
        """
        Test a simple telemetered and recovered case 
        """
        with open(os.path.join(RESOURCE_PATH, 'ctdmo01_20140712_120719.DAT'), 'r') as file_handle:
            parser = CtdmoGhqrImodemParser(file_handle, self.exception_callback, is_telemetered=True)
        
            particles = parser.get_records(4)

            self.assert_particles(particles, "ctdmo01_20140712_120719_telem.yml", RESOURCE_PATH)
        
        with open(os.path.join(RESOURCE_PATH, 'ctdmo01_20140712_120719.DAT'), 'r') as file_handle:
            parser = CtdmoGhqrImodemParser(file_handle, self.exception_callback, is_telemetered=False)

            particles = parser.get_records(4)

            self.assert_particles(particles, "ctdmo01_20140712_120719_recov.yml", RESOURCE_PATH)

        self.assertEqual(self.exception_callback_value, [])
        

    def test_bad_number_samples(self):
        """
        Test bad number of samples. The test file bad_num_samples_ctdmo01_20140712_120719.DAT has 4 samples
        but header reports 62.
        """
        with open(os.path.join(RESOURCE_PATH, 'bad_num_samples_ctdmo01_20140712_120719.DAT'), 'r') as file_handle:
            parser = CtdmoGhqrImodemParser( file_handle, self.exception_callback, is_telemetered=True)

            n_test = 4  # number of particles to test
            particles = parser.get_records(n_test)

            # make sure none of the particles succeeded
            self.assertEqual(len(particles), 0)
            # check that there were 3 recoverable sample exceptions
            self.assertEqual(len(self.exception_callback_value), n_test)
            for i in range(0, n_test):
                self.assert_(isinstance(self.exception_callback_value[i], RecoverableSampleException))

    def test_unexpected(self):
        """
        Test with an unexpected line, confirm we get an exception. The second line of test file
        unexpected_data_ctdmo_ghqr_imodem.DAT contains unxpected data.
        """
        with self.assertRaises(SampleException):
            with open(os.path.join(RESOURCE_PATH, 'unexpected_data_ctdmo_ghqr_imodem.DAT'), 'r') as file_handle:
                parser = WavssADclParser(file_handle, self.exception_callback, is_telemetered=True)

                n_test = 4  
                particles = parser.get_records(n_test)
