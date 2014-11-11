#!/usr/bin/env python

"""
@package mi.dataset.parser.ctdmo_ghqr_imodem
@file mi-dataset/mi/dataset/parser/ctdmo_ghqr_imodem.py
@author Maria Lutz
@brief Parser for the ctdmo_ghqr_imodem dataset driver

Release notes:

Initial release
"""

__author__ = 'Maria Lutz'
__license__ = 'Apache 2.0'
import ntplib
import binascii
import re
import numpy as np
from datetime import datetime
import time

from mi.core.common import BaseEnum
from mi.core.log import get_logger
log = get_logger()
from mi.core.instrument.data_particle import DataParticle
from mi.core.exceptions import SampleException, RecoverableSampleException
from mi.dataset.parser.utilities import zulu_timestamp_to_ntp_time, formatted_timestamp_utc_time


from mi.dataset.dataset_parser import SimpleParser

from mi.dataset.parser.common_regexes import INT_REGEX

# Define regexes
END_OF_LINE_REGEX = r'(?:\r\n|\n)?'  # end of file might be missing terminator so make optional
ZERO_OR_MORE_WHITESPACES = r'\s*'

STATUS_REGEX = r'#MCAT Status' + END_OF_LINE_REGEX
STATUS_MATCHER = re.compile(STATUS_REGEX)
LOG_STATUS_REGEX = r'#logging|#Begin Data|#End Data' + END_OF_LINE_REGEX
LOG_STATUS_MATCHER = re.compile(LOG_STATUS_REGEX)
DATA_FORMAT_REGEX = r'#Data Format:.+'  + END_OF_LINE_REGEX
DATA_FORMAT_MATCHER = re.compile(DATA_FORMAT_REGEX)

DATETIME_REGEX = r'#.+DateTime:' + ZERO_OR_MORE_WHITESPACES + '(\d{8}\s\d{6})' + END_OF_LINE_REGEX
DATETIME_MATCHER = re.compile(DATETIME_REGEX)
SERIAL_NO_REGEX = r'#.+SERIAL NO.' + ZERO_OR_MORE_WHITESPACES + '(\d+) .+' + END_OF_LINE_REGEX
SERIAL_NO_MATCHER = re.compile(SERIAL_NO_REGEX)
VOLTAGE_REGEX = r'#vMain =' + ZERO_OR_MORE_WHITESPACES + '(\d+\.\d{2}).+vLith =' + ZERO_OR_MORE_WHITESPACES + '(\d+\.\d{2})' + END_OF_LINE_REGEX
VOLTAGE_MATCHER = re.compile(VOLTAGE_REGEX)
SAMP_NUM_MEM_REGEX = r'#samplenumber =' + ZERO_OR_MORE_WHITESPACES + '(\d+).+free = (\d+)' + END_OF_LINE_REGEX
SAMP_NUM_MEM_MATCHER = re.compile(SAMP_NUM_MEM_REGEX)
SAMP_INTERVAL_REGEX = r'#sample interval =' + ZERO_OR_MORE_WHITESPACES + '(\d+) seconds' + END_OF_LINE_REGEX
SAMP_INTERVAL_MATCHER = re.compile(SAMP_INTERVAL_REGEX)
PRESSURE_RANGE_REGEX = r'#PressureRange =' + ZERO_OR_MORE_WHITESPACES + '(\d+)' + END_OF_LINE_REGEX
PRESSURE_RANGE_MATCHER = re.compile(PRESSURE_RANGE_REGEX)
SAMPS_RECORDED_REGEX = r'#SamplesRecorded =' + ZERO_OR_MORE_WHITESPACES + '(\d+)' + END_OF_LINE_REGEX
SAMPS_RECORDED_MATCHER = re.compile(SAMPS_RECORDED_REGEX)

# Table contains state (1 bit indicates data received/parsed), matcher, & value.
METADATA_STATE_TABLE = [
    [0x01, DATETIME_MATCHER, 0],
    [0x02, SERIAL_NO_MATCHER, 0],
    [0x04, VOLTAGE_MATCHER, 0],
    [0x08, SAMP_NUM_MEM_MATCHER, 0],
    [0x10, SAMP_INTERVAL_MATCHER, 0],
    [0x20, PRESSURE_RANGE_MATCHER, 0],
    [0x40, SAMPS_RECORDED_MATCHER, 0]
]
ALL_METADATA_RECEIVED = 0x7F
METADATA_VALUE_INDEX = 2

# each line of format tttttcccccppppTTTTTTTT.
# ttttt: temp
# ccccc: conductivity
# pppp: pressure
# TTTTTTTT: time
SCIENCE_DATA_REGEX = b'([0-9a-fA-F]{5})([0-9a-fA-F]{5})([0-9a-fA-F]{4})([0-9a-fA-F]{8})' + END_OF_LINE_REGEX
SCIENCE_DATA_MATCHER = re.compile(SCIENCE_DATA_REGEX)

# particle maps consists of tuples containing the name of parameter, index into regex group, encoding type
# This table is used in the generation of the Instrument data particle.
# Column 1 - particle parameter name
# Column 2 - index into raw_data
# Column 3 - data encoding function (conversion required - int, float, etc)
# each line of format tttttcccccppppTTTTTTTT. 
# is the particle supposed to have that first param 'time'?
CTD_TIME_GROUP = 4
INSTRUMENT_PARTICLE_MAP = [
    ('temperature',     1,  int),
    ('conductivity',    2,  int),
    ('pressure',        3,  int),
    ('ctd_time',        4,  int)
]

# This table is used in the generation of the metadata data particle.
# Column 1 - particle parameter name
# Column 2 - index into raw_data
# Column 3 - data encoding function (conversion required - int, float, etc)
METADATA_DATETIME_GROUP = 1; 
METADATA_PARTICLE_MAP = [
    ('date_time_string',                1,  str),
    ('instrument_serial_number_u32',    1,  int),
    ('battery_voltage_main',            1,  float),
    ('battery_voltage_lithium',         2,  float),
    ('sample_number',                   1,  int),
    ('mem_free',                        2,  int),
    ('sample_interval',                 1,  int),
    ('pressure_range',                  1,  int),
    ('num_samples',                     1,  int)
]

def generate_particle_timestamp(time_2000):
    """
    This function calculates and returns a timestamp in epoch 1900
    based on an ASCII hex time in epoch 2000.
    Parameter:
      time_2000 - number of seconds since Jan 1, 2000
    Returns:
      number of seconds since Jan 1, 1900
    """
    return int(time_2000, 16) + zulu_timestamp_to_ntp_time("2000-01-01T00:00:00.00Z")

class DataParticleType(BaseEnum):
    CTDMO_GHQR_IMODEM_INSTRUMENT_RECOVERED = 'ctdmo_ghqr_imodem_instrument_recovered'
    CTDMO_GHQR_IMODEM_METADATA_RECOVERED = 'ctdmo_ghqr_imodem_metadata_recovered'
    CTDMO_GHQR_IMODEM_INSTRUMENT = 'ctdmo_ghqr_imodem_instrument'
    CTDMO_GHQR_IMODEM_METADATA = 'ctdmo_ghqr_imodem_metadata'
    
class CtdmoGhqrImodemDataParticleKey(BaseEnum):
    # For metadata data particle
    DATE_TIME_STRING = 'date_time_string'
    INSTRUMENT_SERIAL_NUMBER_U32 = 'instrument_serial_number_u32'
    BATTERY_VOLTAGE_MAIN = 'battery_voltage_main'
    BATTERY_VOLTAGE_LITHIUM = 'battery_voltage_lithium'
    SAMPLE_NUMBER = 'sample_number'
    MEM_FREE = 'mem_free'
    SAMPLE_INTERVAL = 'sample_interval'
    PRESSURE_RANGE = 'pressurem_range'
    NUM_SAMPLES = 'num_samples'
    
    # For instrument data particle
    TEMPERATURE = 'temperature'
    CONDUCTIVITY = 'conductivity'
    PRESSURE = 'pressure'
    CTD_TIME = 'ctd_time'

class CtdmoGhqrImodemInstrumentDataParticle(DataParticle):
    """
    Class for generating the CTDMO instrument particle.
    """

    def _build_parsed_values(self):
        """
        Encode parameters from the raw data using the particle maps
        """

        ctd_time = generate_particle_timestamp(self.raw_data[CTD_TIME_GROUP][1])
        self.set_internal_timestamp(timestamp=ctd_time)
     
        # Generate a particle by calling encode_value for each entry
        # in the Instrument Particle Mapping table,
        # where each entry is a tuple containing the particle field name,
        # an index into raw_data and a function to use for data conversion.
        return [self._encode_value(name, self.raw_data[raw_index][1], function)
            for name, raw_index, function in INSTRUMENT_PARTICLE_MAP]

class CtdmoGhqrImodemInstrumentTelemeteredDataParticle(CtdmoGhqrImodemInstrumentDataParticle):
    _data_particle_type = DataParticleType.CTDMO_GHQR_IMODEM_INSTRUMENT

class CtdmoGhqrImodemInstrumentRecoveredDataParticle(CtdmoGhqrImodemInstrumentDataParticle):
    _data_particle_type = DataParticleType.CTDMO_GHQR_IMODEM_INSTRUMENT_RECOVERED


class CtdmoGhqrImodemMetadataDataParticle(DataParticle):
    """
    Class for generating the Metadata particle.
    """
    def _build_parsed_values(self):
        """
        Build parsed values for Recovered and Telemetered Metadata Data Particle.
        """
        # Generate a particle by calling encode_value for each entry
        # in the Metadata Particle Mapping table,
        # where each entry is a tuple containing the particle field name,
        # an index into raw_data and a function to use for data conversion.
        log.debug('self.raw_data[0][1]: %s', self.raw_data[1][1])
        return [self._encode_value(name, self.raw_data[raw_index][1], function)
            for name, raw_index, function in METADATA_PARTICLE_MAP]

class CtdmoGhqrImodemMetadataInstrumentDataParticle(CtdmoGhqrImodemMetadataDataParticle):

    _data_particle_type = DataParticleType.CTDMO_GHQR_IMODEM_METADATA

class CtdmoGhqrImodemMetadataRecoveredDataParticle(CtdmoGhqrImodemMetadataDataParticle):

    _data_particle_type = DataParticleType.CTDMO_GHQR_IMODEM_METADATA_RECOVERED


class CtdmoGhqrImodemParser(SimpleParser):
    def __init__(self,
                 stream_handle,
                 exception_callback,
                 is_telemetered):

        if is_telemetered:
            # this is a telemetered parser
            self.instrument_particle_class = CtdmoGhqrImodemInstrumentTelemeteredDataParticle
            self.metadata_particle_class = CtdmoGhqrImodemMetadataInstrumentDataParticle
            
        else:
            # this is a recovered parser
            self.instrument_particle_class = CtdmoGhqrImodemInstrumentRecoveredDataParticle
            self.metadata_particle_class = CtdmoGhqrImodemMetadataRecoveredDataParticle
           
        # no config for this parser, pass in empty dict
        super(CtdmoGhqrImodemParser, self).__init__({},
                                              stream_handle,
                                              exception_callback)
        
    def parse_file(self):
        
        # initialize
        self._metadata_state = 0
        
        # read the first line in the file
        line = self._stream_handle.readline()

        while line:

            datetime_match = DATETIME_MATCHER.match(line)
            serial_no_match = SERIAL_NO_MATCHER.match(line)
            voltage_match = VOLTAGE_MATCHER.match(line)
            samp_num_mem_match = SAMP_NUM_MEM_MATCHER.match(line)
            samp_interval_match = SAMP_INTERVAL_MATCHER.match(line)
            pressure_range_match = PRESSURE_RANGE_MATCHER.match(line)
            samps_recorded_match = SAMPS_RECORDED_MATCHER.match(line)
            science_data_match = SCIENCE_DATA_MATCHER.match(line)
            num_csv = len(line.split(','))
            
            # Does the line contain data needed for the metadata particle?
            if datetime_match or serial_no_match or voltage_match or samp_num_mem_match or samp_interval_match or \
                pressure_range_match or samps_recorded_match:
                log.debug("Found data for the metadata particle: %s", line)

                # Process the metadata record match
                self._process_metadata_record_part(line)
                
                # Attempt to generate metadata particles
                self._generate_metadata_particles()
                
            # Does the line contain instrument data?
            elif science_data_match:
                # create instrument particle
                particle = self._create_instrument_particle(science_data_match)
                if particle is not None:
                    self._record_buffer.append(particle)

            else:
                # Check for other lines that can be ignored
                log_match = LOG_STATUS_MATCHER.match(line)
                status_match = STATUS_MATCHER.match(line) 
                data_format_match = DATA_FORMAT_MATCHER.match(line)

                # Raise exception if we encounter unexpected data
                if not (log_match or status_match or data_format_match):                   
                    raise SampleException("Ctdmo encountered unexpected data line '%s'" % line)

            # read the next line in the file
            line = self._stream_handle.readline()   
    
    def extract_particle(self, particle_class, match):
        """
        Extract a particle of the specified class and append it to the record buffer
        @param particle_class: particle class to extract
        @param match: regex match to pass in as raw data
        """
   
        particle = self._extract_sample(particle_class, None, match, None)
        self._record_buffer.append(particle)

    def recov_exception(self, error_message):
        """
        Add a warning log message and use the exception callback to pass a recoverable exception
        @param error_message: The error message to use in the log and callback
        """
        log.warn(error_message)
        self._exception_callback(RecoverableSampleException(error_message))
        
    def _process_metadata_record_part(self, line):
        """
        This function checks to see if a metadata record is contained
        in this chunk.
        """

        match_found = False

        for table_data in METADATA_STATE_TABLE:
            state, matcher, value = table_data
            match = matcher.match(line)

            if match is not None:

                match_found = True

                # Update the state to reflect that we've got
                # this particular metadata record.
                self._metadata_state |= state

                # set the value in the metadata table.
                # VOLTAGE_MATCHER and SAMP_NUM_MEM_MATCHER contain 2 groups. All
                # other metadata matchers contain only one group.
                if matcher != VOLTAGE_MATCHER and matcher != SAMP_NUM_MEM_MATCHER:
                    table_data[METADATA_VALUE_INDEX] = match.group(1)
                else:
                    table_data[METADATA_VALUE_INDEX] = [match.group(1), match.group(2)]            

        if match_found is False:
            error_message = 'Unexpected metadata found: ' + line
            log.warn(error_message)
            self._exception_callback(UnexpectedDataException(error_message))
            
    def _generate_metadata_particles(self):
        """
        This function generates a metadata particle.
        """
  
        log.debug("Metadata state: %s", "{0:b}".format(self._metadata_state))

        if self._metadata_state == ALL_METADATA_RECEIVED:

            # Fields for the metadata particle must be
            # in the same order as the RAW_INDEX_META_xxx values.
            meta_fields = [value
                           for state, matcher, value in METADATA_STATE_TABLE]    
            
            utc_time = formatted_timestamp_utc_time(meta_fields[0], "%Y%m%d %H%M%S")
            ntp_time = ntplib.system_to_ntp_time(utc_time)
            
            # For battery_voltage_main and lithium, what's the correct syntax below? same for samp num and mem.
            metadata_tuple = [
                (CtdmoGhqrImodemDataParticleKey.DATE_TIME_STRING,
                 meta_fields[0],
                 str),
                (CtdmoGhqrImodemDataParticleKey.INSTRUMENT_SERIAL_NUMBER_U32,
                 meta_fields[1],
                 int),
                (CtdmoGhqrImodemDataParticleKey.BATTERY_VOLTAGE_MAIN,
                 meta_fields[2][0],
                 float),
                (CtdmoGhqrImodemDataParticleKey.BATTERY_VOLTAGE_LITHIUM,
                 meta_fields[2][1],
                 float),
                (CtdmoGhqrImodemDataParticleKey.SAMPLE_NUMBER,
                 meta_fields[3][0],
                 int),
                (CtdmoGhqrImodemDataParticleKey.MEM_FREE,
                 meta_fields[3][1],
                 int),
                (CtdmoGhqrImodemDataParticleKey.SAMPLE_INTERVAL,
                 meta_fields[4],
                 int),
                (CtdmoGhqrImodemDataParticleKey.PRESSURE_RANGE,
                 meta_fields[5],
                 int),
                (CtdmoGhqrImodemDataParticleKey.NUM_SAMPLES,
                 meta_fields[6],
                 int)]

            # Generate the metadata particle class and add the
            # result to the list of particles to be returned.
            log.debug('about to extract')
            particle = self._extract_sample(self.metadata_particle_class,
                                            None,
                                            metadata_tuple,
                                            ntp_time)

            if particle is not None:
                self._record_buffer.append(particle)
                return particle

        else:
            error_message = 'Incomplete Metadata'
            log.warn(error_message)
            self._exception_callback(RecoverableSampleException(error_message))
  
          
    def _create_instrument_particle(self, inst_match):
        """
        This method will create an instrument particle given
        instrument match data found from parsing an input file.
        """

        # The first parameter 'time' is the ctd_time converted to ntp64
        time_ntp = generate_particle_timestamp(inst_match.group(4))

        # Create the instrument data list of tuples from the instrument match data
        instrument_data_tuple = [
            (CtdmoGhqrImodemDataParticleKey.TIME,
             time_ntp,
             float),
            (CtdmoGhqrImodemDataParticleKey.TEMPERATURE,
             inst_match.group(1),
             int),
            (CtdmoGhqrImodemDataParticleKey.CONDUCTIVITY,
             inst_match.group(2),
             int),
            (CtdmoGhqrImodemDataParticleKey.PRESSURE,
             inst_match.group(3),
             int),
            (CtdmoGhqrImodemDataParticleKey.CTD_TIME,
             inst_match.group(4),
             int)
        ]

        # Extract the instrument particle sample providing the instrument data
        # tuple and ntp timestamp
        particle = self._extract_sample(self.instrument_particle_class,
                                        None,
                                        instrument_data_tuple,
                                        time_ntp)

        return particle
    
    