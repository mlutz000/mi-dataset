#!/usr/bin/env python

"""
@file coi-services/mi/idk/result_set.py
@author Bill French
@brief Read a result set file and use the data to verify
data particles.

Usage:

from mi.core.log import log

rs = ResultSet(result_set_file_path)
if not rs.verify(particles):
    log.info("Particle verified")
else:
    log.error("Particle validate failed")
    log.error(rs.report())

Result Set File Format:
  result files are yml formatted files with a header and data section.
  the data is stored in record elements with the key being the parameter name.
     - two special fields are internal_timestamp and _index.
     - internal timestamp can be input in text string or ntp float format

eg.

# Result data for verifying particles. Comments are ignored.

header:
  particle_object: CtdpfParserDataParticleKey
  particle_type: ctdpf_parsed

data:
  -  _index: 1
     internal_timestamp: 07/26/2013 21:01:03
     temperature: 4.1870
     conductivity: 10.5914
     pressure: 161.06
     oxygen: 2693.0
  -  _index: 2
     internal_timestamp: 07/26/2013 21:01:04
     temperature: 4.1872
     conductivity: 10.5414
     pressure: 161.16
     oxygen: 2693.1

If a driver returns multiple particle types, the particle type must be specified in each particle

header:
  particle_object: 'MULTIPLE'
  particle_type: 'MULTIPLE'

data:
  -  _index: 1
     particle_object: CtdpfParser1DataParticleKey
     particle_type: ctdpf_parsed_1
     internal_timestamp: 07/26/2013 21:01:03
     temperature: 4.1870
     conductivity: 10.5914
     pressure: 161.06
     oxygen: 2693.0
  -  _index: 2
     particle_object: CtdpfParser2DataParticleKey
     particle_type: ctdpf_parsed_2
     internal_timestamp: 07/26/2013 21:01:04
     temperature: 4.1872
     conductivity: 10.5414
     pressure: 161.16
     oxygen: 2693.1


"""

__author__ = 'Bill French'
__license__ = 'Apache 2.0'

from datetime import datetime
import calendar
import re
import yaml
import ntplib
import numpy

from mi.core.instrument.data_particle import DataParticle

from mi.core.log import get_logger ; log = get_logger()

DATE_PATTERN_WITH_MSEC_WITH_Z = r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)Z$'
DATE_PATTERN_WITH_MSEC_WITHOUT_Z = r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)$'
DATE_PATTERN_WITHOUT_MSEC_WITH_Z = r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$'
DATE_PATTERN_WITHOUT_MSEC_WITHOUT_Z = r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}$'
DATE_FORMAT = "%Y-%m-%dT%H:%M:%S.%fZ"

FLOAT_ALLOWED_DIFF = .000001


class ResultSet(object):
    """
    Result Set object
    Read result set files and compare to parsed particles.
    """
    def __init__(self, result_file_path):
        self.yaml = dict()

        log.debug("read result file: %s" % result_file_path)
        stream = file(result_file_path, 'r')
        result_set = yaml.load(stream)

        self._set_result_set(result_set)

        self._clear_report()

    def verify(self, particles):
        """
        Verify particles passed in against result set read
        in the ctor.

        Ensure:
          - Verify particles as a set
          - Verify individual particle data

        store verification result in the object and
        return success or failure.
        @param particles: list of particles to verify.
        @return True if verification successful, False otherwise
        """
        self._clear_report()
        result = True

        if self._verify_set(particles):
            result = self._verify_particles(particles)
        else:
            result = False

        if not result:
            log.error("Failed verification: \n%s", self.report())

        return result

    def report(self):
        """
        Return an ascii formatted verification failure report.
        @return string report
        """
        if len(self._report):
            return "\n".join(self._report)
        else:
            return None

    ###
    #   Helpers
    ###
    def _add_to_report(self, messages, indent = 0):
        """
        Add a message to the report buffer, pass an indent factor to
        indent message in the ascii report.
        """
        if not isinstance(messages, list): messages = [messages]

        for message in messages:
            ind = ""
            for i in range(0, indent):
                ind += "    "
            self._report.append("%s%s" %(ind, message))
            log.warn(message)

    def _clear_report(self):
        """
        Add a message to the report buffer, pass an indent factor to
        indent message in the ascii report.
        """
        self._report = []

    def _set_result_set(self, result_set):
        """
        Take data from yaml file and store it in internal objects for
        verifying data.  Raise an exception on error.
        """
        log.trace("Parsing result set header: %s", result_set)

        self._result_set_header = result_set.get("header")
        if not self._result_set_header: raise IOError("Missing result set header")
        log.trace("Header: %s", self._result_set_header)

        if self._result_set_header.get("particle_object") is None:
            IOError("header.particle_object not defined")

        if self._result_set_header.get("particle_type") is None:
            IOError("header.particle_type not defined")

        self._result_set_data = {}
        data = result_set.get("data")
        if not data: raise IOError("Missing result set data")

        for particle in data:
            index = particle.get("_index")
            if index is None:
                log.error("Particle definition missing _index: %s", particle)
                raise IOError("Particle definition missing _index")

            if self._result_set_data.get(index) is not None:
                log.error("Duplicate particle definition for _index %s: %s", index, particle)
                raise IOError("Duplicate definition found for index: %s"% index)

            self._result_set_data[index] = particle
            log.trace("Result set data: %s", self._result_set_data)

    def _verify_set(self, particles):
        """
        Verify the particles as a set match what we expect.
        - All particles are of the expected type
        - Check particle count
        """
        errors = []

        if len(self._result_set_data) != len(particles):
            errors.append("result set records != particles to verify (%d != %d)" %
                          (len(self._result_set_data), len(particles)))

        # if this driver returns multiple particle classes, type checking happens
        # for each particle in _get_particle_data_errors
        if self._result_set_header.get("particle_object") != 'MULTIPLE' and \
        self._result_set_header.get("particle_type") != 'MULTIPLE':
            for particle in particles:
                if not self._verify_particle_type(particle):
                    log.error("particle type mismatch: %s", particle)
                    errors.append('particle type mismatch')

        if len(errors):
            self._add_to_report("Header verification failure")
            self._add_to_report(errors, 1)
            return False

        return True

    def _verify_particles(self, particles):
        """
        Verify data in the particles individually.
        - Verify order based on _index
        - Verify parameter data values
        - Verify there are extra or missing parameters
        """
        result = True
        index = 1
        for particle in particles:
            particle_def = self._result_set_data.get(index)
            errors = []

            # No particle definition, we fail
            if particle_def is None:
                errors.append("no particle result defined for index %d" % index)

            # Otherwise lets do some validation
            else:
                errors += self._get_particle_header_errors(particle, particle_def)
                errors += self._get_particle_data_errors(particle, particle_def)

            if len(errors):
                self._add_to_report("Failed particle validation for index %d" % index)
                self._add_to_report(errors, 1)
                result = False

            index += 1

        return result

    def _verify_particle_type(self, particle):
        """
        Verify that the object is a DataParticle and is the
        correct type.
        """
        if isinstance(particle, dict):
            return True

        expected = self._result_set_header['particle_object']

        cls = particle.__class__.__name__

        if not issubclass(particle.__class__, DataParticle):
            log.error("type not a data particle")

        if expected != cls:
            log.error("type mismatch: %s != %s", expected, cls)
            return False

        return True

    def _get_particle_header_errors(self, particle, particle_def):
        """
        Verify all parameters defined in the header:
        - Stream type
        - Internal timestamp
        """
        errors = []
        particle_dict = self._particle_as_dict(particle)
        particle_timestamp = particle_dict.get('internal_timestamp')
        expected_time = particle_def.get('internal_timestamp')

        # Verify the timestamp
        if particle_timestamp and not expected_time:
            errors.append("particle_timestamp defined in particle, but not expected")
        elif not particle_timestamp and expected_time:
            errors.append("particle_timestamp expected, but not defined in particle")

        elif particle_timestamp:
            if isinstance(expected_time, str):
                expected = self._string_to_ntp_date_time(expected_time)
            else:
                # if not a string, timestamp should alread be in ntp
                expected = expected_time
            ts_diff = abs(particle_timestamp - expected)
            log.debug("verify timestamp: abs(%s - %s) = %s", expected, particle_timestamp, ts_diff)

            if ts_diff > FLOAT_ALLOWED_DIFF:
                errors.append("expected internal_timestamp mismatch, %.9f != %.9f (%.9f)" %
                              (expected, particle_timestamp, ts_diff))

        # verify the stream name, unless multiple are returned, type checking is done
        # in get_particle_data_errors if so
        particle_stream = particle_dict['stream_name']
        if self._result_set_header['particle_type'] != 'MULTIPLE':
            expected_stream = self._result_set_header['particle_type']
            if particle_stream != expected_stream:
                errors.append("expected stream name mismatch: %s != %s" %
                              (expected_stream, particle_stream))

        return errors

    def _get_particle_data_errors(self, particle, particle_def):
        """
        Verify that all data parameters are present and have the
        expected value
        """
        errors = []
        particle_dict = self._particle_as_dict(particle)
        log.debug("Particle to test: %s", particle_dict)
        log.debug("Particle definition: %s", particle_def)
        particle_values = particle_dict['values']

        # particle object and particle type keys will only be present for drivers
        # returning multiple particle types
        if 'particle_object' in particle_def:
            expected_object = particle_def.get('particle_object')
            expected_type = particle_def.get('particle_type', None)

            # particle is either a class or dictionary, if it is a
            # dictionary there is no class to compare
            if not isinstance(particle, dict):
                # particle is an actual class, check that the class matches
                cls = particle.__class__.__name__
                if not issubclass(particle.__class__, DataParticle):
                    errors.append("Particle class %s is not a subclass of DataParticle" %
                                  particle.__class__)

                if expected_object != cls:
                    errors.append("Class mismatch, expected: %s, received: %s" %
                                  (expected_object, cls))

            particle_stream = particle_dict['stream_name']
            if particle_stream != expected_type:
                log.debug("Stream type mismatch, expected: %s, received: %s" % (expected_type, particle_stream))
                errors.append("Stream type mismatch, expected: %s, received: %s" % (expected_type, particle_stream))

        expected_keys = []
        for (key, value) in particle_def.items():
            if(key not in ['_index', '_new_sequence', 'internal_timestamp', 'particle_object', 'particle_type']):
                expected_keys.append(key)

        particle_keys = []
        pv = {}
        for value in particle_values:
            particle_keys.append(value['value_id'])
            pv[value['value_id']] = value['value']

        if sorted(expected_keys) != sorted(particle_keys):
            errors.append("expected / particle keys mismatch: %s != %s" %
                          (sorted(expected_keys), sorted(particle_keys)))

        else:
            for key in expected_keys:
                expected_value = particle_def[key]
                particle_value = pv[key]
                log.debug("Verify value for '%s'", key)
                e = ResultSet._verify_value(expected_value, particle_value)
                if e:
                    errors.append("'%s' %s"  % (key, e))

        return errors

    @staticmethod
    def _verify_value(expected_value, particle_value):
        """
        Verify a value matches what we expect.  If the expected value (from the yaml)
        is a dict then we expect the value to be in a 'value' field.  Otherwise just
        use the parameter as a raw value.

        when passing a dict you can specify a 'round' factor.
        """

        local_expected_value = expected_value
        round_factor = None

        # Let's first check to see if we have a None for an expected value
        if local_expected_value is None:
            log.debug("No value to compare, ignoring")
            return None

        if isinstance(expected_value, dict):
            local_expected_value = expected_value['value']
            round_factor = expected_value.get('round')

        if isinstance(local_expected_value, float):

            if round_factor is None:
                # unless otherwise specified round all floating point expected values to 5 digits
                round_factor = 5

            if particle_value is None:
                return "value mismatch, expected value is float and particle value is None"
            else:
                if isinstance(particle_value, float):
                    particle_value = round(particle_value, round_factor)

                    log.trace("particle value (rounded) to %s", particle_value)
                    log.trace("expected value %s", local_expected_value)

                    diff = abs(particle_value - local_expected_value)

                    if diff > FLOAT_ALLOWED_DIFF:
                        return "value mismatch, %.6f != %.6f " % \
                               (local_expected_value, particle_value)
                else:
                    return "value mismatch, expected value is float and particle value is %s" \
                           % type(particle_value)

        else:
            if local_expected_value != particle_value:
                # check for nans, two nans will not equal each other but in this
                # case they are considered equal
                if ResultSet._nan_equal_compare(local_expected_value, particle_value):
                    return None

                return "value mismatch, %s != %s" % \
                       (local_expected_value, particle_value)

        return None

    @staticmethod
    def _nan_equal_compare(expected, received):
        """
        Compare the expected and recieved values, considering two NaNs to be equal.
        If the values are equal True is returned, if they are not False is returned
        @param expected: Expected value
        @param received: Received value
        @return: True if expected and received match including NaNs
        """
        # check for a list
        if ResultSet._are_both_lists(expected, received):

            # check for list within a list
            if ResultSet._are_both_lists(expected[0], received[0]):
                # confirm lists are equal size and contain numbers, if not they
                # don't need checking for nans
                if len(expected) == len(received) and \
                        ResultSet._are_both_numbers(expected[0][0], received[0][0]):
                    # loop and compare each list within a list
                    all_lists_equal = True
                    for i in range(0, len(expected)):
                        if not ResultSet._is_equal_nan_list(expected[i], received[i]):
                            all_lists_equal = False
                            break
                    if all_lists_equal:
                        return True

            # check for a list of numbers
            elif ResultSet._are_both_numbers(expected[0], received[0]):
                # this is a single list of float or ints
                if ResultSet._is_equal_nan_list(expected, received):
                    return True

        # check two individual float or ints for nans
        elif ResultSet._are_both_numbers(expected, received):
            if numpy.isnan(expected) and numpy.isnan(received):
                # found two Nans, they match
                return True
        return False

    @staticmethod
    def _are_both_lists(expected, received):
        """
        Compare if the expected and recieved values are both lists
        @param expected: expected value
        @param received: received value
        @return: True if both lists, false if not
        """
        if isinstance(expected, list) and isinstance(received, list):
            return True
        return False

    @staticmethod
    def _are_both_numbers(expected, received):
        """
        Compare if the expected and recieved values are both floats or ints
        @param expected: expected value
        @param received: received value
        @return: True if both lists, false if not
        """
        if isinstance(expected, (float, int)) and isinstance(received, (float, int)):
            return True
        return False

    @staticmethod
    def _is_equal_nan_list(expected_list, received_list):
        """
        Compare two lists that contain nan values, return True if they match,
        False if they do not, considering two NaNs to be equal in values
        @param expected_list: a list of expected values
        @param received_list: the recieved list of values
        @return: True if lists match, False if lists don't match
        """
        # first check list length, if this doesn't match the lists don't match
        if len(expected_list) == len(received_list):
            # get numpy array of True / False matching nan
            isnan_ex_value = numpy.isnan(expected_list)
            isnan_particle_value = numpy.isnan(received_list)
            # make a new list containing only non-nan values
            ex_val_non_nan = []
            particle_non_nan = []
            for i in range(0, len(expected_list)):
                if not isnan_ex_value[i]:
                    ex_val_non_nan.append(expected_list[i])
                if not isnan_particle_value[i]:
                    particle_non_nan.append(received_list[i])
            if ex_val_non_nan == particle_non_nan:
                # all non-nan values match
                return True
        return False

    def _string_to_ntp_date_time(self, datestr):
        """
        Extract an ntp date from a ISO8601 formatted date string.
        @param str an ISO8601 formatted string containing date information
        @retval an ntp date number (seconds since jan 1 1900)
        @throws InstrumentParameterException if datestr cannot be formatted to
        a date.
        """

        if not isinstance(datestr, str):
            raise IOError('Value %s is not a string.' % str(datestr))

        date_match_with_msec_with_z = re.match(DATE_PATTERN_WITH_MSEC_WITH_Z, datestr)
        date_match_with_msec_without_z = re.match(DATE_PATTERN_WITH_MSEC_WITHOUT_Z, datestr)
        date_match_without_msec_with_z = re.match(DATE_PATTERN_WITHOUT_MSEC_WITH_Z, datestr)
        date_match_without_msec_without_z = re.match(DATE_PATTERN_WITHOUT_MSEC_WITHOUT_Z, datestr)

        datestr_to_use = ""

        log.debug("Input datestr: %s", datestr)

        if date_match_without_msec_without_z:

            datestr_to_use = datestr + '.0Z'

        elif date_match_without_msec_with_z:

            datestr_to_use = datestr[:-1] + '.0Z'

        elif date_match_with_msec_without_z:

            datestr_to_use = datestr + 'Z'

        elif date_match_with_msec_with_z:

            datestr_to_use = datestr

        else:
            raise ValueError("NTP date string not in any of the expected formats")

        try:

            log.debug("converting time string '%s'", datestr_to_use)

            dt = datetime.strptime(datestr_to_use, DATE_FORMAT)

            unix_timestamp = calendar.timegm(dt.timetuple()) + (dt.microsecond / 1000000.0)

            timestamp = ntplib.system_to_ntp_time(unix_timestamp)

            log.debug("converted time string '%s', unix_ts: %s ntp: %s", datestr_to_use, unix_timestamp, timestamp)

        except ValueError as e:
            raise ValueError('Value %s could not be formatted to a date. %s' % (str(datestr_to_use), e))

        return timestamp

    def _particle_as_dict(self, particle):
        if isinstance(particle, dict):
            return particle

        return particle.generate_dict()
