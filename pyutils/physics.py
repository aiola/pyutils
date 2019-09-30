""" Physics utilities
"""

import exceptions
import math

class MeasuredQuantity(object):
    """ Class that represents a measured quantiy with errors and units
    """
    def __init__(self, value, error, units=""):
        if error < 0:
            raise exceptions.ValueError()
        if math.isnan(value) or math.isnan(error):
            raise exceptions.ValueError()
        if units == "%":
            value *= 100
            error *= 100
        self.value = value
        self.error = error
        self.units = units
        self.force_exponent = None
        self.precision = float("nan")
        self.exponent = None

    def is_significant(self):
        """ Checks if the uncertainty is smaller than the value
        """
        if self.value == 0:
            return self.error == 0
        return self.error / math.fabs(self.value) < 1.0

    def calculate_precision(self, error=None):
        """ Calculates the number of significant digits from the uncertainty
        """
        if error is None:
            error = self.error
        if error == 0:
            if self.value == 0:
                self.precision = float("nan")
                self.exponent = None
            else:
                v_log_10 = math.log10(abs(self.value))
                if abs(v_log_10) > 3:
                    if v_log_10 > 0:
                        self.exponent = int(math.floor(v_log_10))
                    else:
                        self.exponent = int(math.ceil(v_log_10))
                    self.precision = 6 + self.exponent
                else:
                    self.precision = 6
                    self.exponent = None
        else:
            err_log_10 = math.log10(error)
            # precision is negative for error > 1
            # precision is 0 for error = 1
            # precision is positve for error < 1
            if err_log_10 < 0:
                self.precision = -int(math.floor(err_log_10))
            elif err_log_10 < float('inf'):
                self.precision = -int(math.ceil(err_log_10))
            else:
                self.precision = 0

            if self.precision > 3:
                # error is < 0.001
                self.exponent = -self.precision + 1
            elif self.precision < -3:
                # error is > 1000
                self.exponent = -self.precision - 1
            else:
                self.exponent = None
        if self.force_exponent is not None:
            self.exponent = self.force_exponent

    def to_string(self, plus_minus=r'\pm', times=r'\times'):
        """ String representation
        """
        self.calculate_precision()
        if math.isnan(self.precision):
            return "0"

        if self.exponent is not None:
            prec = abs(self.precision + self.exponent)
            value = 1.0 * self.value / (10**self.exponent)
            error = 1.0 * self.error / (10**self.exponent)
            format_string = "({{value:.{prec}f}} {{pm}} {{error:.{prec}f}}) {{times}} 10^{{{{{{exp}}}}}}".\
                format(prec=prec)
            result = format_string.\
                format(value=value, error=error, exp=self.exponent, pm=plus_minus, times=times)
        else:
            if self.precision < 0:
                prec = 0
            else:
                prec = self.precision
            format_string = "{{value:.{prec}f}} {{pm}} {{error:.{prec}f}}".\
                format(prec=prec)
            result = format_string.\
                format(value=self.value, error=self.error,pm=plus_minus)
        if self.units:
            result += "\, {}".format(self.units)
        return result

    def __str__(self):
        return self.to_string("+/-", "*")

    def __add__(self, other):
        if isinstance(other, MeasuredQuantity):
            if self.units != other.units:
                return NotImplemented
            return MeasuredQuantity(self.value + other.value, math.sqrt(self.error**2 + other.error**2), self.units)
        elif isinstance(other, (int, float)):
            if self.units:
                return NotImplemented
            return MeasuredQuantity(self.value + other, self.error)
        return NotImplemented

    def __neg__(self):
        return MeasuredQuantity(-self.value, self.error, self.units)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, MeasuredQuantity):
            value = self.value * other.value
            error = math.sqrt(self.error**2 / self.value**2 + other.error**2 / other.value**2) * abs(value)
            units = ""
            if units and other.units:
                if units == other.units:
                    if units == "%":
                        value /= 100
                        error /= 100
                    else:
                        units = "({})^2".format(units)
                else:
                    units += " " + other.units
            elif self.units:
                units = self.units
            elif other.units:
                units = other.units
            return MeasuredQuantity(value, error, units)
        elif isinstance(other, (int, float)):
            return MeasuredQuantity(self.value * other, self.error * other, self.units)
        else:
            return NotImplemented

    def __invert__(self):
        if self.units != "%":
            return MeasuredQuantity(1 / self.value, self.error / self.value**2, self.units)
        return NotImplemented

    def __div__(self, other):
        return self * (~other)

    def __lt__(self, other):
        if isinstance(other, MeasuredQuantity):
            return self.value + self.error < other.value - other.error
        elif isinstance(other, (int, float)):
            return self.value + self.error < other
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, MeasuredQuantity):
            return self.value + self.error > other.value - other.error
        elif isinstance(other, (int, float)):
            return self.value + self.error > other
        return NotImplemented

    def ___le__(self, other):
        return not self > other

    def __ge__(self, other):
        return not self < other

    def __ne__(self, other):
        return self > other or self < other

    def __eq__(self, other):
        return not self != other

class MeasuredQuantityAsymmErrors(MeasuredQuantity):
    """ Class that represents a measured quantiy with asymmetric errors and units
    """
    def __init__(self, value, error_up, error_low, units=""):
        if error_up < 0 or error_low < 0:
            raise exceptions.ValueError()
        if math.isnan(value) or math.isnan(error_up) or math.isnan(error_low):
            raise exceptions.ValueError()
        if units == "%":
            value *= 100
            error_up *= 100
            error_low *= 100
        self.value = value
        self.error_up = error_up
        self.error_low = error_low
        self.error = (error_up + error_low) / 2.0
        self.units = units
        self.force_exponent = None
        self.precision = float("nan")
        self.exponent = None
        if max(self.error_low, self.error_up) / min(self.error_low, self.error_up) >= 1.5:
            self.to_string = self.to_string_asymm

    def to_string_asymm(self, _=""):
        """ String representation
        """
        self.calculate_precision(max(self.error_low, self.error_up))
        if math.isnan(self.precision):
            return "0"

        if self.exponent is not None:
            prec = abs(self.precision + self.exponent)
            value = 1.0 * self.value / (10**self.exponent)
            error_up = 1.0 * self.error_up / (10**self.exponent)
            error_low = 1.0 * self.error_low / (10**self.exponent)
            format_string = "({{value:.{prec}f}} + {{error_up:.{prec}f}} - {{error_low:.{prec}f}}) #times 10^{{{{{{exp}}}}}}".\
                format(prec=prec)
            result = format_string.\
                format(value=value, error_up=error_up, error_low=error_low, exp=self.exponent)
        else:
            if self.precision < 0:
                prec = 0
            else:
                prec = self.precision
            format_string = "{{value:.{prec}f}} + {{error_up:.{prec}f}} - {{error_low:.{prec}f}}".\
                format(prec=prec)
            result = format_string.\
                format(value=self.value, error_up=self.error_up, error_low=self.error_low)
        if self.units:
            result += " {}".format(self.units)
        return result
